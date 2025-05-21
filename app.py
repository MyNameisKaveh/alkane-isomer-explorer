import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage
# from PIL import Image # دیگر مستقیما به Image از PIL نیاز نداریم، Gradio خودش مدیریت می‌کند
import gradio as gr
import io # برای کار با بایت‌های تصویر در صورت نیاز، اما Gradio معمولا خودش هندل می‌کند

# تابع برای رسم مولکول از SMILES و بازگرداندن آن به عنوان تصویر PIL
def draw_molecule(smiles_string):
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if mol:
            img = MolToImage(mol, size=(300, 300)) # خروجی یک PIL Image است
            return img
        else:
            print(f"Could not parse SMILES: {smiles_string}")
            return None # یا یک تصویر خالی/پیش‌فرض
    except Exception as e:
        print(f"Error drawing molecule for SMILES {smiles_string}: {e}")
        return None # یا یک تصویر خالی/پیش‌فرض

def find_and_display_isomers(molecule_name):
    """
    تابع اصلی که توسط Gradio فراخوانی می‌شود.
    نام مولکول را گرفته، ایزومرها را پیدا کرده و لیستی از نام‌ها و تصاویر را برمی‌گرداند.
    """
    if not molecule_name.strip():
        return [], "لطفا نام یک مولکول را وارد کنید."

    print(f"Processing request for: {molecule_name}")
    
    isomer_outputs = []
    status_message = ""

    try:
        # 1. جستجوی مولکول اصلی برای یافتن فرمول مولکولی
        print(f"Searching for: {molecule_name} in PubChem...")
        compounds = pcp.get_compounds(molecule_name, 'name', record_type='3d', listkey_count=1)
        if not compounds:
            status_message = f"مولکول '{molecule_name}' در PubChem یافت نشد."
            print(status_message)
            return [], status_message
        
        main_compound = compounds[0]
        molecular_formula = main_compound.molecular_formula
        print(f"Found {molecule_name} with formula: {molecular_formula}")

        if not molecular_formula:
            status_message = f"فرمول مولکولی برای '{molecule_name}' تعیین نشد."
            print(status_message)
            return [], status_message

        # 2. جستجوی ترکیبات با همان فرمول مولکولی (ایزومرها)
        print(f"Searching for isomers with formula: {molecular_formula}...")
        isomers = pcp.get_compounds(molecular_formula, 'formula', listkey_count=20) # محدود کردن اولیه برای جلوگیری از نتایج خیلی زیاد

        if not isomers:
            status_message = f"ایزومری برای فرمول {molecular_formula} یافت نشد."
            print(status_message)
            return [], status_message

        print(f"Found {len(isomers)} potential isomers.")
        
        # برای جلوگیری از نمایش بیش از حد ایزومرها، می‌توانیم یک محدودیت بگذاریم
        # max_isomers_to_display = 10 
        # displayed_isomers = 0

        for isomer in isomers:
            # if displayed_isomers >= max_isomers_to_display:
            #     status_message += f"\n(نمایش به {max_isomers_to_display} ایزومر اول محدود شد)"
            #     break

            smiles = isomer.canonical_smiles
            iupac_name = isomer.iupac_name
            
            if not iupac_name and isomer.synonyms:
                iupac_name = isomer.synonyms[0]
            elif not iupac_name:
                iupac_name = f"N/A (CID: {isomer.cid})"

            if smiles:
                mol_image = draw_molecule(smiles)
                if mol_image:
                    # Gradio می‌تواند لیست‌هایی از تاپل‌ها (یا لیست‌های تو در تو) را برای Gallery یا Outputهای متعدد بپذیرد.
                    # اینجا ما یک لیست از (تصویر، نام) را برای هر ایزومر برمی‌گردانیم.
                    # یا می‌توانیم فرمت خروجی را طوری تنظیم کنیم که Gradio بهتر بفهمد.
                    # برای سادگی، فعلا یک رشته حاوی نام و یک تصویر جدا برمی‌گردانیم.
                    # Gradio با `gr.Gallery` بهتر کار می‌کند اگر لیستی از (تصویر، کپشن) بدهیم.
                    isomer_outputs.append((mol_image, f"{iupac_name}\nSMILES: {smiles}"))
                    # displayed_isomers += 1
            
        if not isomer_outputs:
            status_message = "ایزومرها پیدا شدند اما تصویری برای نمایش وجود ندارد."
            print(status_message)
        elif not status_message: # اگر پیام خطای دیگری از قبل نداشتیم
             status_message = f"{len(isomer_outputs)} ایزومر برای '{molecule_name}' (فرمول: {molecular_formula}) پیدا و نمایش داده شد."
        
        # مرتب‌سازی بر اساس نام (کپشن تصویر)
        isomer_outputs.sort(key=lambda x: x[1])
        
        return isomer_outputs, status_message


    except pcp.PubChemHTTPError as e:
        error_msg = f"خطا در ارتباط با PubChem: {e}"
        print(error_msg)
        return [], error_msg
    except Exception as e:
        error_msg = f"یک خطای غیرمنتظره رخ داد: {e}"
        print(error_msg)
        return [], error_msg

# ساخت رابط کاربری Gradio
# ما می‌خواهیم یک گالری از تصاویر به همراه نامشان نمایش دهیم.
# و یک پیام وضعیت.
iface = gr.Interface(
    fn=find_and_display_isomers,
    inputs=gr.Textbox(label="نام آلکان را وارد کنید (مثلا: pentane, butane)", placeholder="مثال: hexane"),
    outputs=[
        gr.Gallery(label="ایزومرهای یافت شده", પ્રકાર="pil", columns=[3], height="auto", object_fit="contain"),
        gr.Textbox(label="وضعیت")
    ],
    title="یابنده و نمایشگر ایزومرهای آلکان",
    description="نام یک آلکان را وارد کنید تا ایزومرهای آن به همراه ساختارشان نمایش داده شوند. اطلاعات از PubChem و با استفاده از RDKit رسم می‌شوند.",
    examples=[["butane"], ["pentane"], ["hexane"]],
    allow_flagging='never' # برای جلوگیری از نمایش دکمه Flag در هاگینگ فیس
)

# برای اجرا در Hugging Face Spaces، باید متغیر `app` را تعریف کنیم.
# اما چون ما می‌خواهیم مستقیما `iface.launch()` را صدا بزنیم (برای تست محلی و دیپلوی مستقیم Gradio)،
# می‌توانیم این بخش را به این صورت بگذاریم.
# هاگینگ فیس خودش فایل app.py را پیدا کرده و iface.launch() را اجرا می‌کند.

if __name__ == '__main__':
    iface.launch(share=False) # share=True برای ایجاد لینک عمومی موقت (اگر محلی اجرا می‌کنید)
                               # در هاگینگ فیس نیازی به share=True نیست.
