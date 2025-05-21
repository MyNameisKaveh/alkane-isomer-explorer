import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage
import gradio as gr
# import io # دیگر مستقیما به io نیاز نداریم

# تابع برای رسم مولکول از SMILES و بازگرداندن آن به عنوان تصویر PIL
def draw_molecule(smiles_string):
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if mol:
            img = MolToImage(mol, size=(300, 300)) # خروجی یک PIL Image است
            return img
        else:
            print(f"Could not parse SMILES: {smiles_string}")
            # برگرداندن یک تصویر خالی یا یک placeholder در صورت خطا می‌تواند مفید باشد
            # اما فعلا None برمی‌گردانیم
            return None
    except Exception as e:
        print(f"Error drawing molecule for SMILES {smiles_string}: {e}")
        return None

def find_and_display_isomers(molecule_name):
    """
    تابع اصلی که توسط Gradio فراخوانی می‌شود.
    نام مولکول را گرفته، ایزومرها را پیدا کرده و لیستی از نام‌ها و تصاویر را برمی‌گرداند.
    """
    if not molecule_name or not molecule_name.strip():
        return [], "لطفا نام یک مولکول را وارد کنید."

    print(f"Processing request for: {molecule_name}")
    
    isomer_outputs = []
    status_message = ""

    try:
        # 1. جستجوی مولکول اصلی برای یافتن فرمول مولکولی
        print(f"Searching for: {molecule_name} in PubChem...")
        # استفاده از listkey_count=1 برای گرفتن فقط نتیجه اول برای مولکول اصلی
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
        # listkey_count برای کنترل تعداد نتایج (برای آلکان‌های بزرگتر ممکن است نیاز به افزایش باشد)
        isomers = pcp.get_compounds(molecular_formula, 'formula', listkey_count=20) 

        if not isomers:
            status_message = f"ایزومری برای فرمول {molecular_formula} یافت نشد."
            print(status_message)
            return [], status_message

        print(f"Found {len(isomers)} potential isomers.")
        
        processed_cids = set() # برای جلوگیری از نمایش ایزومرهای تکراری (اگر PubChem برگرداند)

        for isomer in isomers:
            if isomer.cid in processed_cids:
                continue # از این ایزومر قبلا پردازش شده، رد شو
            
            smiles = isomer.canonical_smiles
            iupac_name = isomer.iupac_name
            
            if not iupac_name and isomer.synonyms:
                # سعی می‌کنیم نامی را انتخاب کنیم که به نظر آلکان می‌رسد یا ساده‌تر است
                # این بخش می‌تواند هوشمندتر شود
                chosen_synonym = isomer.synonyms[0]
                for syn in isomer.synonyms:
                    if "alkane" in syn.lower() or "iso" in syn.lower() or "neo" in syn.lower() or "ane" in syn.lower()[-3:]:
                        chosen_synonym = syn
                        break
                iupac_name = chosen_synonym
            elif not iupac_name:
                iupac_name = f"N/A (CID: {isomer.cid})"

            if smiles:
                mol_image = draw_molecule(smiles)
                if mol_image: # فقط اگر تصویر با موفقیت رسم شد اضافه کن
                    isomer_outputs.append((mol_image, f"{iupac_name}\nSMILES: {smiles}"))
                    processed_cids.add(isomer.cid)
            
        if not isomer_outputs:
            status_message = "ایزومرها پیدا شدند اما تصویری برای نمایش وجود ندارد (ممکن است خطایی در رسم رخ داده باشد)."
            print(status_message)
        elif not status_message: # اگر پیام خطای دیگری از قبل نداشتیم
             status_message = f"{len(isomer_outputs)} ایزومر برای '{molecule_name}' (فرمول: {molecular_formula}) پیدا و نمایش داده شد."
        
        # مرتب‌سازی بر اساس نام (کپشن تصویر)
        isomer_outputs.sort(key=lambda x: x[1])
        
        return isomer_outputs, status_message

    except pcp.PubChemHTTPError as e:
        error_msg = f"خطا در ارتباط با PubChem: {e}. ممکن است سرویس موقتا در دسترس نباشد یا درخواست شما بیش از حد بزرگ باشد."
        print(error_msg)
        return [], error_msg
    except Exception as e:
        error_msg = f"یک خطای غیرمنتظره در سرور رخ داد: {e}"
        print(error_msg)
        return [], error_msg

# ساخت رابط کاربری Gradio
iface = gr.Interface(
    fn=find_and_display_isomers,
    inputs=gr.Textbox(
        label="نام آلکان را وارد کنید", 
        placeholder="مثال: butane, pentane, hexane",
        info="نام آلکان مورد نظر خود را به انگلیسی وارد کنید."
    ),
    outputs=[
        gr.Gallery(
            label="ایزومرهای یافت شده", 
            columns=[3], # تعداد ستون‌ها برای نمایش تصاویر
            height="auto", 
            object_fit="contain" # نحوه نمایش تصویر در کادر
        ),
        gr.Textbox(label="وضعیت و پیام‌ها")
    ],
    title="یابنده و نمایشگر ایزومرهای آلکان",
    description=(
        "نام یک آلکان (به انگلیسی) را وارد کنید تا ایزومرهای آن به همراه ساختار شیمیایی و نام IUPAC (یا رایج) نمایش داده شوند.\n"
        "اطلاعات از دیتابیس PubChem دریافت شده و ساختارها با استفاده از کتابخانه RDKit رسم می‌شوند."
    ),
    examples=[
        ["butane"], 
        ["pentane"], 
        ["hexane"],
        ["heptane"] # هپتان هم تعداد ایزومر قابل قبولی دارد
    ],
    allow_flagging='never', # غیرفعال کردن دکمه Flag
    theme=gr.themes.Soft() # استفاده از یک تم آماده برای ظاهر بهتر (اختیاری)
)

# برای اجرا در Hugging Face Spaces، این بخش کافی است.
# هاگینگ فیس خودش iface.launch() را صدا می‌زند.
if __name__ == '__main__':
    iface.launch() # share=False به طور پیش‌فرض، برای تست محلی
