import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage
import gradio as gr
import traceback # برای لاگ‌گیری کامل‌تر خطاها

# تابع برای رسم مولکول از SMILES و بازگرداندن آن به عنوان تصویر PIL
def draw_molecule(smiles_string):
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if mol:
            # برای نمایش بهتر اتم‌های هیدروژن در آلکان‌های کوچک (اختیاری)
            # Chem.AddHs(mol) 
            img = MolToImage(mol, size=(300, 300))
            return img
        else:
            print(f"Could not parse SMILES: {smiles_string}")
            return None
    except Exception as e:
        print(f"Error drawing molecule for SMILES {smiles_string}: {e}")
        return None

def find_and_display_isomers(molecule_name_input):
    """
    تابع اصلی که توسط Gradio فراخوانی می‌شود.
    نام مولکول را گرفته، ایزومرها را پیدا کرده و لیستی از نام‌ها و تصاویر را برمی‌گرداند.
    """
    if not molecule_name_input or not molecule_name_input.strip():
        return [], "لطفا نام یک مولکول را وارد کنید."

    molecule_name = molecule_name_input.strip().lower() # تبدیل به حروف کوچک و حذف فضاهای اضافی
    print(f"Processing request for: '{molecule_name}'")
    
    isomer_outputs = []
    status_message = ""

    try:
        # 1. جستجوی مولکول اصلی برای یافتن فرمول مولکولی
        print(f"Searching for exact compound: '{molecule_name}' in PubChem...")
        # استفاده از searchtype و بررسی چندین نتیجه برای یافتن فرمول معتبر
        compounds = pcp.get_compounds(molecule_name, 'name', searchtype='superstructure', listkey_count=5) 
        
        main_compound_obj = None
        molecular_formula = None

        if not compounds:
            status_message = f"مولکول '{molecule_name}' در PubChem یافت نشد."
            print(status_message)
            return [], status_message
        
        print(f"Found {len(compounds)} potential matches for '{molecule_name}'. Checking them...")
        for i, c in enumerate(compounds):
            cid = c.cid
            # سعی در گرفتن نام رایج یا اولین مترادف
            common_name = "N/A"
            if c.synonyms:
                common_name = c.synonyms[0]
            
            formula_attr = hasattr(c, 'molecular_formula')
            actual_formula = c.molecular_formula if formula_attr else None
            
            print(f"  Checking match {i+1}: CID {cid}, Name: '{common_name}', HasFormulaAttr: {formula_attr}, Formula: '{actual_formula}'")

            if actual_formula: # اگر فرمول مولکولی وجود داشت
                # اختیاری: بررسی اینکه آیا این یک هیدروکربن ساده است
                is_simple_hydrocarbon = False
                if c.canonical_smiles:
                    try:
                        mol_obj = Chem.MolFromSmiles(c.canonical_smiles)
                        if mol_obj:
                            elements = {atom.GetSymbol() for atom in mol_obj.GetAtoms()}
                            if elements.issubset({'C', 'H'}):
                                is_simple_hydrocarbon = True
                                print(f"    CID {cid} is a simple hydrocarbon.")
                            else:
                                print(f"    CID {cid} contains elements: {elements}. Not a simple hydrocarbon.")
                        else:
                             print(f"    CID {cid} SMILES '{c.canonical_smiles}' could not be parsed by RDKit.")
                    except Exception as rdkit_ex:
                        print(f"    RDKit error processing SMILES for CID {cid}: {rdkit_ex}")
                else:
                    print(f"    CID {cid} has no SMILES string for element check.")


                # اگر نام ورودی دقیقا یکی از مترادف‌ها باشد و فرمول داشته باشد، اولویت دارد
                if molecule_name in [syn.lower() for syn in c.synonyms]:
                    if is_simple_hydrocarbon or not c.canonical_smiles: # اگر هیدروکربن ساده است یا نمی‌توانیم بررسی کنیم
                        main_compound_obj = c
                        molecular_formula = actual_formula
                        print(f"  Selected main compound (exact name match & hydrocarbon): CID {main_compound_obj.cid}, Name: '{common_name}', Formula: {molecular_formula}")
                        break 
                
                # اگر هنوز main_compound_obj پیدا نشده و این یک هیدروکربن ساده است
                if not main_compound_obj and (is_simple_hydrocarbon or not c.canonical_smiles):
                    main_compound_obj = c
                    molecular_formula = actual_formula
                    print(f"  Tentatively selected main compound (hydrocarbon): CID {main_compound_obj.cid}, Name: '{common_name}', Formula: {molecular_formula}")
                    # ادامه می‌دهیم تا اگر مورد بهتری با تطابق نام پیدا شد، آن را انتخاب کنیم

        # اگر پس از حلقه، هنوز main_compound_obj بر اساس تطابق نام انتخاب نشده بود، اولین هیدروکربن یافت شده را استفاده کن
        if not molecular_formula and main_compound_obj: # این حالت یعنی در بالا انتخاب موقت شده بود
             molecular_formula = main_compound_obj.molecular_formula
             print(f"  Confirmed selection of main compound: CID {main_compound_obj.cid}, Formula: {molecular_formula}")


        if not main_compound_obj or not molecular_formula:
            status_message = f"فرمول مولکولی معتبری برای '{molecule_name}' در نتایج اولیه PubChem یافت نشد. لطفا نام دقیق‌تری وارد کنید یا مطمئن شوید که یک آلکان است."
            print(status_message)
            if compounds: # نمایش اطلاعات بیشتر اگر نتایجی بود اما مناسب نبودند
                 print(f"  PubChem results for '{molecule_name}' that were checked:")
                 for i, comp_item in enumerate(compounds):
                     print(f"    {i+1}. CID: {comp_item.cid}, Synonyms: {comp_item.synonyms[:2] if comp_item.synonyms else 'N/A'}, Formula: {hasattr(comp_item, 'molecular_formula') and comp_item.molecular_formula}")
            return [], status_message
        
        # 2. جستجوی ترکیبات با همان فرمول مولکولی (ایزومرها)
        print(f"Searching for isomers with formula: {molecular_formula}...")
        # افزایش listkey_count برای اطمینان از گرفتن همه ایزومرهای آلکان‌های رایج
        isomers_found_raw = pcp.get_compounds(molecular_formula, 'formula', listkey_count=50) 

        if not isomers_found_raw:
            status_message = f"ایزومری برای فرمول {molecular_formula} (مربوط به '{molecule_name}') یافت نشد."
            print(status_message)
            return [], status_message

        print(f"Found {len(isomers_found_raw)} potential isomer entries from PubChem.")
        
        processed_cids = set() 
        valid_isomers_count = 0

        for isomer_entry in isomers_found_raw:
            if isomer_entry.cid in processed_cids:
                print(f"  Skipping duplicate CID: {isomer_entry.cid}")
                continue 
            
            # اطمینان از اینکه ایزومر هم یک هیدروکربن ساده است (برای آلکان‌ها)
            if isomer_entry.canonical_smiles:
                try:
                    mol_iso = Chem.MolFromSmiles(isomer_entry.canonical_smiles)
                    if mol_iso:
                        elements_iso = {atom.GetSymbol() for atom in mol_iso.GetAtoms()}
                        if not elements_iso.issubset({'C', 'H'}):
                            print(f"  Skipping non-hydrocarbon isomer: CID {isomer_entry.cid}, Name: {isomer_entry.iupac_name or isomer_entry.synonyms[0] if isomer_entry.synonyms else ''}, Elements: {elements_iso}")
                            continue
                    else: # SMILES نامعتبر برای ایزومر
                        print(f"  Skipping isomer with invalid SMILES: CID {isomer_entry.cid}")
                        continue
                except Exception as rdkit_iso_ex:
                    print(f"  RDKit error for isomer SMILES CID {isomer_entry.cid}: {rdkit_iso_ex}")
                    continue
            else: # ایزومر بدون SMILES
                print(f"  Skipping isomer without SMILES: CID {isomer_entry.cid}")
                continue


            smiles = isomer_entry.canonical_smiles
            iupac_name = isomer_entry.iupac_name
            
            if not iupac_name and isomer_entry.synonyms:
                # انتخاب بهترین مترادف
                chosen_synonym = isomer_entry.synonyms[0] # پیش‌فرض
                simple_alkane_names = [s for s in isomer_entry.synonyms if s.lower().endswith("ane") and s.replace('-', '').isalpha()]
                if simple_alkane_names:
                    chosen_synonym = min(simple_alkane_names, key=len) # کوتاه‌ترین نام آلکان مانند
                iupac_name = chosen_synonym.capitalize() # حرف اول بزرگ
            elif not iupac_name:
                iupac_name = f"Unknown Isomer (CID: {isomer_entry.cid})"

            if smiles: # این شرط باید همیشه درست باشد چون در بالا فیلتر کردیم
                mol_image = draw_molecule(smiles)
                if mol_image:
                    isomer_outputs.append((mol_image, f"{iupac_name}\nSMILES: {smiles}"))
                    processed_cids.add(isomer_entry.cid)
                    valid_isomers_count += 1
                else:
                    print(f"  Failed to draw image for isomer: CID {isomer_entry.cid}, SMILES: {smiles}")
            
        print(f"Processed and displayed {valid_isomers_count} valid hydrocarbon isomers.")

        if not isomer_outputs:
            status_message = "ایزومرها پیدا شدند اما هیچکدام معتبر نبودند یا تصویری برای نمایش وجود ندارد."
            print(status_message)
        elif not status_message: # اگر پیام خطای دیگری از قبل نداشتیم
             status_message = f"{len(isomer_outputs)} ایزومر برای '{molecule_name_input}' (فرمول: {molecular_formula}) پیدا و نمایش داده شد."
        
        isomer_outputs.sort(key=lambda x: x[1]) # مرتب‌سازی بر اساس کپشن (نام و SMILES)
        
        return isomer_outputs, status_message

    except pcp.PubChemHTTPError as e:
        error_msg = f"خطا در ارتباط با PubChem: {e}. ممکن است سرویس موقتا در دسترس نباشد یا درخواست شما نامعتبر باشد."
        print(error_msg)
        return [], error_msg
    except Exception as e:
        error_msg = f"یک خطای غیرمنتظره در سرور رخ داد: {e}"
        print(f"FULL TRACEBACK: {traceback.format_exc()}") # نمایش کامل خطا برای عیب‌یابی
        return [], error_msg

# ساخت رابط کاربری Gradio
iface = gr.Interface(
    fn=find_and_display_isomers,
    inputs=gr.Textbox(
        label="نام آلکان را وارد کنید", 
        placeholder="مثال: butane, pentane, hexane",
        info="نام آلکان مورد نظر خود را به انگلیسی و با حروف کوچک وارد کنید."
    ),
    outputs=[
        gr.Gallery(
            label="ایزومرهای یافت شده", 
            columns=[3], 
            height="auto", 
            object_fit="contain"
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
        ["heptane"] 
    ],
    allow_flagging='never',
    theme=gr.themes.Soft() 
)

if __name__ == '__main__':
    iface.launch()
