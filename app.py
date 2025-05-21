import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage
import gradio as gr
import traceback

def draw_molecule(smiles_string):
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if mol:
            img = MolToImage(mol, size=(300, 300))
            return img
        else:
            print(f"Could not parse SMILES: {smiles_string}")
            return None
    except Exception as e:
        print(f"Error drawing molecule for SMILES {smiles_string}: {e}")
        return None

def find_and_display_isomers(molecule_name_input):
    if not molecule_name_input or not molecule_name_input.strip():
        # قبل از هر جستجو یا با ورودی خالی، گالری مخفی و پیام مناسب
        return gr.Gallery.update(value=[], visible=False), gr.Textbox.update(value="لطفاً نام آلکان مورد نظر خود را برای جستجو وارد کنید.")

    molecule_name = molecule_name_input.strip().lower()
    print(f"Processing request for: '{molecule_name}'")
    
    status_message = ""
    isomer_outputs_final = [] # اطمینان از اینکه لیست برای خروجی گالری تعریف شده

    try:
        print(f"Searching for compound: '{molecule_name}' in PubChem...")
        compounds = pcp.get_compounds(molecule_name, 'name', listkey_count=5) 
        
        main_compound_obj = None
        molecular_formula = None

        if not compounds:
            status_message = f"مولکول '{molecule_name}' در PubChem یافت نشد."
            print(status_message)
            # گالری قابل مشاهده اما خالی، به همراه پیام خطا
            return gr.Gallery.update(value=[], visible=True), gr.Textbox.update(value=status_message)
        
        print(f"Found {len(compounds)} potential matches for '{molecule_name}'. Checking them...")
        for i, c in enumerate(compounds):
            cid = c.cid
            common_name = c.synonyms[0] if c.synonyms else "N/A"
            actual_formula = c.molecular_formula if hasattr(c, 'molecular_formula') else None
            
            print(f"  Checking main compound candidate {i+1}: CID {cid}, Name: '{common_name}', Formula: '{actual_formula}'")

            if actual_formula:
                is_standard_hydrocarbon = True 
                if c.canonical_smiles:
                    try:
                        mol_obj = Chem.MolFromSmiles(c.canonical_smiles)
                        if mol_obj:
                            if len(Chem.GetMolFrags(mol_obj)) > 1:
                                print(f"    Main candidate CID {cid} is disconnected.")
                                is_standard_hydrocarbon = False

                            atom_symbols_main = set()
                            for atom in mol_obj.GetAtoms():
                                atom_symbols_main.add(atom.GetSymbol())
                                if atom.GetIsotope() != 0:
                                    is_standard_hydrocarbon = False; break
                            if not atom_symbols_main.issubset({'C', 'H'}):
                                is_standard_hydrocarbon = False
                            
                            if is_standard_hydrocarbon: # فقط اگر شرایط قبلی برقرار بود، پیوندها و حلقه‌ها را چک کن
                                for bond in mol_obj.GetBonds():
                                    if bond.GetBondType() != Chem.BondType.SINGLE:
                                        is_standard_hydrocarbon = False; break
                                if Chem.GetSymmSSSR(mol_obj):
                                    is_standard_hydrocarbon = False
                        else: is_standard_hydrocarbon = False
                    except Exception: is_standard_hydrocarbon = False 
                else: is_standard_hydrocarbon = False 

                if not is_standard_hydrocarbon: continue

                current_compound_name_matches_input = molecule_name in [syn.lower() for syn in c.synonyms]
                if current_compound_name_matches_input: 
                    main_compound_obj = c; molecular_formula = actual_formula; break
                if not main_compound_obj: 
                    main_compound_obj = c; molecular_formula = actual_formula
        
        if not main_compound_obj or not molecular_formula: 
            status_message = f"آلکان استاندارد با نام '{molecule_name}' در PubChem یافت نشد."
            print(status_message)
            return gr.Gallery.update(value=[], visible=True), gr.Textbox.update(value=status_message)
        
        print(f"Proceeding with main compound: CID {main_compound_obj.cid}, Formula: {molecular_formula}")
        isomers_found_raw = pcp.get_compounds(molecular_formula, 'formula', listkey_count=50) 

        if not isomers_found_raw:
            status_message = f"ایزومری برای فرمول {molecular_formula} یافت نشد."
            print(status_message)
            return gr.Gallery.update(value=[], visible=True), gr.Textbox.update(value=status_message)

        print(f"Found {len(isomers_found_raw)} potential isomer entries. Filtering...")
        
        valid_structural_alkanes_entries = [] 
        unique_accepted_smiles = set()

        for isomer_entry in isomers_found_raw:
            smiles = isomer_entry.canonical_smiles
            if not smiles: continue
            try:
                mol_iso = Chem.MolFromSmiles(smiles)
                if not mol_iso: continue
                is_valid_candidate = True
                if len(Chem.GetMolFrags(mol_iso)) > 1: is_valid_candidate = False
                
                if is_valid_candidate:
                    atom_symbols = set()
                    for atom in mol_iso.GetAtoms(): atom_symbols.add(atom.GetSymbol())
                    if not atom_symbols.issubset({'C', 'H'}): is_valid_candidate = False
                
                if is_valid_candidate:
                    for atom in mol_iso.GetAtoms():
                        if atom.GetSymbol() == 'H' and atom.GetDegree() == 0: is_valid_candidate = False; break
                        if atom.GetIsotope() != 0: is_valid_candidate = False; break 
                
                if is_valid_candidate:
                    for bond in mol_iso.GetBonds():
                        if bond.GetBondType() != Chem.BondType.SINGLE: is_valid_candidate = False; break
                    if Chem.GetSymmSSSR(mol_iso): is_valid_candidate = False

                if is_valid_candidate:
                    # استفاده از SMILES کانونی بدون اطلاعات ایزومری برای بررسی یکتایی ساختار
                    canonical_smiles_for_uniqueness = Chem.MolToSmiles(mol_iso, isomericSmiles=False, canonical=True)
                    if canonical_smiles_for_uniqueness not in unique_accepted_smiles:
                        valid_structural_alkanes_entries.append(isomer_entry)
                        unique_accepted_smiles.add(canonical_smiles_for_uniqueness)
            except Exception: continue
        
        print(f"Found {len(valid_structural_alkanes_entries)} unique, valid structural alkane isomers after filtering.")
        
        valid_isomers_count_final = 0

        for final_isomer_entry in valid_structural_alkanes_entries:
            smiles_to_draw = final_isomer_entry.canonical_smiles # SMILES اصلی را برای رسم استفاده کن
            iupac_name = final_isomer_entry.iupac_name
            if not iupac_name and final_isomer_entry.synonyms:
                # ... (منطق انتخاب نام مشابه قبل، ساده شده) ...
                iupac_name = min(final_isomer_entry.synonyms, key=len).capitalize()
            elif not iupac_name:
                iupac_name = f"آلکان (CID: {final_isomer_entry.cid})"

            mol_image = draw_molecule(smiles_to_draw)
            if mol_image:
                isomer_outputs_final.append((mol_image, f"{iupac_name}\nSMILES: {smiles_to_draw}"))
                valid_isomers_count_final += 1
        
        print(f"Displayed {valid_isomers_count_final} isomers.")

        if not isomer_outputs_final:
            status_message = "ایزومر آلکان استاندارد و قابل رسمی برای نمایش پیدا نشد."
        else:
             status_message = f"{len(isomer_outputs_final)} ایزومر ساختاری آلکان برای '{molecule_name_input}' (فرمول: {molecular_formula}) پیدا و نمایش داده شد."
        
        isomer_outputs_final.sort(key=lambda x: x[1])
        # در هر صورت گالری را پس از جستجو نشان بده، حتی اگر خالی باشد
        return gr.Gallery.update(value=isomer_outputs_final, visible=True), gr.Textbox.update(value=status_message)

    except pcp.PubChemHTTPError as e:
        error_msg = f"خطا در ارتباط با PubChem: {e}."
        print(f"FULL TRACEBACK for PubChemHTTPError: {traceback.format_exc()}")
        return gr.Gallery.update(value=[], visible=True), gr.Textbox.update(value=error_msg)
    except Exception as e:
        error_msg = f"یک خطای غیرمنتظره در سرور رخ داد: {e}"
        print(f"FULL TRACEBACK for general Exception: {traceback.format_exc()}")
        return gr.Gallery.update(value=[], visible=True), gr.Textbox.update(value=error_msg)

# تعریف CSS برای جهت‌دهی راست به چپ و فونت
custom_css = """
body { direction: rtl; font-family: 'Tahoma', 'Vazirmatn', sans-serif !important; }
.gradio-container { font-family: 'Tahoma', 'Vazirmatn', sans-serif !important; }
label { font-family: 'Tahoma', 'Vazirmatn', sans-serif !important; }
.gr-button { font-family: 'Tahoma', 'Vazirmatn', sans-serif !important; }
.gr-input { text-align: right !important; }
.output_textClass { text-align: right !important; direction: rtl !important; }
footer { display: none !important; } /* حذف فوتر پیش‌فرض Gradio */
"""

# --- بخش Gradio Interface ---
iface = gr.Interface(
    fn=find_and_display_isomers,
    inputs=gr.Textbox(
        label="نام آلکان را وارد کنید:", # اصلاح لیبل
        placeholder="مثال: butane, pentane, hexane (نام انگلیسی)", # توضیح بیشتر
        info="نام آلکان مورد نظر خود را به انگلیسی و با حروف کوچک وارد کنید."
    ),
    outputs=[
        gr.Gallery(
            label="ایزومرهای یافت شده:", # اصلاح لیبل
            columns=[3], 
            height="auto", 
            object_fit="contain",
            visible=False # <-- گالری در ابتدا مخفی است
        ),
        gr.Textbox(label="وضعیت و پیام‌ها:", elem_classes="output_textClass") # اصلاح لیبل
    ],
    title="یابنده و نمایشگر ایزومرهای آلکان",
    description=(
        "نام یک آلکان را (به انگلیسی) وارد کنید تا ایزومرهای ساختاری آن به همراه ساختار شیمیایی و نام رایج/IUPAC نمایش داده شوند.\n"
        "این برنامه از دیتابیس PubChem برای دریافت اطلاعات و از کتابخانه RDKit برای رسم ساختارها استفاده می‌کند.\n"
        "توجه: فقط ایزومرهای ساختاری استاندارد (غیرحلقوی، بدون پیوند چندگانه، بدون ایزوتوپ غیرمعمول) نمایش داده می‌شوند."
    ),
    examples=[
        ["butane"], 
        ["pentane"], 
        ["hexane"],
        ["heptane"] 
    ],
    allow_flagging='never',
    theme=gr.themes.Soft(), # یا هر تم دیگری که می‌پسندید
    css=custom_css # <-- اعمال CSS سفارشی
)

if __name__ == '__main__':
    iface.launch()
