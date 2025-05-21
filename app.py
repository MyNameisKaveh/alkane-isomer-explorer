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
        return [], "لطفا نام یک مولکول را وارد کنید."

    molecule_name = molecule_name_input.strip().lower()
    print(f"Processing request for: '{molecule_name}'")
    
    isomer_outputs = []
    status_message = ""

    try:
        print(f"Searching for compound: '{molecule_name}' in PubChem (using default name search)...")
        compounds = pcp.get_compounds(molecule_name, 'name', listkey_count=5) 
        
        main_compound_obj = None
        molecular_formula = None

        if not compounds:
            status_message = f"مولکول '{molecule_name}' در PubChem یافت نشد."
            print(status_message)
            return [], status_message
        
        print(f"Found {len(compounds)} potential matches for '{molecule_name}'. Checking them...")
        for i, c in enumerate(compounds):
            cid = c.cid
            common_name = c.synonyms[0] if c.synonyms else "N/A"
            formula_attr = hasattr(c, 'molecular_formula')
            actual_formula = c.molecular_formula if formula_attr else None
            
            print(f"  Checking match {i+1}: CID {cid}, Name: '{common_name}', HasFormulaAttr: {formula_attr}, Formula: '{actual_formula}'")

            if actual_formula:
                is_simple_hydrocarbon = False
                has_non_standard_isotopes_main = False # برای ترکیب اصلی هم چک کنیم
                if c.canonical_smiles:
                    try:
                        mol_obj = Chem.MolFromSmiles(c.canonical_smiles)
                        if mol_obj:
                            elements = set()
                            for atom in mol_obj.GetAtoms():
                                elements.add(atom.GetSymbol())
                                if atom.GetIsotope() != 0: # اتم با ایزوتوپ غیر استاندارد
                                    has_non_standard_isotopes_main = True
                                    break
                            if not has_non_standard_isotopes_main and elements.issubset({'C', 'H'}):
                                is_simple_hydrocarbon = True
                                print(f"    CID {cid} is a simple hydrocarbon with standard isotopes.")
                            elif has_non_standard_isotopes_main:
                                print(f"    CID {cid} has non-standard isotopes.")
                            else:
                                print(f"    CID {cid} contains elements: {elements}. Not a simple hydrocarbon.")
                        else:
                             print(f"    CID {cid} SMILES '{c.canonical_smiles}' could not be parsed by RDKit.")
                    except Exception as rdkit_ex:
                        print(f"    RDKit error processing SMILES for CID {cid}: {rdkit_ex}")
                else:
                    print(f"    CID {cid} has no SMILES string for element/isotope check.")

                current_compound_name_matches_input = molecule_name in [syn.lower() for syn in c.synonyms]

                if current_compound_name_matches_input and is_simple_hydrocarbon and not has_non_standard_isotopes_main:
                    main_compound_obj = c
                    molecular_formula = actual_formula
                    print(f"  Selected main compound (exact name match & standard hydrocarbon): CID {main_compound_obj.cid}, Name: '{common_name}', Formula: {molecular_formula}")
                    break 
                
                if not main_compound_obj and is_simple_hydrocarbon and not has_non_standard_isotopes_main:
                    main_compound_obj = c 
                    print(f"  Tentatively selected main compound (standard hydrocarbon): CID {main_compound_obj.cid}, Name: '{common_name}', Formula: {actual_formula}")
        
        if main_compound_obj and not molecular_formula: 
            molecular_formula = main_compound_obj.molecular_formula 
            print(f"  Confirmed selection of (first found standard hydrocarbon) main compound: CID {main_compound_obj.cid}, Formula: {molecular_formula}")
        elif main_compound_obj and molecular_formula: 
            print(f"  Selection based on name match and standard hydrocarbon confirmed: CID {main_compound_obj.cid}, Formula: {molecular_formula}")

        if not main_compound_obj or not molecular_formula:
            status_message = f"فرمول مولکولی معتبری برای '{molecule_name}' (به عنوان یک آلکان استاندارد) در نتایج اولیه PubChem یافت نشد."
            print(status_message)
            return [], status_message
        
        print(f"Searching for isomers with formula: {molecular_formula}...")
        # listkey_count را کاهش می‌دهیم چون تعداد ایزومرهای واقعی آلکان‌ها معمولا زیاد نیست.
        # اما PubChem ممکن است نتایج زیادی شامل ایزوتوپ‌ها برگرداند، پس فعلا بالا نگه می‌داریم و با فیلتر دقیق جدا می‌کنیم.
        isomers_found_raw = pcp.get_compounds(molecular_formula, 'formula', listkey_count=50) 

        if not isomers_found_raw:
            status_message = f"ایزومری برای فرمول {molecular_formula} (مربوط به '{molecule_name}') یافت نشد."
            print(status_message)
            return [], status_message

        print(f"Found {len(isomers_found_raw)} potential isomer entries from PubChem. Filtering for true structural alkane isomers...")
        
        processed_cids = set() 
        valid_isomers_count = 0

        for isomer_entry in isomers_found_raw:
            if isomer_entry.cid in processed_cids:
                continue 
            
            if isomer_entry.canonical_smiles:
                try:
                    mol_iso = Chem.MolFromSmiles(isomer_entry.canonical_smiles)
                    if mol_iso:
                        is_true_alkane_isomer = True
                        atom_elements = set()
                        for atom in mol_iso.GetAtoms():
                            atom_elements.add(atom.GetSymbol())
                            # *** تغییر کلیدی: بررسی ایزوتوپ اتم ***
                            if atom.GetIsotope() != 0:
                                is_true_alkane_isomer = False
                                print(f"  Skipping isotopically labeled compound: CID {isomer_entry.cid}, Name: {isomer_entry.iupac_name or (isomer_entry.synonyms[0] if isomer_entry.synonyms else '')}, Atom: {atom.GetSymbol()}{atom.GetIdx()+1}, Isotope: {atom.GetIsotope()}")
                                break
                        
                        if not is_true_alkane_isomer: # اگر ایزوتوپ غیر استاندارد داشت
                            continue

                        if not atom_elements.issubset({'C', 'H'}): # اگر شامل عناصری غیر از C و H بود
                            print(f"  Skipping non-CH compound: CID {isomer_entry.cid}, Name: {isomer_entry.iupac_name or (isomer_entry.synonyms[0] if isomer_entry.synonyms else '')}, Elements: {atom_elements}")
                            is_true_alkane_isomer = False # redundant but for clarity
                            continue
                        
                        # اگر به اینجا رسیدیم، یک ایزومر ساختاری آلکان استاندارد است
                        smiles = isomer_entry.canonical_smiles
                        iupac_name = isomer_entry.iupac_name
                        
                        if not iupac_name and isomer_entry.synonyms:
                            chosen_synonym = isomer_entry.synonyms[0]
                            simple_alkane_names = [s for s in isomer_entry.synonyms if s.lower().endswith("ane") and s.replace('-', '').isalpha()]
                            if simple_alkane_names:
                                chosen_synonym = min(simple_alkane_names, key=len)
                            iupac_name = chosen_synonym.capitalize()
                        elif not iupac_name:
                            iupac_name = f"Isomer (CID: {isomer_entry.cid})"

                        mol_image = draw_molecule(smiles)
                        if mol_image:
                            isomer_outputs.append((mol_image, f"{iupac_name}\nSMILES: {smiles}"))
                            processed_cids.add(isomer_entry.cid)
                            valid_isomers_count += 1
                        else:
                            print(f"  Failed to draw image for isomer: CID {isomer_entry.cid}, SMILES: {smiles}")

                    else: # SMILES نامعتبر برای ایزومر
                        print(f"  Skipping isomer with invalid SMILES: CID {isomer_entry.cid}")
                        continue
                except Exception as rdkit_iso_ex:
                    print(f"  RDKit error for isomer SMILES CID {isomer_entry.cid}: {rdkit_iso_ex}")
                    continue
            else: # ایزومر بدون SMILES
                print(f"  Skipping isomer without SMILES: CID {isomer_entry.cid}")
                continue
            
        print(f"Processed and displayed {valid_isomers_count} valid structural alkane isomers.")

        if not isomer_outputs:
            status_message = "ایزومرهای آلکان استاندارد پیدا نشدند یا تصویری برای نمایش وجود ندارد."
            print(status_message)
        elif not status_message:
             status_message = f"{len(isomer_outputs)} ایزومر ساختاری آلکان برای '{molecule_name_input}' (فرمول: {molecular_formula}) پیدا و نمایش داده شد."
        
        isomer_outputs.sort(key=lambda x: x[1])
        
        return isomer_outputs, status_message

    except pcp.PubChemHTTPError as e:
        error_msg = f"خطا در ارتباط با PubChem: {e}."
        print(error_msg)
        print(f"FULL TRACEBACK for PubChemHTTPError: {traceback.format_exc()}")
        return [], error_msg
    except Exception as e:
        error_msg = f"یک خطای غیرمنتظره در سرور رخ داد: {e}"
        print(f"FULL TRACEBACK for general Exception: {traceback.format_exc()}")
        return [], error_msg

# ساخت رابط کاربری Gradio - بدون تغییر نسبت به نسخه قبلی
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
