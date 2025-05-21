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
        return [], "لطفا نام یک مولکول را وارد کنید."

    molecule_name = molecule_name_input.strip().lower()
    print(f"Processing request for: '{molecule_name}'")
    
    # isomer_outputs در انتهای تابع قبل از return ساخته می‌شود
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
            
            print(f"  Checking main compound candidate {i+1}: CID {cid}, Name: '{common_name}', Formula: '{actual_formula}'")

            if actual_formula:
                is_standard_hydrocarbon = True 
                if c.canonical_smiles:
                    try:
                        mol_obj = Chem.MolFromSmiles(c.canonical_smiles)
                        if mol_obj:
                            # Check for disconnected fragments in main compound
                            if len(Chem.GetMolFrags(mol_obj)) > 1:
                                print(f"    Main candidate CID {cid} is disconnected.")
                                is_standard_hydrocarbon = False

                            if is_standard_hydrocarbon:
                                atom_symbols_main = set()
                                for atom in mol_obj.GetAtoms():
                                    atom_symbols_main.add(atom.GetSymbol())
                                    if atom.GetIsotope() != 0:
                                        print(f"    Main candidate CID {cid} has non-standard isotope: {atom.GetSymbol()}{atom.GetIsotope()}")
                                        is_standard_hydrocarbon = False
                                        break
                                if not is_standard_hydrocarbon: continue 

                                if not atom_symbols_main.issubset({'C', 'H'}):
                                    print(f"    Main candidate CID {cid} is not CH only: {atom_symbols_main}")
                                    is_standard_hydrocarbon = False
                        else: 
                             print(f"    Main candidate CID {cid} SMILES '{c.canonical_smiles}' could not be parsed by RDKit.")
                             is_standard_hydrocarbon = False
                    except Exception as rdkit_ex:
                        print(f"    RDKit error processing SMILES for main candidate CID {cid}: {rdkit_ex}")
                        is_standard_hydrocarbon = False 
                else: 
                    print(f"    Main candidate CID {cid} has no SMILES string for detailed check.")
                    is_standard_hydrocarbon = False 

                if not is_standard_hydrocarbon: continue

                current_compound_name_matches_input = molecule_name in [syn.lower() for syn in c.synonyms]

                if current_compound_name_matches_input: 
                    main_compound_obj = c
                    molecular_formula = actual_formula
                    print(f"  SELECTED main compound (name match & standard hydrocarbon): CID {main_compound_obj.cid}, Formula: {molecular_formula}")
                    break 
                
                if not main_compound_obj: 
                    main_compound_obj = c 
                    molecular_formula = actual_formula 
                    print(f"  TENTATIVELY selected main compound (standard hydrocarbon): CID {main_compound_obj.cid}, Formula: {molecular_formula}")
        
        if not main_compound_obj or not molecular_formula: 
            status_message = f"آلکان استاندارد (غیر ایزوتوپی، فقط C و H، متصل) با نام '{molecule_name}' در PubChem یافت نشد."
            print(status_message)
            return [], status_message
        
        print(f"Proceeding with main compound: CID {main_compound_obj.cid}, Formula: {molecular_formula}")
        print(f"Searching for isomers with formula: {molecular_formula}...")
        isomers_found_raw = pcp.get_compounds(molecular_formula, 'formula', listkey_count=50) 

        if not isomers_found_raw:
            status_message = f"ایزومری برای فرمول {molecular_formula} یافت نشد."
            print(status_message)
            return [], status_message

        print(f"Found {len(isomers_found_raw)} potential isomer entries from PubChem. Filtering for true, connected, structural alkane isomers...")
        
        valid_structural_alkanes_entries = [] 
        unique_accepted_smiles = set()

        for isomer_entry in isomers_found_raw:
            smiles = isomer_entry.canonical_smiles
            if not smiles:
                print(f"  Skipping isomer without SMILES: CID {isomer_entry.cid}")
                continue

            try:
                mol_iso = Chem.MolFromSmiles(smiles)
                if not mol_iso:
                    print(f"  FILTERED (Invalid SMILES): CID {isomer_entry.cid}, SMILES: {smiles}")
                    continue

                is_valid_candidate = True
                
                if len(Chem.GetMolFrags(mol_iso)) > 1:
                    print(f"  FILTERED (Disconnected): CID {isomer_entry.cid}, NumFrags: {len(Chem.GetMolFrags(mol_iso))}, SMILES: {smiles}")
                    is_valid_candidate = False
                
                if is_valid_candidate:
                    atom_symbols = set()
                    for atom in mol_iso.GetAtoms():
                        atom_symbols.add(atom.GetSymbol())
                    
                    if not atom_symbols.issubset({'C', 'H'}):
                        print(f"  FILTERED (Non-CH): CID {isomer_entry.cid}, Elements: {atom_symbols}, SMILES: {smiles}")
                        is_valid_candidate = False
                
                if is_valid_candidate:
                    for atom in mol_iso.GetAtoms():
                        if atom.GetSymbol() == 'H' and atom.GetDegree() == 0 and atom.GetMass() > 1.01: 
                             print(f"  FILTERED (Likely [HH] or isolated heavy H): CID {isomer_entry.cid}, Atom: {atom.GetSymbol()}{atom.GetIdx()+1}, SMILES: {smiles}")
                             is_valid_candidate = False
                             break
                        if atom.GetIsotope() != 0:
                            print(f"  FILTERED (Isotope): CID {isomer_entry.cid}, Atom: {atom.GetSymbol()}{atom.GetIdx()+1}, Isotope: {atom.GetIsotope()}, SMILES: {smiles}")
                            is_valid_candidate = False
                            break 
                
                if is_valid_candidate:
                    if smiles not in unique_accepted_smiles:
                        print(f"  ACCEPTED (Unique SMILES): CID {isomer_entry.cid}, SMILES: {smiles}")
                        valid_structural_alkanes_entries.append(isomer_entry)
                        unique_accepted_smiles.add(smiles)
                    else:
                        print(f"  Skipping (Duplicate SMILES): CID {isomer_entry.cid}, SMILES: {smiles}")

            except Exception as rdkit_iso_ex:
                print(f"  RDKit or processing error for isomer SMILES CID {isomer_entry.cid}: {rdkit_iso_ex}, SMILES: {smiles}")
                continue
        
        print(f"Found {len(valid_structural_alkanes_entries)} unique, valid structural alkane isomer entries after filtering.")
        
        isomer_outputs_final = [] # لیست جدید برای خروجی نهایی
        valid_isomers_count_final = 0

        for final_isomer_entry in valid_structural_alkanes_entries:
            smiles_to_draw = final_isomer_entry.canonical_smiles
            iupac_name = final_isomer_entry.iupac_name
            
            if not iupac_name and final_isomer_entry.synonyms:
                chosen_synonym = final_isomer_entry.synonyms[0]
                simple_names = [s for s in final_isomer_entry.synonyms if s.lower().endswith("ane") and not any(char.isdigit() for char in s.split('-')[0]) and '-' not in s.split(' ')[0]]
                if simple_names:
                    chosen_synonym = min(simple_names, key=len)
                else:
                    non_iupac_synonyms = [s for s in final_isomer_entry.synonyms if s != final_isomer_entry.iupac_name]
                    if non_iupac_synonyms:
                        chosen_synonym = min(non_iupac_synonyms, key=len)
                iupac_name = chosen_synonym.capitalize()
            elif not iupac_name:
                iupac_name = f"Alkane (CID: {final_isomer_entry.cid})"

            mol_image = draw_molecule(smiles_to_draw)
            if mol_image:
                isomer_outputs_final.append((mol_image, f"{iupac_name}\nSMILES: {smiles_to_draw}"))
                valid_isomers_count_final += 1
            else:
                print(f"  Failed to draw image for accepted isomer: CID {final_isomer_entry.cid}, SMILES: {smiles_to_draw}")

        print(f"Displayed {valid_isomers_count_final} isomers in the gallery.")

        if not isomer_outputs_final:
            status_message = "ایزومر آلکان استاندارد و قابل رسمی پیدا نشد."
            if len(valid_structural_alkanes_entries) > 0 : # اگر کاندیداهایی بودند اما رسم نشدند
                 status_message += " (برخی در مرحله رسم ناموفق بودند یا کاندیدای معتبری نبودند)."
        else:
             status_message = f"{len(isomer_outputs_final)} ایزومر ساختاری آلکان برای '{molecule_name_input}' (فرمول: {molecular_formula}) پیدا و نمایش داده شد."
        
        isomer_outputs_final.sort(key=lambda x: x[1])
        
        return isomer_outputs_final, status_message # بازگرداندن لیست جدید

    except pcp.PubChemHTTPError as e:
        error_msg = f"خطا در ارتباط با PubChem: {e}."
        print(error_msg)
        print(f"FULL TRACEBACK for PubChemHTTPError: {traceback.format_exc()}")
        return [], error_msg
    except Exception as e:
        error_msg = f"یک خطای غیرمنتظره در سرور رخ داد: {e}"
        print(f"FULL TRACEBACK for general Exception: {traceback.format_exc()}")
        return [], error_msg

# --- بخش Gradio Interface بدون تغییر ---
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
