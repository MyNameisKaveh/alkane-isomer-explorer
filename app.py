import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage
import gradio as gr
import traceback

def draw_molecule(smiles_string):
    """
    Renders a molecule image from a SMILES string.
    """
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
    """
    Finds and displays structural alkane isomers for a given molecule name.
    """
    if not molecule_name_input or not molecule_name_input.strip():
        return [], "لطفا نام یک مولکول را وارد کنید."

    molecule_name = molecule_name_input.strip().lower()
    print(f"Processing request for: '{molecule_name}'")
    
    status_message = ""

    try:
        print(f"Searching for compound: '{molecule_name}' in PubChem (up to 10 candidates)...")
        # Increased listkey_count for initial search
        compounds = pcp.get_compounds(molecule_name, 'name', listkey_count=10) 
        
        main_compound_obj = None
        molecular_formula = None

        if not compounds:
            status_message = f"مولکول '{molecule_name}' در PubChem یافت نشد. لطفا املای آن را بررسی کنید."
            print(status_message)
            return [], status_message
        
        print(f"Found {len(compounds)} potential matches for '{molecule_name}'. Checking them for standard alkane properties...")
        
        # This loop tries to find the most relevant "alkane" match for the input name
        for i, c in enumerate(compounds):
            cid = c.cid
            common_name = c.synonyms[0] if c.synonyms else "N/A"
            actual_formula = c.molecular_formula if hasattr(c, 'molecular_formula') else None
            
            print(f"  Checking main compound candidate {i+1}: CID {cid}, Name: '{common_name}', Formula: '{actual_formula}'")

            # --- Alkane Specificity Check for the Main Compound ---
            is_standard_alkane_candidate = True 

            # Must have a molecular formula
            if not actual_formula:
                is_standard_alkane_candidate = False
                print(f"    Main candidate CID {cid} has no molecular formula.")
                
            # Must have a SMILES string
            if is_standard_alkane_candidate and not c.canonical_smiles:
                is_standard_alkane_candidate = False
                print(f"    Main candidate CID {cid} has no SMILES string for detailed check.")

            if is_standard_alkane_candidate:
                try:
                    mol_obj = Chem.MolFromSmiles(c.canonical_smiles)
                    if mol_obj:
                        # 1. Check for disconnected fragments (should be a single molecule)
                        if len(Chem.GetMolFrags(mol_obj)) > 1:
                            print(f"    Main candidate CID {cid} is disconnected.")
                            is_standard_alkane_candidate = False

                        # 2. Check for non-standard isotopes and atom types (must be C and H only, no isotopes)
                        if is_standard_alkane_candidate:
                            atom_symbols_main = set()
                            for atom in mol_obj.GetAtoms():
                                atom_symbols_main.add(atom.GetSymbol())
                                if atom.GetIsotope() != 0:
                                    print(f"    Main candidate CID {cid} has non-standard isotope: {atom.GetSymbol()}{atom.GetIsotope()}")
                                    is_standard_alkane_candidate = False
                                    break
                            if not atom_symbols_main.issubset({'C', 'H'}):
                                print(f"    Main candidate CID {cid} is not CH only: {atom_symbols_main}")
                                is_standard_alkane_candidate = False

                        # 3. Check for single bonds only and no rings
                        if is_standard_alkane_candidate:
                            for bond in mol_obj.GetBonds():
                                if bond.GetBondType() != Chem.BondType.SINGLE:
                                    print(f"    Main candidate CID {cid} has non-single bond: {bond.GetBondType()}")
                                    is_standard_alkane_candidate = False
                                    break
                            if Chem.GetSymmSSSR(mol_obj): # Smallest Set of Smallest Rings
                                print(f"    Main candidate CID {cid} has rings.")
                                is_standard_alkane_candidate = False

                    else: # RDKit could not parse SMILES
                        print(f"    Main candidate CID {cid} SMILES '{c.canonical_smiles}' could not be parsed by RDKit.")
                        is_standard_alkane_candidate = False
                except Exception as rdkit_ex:
                    print(f"    RDKit error processing SMILES for main candidate CID {cid}: {rdkit_ex}")
                    is_standard_alkane_candidate = False 
            
            # If not a standard alkane, skip to the next candidate
            if not is_standard_alkane_candidate: 
                continue

            # Prioritize compounds whose synonyms match the input name
            current_compound_name_matches_input = molecule_name in [syn.lower() for syn in c.synonyms]
            if current_compound_name_matches_input: 
                main_compound_obj = c
                molecular_formula = actual_formula
                print(f"  SELECTED main compound based on name match: CID {main_compound_obj.cid}, Formula: {molecular_formula}")
                break 
            
            # If no direct name match, tentatively select the first valid alkane candidate
            if not main_compound_obj: 
                main_compound_obj = c 
                molecular_formula = actual_formula 
                print(f"  TENTATIVELY selected first valid alkane candidate: CID {main_compound_obj.cid}, Formula: {molecular_formula}")
        
        if not main_compound_obj or not molecular_formula: 
            status_message = f"آلکان ساختاری استاندارد با نام '{molecule_name}' در PubChem یافت نشد. (این ابزار تنها آلکان‌های شامل کربن و هیدروژن، بدون حلقه، بدون پیوند دوگانه/سه‌گانه و بدون ایزوتوپ را جستجو می‌کند.)"
            print(status_message)
            return [], status_message
        
        print(f"Proceeding with main compound: CID {main_compound_obj.cid}, Formula: {molecular_formula}")
        print(f"Searching for isomers with formula: {molecular_formula} (up to 200 candidates)...")
        # Increased listkey_count significantly for isomers
        isomers_found_raw = pcp.get_compounds(molecular_formula, 'formula', listkey_count=200) 

        if not isomers_found_raw:
            status_message = f"ایزومری برای فرمول {molecular_formula} یافت نشد."
            print(status_message)
            return [], status_message

        print(f"Found {len(isomers_found_raw)} potential isomer entries from PubChem. Filtering for true structural alkane isomers...")
        
        valid_structural_alkanes_entries = [] 
        unique_accepted_smiles = set() # To ensure uniqueness based on canonical SMILES

        for isomer_entry in isomers_found_raw:
            smiles = isomer_entry.canonical_smiles
            if not smiles:
                print(f"  Skipping isomer without SMILES: CID {isomer_entry.cid}")
                continue

            try:
                mol_iso = Chem.MolFromSmiles(smiles)
                if not mol_iso: # RDKit failed to parse SMILES
                    print(f"  FILTERED (Invalid SMILES parse): CID {isomer_entry.cid}, SMILES: {smiles}")
                    continue

                is_valid_alkane_isomer = True
                
                # 1. Check for disconnected fragments
                if len(Chem.GetMolFrags(mol_iso)) > 1:
                    print(f"  FILTERED (Disconnected): CID {isomer_entry.cid}, NumFrags: {len(Chem.GetMolFrags(mol_iso))}, SMILES: {smiles}")
                    is_valid_alkane_isomer = False
                
                # 2. Check for C and H only (no other atoms, no isotopes)
                if is_valid_alkane_isomer:
                    atom_symbols = set()
                    for atom in mol_iso.GetAtoms():
                        atom_symbols.add(atom.GetSymbol())
                        if atom.GetIsotope() != 0: # Check for isotopes
                            print(f"  FILTERED (Isotope present): CID {isomer_entry.cid}, Atom: {atom.GetSymbol()}{atom.GetIdx()+1}, Isotope: {atom.GetIsotope()}, SMILES: {smiles}")
                            is_valid_alkane_isomer = False
                            break
                    if not atom_symbols.issubset({'C', 'H'}): # Check for non-C/H atoms
                        print(f"  FILTERED (Non-CH elements): CID {isomer_entry.cid}, Elements: {atom_symbols}, SMILES: {smiles}")
                        is_valid_alkane_isomer = False
                
                # 3. Check for single bonds only and no rings
                if is_valid_alkane_isomer:
                    for bond in mol_iso.GetBonds():
                        if bond.GetBondType() != Chem.BondType.SINGLE:
                            print(f"  FILTERED (Non-single bond): CID {isomer_entry.cid}, BondType: {bond.GetBondType()}, SMILES: {smiles}")
                            is_valid_alkane_isomer = False
                            break
                    if Chem.GetSymmSSSR(mol_iso): # Check for rings
                        print(f"  FILTERED (Has rings): CID {isomer_entry.cid}, SMILES: {smiles}")
                        is_valid_alkane_isomer = False

                if is_valid_alkane_isomer:
                    # Use non-isomeric canonical SMILES for uniqueness to avoid stereo-isomers
                    canonical_smiles_for_uniqueness = Chem.MolToSmiles(mol_iso, isomericSmiles=False)
                    if canonical_smiles_for_uniqueness not in unique_accepted_smiles:
                        print(f"  ACCEPTED: CID {isomer_entry.cid}, SMILES: {smiles}")
                        valid_structural_alkanes_entries.append(isomer_entry)
                        unique_accepted_smiles.add(canonical_smiles_for_uniqueness)
                    else:
                        print(f"  Skipping (Duplicate SMILES): CID {isomer_entry.cid}, SMILES: {smiles}")

            except Exception as rdkit_iso_ex:
                print(f"  RDKit or processing error for isomer SMILES CID {isomer_entry.cid}: {rdkit_iso_ex}, SMILES: {smiles}")
                continue
        
        print(f"Found {len(valid_structural_alkanes_entries)} unique, valid structural alkane isomers after filtering.")
        
        isomer_outputs_final = []
        valid_isomers_count_final = 0

        # Prepare results for display
        for final_isomer_entry in valid_structural_alkanes_entries:
            smiles_to_draw = final_isomer_entry.canonical_smiles
            
            # --- Improved Naming Logic ---
            isomer_display_name = final_isomer_entry.iupac_name
            if not isomer_display_name and final_isomer_entry.synonyms:
                # Prioritize simple alkane names (e.g., "Butane" over "n-Butane")
                simple_alkane_names = [
                    s for s in final_isomer_entry.synonyms 
                    if s.lower().endswith("ane") and not any(char.isdigit() for char in s) and '-' not in s
                ]
                if simple_alkane_names:
                    isomer_display_name = min(simple_alkane_names, key=len) # Choose the shortest simple alkane name
                else:
                    isomer_display_name = final_isomer_entry.synonyms[0] # Fallback to the first synonym
            
            if not isomer_display_name: # Final fallback if no IUPAC or synonyms found
                isomer_display_name = f"Alkane (CID: {final_isomer_entry.cid})"
            
            isomer_display_name = isomer_display_name.capitalize() # Capitalize for consistent display

            mol_image = draw_molecule(smiles_to_draw)
            if mol_image:
                isomer_outputs_final.append((mol_image, f"{isomer_display_name}\nSMILES: {smiles_to_draw}"))
                valid_isomers_count_final += 1
            else:
                print(f"  Failed to draw image for accepted isomer: CID {final_isomer_entry.cid}, SMILES: {smiles_to_draw}")

        print(f"Displayed {valid_isomers_count_final} isomers in the gallery.")

        if not isomer_outputs_final:
            status_message = "ایزومر ساختاری آلکان استاندارد و قابل رسمی برای مولکول وارد شده پیدا نشد."
            if len(valid_structural_alkanes_entries) > 0:
                status_message += " (برخی ایزومرهای شناسایی شده در مرحله رسم ناموفق بودند یا در فیلترهای نهایی رد شدند.)"
        else:
            status_message = (
                f"{len(isomer_outputs_final)} ایزومر ساختاری آلکان برای '{molecule_name_input}' "
                f"(فرمول: {molecular_formula}) پیدا و نمایش داده شد. "
                f"توجه: این ابزار تنها ایزومرهای شامل کربن و هیدروژن، بدون حلقه، بدون پیوند چندگانه و بدون ایزوتوپ را شناسایی می‌کند. "
                f"(ممکن است ایزومرهای بیشتری نیز در PubChem وجود داشته باشند که در این جستجو دریافت نشده‌اند.)"
            )
        
        # Sort by name for better presentation
        isomer_outputs_final.sort(key=lambda x: x[1])
        
        return isomer_outputs_final, status_message

    except pcp.PubChemHTTPError as e:
        error_msg = f"خطا در ارتباط با PubChem: {e}. لطفا اتصال اینترنت خود را بررسی کنید یا بعداً امتحان کنید."
        print(error_msg)
        print(f"FULL TRACEBACK for PubChemHTTPError: {traceback.format_exc()}")
        return [], error_msg
    except Exception as e:
        error_msg = f"یک خطای غیرمنتظره در سرور رخ داد: {e}"
        print(f"FULL TRACEBACK for general Exception: {traceback.format_exc()}")
        return [], error_msg

# --- Gradio Interface ---
iface = gr.Interface(
    fn=find_and_display_isomers,
    inputs=gr.Textbox(
        label="نام آلکان را وارد کنید", 
        placeholder="مثال: butane, pentane, hexane",
        info="نام آلکان مورد نظر خود را به انگلیسی و با حروف کوچک وارد کنید. (مانند: pentane)"
    ),
    outputs=[
        gr.Gallery(
            label="ایزومرهای یافت شده", 
            columns=[3], # 3 columns for better layout
            height="auto", # Adjust height automatically based on content
            object_fit="contain" # Ensures images fit within their containers
        ),
        gr.Textbox(label="وضعیت و پیام‌ها")
    ],
    title="یابنده و نمایشگر ایزومرهای ساختاری آلکان",
    description=(
        "نام یک آلکان (به انگلیسی) را وارد کنید تا ایزومرهای ساختاری آن (تنها شامل کربن و هیدروژن، "
        "بدون حلقه، بدون پیوند چندگانه، بدون ایزوتوپ و بدون قطعات جدا شده) به همراه ساختار شیمیایی و نام IUPAC (یا رایج) نمایش داده شوند.\n"
        "اطلاعات از دیتابیس PubChem دریافت شده و ساختارها با استفاده از کتابخانه RDKit رسم می‌شوند."
    ),
    examples=[
        ["butane"], 
        ["pentane"], 
        ["hexane"],
        ["heptane"], # Has 9 isomers, a good test case
        ["octane"] # Has 18 isomers
    ],
    allow_flagging='never', # Disables the "Flag" button
    theme=gr.themes.Soft(), # A pleasant, soft theme
    live=False # Only run function on submit, not on every keystroke
)

if __name__ == '__main__':
    iface.launch()
