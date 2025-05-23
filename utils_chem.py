# utils_chem.py

import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage
from rdkit.Chem import rdDepictor
import traceback
from PIL import Image
import os # برای get_sdf_content

# --- Drawing and SDF ---
def draw_molecule_pil(smiles_string, size=(350, 300)):
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if mol:
            rdDepictor.Compute2DCoords(mol)
            img = MolToImage(mol, size=size)
            return img
        else: return None
    except Exception as e:
        print(f"Error in draw_molecule_pil for SMILES {smiles_string}: {e}")
        return None

def get_sdf_content(cid):
    if cid is None: return None, "CID not provided for SDF content."
    print(f"Fetching SDF content for CID: {cid}...")
    temp_sdf_file_dir = "/tmp" 
    if not os.path.exists(temp_sdf_file_dir):
        try: os.makedirs(temp_sdf_file_dir, exist_ok=True)
        except OSError: temp_sdf_file_dir = "."
    temp_sdf_file = os.path.join(temp_sdf_file_dir, f'temp_3d_structure_{cid}.sdf')
    sdf_data, error_message = None, None
    try:
        pcp.download('SDF', temp_sdf_file, str(cid), 'cid', record_type='3d', overwrite=True)
        with open(temp_sdf_file, 'r') as f: sdf_data = f.read()
        if not sdf_data or sdf_data.strip() == "$$$$\n" or sdf_data.strip() == "":
            error_message, sdf_data = f"SDF file for CID {cid} is empty or invalid.", None
    except pcp.NotFoundError: error_message = f"3D structure (SDF) for CID {cid} not found on PubChem."
    except FileNotFoundError: error_message = f"Temporary SDF file for CID {cid} not found."
    except Exception as e: error_message = f"Error downloading/reading SDF for CID {cid}: {e}\n{traceback.format_exc()}"
    finally:
        if os.path.exists(temp_sdf_file):
            try: os.remove(temp_sdf_file)
            except Exception: pass
    # هشدار در Streamlit UI باید توسط تابع فراخواننده مدیریت شود
    return sdf_data, error_message

# --- Property Extraction ---
def get_compound_properties(compound_obj):
    if not compound_obj: return None
    properties = {
        "IUPAC Name": getattr(compound_obj, 'iupac_name', 'N/A'),
        "Molecular Formula": getattr(compound_obj, 'molecular_formula', 'N/A'),
        "Molecular Weight": f"{getattr(compound_obj, 'molecular_weight', 'N/A')} g/mol",
        "Canonical SMILES": getattr(compound_obj, 'canonical_smiles', 'N/A'),
        "Isomeric SMILES": getattr(compound_obj, 'isomeric_smiles', 'N/A'),
        "PubChem CID": getattr(compound_obj, 'cid', 'N/A'),
        "Charge": getattr(compound_obj, 'charge', 'N/A'),
        "Exact Mass": f"{getattr(compound_obj, 'exact_mass', 'N/A')} amu",
        "Monoisotopic Mass": f"{getattr(compound_obj, 'monoisotopic_mass', 'N/A')} amu",
        "TPSA": f"{getattr(compound_obj, 'tpsa', 'N/A')} Å²" if hasattr(compound_obj, 'tpsa') and compound_obj.tpsa is not None else 'N/A',
        "XLogP": getattr(compound_obj, 'xlogp', 'N/A'),
        "Heavy Atom Count": getattr(compound_obj, 'heavy_atom_count', 'N/A'),
        "H-Bond Donor Count": getattr(compound_obj, 'hydrogen_bond_donor_count', 'N/A'),
        "H-Bond Acceptor Count": getattr(compound_obj, 'hydrogen_bond_acceptor_count', 'N/A'),
        "Rotatable Bond Count": getattr(compound_obj, 'rotatable_bond_count', 'N/A'),
    }
    return properties

# --- Processing Logic ---
def process_alkane_request(molecule_name_input):
    # ... (کد کامل process_alkane_request از پاسخ قبلی اینجا قرار می‌گیرد) ...
    # فقط مطمئن شوید که import های لازم (مانند Chem, pcp, rdDepictor) در بالای این فایل موجود باشند.
    if not molecule_name_input or not molecule_name_input.strip():
        return [], "Please enter a molecule name.", None 
    molecule_name = molecule_name_input.strip().lower()
    status_message, isomer_details_list = f"Searching for '{molecule_name}'...", []
    main_molecule_properties = None
    try:
        compounds = pcp.get_compounds(molecule_name, 'name', listkey_count=5) 
        main_compound_obj, molecular_formula = None, None
        if not compounds: return [], f"Molecule '{molecule_name}' not found on PubChem.", None
        for c in compounds:
            cid = c.cid; actual_formula = c.molecular_formula if hasattr(c, 'molecular_formula') else None
            if actual_formula:
                is_standard_hydrocarbon = True 
                if c.canonical_smiles:
                    try:
                        mol_obj = Chem.MolFromSmiles(c.canonical_smiles)
                        if mol_obj:
                            rdDepictor.Compute2DCoords(mol_obj)
                            if len(Chem.GetMolFrags(mol_obj)) > 1: is_standard_hydrocarbon = False
                            if is_standard_hydrocarbon:
                                atom_symbols_main = set(atom.GetSymbol() for atom in mol_obj.GetAtoms())
                                if not atom_symbols_main.issubset({'C', 'H'}): is_standard_hydrocarbon = False
                                for atom in mol_obj.GetAtoms():
                                    if atom.GetIsotope() != 0: is_standard_hydrocarbon = False; break
                            if is_standard_hydrocarbon:
                                for bond in mol_obj.GetBonds():
                                    if bond.GetBondType() != Chem.BondType.SINGLE: is_standard_hydrocarbon = False; break
                                if Chem.GetSymmSSSR(mol_obj): is_standard_hydrocarbon = False
                        else: is_standard_hydrocarbon = False
                    except Exception: is_standard_hydrocarbon = False 
                else: is_standard_hydrocarbon = False 
                if not is_standard_hydrocarbon: continue
                if molecule_name in [syn.lower() for syn in c.synonyms]: 
                    main_compound_obj = c; molecular_formula = actual_formula; break
                if not main_compound_obj: main_compound_obj = c; molecular_formula = actual_formula
        if not main_compound_obj or not molecular_formula:
            return [], f"Standard alkane '{molecule_name}' not found or does not meet criteria.", None
        main_molecule_properties = get_compound_properties(main_compound_obj)
        main_name = main_molecule_properties.get("IUPAC Name", molecule_name.capitalize()) if main_molecule_properties else molecule_name.capitalize()
        status_message = f"Finding isomers for {main_name} (Formula: {molecular_formula})..."
        isomers_found_raw = pcp.get_compounds(molecular_formula, 'formula', listkey_count=50)
        if not isomers_found_raw: return [], f"No isomers found for formula {molecular_formula}.", main_molecule_properties
        valid_structural_alkanes_entries, unique_accepted_smiles = [], set()
        for isomer_entry in isomers_found_raw:
            smiles = isomer_entry.canonical_smiles;
            if not smiles: continue
            try:
                mol_iso = Chem.MolFromSmiles(smiles)
                if not mol_iso: continue
                is_valid_candidate = True
                if len(Chem.GetMolFrags(mol_iso)) > 1: is_valid_candidate = False
                if is_valid_candidate:
                    atom_symbols = set(atom.GetSymbol() for atom in mol_iso.GetAtoms())
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
                    canonical_smiles_for_uniqueness = Chem.MolToSmiles(mol_iso, isomericSmiles=False)
                    if canonical_smiles_for_uniqueness not in unique_accepted_smiles:
                        valid_structural_alkanes_entries.append(isomer_entry); unique_accepted_smiles.add(canonical_smiles_for_uniqueness)
            except Exception: continue
        if not valid_structural_alkanes_entries:
             return [], f"No valid alkane isomers found for {molecular_formula} after filtering.", main_molecule_properties
        for entry in sorted(valid_structural_alkanes_entries, key=lambda x: (len(x.canonical_smiles), x.cid)):
            pil_image = draw_molecule_pil(entry.canonical_smiles)
            if pil_image:
                name = entry.iupac_name
                if not name and entry.synonyms:
                    simple_names = [s for s in entry.synonyms if s.lower().endswith("ane") and not any(char.isdigit() for char in s.split('-')[0]) and '-' not in s.split(' ')[0]]
                    if simple_names: name = min(simple_names, key=len)
                    else: name = entry.synonyms[0]
                elif not name: name = f"Alkane (CID: {entry.cid})"
                isomer_details_list.append({"cid": entry.cid, "name": name.capitalize(), "smiles": entry.canonical_smiles, "image": pil_image})
        if isomer_details_list: status_message = f"{len(isomer_details_list)} isomers found for '{molecule_name}' ({molecular_formula})."
        else: status_message = f"No displayable isomers found for '{molecule_name}'."
        return isomer_details_list, status_message, main_molecule_properties
    except pcp.PubChemHTTPError as e: return [], f"PubChem API Error: {e}", None
    except Exception as e: return [], f"An unexpected error occurred: {e}\n{traceback.format_exc()}", None


def process_general_molecule_search(search_term):
    # ... (کد کامل process_general_molecule_search از پاسخ قبلی اینجا قرار می‌گیرد) ...
    if not search_term or not search_term.strip():
        return None, "Please enter a chemical name or identifier."
    term = search_term.strip()
    status_message = f"Searching for '{term}'..."
    molecule_details = None
    try:
        compounds = pcp.get_compounds(term, 'name', listkey_count=1)
        if compounds:
            compound_obj = compounds[0]
            props = get_compound_properties(compound_obj)
            name_to_display = props.get("IUPAC Name", term.capitalize()) if props else term.capitalize()
            img_2d = None
            if compound_obj.canonical_smiles:
                img_2d = draw_molecule_pil(compound_obj.canonical_smiles, size=(400,350))
            molecule_details = {"cid": compound_obj.cid, "name": name_to_display, "properties": props, "image_2d": img_2d, "smiles": compound_obj.canonical_smiles}
            status_message = f"Information found for: {name_to_display} (CID: {compound_obj.cid})"
        else: status_message = f"Molecule '{term}' not found on PubChem."
        return molecule_details, status_message
    except pcp.PubChemHTTPError as e: return None, f"PubChem API Error during general search: {e}"
    except Exception as e: return None, f"An unexpected error occurred during general search: {e}\n{traceback.format_exc()}"
