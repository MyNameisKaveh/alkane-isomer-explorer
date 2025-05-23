import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage
import py3Dmol
import os
import traceback
from PIL import Image
from rdkit.Chem import rdDepictor
from rdkit import Chem
from rdkit.Chem import AllChem
from io import StringIO

# --- Helper Functions ---
def draw_molecule_pil(smiles_string, size=(350, 300)):
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if mol:
            rdDepictor.Compute2DCoords(mol)
            img = MolToImage(mol, size=size)
            return img
        else:
            return None
    except Exception as e:
        print(f"Error in draw_molecule_pil for SMILES {smiles_string}: {e}")
        return None

def get_sdf_content(cid):
    if cid is None:
        return None, "CID not provided for SDF content."
    print(f"Fetching SDF content for CID: {cid}...")
    temp_sdf_file_dir = "/tmp"
    if not os.path.exists(temp_sdf_file_dir):
        try:
            os.makedirs(temp_sdf_file_dir, exist_ok=True)
        except OSError:
            temp_sdf_file_dir = "."
    temp_sdf_file = os.path.join(temp_sdf_file_dir, f'temp_3d_structure_{cid}.sdf')
    sdf_data, error_message = None, None
    try:
        pcp.download('SDF', temp_sdf_file, str(cid), 'cid', record_type='3d', overwrite=True)
        with open(temp_sdf_file, 'r') as f:
            sdf_data = f.read()
        if not sdf_data or sdf_data.strip() == "$$$$\n" or sdf_data.strip() == "":
            error_message, sdf_data = f"SDF file for CID {cid} is empty or invalid.", None
    except pcp.NotFoundError:
        error_message = f"3D structure (SDF) for CID {cid} not found on PubChem."
    except FileNotFoundError:
        error_message = f"Temporary SDF file for CID {cid} not found."
    except Exception as e:
        error_message = f"Error downloading/reading SDF for CID {cid}: {e}\n{traceback.format_exc()}"
    finally:
        if os.path.exists(temp_sdf_file):
            try:
                os.remove(temp_sdf_file)
            except Exception:
                pass
    if error_message and 'st' in globals() and hasattr(st, 'warning'):
        st.warning(error_message)
    return sdf_data, error_message

def generate_3d_viewer_html(sdf_data, molecule_name, display_style='stick', component_width=600, component_height=450):
    if not sdf_data:
        return "<p style='color:orange; text-align:center;'>SDF data not available for 3D view.</p>"
    try:
        # Step 1: Load molecule from SDF data using RDKit
        suppl = Chem.SDMolSupplier()
        suppl.SetData(sdf_data, removeHs=False)
        mol = next(suppl)
        if mol:
            # Step 2: Compute the centroid of the molecule
            conf = mol.GetConformer()
            center = AllChem.ComputeCentroid(conf)
            
            # Step 3: Translate all atoms to center the molecule at origin
            for i in range(mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(i)
                new_pos = (pos.x - center.x, pos.y - center.y, pos.z - center.z)
                conf.SetAtomPosition(i, new_pos)
            
            # Step 4: Write the translated molecule back to SDF format
            sio = StringIO()
            writer = Chem.SDWriter(sio)
            writer.write(mol)
            writer.flush()
            centered_sdf_data = sio.getvalue()
        else:
            centered_sdf_data = sdf_data  # Fallback to original if loading fails

        # Step 5: Render the centered molecule in py3Dmol
        viewer = py3Dmol.view(width=component_width, height=component_height)
        viewer.addModel(centered_sdf_data, 'sdf')
        
        # Apply display style
        if display_style == 'stick':
            viewer.setStyle({'stick': {}})
        elif display_style == 'sphere':
            viewer.setStyle({'sphere': {'scale': 0.35}})
        elif display_style == 'line':
            viewer.setStyle({'line': {'linewidth': 1.5, 'colorscheme': 'blackCarbon'}})
        elif display_style == 'ball_and_stick':
            viewer.setStyle({'stick': {'radius': 0.08}, 'sphere': {'scale': 0.25}})
        else:
            viewer.setStyle({'stick': {}})
            
        # Set orthographic projection for better centering
        viewer.setProjection('orthographic')
        viewer.center()  # Center the view on the molecule (now at origin)
        viewer.zoomTo()  # Adjust zoom to fit the molecule
        
        # Manual translation to adjust position (tweak these values as needed)
        viewer.translate(-10, -10, 0)  # Move left and down by 10 units
        
        return viewer._make_html()
    except Exception as e:
        error_msg_html = f"<p style='color:red; text-align:center;'>Error rendering 3D view for {molecule_name}: {str(e)}</p>"
        if 'st' in globals() and hasattr(st, 'error'):
            st.error(f"Error creating 3D viewer for {molecule_name} with style {display_style}: {e}")
        return error_msg_html

def get_compound_properties(compound_obj):
    if not compound_obj:
        return None, None
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
        "Synonyms": ", ".join(compound_obj.synonyms[:5]) + ("..." if len(compound_obj.synonyms) > 5 else "") if compound_obj.synonyms else 'N/A',
    }
    pubchem_link = f"https://pubchem.ncbi.nlm.nih.gov/compound/{compound_obj.cid}" if compound_obj.cid else None
    return properties, pubchem_link

def process_alkane_request(molecule_name_input):
    if not molecule_name_input or not molecule_name_input.strip():
        return [], "Please enter a molecule name.", None
    molecule_name = molecule_name_input.strip().lower()
    status_message, isomer_details_list = f"Searching for '{molecule_name}'...", []
    main_molecule_props_tuple = None
    try:
        compounds = pcp.get_compounds(molecule_name, 'name', listkey_count=5)
        main_compound_obj, molecular_formula = None, None
        if not compounds:
            return [], f"Molecule '{molecule_name}' not found on PubChem.", None
        for c_alk in compounds:
            cid_alk = c_alk.cid
            actual_formula_alk = c_alk.molecular_formula if hasattr(c_alk, 'molecular_formula') else None
            if actual_formula_alk:
                is_standard_hydrocarbon = True
                if c_alk.canonical_smiles:
                    try:
                        mol_obj_alk = Chem.MolFromSmiles(c_alk.canonical_smiles)
                        if mol_obj_alk:
                            rdDepictor.Compute2DCoords(mol_obj_alk)
                            if len(Chem.GetMolFrags(mol_obj_alk)) > 1:
                                is_standard_hydrocarbon = False
                            if is_standard_hydrocarbon:
                                atom_symbols_main_alk = set(atom.GetSymbol() for atom in mol_obj_alk.GetAtoms())
                                if not atom_symbols_main_alk.issubset({'C', 'H'}):
                                    is_standard_hydrocarbon = False
                                for atom_alk in mol_obj_alk.GetAtoms():
                                    if atom_alk.GetIsotope() != 0:
                                        is_standard_hydrocarbon = False
                                        break
                            if is_standard_hydrocarbon:
                                for bond_alk in mol_obj_alk.GetBonds():
                                    if bond_alk.GetBondType() != Chem.BondType.SINGLE:
                                        is_standard_hydrocarbon = False
                                        break
                                if Chem.GetSymmSSSR(mol_obj_alk):
                                    is_standard_hydrocarbon = False
                        else:
                            is_standard_hydrocarbon = False
                    except Exception:
                        is_standard_hydrocarbon = False
                else:
                    is_standard_hydrocarbon = False
                if not is_standard_hydrocarbon:
                    continue
                if molecule_name in [syn.lower() for syn in c_alk.synonyms]:
                    main_compound_obj = c_alk
                    molecular_formula = actual_formula_alk
                    break
                if not main_compound_obj:
                    main_compound_obj = c_alk
                    molecular_formula = actual_formula_alk
        if not main_compound_obj or not molecular_formula:
            return [], f"Standard alkane '{molecule_name}' not found or does not meet criteria.", None
        main_molecule_properties_dict, main_pubchem_link = get_compound_properties(main_compound_obj)
        main_molecule_props_tuple = (main_molecule_properties_dict, main_pubchem_link)
        main_name_alk = main_molecule_properties_dict.get("IUPAC Name", molecule_name.capitalize()) if main_molecule_properties_dict else molecule_name.capitalize()
        status_message = f"Finding isomers for {main_name_alk} (Formula: {molecular_formula})..."
        isomers_found_raw = pcp.get_compounds(molecular_formula, 'formula', listkey_count=50)
        if not isomers_found_raw:
            return [], f"No isomers found for formula {molecular_formula}.", main_molecule_props_tuple
        valid_structural_alkanes_entries, unique_accepted_smiles = [], set()
        for isomer_entry in isomers_found_raw:
            smiles_iso = isomer_entry.canonical_smiles
            if not smiles_iso:
                continue
            try:
                mol_iso_obj = Chem.MolFromSmiles(smiles_iso)
                if not mol_iso_obj:
                    continue
                is_valid_candidate = True
                if len(Chem.GetMolFrags(mol_iso_obj)) > 1:
                    is_valid_candidate = False
                if is_valid_candidate:
                    atom_symbols_iso = set(atom.GetSymbol() for atom in mol_iso_obj.GetAtoms())
                    if not atom_symbols_iso.issubset({'C', 'H'}):
                        is_valid_candidate = False
                if is_valid_candidate:
                    for atom_iso in mol_iso_obj.GetAtoms():
                        if atom_iso.GetSymbol() == 'H' and atom_iso.GetDegree() == 0:
                            is_valid_candidate = False
                            break
                        if atom_iso.GetIsotope() != 0:
                            is_valid_candidate = False
                            break
                if is_valid_candidate:
                    for bond_iso in mol_iso_obj.GetBonds():
                        if bond_iso.GetBondType() != Chem.BondType.SINGLE:
                            is_valid_candidate = False
                            break
                    if Chem.GetSymmSSSR(mol_iso_obj):
                        is_valid_candidate = False
                if is_valid_candidate:
                    canonical_smiles_for_uniqueness = Chem.MolToSmiles(mol_iso_obj, isomericSmiles=False)
                    if canonical_smiles_for_uniqueness not in unique_accepted_smiles:
                        valid_structural_alkanes_entries.append(isomer_entry)
                        unique_accepted_smiles.add(canonical_smiles_for_uniqueness)
            except Exception:
                continue
        if not valid_structural_alkanes_entries:
            return [], f"No valid alkane isomers found for {molecular_formula} after filtering.", main_molecule_props_tuple
        for entry_iso in valid_structural_alkanes_entries:
            pil_image = draw_molecule_pil(entry_iso.canonical_smiles)
            if pil_image:
                name_iso = entry_iso.iupac_name
                if not name_iso and entry_iso.synonyms:
                    simple_names_iso = [s for s in entry_iso.synonyms if s.lower().endswith("ane") and not any(char.isdigit() for char in s.split('-')[0]) and '-' not in s.split(' ')[0]]
                    if simple_names_iso:
                        name_iso = min(simple_names_iso, key=len)
                    else:
                        name_iso = entry_iso.synonyms[0]
                elif not name_iso:
                    name_iso = f"Alkane (CID: {entry_iso.cid})"
                isomer_details_list.append({"cid": entry_iso.cid, "name": name_iso.capitalize(), "smiles": entry_iso.canonical_smiles, "image": pil_image})
        if isomer_details_list:
            status_message = f"{len(isomer_details_list)} isomers found for '{molecule_name}' ({molecular_formula})."
        else:
            status_message = f"No displayable isomers found for '{molecule_name}'."
        return isomer_details_list, status_message, main_molecule_props_tuple
    except pcp.PubChemHTTPError as e:
        return [], f"PubChem API Error: {e}", None
    except Exception as e:
        return [], f"An unexpected error occurred: {e}\n{traceback.format_exc()}", None

def process_general_molecule_search(search_term):
    if not search_term or not search_term.strip():
        return None, "Please enter a chemical name or identifier."
    term = search_term.strip()
    status_message = f"Searching for '{term}'..."
    molecule_details_dict = None
    try:
        compounds_gen = pcp.get_compounds(term, 'name', listkey_count=1)
        if compounds_gen:
            compound_obj_gen = compounds_gen[0]
            props_dict_gen, pubchem_link_gen = get_compound_properties(compound_obj_gen)
            name_to_display_gen = props_dict_gen.get("IUPAC Name", term.capitalize()) if props_dict_gen else term.capitalize()
            img_2d_gen = None
            if compound_obj_gen.canonical_smiles:
                img_2d_gen = draw_molecule_pil(compound_obj_gen.canonical_smiles, size=(400,350))
            molecule_details_dict = {"cid": compound_obj_gen.cid, "name": name_to_display_gen, "properties_dict": props_dict_gen,"pubchem_link": pubchem_link_gen, "image_2d": img_2d_gen, "smiles": compound_obj_gen.canonical_smiles}
            status_message = f"Information found for: {name_to_display_gen} (CID: {compound_obj_gen.cid})"
        else:
            status_message = f"Molecule '{term}' not found on PubChem."
        return molecule_details_dict, status_message
    except pcp.PubChemHTTPError as e:
        return None, f"PubChem API Error during general search: {e}"
    except Exception as e:
        return None, f"An unexpected error occurred during general search: {e}\n{traceback.format_exc()}"

# --- Streamlit UI ---
st.set_page_config(page_title="Chemical Compound Explorer", layout="wide", initial_sidebar_state="expanded")
st.title("Chemical Compound Explorer")

if 'alkane_isomer_data' not in st.session_state:
    st.session_state.alkane_isomer_data = []
if 'alkane_main_molecule_props_tuple' not in st.session_state:
    st.session_state.alkane_main_molecule_props_tuple = None
if 'selected_isomer_cid_for_3d' not in st.session_state:
    st.session_state.selected_isomer_cid_for_3d = None
if 'selected_isomer_name_for_3d' not in st.session_state:
    st.session_state.selected_isomer_name_for_3d = ""
if 'current_isomer_3d_index' not in st.session_state:
    st.session_state.current_isomer_3d_index = -1
if 'alkane_molecule_searched' not in st.session_state:
    st.session_state.alkane_molecule_searched = ""
if 'alkane_name_input' not in st.session_state:
    st.session_state.alkane_name_input = ""
if 'run_alkane_search_after_example' not in st.session_state:
    st.session_state.run_alkane_search_after_example = False
if 'general_molecule_data_dict' not in st.session_state:
    st.session_state.general_molecule_data_dict = None
if 'general_molecule_name_input' not in st.session_state:
    st.session_state.general_molecule_name_input = ""
if 'status_message' not in st.session_state:
    st.session_state.status_message = ""
if 'selected_3d_style' not in st.session_state:
    st.session_state.selected_3d_style = 'stick'
if 'last_search_type' not in st.session_state:
    st.session_state.last_search_type = None

st.sidebar.header("Search Options")
st.sidebar.subheader("1. Alkane Isomer Search")
current_alkane_input = st.sidebar.text_input(label="Enter alkane name (for isomers):",value=st.session_state.alkane_name_input,placeholder="e.g., butane, pentane",key="sidebar_alkane_input",help="Finds structural isomers for standard alkanes.")
if current_alkane_input != st.session_state.alkane_name_input:
    st.session_state.alkane_name_input = current_alkane_input
    st.session_state.run_alkane_search_after_example = False
example_alkanes = ["butane", "pentane", "hexane"]
st.sidebar.caption("Alkane Examples:")
cols_examples_alkane = st.sidebar.columns(len(example_alkanes))
for i_ex_alk, example_alk in enumerate(example_alkanes):
    if cols_examples_alkane[i_ex_alk].button(example_alk.capitalize(), key=f"example_alkane_{example_alk}", use_container_width=True):
        st.session_state.alkane_name_input = example_alk
        st.session_state.run_alkane_search_after_example = True
        st.rerun()
if st.sidebar.button("Search Alkane Isomers", type="primary", key="search_alkane_button", use_container_width=True):
    if st.session_state.alkane_name_input:
        st.session_state.last_search_type = "alkane"
        with st.spinner(f"Searching isomers for {st.session_state.alkane_name_input}..."):
            isomers_res, status_msg_res, main_props_tuple_res = process_alkane_request(st.session_state.alkane_name_input)
            st.session_state.alkane_isomer_data, st.session_state.status_message, st.session_state.alkane_main_molecule_props_tuple = isomers_res, status_msg_res, main_props_tuple_res
            st.session_state.selected_isomer_cid_for_3d, st.session_state.selected_isomer_name_for_3d, st.session_state.current_isomer_3d_index = None, "", -1
            st.session_state.alkane_molecule_searched = st.session_state.alkane_name_input
            st.session_state.general_molecule_data_dict = None
    else:
        st.session_state.status_message = "Please enter an alkane name for isomer search."
if st.session_state.run_alkane_search_after_example and st.session_state.alkane_name_input:
    st.session_state.last_search_type = "alkane"
    st.session_state.run_alkane_search_after_example = False
    with st.spinner(f"Searching isomers for {st.session_state.alkane_name_input}..."):
        isomers_res_ex, status_msg_res_ex, main_props_tuple_res_ex = process_alkane_request(st.session_state.alkane_name_input)
        st.session_state.alkane_isomer_data, st.session_state.status_message, st.session_state.alkane_main_molecule_props_tuple = isomers_res_ex, status_msg_res_ex, main_props_tuple_res_ex
        st.session_state.selected_isomer_cid_for_3d, st.session_state.selected_isomer_name_for_3d, st.session_state.current_isomer_3d_index = None, "", -1
        st.session_state.alkane_molecule_searched = st.session_state.alkane_name_input
        st.session_state.general_molecule_data_dict = None
    st.rerun()
st.sidebar.markdown("---")
st.sidebar.subheader("2. General Molecule Information")
general_molecule_search_term = st.sidebar.text_input(label="Enter any chemical name/identifier:",value=st.session_state.general_molecule_name_input,placeholder="e.g., aspirin, C6H6",key="sidebar_general_input",help="Get info, 2D, and 3D structure.")
if general_molecule_search_term != st.session_state.general_molecule_name_input:
    st.session_state.general_molecule_name_input = general_molecule_search_term
if st.sidebar.button("Search Molecule Info", type="primary", key="search_general_button", use_container_width=True):
    if st.session_state.general_molecule_name_input:
        st.session_state.last_search_type = "general"
        with st.spinner(f"Searching for {st.session_state.general_molecule_name_input}..."):
            mol_data_res, status_msg_g_res = process_general_molecule_search(st.session_state.general_molecule_name_input)
            st.session_state.general_molecule_data_dict, st.session_state.status_message = mol_data_res, status_msg_g_res
            st.session_state.alkane_isomer_data, st.session_state.alkane_main_molecule_props_tuple = [], None
            st.session_state.selected_isomer_cid_for_3d, st.session_state.selected_isomer_name_for_3d, st.session_state.current_isomer_3d_index = None, "", -1
    else:
        st.session_state.status_message = "Please enter a chemical name for general search."
if st.session_state.status_message:
    is_error = any(keyword in st.session_state.status_message.lower() for keyword in ["error", "not found", "empty", "invalid"])
    if is_error:
        st.sidebar.error(st.session_state.status_message)
    else:
        st.sidebar.info(st.session_state.status_message)

gallery_title, props_title, view2d_title, view3d_title = "Isomer Gallery", "Properties", "2D Structure", "3D View"
tab_objects = []
if st.session_state.last_search_type == "alkane":
    main_props_dict_tuple = st.session_state.alkane_main_molecule_props_tuple if st.session_state.alkane_main_molecule_props_tuple else (None, None)
    main_props_dict_alk, _ = main_props_dict_tuple
    if st.session_state.alkane_isomer_data:
        gallery_title = f"Alkane Isomers ({len(st.session_state.alkane_isomer_data)})"
    if main_props_dict_alk:
        props_title = f"Alkane: {main_props_dict_alk.get('IUPAC Name', st.session_state.alkane_molecule_searched.capitalize())[:20]}"
    if st.session_state.selected_isomer_cid_for_3d:
        view3d_title = f"Isomer 3D: {st.session_state.selected_isomer_name_for_3d[:20]}"
    if st.session_state.alkane_isomer_data or main_props_dict_alk or st.session_state.selected_isomer_cid_for_3d:
        tab_objects = st.tabs([gallery_title, props_title, view3d_title])
    else:
        if st.session_state.alkane_molecule_searched:
            st.info(f"No results found for '{st.session_state.alkane_molecule_searched}'. Try another search.")
        else:
            st.info("Use the sidebar to search for an alkane or a general chemical compound.")
elif st.session_state.last_search_type == "general" and st.session_state.general_molecule_data_dict:
    g_data_dict = st.session_state.general_molecule_data_dict
    props_title = f"Properties: {g_data_dict['name'][:20]}"
    view2d_title = f"2D: {g_data_dict['name'][:20]}" if g_data_dict.get("image_2d") else None
    view3d_title = f"3D: {g_data_dict['name'][:20]}"
    tabs_to_create = [props_title]
    if view2d_title:
        tabs_to_create.append(view2d_title)
    tabs_to_create.append(view3d_title)
    tab_objects = st.tabs(tabs_to_create)
else:
    st.info("To begin, use the search options in the sidebar.")

if st.session_state.last_search_type == "alkane" and tab_objects and len(tab_objects) == 3:
    tab_gallery, tab_main_props, tab_3d_isomer = tab_objects[0], tab_objects[1], tab_objects[2]
    with tab_gallery:
        if st.session_state.alkane_isomer_data:
            st.subheader(f"Isomers found for: {st.session_state.alkane_molecule_searched.capitalize()}")
            num_columns_gallery = 3
            gallery_cols = st.columns(num_columns_gallery)
            for i_iso, isomer_item in enumerate(st.session_state.alkane_isomer_data):
                with gallery_cols[i_iso % num_columns_gallery]:
                    container = st.container(border=True)
                    container.image(isomer_item["image"], caption=f"{isomer_item['name']}", use_container_width=True)
                    container.markdown(f"<small>SMILES: {isomer_item['smiles']}<br>CID: {isomer_item['cid']}</small>", unsafe_allow_html=True)
                    if container.button(f"View 3D", key=f"btn_3d_isomer_{isomer_item['cid']}"):
                        st.session_state.selected_isomer_cid_for_3d, st.session_state.selected_isomer_name_for_3d, st.session_state.current_isomer_3d_index = isomer_item['cid'], isomer_item['name'], i_iso
                        st.rerun()
    if st.session_state.alkane_main_molecule_props_tuple:
        main_props_dict_tuple_val, main_pubchem_link_val = st.session_state.alkane_main_molecule_props_tuple
        if tab_main_props and main_props_dict_tuple_val:
            with tab_main_props:
                main_mol_name_val = main_props_dict_tuple_val.get("IUPAC Name", st.session_state.alkane_molecule_searched.capitalize())
                st.subheader(f"Properties for Searched Alkane: {main_mol_name_val}")
                if main_pubchem_link_val:
                    st.markdown(f"**View on PubChem:** [{main_mol_name_val}]({main_pubchem_link_val})", unsafe_allow_html=True)
                    st.markdown("---")
                prop_cols = st.columns(2)
                prop_list_val = list(main_props_dict_tuple_val.items())
                for i_prop, (key_prop, value_prop) in enumerate(prop_list_val):
                    if value_prop != 'N/A' and value_prop is not None:
                        prop_cols[i_prop % 2].markdown(f"**{key_prop}:** {value_prop}")
                st.markdown("---")
    if tab_3d_isomer:
        with tab_3d_isomer:
            if st.session_state.selected_isomer_cid_for_3d:
                st.subheader(f"3D Structure for Isomer: {st.session_state.selected_isomer_name_for_3d}")
                style_options_map = {'Stick': 'stick', 'Line': 'line', 'Ball and Stick': 'ball_and_stick'}
                style_labels = list(style_options_map.keys())
                try:
                    current_style_label = [k_style for k_style, v_style in style_options_map.items() if v_style == st.session_state.selected_3d_style][0]
                    current_style_index = style_labels.index(current_style_label)
                except IndexError:
                    current_style_index = 0
                    if style_labels:
                        st.session_state.selected_3d_style = style_options_map[style_labels[0]]
                    else:
                        st.session_state.selected_3d_style = 'stick'
                selected_style_label = st.radio("Select Display Style:", options=style_labels, key="radio_3d_style_isomer", horizontal=True, index=current_style_index)
                if style_options_map.get(selected_style_label) != st.session_state.selected_3d_style:
                    st.session_state.selected_3d_style = style_options_map[selected_style_label]
                    st.rerun()
                nav_col_layout = [1, 0.1, 1.5, 0.1, 1]
                prev_col, _, clear_col, _, next_col = st.columns(nav_col_layout)
                with prev_col:
                    if st.button("⬅️ Previous", key="prev_3d_isomer", help="View Previous Isomer", use_container_width=True, disabled=(st.session_state.current_isomer_3d_index <= 0)):
                        st.session_state.current_isomer_3d_index -= 1
                        prev_isomer = st.session_state.alkane_isomer_data[st.session_state.current_isomer_3d_index]
                        st.session_state.selected_isomer_cid_for_3d, st.session_state.selected_isomer_name_for_3d = prev_isomer['cid'], prev_isomer['name']
                        st.rerun()
                with clear_col:
                    if st.button("Clear Isomer 3D View", key="clear_3d_isomer", use_container_width=True):
                        st.session_state.selected_isomer_cid_for_3d, st.session_state.selected_isomer_name_for_3d, st.session_state.current_isomer_3d_index = None, "", -1
                        st.rerun()
                with next_col:
                    if st.button("Next ➡️", key="next_3d_isomer", help="View Next Isomer", use_container_width=True, disabled=(st.session_state.current_isomer_3d_index >= len(st.session_state.alkane_isomer_data) - 1)):
                        st.session_state.current_isomer_3d_index += 1
                        next_isomer = st.session_state.alkane_isomer_data[st.session_state.current_isomer_3d_index]
                        st.session_state.selected_isomer_cid_for_3d, st.session_state.selected_isomer_name_for_3d = next_isomer['cid'], next_isomer['name']
                        st.rerun()
                st.markdown("---")
                with st.spinner(f"Loading 3D structure for {st.session_state.selected_isomer_name_for_3d}..."):
                    sdf_data_iso, error_iso = get_sdf_content(st.session_state.selected_isomer_cid_for_3d)
                    if sdf_data_iso:
                        html_3d_iso = generate_3d_viewer_html(sdf_data_iso, st.session_state.selected_isomer_name_for_3d, display_style=st.session_state.selected_3d_style, component_width=600, component_height=450)
                        st.components.v1.html(html_3d_iso, height=470, width=600, scrolling=False)
            else:
                if st.session_state.alkane_molecule_searched and st.session_state.alkane_isomer_data:
                    st.info("Select an isomer from the 'Isomer Gallery' tab to view its 3D structure.")
                elif st.session_state.alkane_molecule_searched:
                    st.info("No isomers were found for the searched alkane to display in 3D.")

if st.session_state.last_search_type == "general" and st.session_state.general_molecule_data_dict and tab_objects:
    g_data_dict_val = st.session_state.general_molecule_data_dict
    num_general_tabs_created = len(tab_objects)
    props_tab_index = 0
    view2d_tab_index = 1 if g_data_dict_val.get("image_2d") and num_general_tabs_created > (props_tab_index + 1) else -1
    view3d_tab_index = (view2d_tab_index + 1) if view2d_tab_index != -1 and num_general_tabs_created > view2d_tab_index else (props_tab_index + 1)
    
    if num_general_tabs_created > props_tab_index:
        with tab_objects[props_tab_index]:
            st.subheader(f"Properties for: {g_data_dict_val['name']}")
            if g_data_dict_val["properties_dict"]:
                props_dict_g = g_data_dict_val["properties_dict"]
                if g_data_dict_val.get("pubchem_link"):
                    st.markdown(f"**View on PubChem:** [{g_data_dict_val['name']}]({g_data_dict_val['pubchem_link']})", unsafe_allow_html=True)
                    st.markdown("---")
                prop_cols = st.columns(2)
                prop_list_g = list(props_dict_g.items())
                for i_prop_g, (key_prop_g, value_prop_g) in enumerate(prop_list_g):
                    if value_prop_g != 'N/A' and value_prop_g is not None:
                        if key_prop_g == "Synonyms" and isinstance(value_prop_g, str) and len(value_prop_g) > 100:
                            prop_cols[i_prop_g % 2].markdown(f"**{key_prop_g}:** {value_prop_g[:100]}...")
                        else:
                            prop_cols[i_prop_g % 2].markdown(f"**{key_prop_g}:** {value_prop_g}")
            else:
                st.info("No detailed properties found.")
            st.markdown("---")
    if view2d_tab_index != -1 and num_general_tabs_created > view2d_tab_index:
        with tab_objects[view2d_tab_index]:
            if g_data_dict_val.get("image_2d"):
                st.subheader(f"2D Structure for: {g_data_dict_val['name']}")
                st.image(g_data_dict_val["image_2d"], use_container_width=True)
            else:
                st.info("2D image not available for this compound.")
    if num_general_tabs_created > view3d_tab_index:
        with tab_objects[view3d_tab_index]:
            st.subheader(f"3D Structure for: {g_data_dict_val['name']}")
            style_options_map_g = {'Stick': 'stick', 'Line': 'line', 'Ball and Stick': 'ball_and_stick'}
            style_labels_g = list(style_options_map_g.keys())
            try:
                current_style_label_g = [k_style_g for k_style_g, v_style_g in style_options_map_g.items() if v_style_g == st.session_state.selected_3d_style][0]
                current_style_index_g = style_labels_g.index(current_style_label_g)
            except IndexError:
                current_style_index_g = 0
                if style_labels_g:
                    st.session_state.selected_3d_style = style_options_map_g[style_labels_g[0]]
                else:
                    st.session_state.selected_3d_style = 'stick'
            selected_style_label_g = st.radio("Select Display Style:", options=style_labels_g, key="radio_3d_style_general", horizontal=True, index=current_style_index_g)
            if style_options_map_g.get(selected_style_label_g) != st.session_state.selected_3d_style:
                st.session_state.selected_3d_style = style_options_map_g[selected_style_label_g]
                st.rerun()
            st.markdown("---")
            with st.spinner(f"Loading 3D structure for {g_data_dict_val['name']}..."):
                sdf_data_g, error_g = get_sdf_content(g_data_dict_val['cid'])
                if sdf_data_g:
                    html_3d_g = generate_3d_viewer_html(sdf_data_g, g_data_dict_val['name'], display_style=st.session_state.selected_3d_style, component_width=600, component_height=450)
                    st.components.v1.html(html_3d_g, height=470, width=600, scrolling=False)
                else:
                    st.info(f"3D structure could not be loaded. {error_g if error_g else ''}")

st.sidebar.markdown("---")
st.sidebar.markdown("Built with Streamlit, RDKit, PubChemPy, and py3Dmol.")
