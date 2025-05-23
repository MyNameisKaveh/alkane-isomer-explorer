import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage
import py3Dmol
import os
import traceback
from PIL import Image
from rdkit.Chem import rdDepictor

# --- Helper Functions (توابع کمکی مانند قبل، فقط process_alkane_request تغییر می‌کند) ---
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
    if error_message and 'st' in globals() and hasattr(st, 'warning'): 
        st.warning(error_message)
    return sdf_data, error_message

def generate_3d_viewer_html(sdf_data, molecule_name, display_style='stick', width=500, height=400):
    if not sdf_data: return "<p style='color:orange; text-align:center;'>SDF data not available for 3D view.</p>"
    try:
        viewer = py3Dmol.view(width=width, height=height)
        viewer.addModel(sdf_data, 'sdf')
        if display_style == 'stick': viewer.setStyle({'stick': {}})
        elif display_style == 'sphere': viewer.setStyle({'sphere': {'scale': 0.35}})
        elif display_style == 'line': viewer.setStyle({'line': {'linewidth': 2.0, 'colorscheme': 'blackCarbon'}})
        elif display_style == 'ball_and_stick': viewer.setStyle({'stick': {'radius': 0.08}, 'sphere': {'scale': 0.25}})
        else: viewer.setStyle({'stick': {}})
        viewer.setBackgroundColor('0xeeeeee')
        viewer.zoomTo()
        return viewer._make_html()
    except Exception as e:
        error_msg_html = f"<p style='color:red; text-align:center;'>Error rendering 3D view for {molecule_name}: {str(e)}</p>"
        if 'st' in globals() and hasattr(st, 'error'): 
            st.error(f"Error creating 3D viewer for {molecule_name} with style {display_style}: {e}")
        return error_msg_html

def get_compound_properties(compound_obj):
    """Extracts selected properties from a PubChemPy Compound object."""
    if not compound_obj:
        return None
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
    # Filter out N/A values if desired, or keep them
    # return {k: v for k, v in properties.items() if v != 'N/A'}
    return properties

def process_alkane_request(molecule_name_input):
    if not molecule_name_input or not molecule_name_input.strip():
        return [], "Please enter a molecule name.", None # Isomers, Status, Main Molecule Properties
    molecule_name = molecule_name_input.strip().lower()
    status_message, isomer_details_list = f"Searching for '{molecule_name}'...", []
    main_molecule_properties = None
    try:
        compounds = pcp.get_compounds(molecule_name, 'name', listkey_count=5) # Get a few matches
        main_compound_obj, molecular_formula = None, None
        if not compounds: return [], f"Molecule '{molecule_name}' not found on PubChem.", None
        
        # Select the most relevant main compound (e.g., first exact match or first valid alkane)
        for c in compounds:
            # ... (منطق قبلی برای فیلتر کردن و انتخاب main_compound_obj) ...
            cid = c.cid
            actual_formula = c.molecular_formula if hasattr(c, 'molecular_formula') else None
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
                # Prioritize exact name match if possible
                if molecule_name in [syn.lower() for syn in c.synonyms]: 
                    main_compound_obj = c; molecular_formula = actual_formula; break
                if not main_compound_obj: # Tentatively select the first valid one
                    main_compound_obj = c; molecular_formula = actual_formula
        
        if not main_compound_obj or not molecular_formula:
            return [], f"Standard alkane '{molecule_name}' not found or does not meet criteria.", None
        
        # Get properties of the main identified compound
        main_molecule_properties = get_compound_properties(main_compound_obj)
        
        main_name = main_molecule_properties.get("IUPAC Name", molecule_name.capitalize()) if main_molecule_properties else molecule_name.capitalize()
        status_message = f"Finding isomers for {main_name} (Formula: {molecular_formula})..."
        isomers_found_raw = pcp.get_compounds(molecular_formula, 'formula', listkey_count=50)
        if not isomers_found_raw: return [], f"No isomers found for formula {molecular_formula}.", main_molecule_properties
        
        valid_structural_alkanes_entries, unique_accepted_smiles = [], set()
        # ... (منطق فیلتر کردن ایزومرها مانند قبل) ...
        for isomer_entry in isomers_found_raw:
            smiles = isomer_entry.canonical_smiles
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

# --- Streamlit UI ---
st.set_page_config(page_title="Alkane Isomer Viewer", layout="wide", initial_sidebar_state="expanded")
st.title("Alkane Isomer Finder and Viewer")

# Initialize session state
if 'isomer_data' not in st.session_state: st.session_state.isomer_data = []
if 'main_molecule_props' not in st.session_state: st.session_state.main_molecule_props = None # For main molecule properties
if 'status_message' not in st.session_state: st.session_state.status_message = ""
if 'selected_cid_for_3d' not in st.session_state: st.session_state.selected_cid_for_3d = None
if 'selected_name_for_3d' not in st.session_state: st.session_state.selected_name_for_3d = ""
if 'current_3d_index' not in st.session_state: st.session_state.current_3d_index = -1
if 'selected_3d_style' not in st.session_state: st.session_state.selected_3d_style = 'stick' 
if 'molecule_searched' not in st.session_state: st.session_state.molecule_searched = ""
if 'molecule_name_input' not in st.session_state: st.session_state.molecule_name_input = "" 
if 'run_search_after_example' not in st.session_state: st.session_state.run_search_after_example = False

# --- Sidebar ---
st.sidebar.header("Search and Examples")
current_molecule_input = st.sidebar.text_input(
    label="Enter alkane name:", value=st.session_state.molecule_name_input,
    placeholder="e.g., butane, pentane", key="sidebar_molecule_input",
    help="Enter the name of the alkane (e.g., hexane) in lowercase English."
)
if current_molecule_input != st.session_state.molecule_name_input:
    st.session_state.molecule_name_input = current_molecule_input
    st.session_state.run_search_after_example = False 

example_alkanes = ["butane", "pentane", "hexane", "heptane", "octane"]
st.sidebar.subheader("Or select an example:")
cols_examples = st.sidebar.columns(len(example_alkanes) if len(example_alkanes) <=3 else 3) 
for i, example in enumerate(example_alkanes):
    if cols_examples[i % len(cols_examples)].button(example.capitalize(), key=f"example_{example}"):
        st.session_state.molecule_name_input = example 
        st.session_state.selected_cid_for_3d, st.session_state.selected_name_for_3d, st.session_state.current_3d_index = None, "", -1
        st.session_state.main_molecule_props = None # Clear old main molecule props
        st.session_state.run_search_after_example = True 
        st.rerun() 

search_button_sidebar = st.sidebar.button("Search Isomers", type="primary", key="sidebar_search_button")

should_run_search = False
if search_button_sidebar and st.session_state.molecule_name_input: should_run_search = True
elif st.session_state.run_search_after_example and st.session_state.molecule_name_input:
    should_run_search = True
    st.session_state.run_search_after_example = False 

if should_run_search:
    with st.spinner(f"Processing for {st.session_state.molecule_name_input}..."):
        isomers, status_msg, main_props = process_alkane_request(st.session_state.molecule_name_input)
        st.session_state.isomer_data = isomers
        st.session_state.status_message = status_msg
        st.session_state.main_molecule_props = main_props # Store main molecule properties
        st.session_state.selected_cid_for_3d, st.session_state.selected_name_for_3d, st.session_state.current_3d_index = None, "", -1
        st.session_state.molecule_searched = st.session_state.molecule_name_input
elif search_button_sidebar and not st.session_state.molecule_name_input: 
    st.session_state.status_message = "Please enter a molecule name or select an example."
    st.session_state.isomer_data, st.session_state.main_molecule_props = [], None
    st.session_state.selected_cid_for_3d, st.session_state.selected_name_for_3d, st.session_state.current_3d_index = None, "", -1
    st.session_state.molecule_searched = ""

if st.session_state.status_message:
    if "Error" in st.session_state.status_message or "not found" in st.session_state.status_message or "empty" in st.session_state.status_message:
        st.sidebar.error(st.session_state.status_message)
    else: st.sidebar.info(st.session_state.status_message)

# --- Main Page Tabs ---
if st.session_state.isomer_data or st.session_state.selected_cid_for_3d or st.session_state.main_molecule_props :
    gallery_tab_title = f"Isomer Gallery ({len(st.session_state.isomer_data)} found)" if st.session_state.isomer_data else "Isomer Gallery"
    viewer_3d_tab_title = f"3D View: {st.session_state.selected_name_for_3d}" if st.session_state.selected_cid_for_3d and st.session_state.selected_name_for_3d else "3D View"
    
    tab_gallery, tab_3d_viewer = st.tabs([gallery_tab_title, viewer_3d_tab_title])

    with tab_gallery:
        # Display Main Molecule Properties if available
        if st.session_state.main_molecule_props:
            props = st.session_state.main_molecule_props
            main_mol_name = props.get("IUPAC Name", st.session_state.molecule_searched.capitalize())
            with st.expander(f"Properties for Searched Molecule: {main_mol_name}", expanded=True):
                # st.write(props) # نمایش خام دیکشنری
                prop_cols = st.columns(2)
                prop_list = list(props.items())
                for i, (key, value) in enumerate(prop_list):
                    if value != 'N/A' and value is not None: # فقط نمایش ویژگی‌های موجود
                        prop_cols[i % 2].markdown(f"**{key}:** {value}")
                st.markdown("---")


        if st.session_state.isomer_data:
            st.subheader(f"Isomers found for: {st.session_state.molecule_searched.capitalize()}")
            num_columns_gallery = 3 
            gallery_cols = st.columns(num_columns_gallery)
            for i, isomer in enumerate(st.session_state.isomer_data):
                with gallery_cols[i % num_columns_gallery]:
                    container = st.container(border=True) 
                    container.image(isomer["image"], caption=f"{isomer['name']}", use_container_width=True) 
                    container.markdown(f"<small>SMILES: {isomer['smiles']}<br>CID: {isomer['cid']}</small>", unsafe_allow_html=True)
                    if container.button(f"View 3D", key=f"btn_3d_tab_{isomer['cid']}"):
                        st.session_state.selected_cid_for_3d = isomer['cid']
                        st.session_state.selected_name_for_3d = isomer['name'] 
                        st.session_state.current_3d_index = i 
                        st.rerun() 
        elif st.session_state.molecule_searched:
            st.info("No isomers to display in the gallery for the current search.")
        else:
            st.info("Search for an alkane or select an example to see the isomer gallery.")

    with tab_3d_viewer:
        if st.session_state.selected_cid_for_3d:
            st.subheader(f"3D Structure for: {st.session_state.selected_name_for_3d}")
            style_options_map = {'Stick': 'stick', 'Line': 'line', 'Ball and Stick': 'ball_and_stick'}
            style_labels = list(style_options_map.keys())
            try:
                current_style_label = [k for k, v in style_options_map.items() if v == st.session_state.selected_3d_style][0]
                current_style_index = style_labels.index(current_style_label)
            except IndexError: 
                current_style_index = 0 
                if style_labels: st.session_state.selected_3d_style = style_options_map[style_labels[0]]
                else: st.session_state.selected_3d_style = 'stick'
            selected_style_label = st.radio("Select Display Style:", options=style_labels, key="radio_3d_style", horizontal=True, index=current_style_index)
            if style_options_map.get(selected_style_label) != st.session_state.selected_3d_style:
                st.session_state.selected_3d_style = style_options_map[selected_style_label]
                st.rerun() 
            nav_col_layout = [1, 0.1, 1.5, 0.1, 1]; prev_col, _, clear_col, _, next_col = st.columns(nav_col_layout)
            with prev_col:
                if st.button("⬅️ Previous", key="prev_3d_icon", help="View Previous Isomer", use_container_width=True, disabled=(st.session_state.current_3d_index <= 0)):
                    st.session_state.current_3d_index -= 1; prev_isomer = st.session_state.isomer_data[st.session_state.current_3d_index]
                    st.session_state.selected_cid_for_3d, st.session_state.selected_name_for_3d = prev_isomer['cid'], prev_isomer['name']; st.rerun()
            with clear_col:
                if st.button("Clear 3D View", key="clear_3d_view_tab_main", use_container_width=True):
                    st.session_state.selected_cid_for_3d, st.session_state.selected_name_for_3d, st.session_state.current_3d_index = None, "", -1; st.rerun()
            with next_col:
                if st.button("Next ➡️", key="next_3d_icon", help="View Next Isomer", use_container_width=True, disabled=(st.session_state.current_3d_index >= len(st.session_state.isomer_data) - 1)):
                    st.session_state.current_3d_index += 1; next_isomer = st.session_state.isomer_data[st.session_state.current_3d_index]
                    st.session_state.selected_cid_for_3d, st.session_state.selected_name_for_3d = next_isomer['cid'], next_isomer['name']; st.rerun()
            st.markdown("---") 
            with st.spinner(f"Loading 3D structure for {st.session_state.selected_name_for_3d} ({st.session_state.selected_3d_style} style)..."):
                sdf_data, error = get_sdf_content(st.session_state.selected_cid_for_3d)
                if sdf_data:
                    html_3d = generate_3d_viewer_html(sdf_data, st.session_state.selected_name_for_3d, display_style=st.session_state.selected_3d_style, width=600, height=450)
                    st.components.v1.html(html_3d, height=470, width=620, scrolling=False)
        else: st.info("Select an isomer from the 'Isomer Gallery' tab to view its 3D structure.")
else: 
    if st.session_state.molecule_searched and not st.session_state.isomer_data: pass 
    else: st.info("To begin, enter an alkane name in the sidebar and click 'Search Isomers', or select an example.")

st.sidebar.markdown("---")
st.sidebar.markdown("Built with Streamlit, RDKit, PubChemPy, and py3Dmol.")
