import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage
import py3Dmol
import os
import traceback
from PIL import Image
from rdkit.Chem import rdDepictor # For better 2D coordinates
# from rdkit.Chem.Draw import rdMolDraw2D # Would be needed for SVG or more advanced drawing options

# --- Helper Functions ---
def draw_molecule_pil(smiles_string, size=(350, 300)): # Increased size for better quality
    """Draws a molecule from a SMILES string using RDKit and returns a PIL Image object."""
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if mol:
            # Generate 2D coordinates for better drawing
            rdDepictor.Compute2DCoords(mol)
            img = MolToImage(mol, size=size)
            return img
        else:
            # Error logging can be done in the main processing function
            return None
    except Exception as e:
        print(f"Error in draw_molecule_pil for SMILES {smiles_string}: {e}") # Server-side log
        return None

def get_sdf_content(cid):
    """Fetches the 3D SDF content for a given CID from PubChem."""
    if cid is None:
        return None, "CID not provided for SDF content."
    print(f"Fetching SDF content for CID: {cid}...") # Server-side log
    temp_sdf_file_dir = "/tmp" 
    if not os.path.exists(temp_sdf_file_dir):
        try: os.makedirs(temp_sdf_file_dir, exist_ok=True)
        except OSError: temp_sdf_file_dir = "." # Fallback
    temp_sdf_file = os.path.join(temp_sdf_file_dir, f'temp_3d_structure_{cid}.sdf')
    sdf_data, error_message = None, None
    try:
        pcp.download('SDF', temp_sdf_file, str(cid), 'cid', record_type='3d', overwrite=True)
        with open(temp_sdf_file, 'r') as f: sdf_data = f.read()
        if not sdf_data or sdf_data.strip() == "$$$$\n" or sdf_data.strip() == "":
            error_message, sdf_data = f"SDF file for CID {cid} is empty or invalid.", None
    except pcp.NotFoundError: error_message = f"3D structure (SDF) for CID {cid} not found on PubChem."
    except FileNotFoundError: error_message = f"Temporary SDF file for CID {cid} not found after download attempt."
    except Exception as e: error_message = f"Error downloading/reading SDF for CID {cid}: {e}\n{traceback.format_exc()}"
    finally:
        if os.path.exists(temp_sdf_file):
            try: os.remove(temp_sdf_file)
            except Exception: pass
    
    # Display warning in Streamlit UI if an error occurred
    if error_message and 'st' in globals() and hasattr(st, 'warning'): 
        st.warning(error_message)
    return sdf_data, error_message

def generate_3d_viewer_html(sdf_data, molecule_name, width=500, height=400):
    """Generates the HTML for the py3Dmol viewer from SDF data."""
    if not sdf_data: return "<p style='color:orange;'>SDF data not available for 3D view.</p>"
    try:
        viewer = py3Dmol.view(width=width, height=height)
        viewer.addModel(sdf_data, 'sdf')
        viewer.setStyle({'stick': {}})
        viewer.setBackgroundColor('0xeeeeee')
        viewer.zoomTo()
        return viewer._make_html()
    except Exception as e:
        if 'st' in globals() and hasattr(st, 'error'): 
            st.error(f"Error creating 3D viewer for {molecule_name}: {e}")
        return f"<p style='color:red;'>Error rendering 3D view for {molecule_name}: {e}</p>"

def process_alkane_request(molecule_name_input):
    """
    Processes an alkane name, finds its isomers, and returns a list of dictionaries
    containing isomer details (PIL image, name, CID, SMILES) and a status message.
    """
    if not molecule_name_input or not molecule_name_input.strip():
        return [], "Please enter a molecule name."
    molecule_name = molecule_name_input.strip().lower()
    status_message = f"Searching for '{molecule_name}'..."
    isomer_details_list = []
    try:
        compounds = pcp.get_compounds(molecule_name, 'name', listkey_count=5)
        main_compound_obj, molecular_formula = None, None
        if not compounds: return [], f"Molecule '{molecule_name}' not found on PubChem."
        
        # Logic to select the main compound and its molecular formula
        for i, c in enumerate(compounds):
            cid = c.cid
            actual_formula = c.molecular_formula if hasattr(c, 'molecular_formula') else None
            if actual_formula:
                is_standard_hydrocarbon = True 
                if c.canonical_smiles:
                    try:
                        mol_obj = Chem.MolFromSmiles(c.canonical_smiles)
                        if mol_obj:
                            rdDepictor.Compute2DCoords(mol_obj) # For better drawing
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
                current_compound_name_matches_input = molecule_name in [syn.lower() for syn in c.synonyms]
                if current_compound_name_matches_input: 
                    main_compound_obj = c; molecular_formula = actual_formula; break
                if not main_compound_obj: main_compound_obj = c; molecular_formula = actual_formula
        
        if not main_compound_obj or not molecular_formula:
            return [], f"Standard alkane '{molecule_name}' not found or does not meet criteria."
        
        main_name = main_compound_obj.iupac_name or molecule_name
        status_message = f"Finding isomers for {main_name} (Formula: {molecular_formula})..."
        isomers_found_raw = pcp.get_compounds(molecular_formula, 'formula', listkey_count=50)
        if not isomers_found_raw: return [], f"No isomers found for formula {molecular_formula}."
        
        valid_structural_alkanes_entries, unique_accepted_smiles = [], set()
        # Logic to filter for valid structural alkane isomers
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
             return [], f"No valid alkane isomers found for {molecular_formula} after filtering."
        
        # Prepare details for display
        for entry in sorted(valid_structural_alkanes_entries, key=lambda x: (len(x.canonical_smiles), x.cid)):
            pil_image = draw_molecule_pil(entry.canonical_smiles)
            if pil_image:
                name = entry.iupac_name
                if not name and entry.synonyms:
                    # Simplified name selection logic
                    simple_names = [s for s in entry.synonyms if s.lower().endswith("ane") and not any(char.isdigit() for char in s.split('-')[0]) and '-' not in s.split(' ')[0]]
                    if simple_names: name = min(simple_names, key=len)
                    else: name = entry.synonyms[0]
                elif not name: name = f"Alkane (CID: {entry.cid})" # Fallback if no name
                isomer_details_list.append({"cid": entry.cid, "name": name.capitalize(), "smiles": entry.canonical_smiles, "image": pil_image})
        
        if isomer_details_list: status_message = f"{len(isomer_details_list)} structural alkane isomers found for '{molecule_name}' (Formula: {molecular_formula})."
        else: status_message = f"No displayable isomers found for '{molecule_name}'."
        return isomer_details_list, status_message
    except pcp.PubChemHTTPError as e: return [], f"PubChem API Error: {e}"
    except Exception as e: return [], f"An unexpected error occurred: {e}\n{traceback.format_exc()}"

# --- Streamlit UI ---
st.set_page_config(page_title="Alkane Isomer Viewer", layout="wide")
st.title("Alkane Isomer Finder and Viewer")

# Initialize session state variables
if 'isomer_data' not in st.session_state: st.session_state.isomer_data = [] # List of isomer dicts
if 'status_message' not in st.session_state: st.session_state.status_message = ""
if 'selected_cid_for_3d' not in st.session_state: st.session_state.selected_cid_for_3d = None
if 'selected_name_for_3d' not in st.session_state: st.session_state.selected_name_for_3d = ""
if 'current_3d_index' not in st.session_state: st.session_state.current_3d_index = -1 # Index of the currently shown 3D model in isomer_data
if 'molecule_searched' not in st.session_state: st.session_state.molecule_searched = ""
if 'molecule_name_input' not in st.session_state: st.session_state.molecule_name_input = "" 
if 'run_search_after_example' not in st.session_state: st.session_state.run_search_after_example = False


# --- Sidebar for Input and Examples ---
st.sidebar.header("Search and Examples")
current_molecule_input = st.sidebar.text_input(
    label="Enter alkane name:",
    value=st.session_state.molecule_name_input,
    placeholder="e.g., butane, pentane",
    key="sidebar_molecule_input",
    help="Enter the name of the alkane (e.g., hexane) in lowercase English."
)
if current_molecule_input != st.session_state.molecule_name_input: # User typed in the box
    st.session_state.molecule_name_input = current_molecule_input
    st.session_state.run_search_after_example = False 

example_alkanes = ["butane", "pentane", "hexane", "heptane"]
st.sidebar.subheader("Or select an example:")
cols_examples = st.sidebar.columns(2) 
for i, example in enumerate(example_alkanes):
    if cols_examples[i % 2].button(example.capitalize(), key=f"example_{example}"):
        st.session_state.molecule_name_input = example 
        st.session_state.selected_cid_for_3d = None 
        st.session_state.selected_name_for_3d = ""
        st.session_state.current_3d_index = -1
        st.session_state.run_search_after_example = True 
        st.rerun() 

search_button_sidebar = st.sidebar.button("Search Isomers", type="primary", key="sidebar_search_button")

# Logic for running the search
should_run_search = False
if search_button_sidebar and st.session_state.molecule_name_input:
    should_run_search = True
elif st.session_state.run_search_after_example and st.session_state.molecule_name_input:
    should_run_search = True
    st.session_state.run_search_after_example = False # Reset flag

if should_run_search:
    with st.spinner(f"Processing for {st.session_state.molecule_name_input}..."):
        isomers, status_msg = process_alkane_request(st.session_state.molecule_name_input)
        st.session_state.isomer_data = isomers
        st.session_state.status_message = status_msg
        st.session_state.selected_cid_for_3d = None 
        st.session_state.selected_name_for_3d = ""
        st.session_state.current_3d_index = -1
        st.session_state.molecule_searched = st.session_state.molecule_name_input
elif search_button_sidebar and not st.session_state.molecule_name_input: 
    st.session_state.status_message = "Please enter a molecule name or select an example."
    st.session_state.isomer_data = []
    st.session_state.selected_cid_for_3d = None
    st.session_state.selected_name_for_3d = ""
    st.session_state.current_3d_index = -1
    st.session_state.molecule_searched = ""

# Display status message in the sidebar
if st.session_state.status_message:
    if "Error" in st.session_state.status_message or "not found" in st.session_state.status_message :
        st.sidebar.error(st.session_state.status_message)
    else:
        st.sidebar.info(st.session_state.status_message)


# --- Tabs for Displaying Results ---
if st.session_state.isomer_data or st.session_state.selected_cid_for_3d :
    gallery_tab_title = f"Isomer Gallery ({len(st.session_state.isomer_data)} found)" if st.session_state.isomer_data else "Isomer Gallery"
    viewer_3d_tab_title = f"3D View ({st.session_state.selected_name_for_3d})" if st.session_state.selected_cid_for_3d and st.session_state.selected_name_for_3d else "3D View"
    
    tab_gallery, tab_3d_viewer = st.tabs([gallery_tab_title, viewer_3d_tab_title])

    with tab_gallery:
        if st.session_state.isomer_data:
            st.subheader(f"Isomers found for: {st.session_state.molecule_searched.capitalize()}")
            num_columns_gallery = 3 # Adjust as needed
            gallery_cols = st.columns(num_columns_gallery)
            for i, isomer in enumerate(st.session_state.isomer_data):
                with gallery_cols[i % num_columns_gallery]:
                    container = st.container(border=True) 
                    container.image(isomer["image"], 
                                    caption=f"{isomer['name']}", 
                                    use_container_width=True) 
                    container.markdown(f"<small>SMILES: {isomer['smiles']}<br>CID: {isomer['cid']}</small>", unsafe_allow_html=True)
                    if container.button(f"View 3D", key=f"btn_3d_tab_{isomer['cid']}"):
                        st.session_state.selected_cid_for_3d = isomer['cid']
                        st.session_state.selected_name_for_3d = isomer['name'] 
                        st.session_state.current_3d_index = i # Store index for Next/Previous
                        st.rerun() 
        elif st.session_state.molecule_searched:
            st.info("No isomers to display in the gallery for the current search.")
        else:
            st.info("Search for an alkane to see the isomer gallery.")

    with tab_3d_viewer:
        if st.session_state.selected_cid_for_3d:
            st.subheader(f"3D Structure for: {st.session_state.selected_name_for_3d}")
            
            # Navigation buttons
            nav_cols = st.columns([1, 1, 5, 1, 1]) # Adjust ratios for spacing
            
            # Previous button
            if nav_cols[0].button("⬅️ Previous", key="prev_3d", disabled=(st.session_state.current_3d_index <= 0)):
                st.session_state.current_3d_index -= 1
                prev_isomer = st.session_state.isomer_data[st.session_state.current_3d_index]
                st.session_state.selected_cid_for_3d = prev_isomer['cid']
                st.session_state.selected_name_for_3d = prev_isomer['name']
                st.rerun()

            # Clear button in the middle
            if nav_cols[2].button("Clear 3D View", key="clear_3d_view_tab"):
                st.session_state.selected_cid_for_3d = None
                st.session_state.selected_name_for_3d = ""
                st.session_state.current_3d_index = -1
                st.rerun()

            # Next button
            if nav_cols[4].button("Next ➡️", key="next_3d", disabled=(st.session_state.current_3d_index >= len(st.session_state.isomer_data) - 1)):
                st.session_state.current_3d_index += 1
                next_isomer = st.session_state.isomer_data[st.session_state.current_3d_index]
                st.session_state.selected_cid_for_3d = next_isomer['cid']
                st.session_state.selected_name_for_3d = next_isomer['name']
                st.rerun()

            # Display 3D viewer
            with st.spinner(f"Loading 3D structure for {st.session_state.selected_name_for_3d} (CID: {st.session_state.selected_cid_for_3d})..."):
                sdf_data, error = get_sdf_content(st.session_state.selected_cid_for_3d)
                if sdf_data:
                    html_3d = generate_3d_viewer_html(sdf_data, st.session_state.selected_name_for_3d, width=600, height=450)
                    st.components.v1.html(html_3d, height=470, width=620, scrolling=False)
                # Error is handled by get_sdf_content or generate_3d_viewer_html by displaying st.warning/st.error
            
        else:
            st.info("Select an isomer from the 'Isomer Gallery' tab to view its 3D structure.")
else: # No search done yet, or search yielded no results and nothing is selected for 3D
    if st.session_state.molecule_searched and not st.session_state.isomer_data:
        # Status message already shown in sidebar if search yielded no results
        pass 
    else:
        st.info("To begin, enter an alkane name in the sidebar and click 'Search Isomers', or select an example.")

st.sidebar.markdown("---")
st.sidebar.markdown("Built with Streamlit, RDKit, PubChemPy, and py3Dmol.")
