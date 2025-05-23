import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage, MolDrawOptions
import py3Dmol
import os
import traceback
from PIL import Image
from rdkit.Chem import rdDepictor
from io import BytesIO

# --- Helper Functions ---
def draw_molecule_pil(smiles_string, size=(400, 350), legend=""):
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if mol:
            rdDepictor.Compute2DCoords(mol)
            draw_options = MolDrawOptions()
            draw_options.bondLineWidth = 2
            draw_options.padding = 0.05     
            img = MolToImage(mol, size=size, legend=legend if legend else "", options=draw_options)
            return img
        else: return None
    except Exception as e:
        print(f"Error in draw_molecule_pil for SMILES {smiles_string} with legend '{legend}': {e}")
        return None

def image_to_bytes(pil_image, format="PNG"):
    if pil_image is None: return None
    img_byte_arr = BytesIO()
    pil_image.save(img_byte_arr, format=format)
    return img_byte_arr.getvalue()

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
        elif display_style == 'sphere': viewer.setStyle({'sphere': {'scale': 0.35}}) # CPK-like
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
    if not compound_obj: return None
    properties = {
        "IUPAC Name": getattr(compound_obj, 'iupac_name', 'N/A'),
        "Molecular Formula": getattr(compound_obj, 'molecular_formula', 'N/A'),
        "Molecular Weight": f"{getattr(compound_obj, 'molecular_weight', 'N/A')} g/mol",
        "Canonical SMILES": getattr(compound_obj, 'canonical_smiles', 'N/A'),
        # ... (سایر خواص مانند قبل)
        "PubChem CID": getattr(compound_obj, 'cid', 'N/A'),
        "Charge": getattr(compound_obj, 'charge', 'N/A'),
        "Heavy Atom Count": getattr(compound_obj, 'heavy_atom_count', 'N/A'),
        "H-Bond Donor Count": getattr(compound_obj, 'hydrogen_bond_donor_count', 'N/A'),
        "H-Bond Acceptor Count": getattr(compound_obj, 'hydrogen_bond_acceptor_count', 'N/A'),
    }
    return properties

def process_alkane_request(molecule_name_input):
    # ... (کد کامل process_alkane_request مانند قبل) ...
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
            smiles = isomer_entry.canonical_smiles;_ = smiles
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
            name = entry.iupac_name
            if not name and entry.synonyms:
                simple_names = [s for s in entry.synonyms if s.lower().endswith("ane") and not any(char.isdigit() for char in s.split('-')[0]) and '-' not in s.split(' ')[0]]
                if simple_names: name = min(simple_names, key=len)
                else: name = entry.synonyms[0]
            elif not name: name = f"Alkane (CID: {entry.cid})"
            pil_image = draw_molecule_pil(entry.canonical_smiles, legend=name.capitalize()) 
            if pil_image:
                isomer_details_list.append({"cid": entry.cid, "name": name.capitalize(), "smiles": entry.canonical_smiles, "image": pil_image})
        if isomer_details_list: status_message = f"{len(isomer_details_list)} isomers found for '{molecule_name}' ({molecular_formula})."
        else: status_message = f"No displayable isomers found for '{molecule_name}'."
        return isomer_details_list, status_message, main_molecule_properties
    except pcp.PubChemHTTPError as e: return [], f"PubChem API Error: {e}", None
    except Exception as e: return [], f"An unexpected error occurred: {e}\n{traceback.format_exc()}", None

def process_general_molecule_search(search_term):
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
                img_2d = draw_molecule_pil(compound_obj.canonical_smiles, size=(450,400), legend=name_to_display)
            molecule_details = {
                "cid": compound_obj.cid, "name": name_to_display, "properties": props,
                "image_2d": img_2d, "smiles": compound_obj.canonical_smiles 
            }
            status_message = f"Information found for: {name_to_display} (CID: {compound_obj.cid})"
        else:
            status_message = f"Molecule '{term}' not found on PubChem."
        return molecule_details, status_message
    except pcp.PubChemHTTPError as e: return None, f"PubChem API Error during general search: {e}"
    except Exception as e: return None, f"An unexpected error occurred during general search: {e}\n{traceback.format_exc()}"

# --- New Helper Function for Element Legend ---
def generate_element_legend_html(smiles_or_elements_list=None):
    """
    Generates an HTML string for an element color legend.
    If smiles_or_elements_list is a SMILES string, it tries to find elements.
    If it's a list of element symbols, it uses that.
    Otherwise, shows a default legend.
    """
    # Jmol/py3Dmol default-like colors (approximate hex)
    element_colors = {
        'C': '#888888', 'H': '#FFFFFF', 'O': '#FF kırmızı', 'N': '#8F8FFF', # Adjusted N to be lighter blue
        'S': '#FFC832', 'P': '#FF8000', 'F': '#DAFF7D', 'Cl': '#50D250',
        'Br': '#C82828', 'I': '#A020F0', 'He': '#D9FFFF', 'Li': '#CC80FF',
        'B': '#FFB5B5', 'Na': '#AB5CF2', 'Mg': '#8AFF00', 'Al': '#BFA6A6',
        'Si': '#F0C8A0', 'K': '#8F40D4', 'Ca': '#3CFF00', 'Fe': '#E06633',
        # Add more as needed, or use a more comprehensive mapping
        'OTHER': '#BEBEBE' # Default for unlisted elements
    }
    
    elements_in_molecule = set()
    if isinstance(smiles_or_elements_list, str): # Assume SMILES
        try:
            mol = Chem.MolFromSmiles(smiles_or_elements_list)
            if mol:
                for atom in mol.GetAtoms():
                    elements_in_molecule.add(atom.GetSymbol())
        except: # If SMILES parsing fails, fallback to default
            elements_in_molecule.clear() 
            
    elif isinstance(smiles_or_elements_list, list):
        elements_in_molecule = set(smiles_or_elements_list)

    # If no specific elements found/provided, show a common default set
    if not elements_in_molecule:
        elements_to_show = ['C', 'H', 'O', 'N', 'S', 'P', 'Cl', 'F']
    else:
        # Sort for consistent order, maybe prioritize common ones
        common_first = [el for el in ['C', 'H', 'O', 'N'] if el in elements_in_molecule]
        other_elements = sorted(list(elements_in_molecule - set(common_first)))
        elements_to_show = common_first + other_elements

    legend_html = "<div style='padding: 5px; border: 1px solid #444; border-radius: 5px; font-size: 0.75em; background-color: #222;'>"
    legend_html += "<b style='color: #eee;'>Element Legend:</b><br>"
    
    for symbol in elements_to_show:
        color = element_colors.get(symbol, element_colors['OTHER'])
        legend_html += f"<span style='display:inline-block; width:10px; height:10px; background-color:{color}; border: 1px solid #555; margin-right:5px; vertical-align: middle;'></span><span style='color: #ddd; vertical-align: middle;'>{symbol}</span><br>"
            
    legend_html += "</div>"
    return legend_html


# --- Streamlit UI ---
st.set_page_config(page_title="Chemical Compound Explorer", layout="centered", initial_sidebar_state="expanded")
st.title("Chemical Compound Explorer")

# ... (Session state initialization like before) ...
if 'alkane_isomer_data' not in st.session_state: st.session_state.alkane_isomer_data = []
# ... (all other session state vars)
if 'alkane_main_molecule_props' not in st.session_state: st.session_state.alkane_main_molecule_props = None
if 'selected_isomer_cid_for_3d' not in st.session_state: st.session_state.selected_isomer_cid_for_3d = None
if 'selected_isomer_name_for_3d' not in st.session_state: st.session_state.selected_isomer_name_for_3d = ""
if 'current_isomer_3d_index' not in st.session_state: st.session_state.current_isomer_3d_index = -1
if 'general_molecule_data' not in st.session_state: st.session_state.general_molecule_data = None
if 'general_molecule_name_input' not in st.session_state: st.session_state.general_molecule_name_input = ""
if 'status_message' not in st.session_state: st.session_state.status_message = ""
if 'selected_3d_style' not in st.session_state: st.session_state.selected_3d_style = 'stick' 
if 'last_search_type' not in st.session_state: st.session_state.last_search_type = None 
if 'alkane_name_input' not in st.session_state: st.session_state.alkane_name_input = "" 
if 'run_alkane_search_after_example' not in st.session_state: st.session_state.run_alkane_search_after_example = False
if 'alkane_molecule_searched' not in st.session_state: st.session_state.alkane_molecule_searched = ""


# --- Sidebar (like before) ---
# ... (کد کامل سایدبار از پاسخ قبلی) ...
st.sidebar.header("Search Options")
st.sidebar.subheader("1. Alkane Isomer Search")
current_alkane_input = st.sidebar.text_input(
    label="Enter alkane name (for isomers):", value=st.session_state.alkane_name_input,
    placeholder="e.g., butane, pentane", key="sidebar_alkane_input",
    help="Finds structural isomers for standard alkanes."
)
if current_alkane_input != st.session_state.alkane_name_input:
    st.session_state.alkane_name_input = current_alkane_input
    st.session_state.run_alkane_search_after_example = False 
example_alkanes = ["butane", "pentane", "hexane", "heptane", "octane"]
st.sidebar.caption("Alkane Examples:")
num_example_cols = 3 
cols_examples_alkane = st.sidebar.columns(num_example_cols) 
for i, example in enumerate(example_alkanes):
    if cols_examples_alkane[i % num_example_cols].button(example.capitalize(), key=f"example_alkane_{example}", use_container_width=True):
        st.session_state.alkane_name_input = example 
        st.session_state.run_alkane_search_after_example = True 
        st.rerun() 
if st.sidebar.button("Search Alkane Isomers", type="primary", key="search_alkane_button", use_container_width=True):
    if st.session_state.alkane_name_input:
        st.session_state.last_search_type = "alkane"
        with st.spinner(f"Searching isomers for {st.session_state.alkane_name_input}..."):
            isomers, status_msg, main_props = process_alkane_request(st.session_state.alkane_name_input)
            st.session_state.alkane_isomer_data, st.session_state.status_message, st.session_state.alkane_main_molecule_props = isomers, status_msg, main_props
            st.session_state.selected_isomer_cid_for_3d, st.session_state.selected_isomer_name_for_3d, st.session_state.current_isomer_3d_index = None, "", -1
            st.session_state.alkane_molecule_searched = st.session_state.alkane_name_input
            st.session_state.general_molecule_data = None 
    else: st.session_state.status_message = "Please enter an alkane name for isomer search."
if st.session_state.run_alkane_search_after_example and st.session_state.alkane_name_input:
    st.session_state.last_search_type = "alkane"
    st.session_state.run_alkane_search_after_example = False 
    with st.spinner(f"Searching isomers for {st.session_state.alkane_name_input}..."):
        isomers, status_msg, main_props = process_alkane_request(st.session_state.alkane_name_input)
        st.session_state.alkane_isomer_data, st.session_state.status_message, st.session_state.alkane_main_molecule_props = isomers, status_msg, main_props
        st.session_state.selected_isomer_cid_for_3d, st.session_state.selected_isomer_name_for_3d, st.session_state.current_isomer_3d_index = None, "", -1
        st.session_state.alkane_molecule_searched = st.session_state.alkane_name_input
        st.session_state.general_molecule_data = None 
    st.rerun() 
st.sidebar.markdown("---")
st.sidebar.subheader("2. General Molecule Information")
general_molecule_search_term = st.sidebar.text_input(
    label="Enter any chemical name/identifier:", value=st.session_state.general_molecule_name_input,
    placeholder="e.g., aspirin, caffeine, C6H6", key="sidebar_general_input",
    help="Get information, 2D, and 3D structure for any chemical."
)
if general_molecule_search_term != st.session_state.general_molecule_name_input:
    st.session_state.general_molecule_name_input = general_molecule_search_term
if st.sidebar.button("Search Molecule Info", type="primary", key="search_general_button", use_container_width=True):
    if st.session_state.general_molecule_name_input:
        st.session_state.last_search_type = "general"
        with st.spinner(f"Searching for {st.session_state.general_molecule_name_input}..."):
            mol_data, status_msg = process_general_molecule_search(st.session_state.general_molecule_name_input)
            st.session_state.general_molecule_data, st.session_state.status_message = mol_data, status_msg
            st.session_state.alkane_isomer_data, st.session_state.alkane_main_molecule_props = [], None
            st.session_state.selected_isomer_cid_for_3d, st.session_state.selected_isomer_name_for_3d, st.session_state.current_isomer_3d_index = None, "", -1
    else: st.session_state.status_message = "Please enter a chemical name for general search."
if st.session_state.status_message:
    is_error = any(keyword in st.session_state.status_message.lower() for keyword in ["error", "not found", "empty", "invalid"])
    if is_error: st.sidebar.error(st.session_state.status_message)
    else: st.sidebar.info(st.session_state.status_message)

# --- Main Page Tabs (منطق تب‌ها مانند قبل) ---
# ... (کپی از پاسخ قبلی) ...
gallery_title, props_title, view2d_title, view3d_title = "Isomer Gallery", "Properties", "2D Structure", "3D View"
if st.session_state.last_search_type == "alkane":
    if st.session_state.alkane_isomer_data: gallery_title = f"Alkane Isomers ({len(st.session_state.alkane_isomer_data)})"
    if st.session_state.alkane_main_molecule_props: props_title = f"Alkane Properties: {st.session_state.alkane_main_molecule_props.get('IUPAC Name', '')[:20]}"
    if st.session_state.selected_isomer_cid_for_3d: view3d_title = f"Isomer 3D: {st.session_state.selected_isomer_name_for_3d[:20]}"
    tabs_to_show = [gallery_title, props_title, view3d_title]
    tab_gallery, tab_main_props, tab_3d_isomer = st.tabs(tabs_to_show)
elif st.session_state.last_search_type == "general" and st.session_state.general_molecule_data:
    g_data = st.session_state.general_molecule_data
    props_title = f"Properties: {g_data['name'][:20]}"
    if g_data.get("image_2d"): view2d_title = f"2D: {g_data['name'][:20]}"
    view3d_title = f"3D: {g_data['name'][:20]}"
    tabs_to_show = [props_title, view2d_title, view3d_title] if g_data.get("image_2d") else [props_title, view3d_title]
    if g_data.get("image_2d"): tab_g_props, tab_g_2d, tab_g_3d = st.tabs(tabs_to_show)
    else: tab_g_props, tab_g_3d = st.tabs(tabs_to_show); tab_g_2d = None 
else:
    st.info("To begin, use the search options in the sidebar.")
    tab_gallery, tab_main_props, tab_3d_isomer = (None, None, None)
    tab_g_props, tab_g_2d, tab_g_3d = (None, None, None)


# Content for Alkane Isomer Search
if st.session_state.last_search_type == "alkane":
    if tab_gallery and st.session_state.alkane_isomer_data:
        # ... (محتوای تب گالری ایزومر مانند قبل با دکمه دانلود) ...
        with tab_gallery:
            st.subheader(f"Isomers found for: {st.session_state.alkane_molecule_searched.capitalize()}")
            num_columns_gallery = 3 
            gallery_cols = st.columns(num_columns_gallery)
            for i, isomer in enumerate(st.session_state.alkane_isomer_data):
                with gallery_cols[i % num_columns_gallery]:
                    container = st.container(border=True) 
                    pil_image_isomer = isomer["image"]
                    container.image(pil_image_isomer, caption=f"CID: {isomer['cid']}", use_container_width=True) 
                    if pil_image_isomer:
                        img_bytes_isomer = image_to_bytes(pil_image_isomer)
                        container.download_button(label="Download 2D", data=img_bytes_isomer, file_name=f"{isomer['name'].replace(' ', '_')}_CID_{isomer['cid']}_2D.png", mime="image/png", key=f"download_iso_{isomer['cid']}")
                    container.markdown(f"<small>SMILES: {isomer['smiles']}</small>", unsafe_allow_html=True)
                    if container.button(f"View 3D", key=f"btn_3d_isomer_{isomer['cid']}"):
                        st.session_state.selected_isomer_cid_for_3d, st.session_state.selected_isomer_name_for_3d, st.session_state.current_isomer_3d_index = isomer['cid'], isomer['name'], i; st.rerun()
    if tab_main_props and st.session_state.alkane_main_molecule_props:
        # ... (محتوای تب خواص مولکول اصلی مانند قبل) ...
        with tab_main_props:
            props = st.session_state.alkane_main_molecule_props
            main_mol_name = props.get("IUPAC Name", st.session_state.alkane_molecule_searched.capitalize())
            st.subheader(f"Properties for Searched Alkane: {main_mol_name}")
            prop_cols = st.columns(2)
            prop_list = list(props.items())
            for i_prop, (key, value) in enumerate(prop_list):
                if value != 'N/A' and value is not None: prop_cols[i_prop % 2].markdown(f"**{key}:** {value}")
            st.markdown("---")
    if tab_3d_isomer and st.session_state.selected_isomer_cid_for_3d:
        with tab_3d_isomer:
            # ... (محتوای تب سه‌بعدی ایزومر با انتخابگر سبک و دکمه‌های ناوبری مانند قبل) ...
            st.subheader(f"3D Structure for Isomer: {st.session_state.selected_isomer_name_for_3d}")
            style_options_map = {'Stick': 'stick', 'Line': 'line', 'Ball and Stick': 'ball_and_stick'}
            style_labels = list(style_options_map.keys())
            try:
                current_style_label = [k for k, v in style_options_map.items() if v == st.session_state.selected_3d_style][0]
                current_style_index = style_labels.index(current_style_label)
            except IndexError: 
                current_style_index = 0 
                if style_labels: st.session_state.selected_3d_style = style_options_map[style_labels[0]]
                else: st.session_state.selected_3d_style = 'stick'
            selected_style_label = st.radio("Select Display Style:", options=style_labels, key="radio_3d_style_isomer", horizontal=True, index=current_style_index)
            if style_options_map.get(selected_style_label) != st.session_state.selected_3d_style:
                st.session_state.selected_3d_style = style_options_map[selected_style_label]; st.rerun() 
            nav_col_layout = [1, 0.1, 1.5, 0.1, 1]; prev_col, _, clear_col, _, next_col = st.columns(nav_col_layout)
            with prev_col:
                if st.button("⬅️ Previous", key="prev_3d_isomer", help="View Previous Isomer", use_container_width=True, disabled=(st.session_state.current_isomer_3d_index <= 0)):
                    st.session_state.current_isomer_3d_index -= 1; prev_isomer = st.session_state.alkane_isomer_data[st.session_state.current_isomer_3d_index]
                    st.session_state.selected_isomer_cid_for_3d, st.session_state.selected_isomer_name_for_3d = prev_isomer['cid'], prev_isomer['name']; st.rerun()
            with clear_col:
                if st.button("Clear Isomer 3D View", key="clear_3d_isomer", use_container_width=True):
                    st.session_state.selected_isomer_cid_for_3d, st.session_state.selected_isomer_name_for_3d, st.session_state.current_isomer_3d_index = None, "", -1; st.rerun()
            with next_col:
                if st.button("Next ➡️", key="next_3d_isomer", help="View Next Isomer", use_container_width=True, disabled=(st.session_state.current_isomer_3d_index >= len(st.session_state.alkane_isomer_data) - 1)):
                    st.session_state.current_isomer_3d_index += 1; next_isomer = st.session_state.alkane_isomer_data[st.session_state.current_isomer_3d_index]
                    st.session_state.selected_isomer_cid_for_3d, st.session_state.selected_isomer_name_for_3d = next_isomer['cid'], next_isomer['name']; st.rerun()
            st.markdown("---") 
            with st.spinner(f"Loading 3D structure for {st.session_state.selected_isomer_name_for_3d}..."):
                sdf_data, error = get_sdf_content(st.session_state.selected_isomer_cid_for_3d)
                if sdf_data:
                    html_3d = generate_3d_viewer_html(sdf_data, st.session_state.selected_isomer_name_for_3d, display_style=st.session_state.selected_3d_style, width=600, height=450)
                    st.components.v1.html(html_3d, height=470, width=620, scrolling=False)


# Content for General Molecule Search
if st.session_state.last_search_type == "general" and st.session_state.general_molecule_data:
    g_data = st.session_state.general_molecule_data
    if tab_g_props:
        with tab_g_props:
            st.subheader(f"Properties for: {g_data['name']}")
            if g_data["properties"]:
                prop_cols = st.columns(2)
                prop_list = list(g_data["properties"].items())
                for i_prop_g, (key_g, value_g) in enumerate(prop_list):
                    if value_g != 'N/A' and value_g is not None: prop_cols[i_prop_g % 2].markdown(f"**{key_g}:** {value_g}")
            else: st.info("No detailed properties found.")
            st.markdown("---")

    if tab_g_2d and g_data.get("image_2d"):
        with tab_g_2d:
            st.subheader(f"2D Structure for: {g_data['name']}")
            pil_image_general = g_data["image_2d"]
            st.image(pil_image_general, use_container_width=True) 
            if pil_image_general:
                img_bytes_general = image_to_bytes(pil_image_general)
                st.download_button(label="Download 2D Structure", data=img_bytes_general, file_name=f"{g_data['name'].replace(' ', '_')}_CID_{g_data['cid']}_2D.png", mime="image/png", key=f"download_general_{g_data['cid']}")

    if tab_g_3d:
        with tab_g_3d:
            st.subheader(f"3D Structure for: {g_data['name']}")
            # Layout for 3D viewer and legend
            viewer_col, legend_col = st.columns([3, 1]) # Viewer gets 3/4, legend gets 1/4 width

            with viewer_col:
                style_options_map_g = {'Stick': 'stick', 'Line': 'line', 'Ball and Stick': 'ball_and_stick'}
                style_labels_g = list(style_options_map_g.keys())
                try:
                    current_style_label_g = [k for k, v in style_options_map_g.items() if v == st.session_state.selected_3d_style][0]
                    current_style_index_g = style_labels_g.index(current_style_label_g)
                except IndexError: 
                    current_style_index_g = 0 
                    if style_labels_g: st.session_state.selected_3d_style = style_options_map_g[style_labels_g[0]]
                    else: st.session_state.selected_3d_style = 'stick'
                
                selected_style_label_g = st.radio(
                    "Select Display Style:", 
                    options=style_labels_g, 
                    key="radio_3d_style_general", 
                    horizontal=True, 
                    index=current_style_index_g
                )
                if style_options_map_g.get(selected_style_label_g) != st.session_state.selected_3d_style:
                    st.session_state.selected_3d_style = style_options_map_g[selected_style_label_g]
                    st.rerun() 

                with st.spinner(f"Loading 3D structure for {g_data['name']}..."):
                    sdf_data, error = get_sdf_content(g_data['cid'])
                    if sdf_data:
                        html_3d = generate_3d_viewer_html(
                            sdf_data, 
                            g_data['name'], 
                            display_style=st.session_state.selected_3d_style, 
                            width=450, # Adjust width for column
                            height=400
                        )
                        st.components.v1.html(html_3d, height=420, width=470, scrolling=False)
                    else: 
                        st.info(f"3D structure could not be loaded. {error if error else ''}")
            
            with legend_col:
                st.markdown("##### ") # Add some space above legend
                if g_data.get('smiles'):
                    legend_html_content = generate_element_legend_html(g_data['smiles'])
                    st.markdown(legend_html_content, unsafe_allow_html=True)
                else: # Fallback to default legend if no SMILES
                    legend_html_content = generate_element_legend_html() 
                    st.markdown(legend_html_content, unsafe_allow_html=True)


st.sidebar.markdown("---")
st.sidebar.markdown("Built with Streamlit, RDKit, PubChemPy, and py3Dmol.")
