import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage, MolDrawOptions, rdMolDraw2D # Import rdMolDraw2D
import py3Dmol
import os
import traceback
from PIL import Image
from rdkit.Chem import rdDepictor, Descriptors, Lipinski
from io import BytesIO

# --- Helper Functions ---
def draw_molecule_pil(smiles_string, size=(400, 350), legend="", highlight_atoms=None):
    """
    Draws a molecule with enhanced styling (CPK colors, stereo annotations) 
    and returns a PIL Image object.
    """
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if not mol:
            return None

        # Generate 2D coordinates
        rdDepictor.Compute2DCoords(mol)
        Chem.WedgeMolBonds(mol, mol.GetConformer()) # For better stereo bond display

        # Use MolDraw2DCairo for higher quality PNG output
        drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
        opts = drawer.drawOptions()

        # --- Enhanced Drawing Options ---
        # 1. CPK Coloring for atoms
        # opts.useBWAtomPalette() # Uncomment for black & white
        # For CPK colors, RDKit usually applies them by default if not B&W.
        # We can ensure it by clearing any specific palette or setting a default one.
        # opts. ομάδαColourPalette.clear() # Or use a predefined CPK-like palette if needed
        # Default RDKit colors are generally good for CPK.

        # 2. Bond line width
        opts.bondLineWidth = 2
        opts.padding = 0.1 # Increased padding a bit for legend and clarity

        # 3. Atom Label Font Size (optional)
        # opts.atomLabelFontSize = 18 

        # 4. Add Stereo Annotation (R/S, E/Z)
        opts.addStereoAnnotation = True
        opts.includeChiralFlagLabel = True # Shows R/S label on chiral centers if specified

        # 5. Legend (Molecule Name)
        if legend:
            opts.legendFraction = 0.10 # Space for legend
            opts.legendFontSize = 18   # Font size for legend

        # Highlight atoms if specified (example of advanced feature)
        # highlight_atoms should be a list of atom indices
        # highlight_bonds can also be used
        # Default to an empty list if None
        atom_highlights = highlight_atoms if highlight_atoms is not None else []


        # Draw the molecule
        # For legend, MolToImage handles it well. If drawing directly with drawer, you'd add it manually.
        # drawer.DrawMolecule(mol, legend=legend if legend else "", highlightAtoms=atom_highlights, options=opts)
        # drawer.FinishDrawing()
        # png_data = drawer.GetDrawingText() # This is PNG data as bytes
        # return Image.open(BytesIO(png_data))

        # Using MolToImage is simpler if it supports the options well enough
        # MolToImage internally can use MolDraw2DCairo or MolDraw2DSVG
        # We pass the options object to MolToImage
        img = MolToImage(mol, size=size, legend=legend if legend else "", options=opts, highlightAtoms=atom_highlights)
        return img

    except Exception as e:
        print(f"Error in draw_molecule_pil for SMILES {smiles_string} with legend '{legend}': {e}")
        print(traceback.format_exc())
        return None

# ... (بقیه توابع image_to_bytes, get_sdf_content, generate_3d_viewer_html, get_rdkit_properties, get_compound_properties, process_alkane_request, process_general_molecule_search مانند قبل خواهند بود، فقط فراخوانی draw_molecule_pil در آنها از legend استفاده می‌کند)

# --- توابع دیگر (بدون تغییر از نسخه قبلی که کش اضافه شد) ---
@st.cache_data(ttl=3600) 
def cached_process_alkane_request(molecule_name_input):
    return process_alkane_request(molecule_name_input)
@st.cache_data(ttl=3600)
def cached_process_general_molecule_search(search_term):
    return process_general_molecule_search(search_term)
@st.cache_data(ttl=3600) 
def cached_get_sdf_content(cid):
    return get_sdf_content(cid)

# ... (کد کامل get_sdf_content, generate_3d_viewer_html, get_compound_properties, process_alkane_request, process_general_molecule_search از پاسخ قبلی کپی شود. فقط مطمئن شوید که در process_alkane_request و process_general_molecule_search، هنگام فراخوانی draw_molecule_pil، مقدار legend (نام مولکول) پاس داده می‌شود.)
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

def get_rdkit_properties(mol_obj):
    if not mol_obj: return {}
    return {
        "Number of Rings": Lipinski.RingCount(mol_obj),
        "Number of Rotatable Bonds (RDKit)": Lipinski.NumRotatableBonds(mol_obj),
        "Fraction Csp3": Lipinski.FractionCSP3(mol_obj),
        "LogP (MolLogP)": Descriptors.MolLogP(mol_obj),
    }

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
        "TPSA": f"{getattr(compound_obj, 'tpsa', 'N/A')} Å²" if hasattr(compound_obj, 'tpsa') and compound_obj.tpsa is not None else 'N/A',
        "XLogP (PubChem)": getattr(compound_obj, 'xlogp', 'N/A'),
        "Heavy Atom Count": getattr(compound_obj, 'heavy_atom_count', 'N/A'),
        "H-Bond Donor Count": getattr(compound_obj, 'hydrogen_bond_donor_count', 'N/A'),
        "H-Bond Acceptor Count": getattr(compound_obj, 'hydrogen_bond_acceptor_count', 'N/A'),
        "Rotatable Bond Count (PubChem)": getattr(compound_obj, 'rotatable_bond_count', 'N/A'),
    }
    if compound_obj.canonical_smiles:
        mol = Chem.MolFromSmiles(compound_obj.canonical_smiles)
        if mol:
            rdkit_props = get_rdkit_properties(mol)
            for key, value in rdkit_props.items():
                if isinstance(value, float): properties[f"{key} (RDKit)"] = f"{value:.2f}"
                else: properties[f"{key} (RDKit)"] = value
    return properties

def process_alkane_request(molecule_name_input):
    if not molecule_name_input or not molecule_name_input.strip(): return [], "Please enter a molecule name.", None 
    molecule_name = molecule_name_input.strip().lower()
    status_message, isomer_details_list, main_molecule_properties = f"Searching for '{molecule_name}'...", [], None
    try:
        compounds = pcp.get_compounds(molecule_name, 'name', listkey_count=5) 
        main_compound_obj, molecular_formula = None, None
        if not compounds: return [], f"Molecule '{molecule_name}' not found on PubChem.", None
        for c in compounds:
            cid, actual_formula = c.cid, (c.molecular_formula if hasattr(c, 'molecular_formula') else None)
            if actual_formula:
                is_standard_hydrocarbon = True 
                if c.canonical_smiles:
                    try:
                        mol_obj = Chem.MolFromSmiles(c.canonical_smiles)
                        if mol_obj:
                            rdDepictor.Compute2DCoords(mol_obj)
                            if len(Chem.GetMolFrags(mol_obj)) > 1: is_standard_hydrocarbon = False
                            if is_standard_hydrocarbon: # Simplified checks, assuming alkane focus
                                atom_symbols_main = set(atom.GetSymbol() for atom in mol_obj.GetAtoms())
                                if not atom_symbols_main.issubset({'C', 'H'}): is_standard_hydrocarbon = False
                                for atom in mol_obj.GetAtoms():
                                    if atom.GetIsotope() != 0: is_standard_hydrocarbon = False; break
                            if is_standard_hydrocarbon:
                                for bond in mol_obj.GetBonds():
                                    if bond.GetBondType() != Chem.BondType.SINGLE: is_standard_hydrocarbon = False; break
                                if Chem.GetSymmSSSR(mol_obj): is_standard_hydrocarbon = False # Alkanes are acyclic
                        else: is_standard_hydrocarbon = False
                    except Exception: is_standard_hydrocarbon = False 
                else: is_standard_hydrocarbon = False 
                if not is_standard_hydrocarbon: continue
                if molecule_name in [syn.lower() for syn in c.synonyms]: 
                    main_compound_obj, molecular_formula = c, actual_formula; break
                if not main_compound_obj: main_compound_obj, molecular_formula = c, actual_formula
        if not main_compound_obj or not molecular_formula: return [], f"Standard alkane '{molecule_name}' not found or does not meet criteria.", None
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
                mol_iso = Chem.MolFromSmiles(smiles);
                if not mol_iso: continue
                is_valid_candidate = True # Simplified filtering for structural alkane isomers
                if len(Chem.GetMolFrags(mol_iso)) > 1: is_valid_candidate = False
                if is_valid_candidate:
                    atom_symbols = set(atom.GetSymbol() for atom in mol_iso.GetAtoms())
                    if not atom_symbols.issubset({'C', 'H'}): is_valid_candidate = False
                if is_valid_candidate: # No isotopes, single bonds, acyclic
                    for atom in mol_iso.GetAtoms():
                        if atom.GetIsotope() != 0: is_valid_candidate = False; break
                    if not is_valid_candidate: continue
                    for bond in mol_iso.GetBonds():
                        if bond.GetBondType() != Chem.BondType.SINGLE: is_valid_candidate = False; break
                    if not is_valid_candidate: continue
                    if Chem.GetSymmSSSR(mol_iso): is_valid_candidate = False
                if is_valid_candidate:
                    canonical_smiles_for_uniqueness = Chem.MolToSmiles(mol_iso, isomericSmiles=False)
                    if canonical_smiles_for_uniqueness not in unique_accepted_smiles:
                        valid_structural_alkanes_entries.append(isomer_entry); unique_accepted_smiles.add(canonical_smiles_for_uniqueness)
            except Exception: continue
        if not valid_structural_alkanes_entries: return [], f"No valid alkane isomers for {molecular_formula} after filtering.", main_molecule_properties
        for entry in sorted(valid_structural_alkanes_entries, key=lambda x: (len(x.canonical_smiles), x.cid)):
            name = entry.iupac_name
            if not name and entry.synonyms:
                simple_names = [s for s in entry.synonyms if s.lower().endswith("ane") and not any(char.isdigit() for char in s.split('-')[0]) and '-' not in s.split(' ')[0]]
                if simple_names: name = min(simple_names, key=len)
                else: name = entry.synonyms[0]
            elif not name: name = f"Alkane (CID: {entry.cid})"
            pil_image = draw_molecule_pil(entry.canonical_smiles, legend=name.capitalize()) # Pass name as legend
            if pil_image: isomer_details_list.append({"cid": entry.cid, "name": name.capitalize(), "smiles": entry.canonical_smiles, "image": pil_image})
        status_message = f"{len(isomer_details_list)} isomers found." if isomer_details_list else f"No displayable isomers for '{molecule_name}'."
        return isomer_details_list, status_message, main_molecule_properties
    except pcp.PubChemHTTPError as e: return [], f"PubChem API Error: {e}", None
    except Exception as e: return [], f"An unexpected error: {e}\n{traceback.format_exc()}", None

def process_general_molecule_search(search_term):
    if not search_term or not search_term.strip(): return None, "Please enter a chemical name or identifier."
    term, status_message, molecule_details = search_term.strip(), f"Searching for '{term}'...", None
    try:
        compounds = pcp.get_compounds(term, 'name', listkey_count=1) 
        if compounds:
            compound_obj = compounds[0]
            props = get_compound_properties(compound_obj)
            name_to_display = props.get("IUPAC Name", term.capitalize()) if props else term.capitalize()
            img_2d = None
            if compound_obj.canonical_smiles:
                img_2d = draw_molecule_pil(compound_obj.canonical_smiles, size=(450,400), legend=name_to_display)
            molecule_details = {"cid": compound_obj.cid, "name": name_to_display, "properties": props, "image_2d": img_2d, "smiles": compound_obj.canonical_smiles}
            status_message = f"Info for: {name_to_display} (CID: {compound_obj.cid})"
        else: status_message = f"Molecule '{term}' not found."
        return molecule_details, status_message
    except pcp.PubChemHTTPError as e: return None, f"PubChem API Error: {e}"
    except Exception as e: return None, f"An unexpected error: {e}\n{traceback.format_exc()}"

# --- Streamlit UI (بخش UI با تغییرات جزئی برای نمایش legend در کپشن تصویر) ---
# ... (کد کامل UI از پاسخ قبلی، فقط در بخش گالری، caption تصویر را به صورت زیر تغییر دهید
#    چون legend حالا بخشی از خود تصویر است، می‌توانید CID را به تنهایی در کپشن بگذارید
#    یا اصلاً کپشن نگذارید اگر لجند کافی است.) ...
# ... (و در بخش تصویر دوبعدی مولکول عمومی نیز مشابه)

# --- Streamlit UI (کامل) ---
st.set_page_config(page_title="Chemical Compound Explorer", layout="centered", initial_sidebar_state="expanded")
st.title("Chemical Compound Explorer")

# Initialize session state (کامل)
if 'alkane_isomer_data' not in st.session_state: st.session_state.alkane_isomer_data = []
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


# Sidebar (کامل)
st.sidebar.header("Search Options")
st.sidebar.info("Enter a chemical name or select an example to get started.")
st.sidebar.divider()
st.sidebar.subheader("1. Alkane Isomer Search")
current_alkane_input = st.sidebar.text_input(label="Enter alkane name (for isomers):", value=st.session_state.alkane_name_input, placeholder="e.g., butane, pentane", key="sidebar_alkane_input_v3", help="Finds structural isomers for standard alkanes.")
if current_alkane_input != st.session_state.alkane_name_input:
    st.session_state.alkane_name_input = current_alkane_input
    st.session_state.run_alkane_search_after_example = False 
example_alkanes = ["butane", "pentane", "hexane", "heptane", "octane"]
st.sidebar.caption("Alkane Examples:")
num_example_cols = 3 
cols_examples_alkane = st.sidebar.columns(num_example_cols) 
for i, example in enumerate(example_alkanes):
    if cols_examples_alkane[i % num_example_cols].button(example.capitalize(), key=f"example_alkane_{example}_v3", use_container_width=True):
        st.session_state.alkane_name_input = example 
        st.session_state.run_alkane_search_after_example = True 
        st.rerun() 
if st.sidebar.button("Search Alkane Isomers", type="primary", key="search_alkane_button_v3", use_container_width=True):
    if st.session_state.alkane_name_input:
        st.session_state.last_search_type = "alkane"
        isomers, status_msg, main_props = cached_process_alkane_request(st.session_state.alkane_name_input)
        st.session_state.alkane_isomer_data, st.session_state.status_message, st.session_state.alkane_main_molecule_props = isomers, status_msg, main_props
        st.session_state.selected_isomer_cid_for_3d, st.session_state.selected_isomer_name_for_3d, st.session_state.current_isomer_3d_index = None, "", -1
        st.session_state.alkane_molecule_searched = st.session_state.alkane_name_input
        st.session_state.general_molecule_data = None 
    else: st.session_state.status_message = "Please enter an alkane name for isomer search."
if st.session_state.run_alkane_search_after_example and st.session_state.alkane_name_input:
    st.session_state.last_search_type = "alkane"
    st.session_state.run_alkane_search_after_example = False 
    isomers, status_msg, main_props = cached_process_alkane_request(st.session_state.alkane_name_input)
    st.session_state.alkane_isomer_data, st.session_state.status_message, st.session_state.alkane_main_molecule_props = isomers, status_msg, main_props
    st.session_state.selected_isomer_cid_for_3d, st.session_state.selected_isomer_name_for_3d, st.session_state.current_isomer_3d_index = None, "", -1
    st.session_state.alkane_molecule_searched = st.session_state.alkane_name_input
    st.session_state.general_molecule_data = None 
    st.rerun() 
st.sidebar.divider()
st.sidebar.subheader("2. General Molecule Information")
general_molecule_search_term = st.sidebar.text_input(label="Enter any chemical name/identifier:", value=st.session_state.general_molecule_name_input, placeholder="e.g., aspirin, caffeine, C6H6", key="sidebar_general_input_v3", help="Get information, 2D, and 3D structure for any chemical.")
if general_molecule_search_term != st.session_state.general_molecule_name_input: st.session_state.general_molecule_name_input = general_molecule_search_term
if st.sidebar.button("Search Molecule Info", type="primary", key="search_general_button_v3", use_container_width=True):
    if st.session_state.general_molecule_name_input:
        st.session_state.last_search_type = "general"
        mol_data, status_msg = cached_process_general_molecule_search(st.session_state.general_molecule_name_input)
        st.session_state.general_molecule_data, st.session_state.status_message = mol_data, status_msg
        st.session_state.alkane_isomer_data, st.session_state.alkane_main_molecule_props = [], None
        st.session_state.selected_isomer_cid_for_3d, st.session_state.selected_isomer_name_for_3d, st.session_state.current_isomer_3d_index = None, "", -1
    else: st.session_state.status_message = "Please enter a chemical name for general search."
st.sidebar.divider()
if st.sidebar.button("Clear Current Search Results", key="clear_all_search_v3", use_container_width=True):
    st.session_state.alkane_isomer_data, st.session_state.alkane_main_molecule_props = [], None
    st.session_state.selected_isomer_cid_for_3d, st.session_state.selected_isomer_name_for_3d, st.session_state.current_isomer_3d_index = None, "", -1
    st.session_state.general_molecule_data = None
    st.session_state.status_message, st.session_state.last_search_type, st.session_state.alkane_molecule_searched = "Search results cleared.", None, ""
    st.rerun()
if st.session_state.status_message:
    is_error = any(keyword in st.session_state.status_message.lower() for keyword in ["error", "not found", "empty", "invalid", "unexpected"])
    if is_error: st.sidebar.error(st.session_state.status_message)
    else: st.sidebar.info(st.session_state.status_message)

# Main Page Tabs (کامل)
gallery_title, props_title, view2d_title, view3d_title = "Isomer Gallery", "Properties", "2D Structure", "3D View"
active_tabs, tabs_list = {}, [] # Define before conditional assignment
if st.session_state.last_search_type == "alkane":
    if st.session_state.alkane_isomer_data or st.session_state.alkane_main_molecule_props: tabs_list.append(gallery_title)
    if st.session_state.alkane_main_molecule_props and gallery_title not in tabs_list: tabs_list.append(props_title) # Add if not part of gallery implicitly
    if st.session_state.selected_isomer_cid_for_3d or st.session_state.alkane_isomer_data: tabs_list.append(view3d_title)
    if tabs_list: created_tabs = st.tabs(tabs_list); active_tabs = {title: tab for title, tab in zip(tabs_list, created_tabs)}
elif st.session_state.last_search_type == "general" and st.session_state.general_molecule_data:
    g_data = st.session_state.general_molecule_data
    if g_data.get("properties"): tabs_list.append(props_title)
    if g_data.get("image_2d"): tabs_list.append(view2d_title)
    if g_data.get("cid"): tabs_list.append(view3d_title)
    if tabs_list: created_tabs = st.tabs(tabs_list); active_tabs = {title: tab for title, tab in zip(tabs_list, created_tabs)}
else: 
    if not (st.session_state.last_search_type == "alkane" or st.session_state.last_search_type == "general"):
        st.info("To begin, use the search options in the sidebar.")

# Content for Alkane Isomer Search
if st.session_state.last_search_type == "alkane":
    if gallery_title in active_tabs and (st.session_state.alkane_isomer_data or st.session_state.alkane_main_molecule_props):
        with active_tabs[gallery_title]:
            if st.session_state.alkane_main_molecule_props:
                main_props = st.session_state.alkane_main_molecule_props
                main_mol_name_display = main_props.get("IUPAC Name", st.session_state.alkane_molecule_searched.capitalize())
                with st.expander(f"Properties for Searched Alkane: {main_mol_name_display}", expanded=True):
                    display_properties(main_props)
                st.divider()
            if st.session_state.alkane_isomer_data:
                st.subheader(f"Isomers found for: {st.session_state.alkane_molecule_searched.capitalize()}")
                num_columns_gallery = 3 
                gallery_cols = st.columns(num_columns_gallery)
                for i, isomer in enumerate(st.session_state.alkane_isomer_data):
                    with gallery_cols[i % num_columns_gallery]:
                        container = st.container(border=True) 
                        pil_image_isomer = isomer["image"]
                        # Legend is now part of the image, caption can show CID or be removed if legend is sufficient
                        container.image(pil_image_isomer, caption=f"CID: {isomer['cid']}", use_container_width=True) 
                        if pil_image_isomer:
                            img_bytes_isomer = image_to_bytes(pil_image_isomer)
                            container.download_button(label="Download 2D", data=img_bytes_isomer, file_name=f"{isomer['name'].replace(' ', '_')}_CID_{isomer['cid']}_2D.png", mime="image/png", key=f"download_iso_{isomer['cid']}_v3")
                        container.markdown(f"<small>SMILES: {isomer['smiles']}</small>", unsafe_allow_html=True)
                        if container.button(f"View 3D", key=f"btn_3d_isomer_{isomer['cid']}_v3"):
                            st.session_state.selected_isomer_cid_for_3d, st.session_state.selected_isomer_name_for_3d, st.session_state.current_isomer_3d_index = isomer['cid'], isomer['name'], i; st.rerun()
            elif st.session_state.alkane_molecule_searched: st.info("No structural isomers (matching criteria) found to display in the gallery.")
    
    # If props_title was not part of gallery_title implicitly (e.g., no isomers but main_props exist)
    if props_title in active_tabs and st.session_state.alkane_main_molecule_props and gallery_title not in active_tabs:
         with active_tabs[props_title]:
            main_props_standalone = st.session_state.alkane_main_molecule_props
            main_mol_name_standalone = main_props_standalone.get("IUPAC Name", st.session_state.alkane_molecule_searched.capitalize())
            st.subheader(f"Properties for Searched Alkane: {main_mol_name_standalone}")
            display_properties(main_props_standalone)

    if view3d_title in active_tabs and st.session_state.selected_isomer_cid_for_3d:
        with active_tabs[view3d_title]:
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
            selected_style_label = st.radio("Select Display Style:", options=style_labels, key="radio_3d_style_isomer_v3", horizontal=True, index=current_style_index)
            if style_options_map.get(selected_style_label) != st.session_state.selected_3d_style:
                st.session_state.selected_3d_style = style_options_map[selected_style_label]; st.rerun() 
            nav_col_layout = [1, 0.1, 1.5, 0.1, 1]; prev_col, _, clear_col, _, next_col = st.columns(nav_col_layout)
            with prev_col:
                if st.button("⬅️ Previous", key="prev_3d_isomer_v3", help="View Previous Isomer", use_container_width=True, disabled=(st.session_state.current_isomer_3d_index <= 0)):
                    st.session_state.current_isomer_3d_index -= 1; prev_isomer = st.session_state.alkane_isomer_data[st.session_state.current_isomer_3d_index]
                    st.session_state.selected_isomer_cid_for_3d, st.session_state.selected_isomer_name_for_3d = prev_isomer['cid'], prev_isomer['name']; st.rerun()
            with clear_col:
                if st.button("Clear Isomer 3D View", key="clear_3d_isomer_v3", use_container_width=True):
                    st.session_state.selected_isomer_cid_for_3d, st.session_state.selected_isomer_name_for_3d, st.session_state.current_isomer_3d_index = None, "", -1; st.rerun()
            with next_col:
                if st.button("Next ➡️", key="next_3d_isomer_v3", help="View Next Isomer", use_container_width=True, disabled=(st.session_state.current_isomer_3d_index >= len(st.session_state.alkane_isomer_data) - 1)):
                    st.session_state.current_isomer_3d_index += 1; next_isomer = st.session_state.alkane_isomer_data[st.session_state.current_isomer_3d_index]
                    st.session_state.selected_isomer_cid_for_3d, st.session_state.selected_isomer_name_for_3d = next_isomer['cid'], next_isomer['name']; st.rerun()
            st.divider()
            with st.spinner(f"Loading 3D structure for {st.session_state.selected_isomer_name_for_3d}..."):
                sdf_data, error = cached_get_sdf_content(st.session_state.selected_isomer_cid_for_3d)
                if sdf_data:
                    html_3d = generate_3d_viewer_html(sdf_data, st.session_state.selected_isomer_name_for_3d, display_style=st.session_state.selected_3d_style, width=600, height=450)
                    st.components.v1.html(html_3d, height=470, width=620, scrolling=False)

# Content for General Molecule Search
if st.session_state.last_search_type == "general" and st.session_state.general_molecule_data:
    g_data = st.session_state.general_molecule_data
    if props_title in active_tabs and g_data.get("properties"): # Check if tab was created
        with active_tabs[props_title]:
            st.subheader(f"Properties for: {g_data['name']}")
            display_properties(g_data["properties"])
            st.markdown("---")
    if view2d_title in active_tabs and g_data.get("image_2d"):
        with active_tabs[view2d_title]:
            st.subheader(f"2D Structure for: {g_data['name']}")
            pil_image_general = g_data["image_2d"]
            st.image(pil_image_general, use_container_width=True) # Legend is now part of the image
            if pil_image_general:
                img_bytes_general = image_to_bytes(pil_image_general)
                st.download_button(label="Download 2D Structure", data=img_bytes_general, file_name=f"{g_data['name'].replace(' ', '_')}_CID_{g_data['cid']}_2D.png", mime="image/png", key=f"download_general_{g_data['cid']}_v3")
    if view3d_title in active_tabs and g_data.get("cid"):
        with active_tabs[view3d_title]:
            st.subheader(f"3D Structure for: {g_data['name']}")
            style_options_map_g = {'Stick': 'stick', 'Line': 'line', 'Ball and Stick': 'ball_and_stick'}
            style_labels_g = list(style_options_map_g.keys())
            try:
                current_style_label_g = [k for k, v in style_options_map_g.items() if v == st.session_state.selected_3d_style][0]
                current_style_index_g = style_labels_g.index(current_style_label_g)
            except IndexError: 
                current_style_index_g = 0 
                if style_labels_g: st.session_state.selected_3d_style = style_options_map_g[style_labels_g[0]]
                else: st.session_state.selected_3d_style = 'stick'
            selected_style_label_g = st.radio("Select Display Style:", options=style_labels_g, key="radio_3d_style_general_v3", horizontal=True, index=current_style_index_g)
            if style_options_map_g.get(selected_style_label_g) != st.session_state.selected_3d_style:
                st.session_state.selected_3d_style = style_options_map_g[selected_style_label_g]; st.rerun() 
            st.divider()
            with st.spinner(f"Loading 3D structure for {g_data['name']}..."):
                sdf_data, error = cached_get_sdf_content(g_data['cid'])
                if sdf_data:
                    html_3d = generate_3d_viewer_html(sdf_data, g_data['name'], display_style=st.session_state.selected_3d_style, width=600, height=450)
                    st.components.v1.html(html_3d, height=470, width=620, scrolling=False)
                else: st.info(f"3D structure could not be loaded. {error if error else ''}")

st.sidebar.markdown("---")
st.sidebar.caption("This application helps visualize alkane isomers and general chemical compounds. "
                 "Data is sourced from PubChem via PubChemPy. "
                 "2D structures are rendered by RDKit, and 3D structures by py3Dmol.")
st.sidebar.markdown("Built with Streamlit.")
