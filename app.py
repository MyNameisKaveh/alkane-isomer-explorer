import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage
from rdkit.Chem import AllChem
import gradio as gr
import traceback
import py3Dmol

def draw_molecule(smiles_string):
    """
    Renders a 2D molecule image from a SMILES string.
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

def render_3d_molecule(smiles_string):
    """
    Renders a 3D molecule viewer using py3Dmol for the given SMILES string,
    generating a standalone HTML content compatible with Gradio.
    """
    if not smiles_string:
        return "Please select a molecule to display its 3D structure.", gr.update(value="", visible=False)

    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if not mol:
            print(f"DEBUG: RDKit could not parse SMILES for 3D: {smiles_string}")
            return f"Error: SMILES '{smiles_string}' cannot be processed.", gr.update(value="", visible=False)
        
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv2())
        AllChem.MMFFOptimizeMolecule(mol)
        sdf_string = Chem.MolToMolBlock(mol)
        
        # Debugging: Print SDF string to verify
        print(f"DEBUG: SDF String for {smiles_string}:\n{sdf_string}")
        
        # Escape SDF string for JavaScript
        sdf_escaped = sdf_string.replace("'", "\\'").replace("\n", "\\n")
        full_html = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>3D Molecule Viewer</title>
            <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
            <style>
                body {{ margin: 0; padding: 0; }}
                #viewer {{ width: 400px; height: 400px; border: 1px solid #ccc; }}
            </style>
        </head>
        <body>
            <div id="viewer"></div>
            <script>
                // Check if 3Dmol is loaded
                if (typeof $3Dmol === 'undefined') {{
                    document.getElementById('viewer').innerText = 'Error: 3Dmol.js failed to load. Check your internet connection or browser settings.';
                }} else {{
                    var viewer = $3Dmol.createViewer(document.getElementById('viewer'), {{ backgroundColor: "white" }});
                    try {{
                        viewer.addModel('{sdf_escaped}', 'sdf');
                        viewer.setStyle({{ stick: {{}} }});
                        viewer.zoomTo();
                        viewer.render();
                    }} catch (e) {{
                        document.getElementById('viewer').innerText = 'Error rendering 3D model: ' + e.message;
                    }}
                }}
            </script>
        </body>
        </html>
        """
        
        print(f"DEBUG: Generated 3D HTML content length for '{smiles_string}': {len(full_html)}")
        return "3D structure of the molecule:", gr.update(value=full_html, visible=True)
    except Exception as e:
        print(f"Error rendering 3D molecule for SMILES {smiles_string}: {e}")
        return f"Error in 3D display: {e}", gr.update(value="", visible=False)

def on_3d_dropdown_select(selected_name, all_isomers_data):
    """
    Called when an item is selected from the 3D dropdown.
    Finds the SMILES for the selected name and renders the 3D view.
    """
    if not selected_name:
        return "", gr.update(value="", visible=False)
    
    selected_smiles = None
    for name, smiles in all_isomers_data:
        if name == selected_name:
            selected_smiles = smiles
            break
    
    if selected_smiles:
        status_text, html_content_update = render_3d_molecule(selected_smiles)
        return status_text, html_content_update
    else:
        return "Error: Selected isomer not found.", gr.update(value="", visible=False)

def find_and_display_isomers(molecule_name_input):
    """
    Finds and displays structural alkane isomers for a given molecule name.
    Also prepares data for the 3D viewer.
    """
    empty_gallery = gr.update(value=[], visible=False)
    initial_status_text = ""
    initial_status_text_visible = False
    empty_isomer_data_state = []
    empty_3d_dropdown_choices = gr.update(choices=[], value=None, visible=False)
    empty_3d_status_text = ""
    empty_3d_html_viewer = gr.update(value="", visible=False)

    if not molecule_name_input or not molecule_name_input.strip():
        initial_status_text = "Please enter a molecule name."
        initial_status_text_visible = True
        return empty_gallery, gr.update(value=initial_status_text, visible=initial_status_text_visible), \
               empty_isomer_data_state, empty_3d_dropdown_choices, empty_3d_status_text, empty_3d_html_viewer

    molecule_name = molecule_name_input.strip().lower()
    print(f"Processing request for: '{molecule_name}'")
    
    try:
        print(f"Searching for compound: '{molecule_name}' in PubChem (up to 10 candidates)...")
        compounds = pcp.get_compounds(molecule_name, 'name', listkey_count=10)
        
        main_compound_obj = None
        molecular_formula = None

        if not compounds:
            status_message = f"Molecule '{molecule_name}' not found in PubChem. Please check its spelling."
            print(status_message)
            return empty_gallery, gr.update(value=status_message, visible=True), \
                   empty_isomer_data_state, empty_3d_dropdown_choices, empty_3d_status_text, empty_3d_html_viewer
        
        print(f"Found {len(compounds)} potential matches for '{molecule_name}'. Checking them for standard alkane properties...")
        
        for i, c in enumerate(compounds):
            cid = c.cid
            common_name = c.synonyms[0] if c.synonyms else "N/A"
            actual_formula = c.molecular_formula if hasattr(c, 'molecular_formula') else None
            
            print(f"  Checking main compound candidate {i+1}: CID {cid}, Name: '{common_name}', Formula: '{actual_formula}'")

            is_standard_alkane_candidate = True

            if not actual_formula:
                is_standard_alkane_candidate = False
                print(f"    Main candidate CID {cid} has no molecular_formula.")
                
            if is_standard_alkane_candidate and not c.canonical_smiles:
                is_standard_alkane_candidate = False
                print(f"    Main candidate CID {cid} has no SMILES string for detailed check.")

            if is_standard_alkane_candidate:
                try:
                    mol_obj = Chem.MolFromSmiles(c.canonical_smiles)
                    if mol_obj:
                        if len(Chem.GetMolFrags(mol_obj)) > 1:
                            print(f"    Main candidate CID {cid} is disconnected.")
                            is_standard_alkane_candidate = False

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

                        if is_standard_alkane_candidate:
                            for bond in mol_obj.GetBonds():
                                if bond.GetBondType() != Chem.BondType.SINGLE:
                                    print(f"    Main candidate CID {cid} has non-single bond: {bond.GetBondType()}")
                                    is_standard_alkane_candidate = False
                                    break
                            if Chem.GetSymmSSSR(mol_obj):
                                print(f"    Main candidate CID {cid} has rings.")
                                is_standard_alkane_candidate = False

                    else:
                        print(f"    Main candidate CID {cid} SMILES '{c.canonical_smiles}' could not be parsed by RDKit.")
                        is_standard_alkane_candidate = False
                except Exception as rdkit_ex:
                    print(f"    RDKit error processing SMILES for main candidate CID {cid}: {rdkit_ex}")
                    is_standard_alkane_candidate = False
            
            if not is_standard_alkane_candidate:
                continue

            current_compound_name_matches_input = molecule_name in [syn.lower() for syn in c.synonyms]
            if current_compound_name_matches_input:
                main_compound_obj = c
                molecular_formula = actual_formula
                print(f"  SELECTED main compound based on name match: CID {main_compound_obj.cid}, Formula: {molecular_formula}")
                break
            
            if not main_compound_obj:
                main_compound_obj = c
                molecular_formula = actual_formula
                print(f"  TENTATIVELY selected first valid alkane candidate: CID {main_compound_obj.cid}, Formula: {molecular_formula}")
        
        if not main_compound_obj or not molecular_formula:
            status_message = f"Standard structural alkane with name '{molecule_name}' not found in PubChem. (This tool only searches for alkanes containing carbon and hydrogen, without rings, without double/triple bonds, and without isotopes.)"
            print(status_message)
            return empty_gallery, gr.update(value=status_message, visible=True), \
                   empty_isomer_data_state, empty_3d_dropdown_choices, empty_3d_status_text, empty_3d_html_viewer
        
        print(f"Proceeding with main compound: CID {main_compound_obj.cid}, Formula: {molecular_formula}")
        print(f"Searching for isomers with formula: {molecular_formula} (up to 200 candidates)...")
        isomers_found_raw = pcp.get_compounds(molecular_formula, 'formula', listkey_count=200)

        if not isomers_found_raw:
            status_message = f"No isomer found for formula {molecular_formula}."
            print(status_message)
            return empty_gallery, gr.update(value=status_message, visible=True), \
                   empty_isomer_data_state, empty_3d_dropdown_choices, empty_3d_status_text, empty_3d_html_viewer

        print(f"Found {len(isomers_found_raw)} potential isomer entries from PubChem. Filtering for true structural alkane isomers...")
        
        valid_structural_alkanes_entries = []
        unique_accepted_smiles = set()
        all_isomers_for_3d_data = []

        for isomer_entry in isomers_found_raw:
            smiles = isomer_entry.canonical_smiles
            if not smiles:
                print(f"  Skipping isomer without SMILES: CID {isomer_entry.cid}")
                continue

            try:
                mol_iso = Chem.MolFromSmiles(smiles)
                if not mol_iso:
                    print(f"  FILTERED (Invalid SMILES parse): CID {isomer_entry.cid}, SMILES: {smiles}")
                    continue

                is_valid_alkane_isomer = True
                
                if len(Chem.GetMolFrags(mol_iso)) > 1:
                    print(f"  FILTERED (Disconnected): CID {isomer_entry.cid}, NumFrags: {len(Chem.GetMolFrags(mol_iso))}, SMILES: {smiles}")
                    is_valid_alkane_isomer = False
                
                if is_valid_alkane_isomer:
                    atom_symbols = set()
                    for atom in mol_iso.GetAtoms():
                        atom_symbols.add(atom.GetSymbol())
                        if atom.GetIsotope() != 0:
                            print(f"  FILTERED (Isotope present): CID {isomer_entry.cid}, Atom: {atom.GetSymbol()}{atom.GetIdx()+1}, Isotope: {atom.GetIsotope()}, SMILES: {smiles}")
                            is_valid_alkane_isomer = False
                            break
                    if not atom_symbols.issubset({'C', 'H'}):
                        print(f"  FILTERED (Non-CH elements): CID {isomer_entry.cid}, Elements: {atom_symbols}, SMILES: {smiles}")
                        is_valid_alkane_isomer = False
                
                if is_valid_alkane_isomer:
                    for bond in mol_iso.GetBonds():
                        if bond.GetBondType() != Chem.BondType.SINGLE:
                            print(f"  FILTERED (Non-single bond): CID {isomer_entry.cid}, BondType: {bond.GetBondType()}, SMILES: {smiles}")
                            is_valid_alkane_isomer = False
                            break
                    if Chem.GetSymmSSSR(mol_iso):
                        print(f"  FILTERED (Has rings): CID {isomer_entry.cid}, SMILES: {smiles}")
                        is_valid_alkane_isomer = False

                if is_valid_alkane_isomer:
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
        
        isomer_outputs_final_2d = []
        
        for final_isomer_entry in valid_structural_alkanes_entries:
            smiles_to_draw = final_isomer_entry.canonical_smiles
            
            isomer_display_name = final_isomer_entry.iupac_name
            if not isomer_display_name and final_isomer_entry.synonyms:
                simple_alkane_names = [
                    s for s in final_isomer_entry.synonyms 
                    if s.lower().endswith("ane") and not any(char.isdigit() for char in s) and '-' not in s
                ]
                if simple_alkane_names:
                    isomer_display_name = min(simple_alkane_names, key=len)
                else:
                    isomer_display_name = final_isomer_entry.synonyms[0]
            
            if not isomer_display_name:
                isomer_display_name = f"Alkane (CID: {final_isomer_entry.cid})"
            
            isomer_display_name = isomer_display_name.capitalize()

            mol_image = draw_molecule(smiles_to_draw)
            if mol_image:
                isomer_outputs_final_2d.append((mol_image, f"{isomer_display_name}\nSMILES: {smiles_to_draw}"))
                all_isomers_for_3d_data.append((isomer_display_name, smiles_to_draw))
            else:
                print(f"  Failed to draw image for accepted isomer: CID {final_isomer_entry.cid}, SMILES: {smiles_to_draw}")

        print(f"Displayed {len(isomer_outputs_final_2d)} isomers in the 2D gallery.")
        print(f"Prepared {len(all_isomers_for_3d_data)} isomers for 3D selection.")

        if not isomer_outputs_final_2d:
            status_message = "No standard and drawable structural alkane isomer found for the entered molecule."
            if len(valid_structural_alkanes_entries) > 0:
                status_message += " (Some identified isomers failed to render or were rejected in final filters.)"
        else:
            status_message = (
                f"{len(isomer_outputs_final_2d)} structural alkane isomer(s) found and displayed for '{molecule_name_input}' "
                f"(Formula: {molecular_formula}). "
                f"Note: This tool only identifies isomers containing carbon and hydrogen, without rings, without multiple bonds, and without isotopes. "
                f"(More isomers might exist in PubChem that were not retrieved in this search.)"
            )
        
        isomer_outputs_final_2d.sort(key=lambda x: x[1])
        all_isomers_for_3d_data.sort(key=lambda x: x[0])
        
        dropdown_choices = [name for name, smiles in all_isomers_for_3d_data]
        
        return gr.update(value=isomer_outputs_final_2d, visible=True), \
               gr.update(value=status_message, visible=True), \
               all_isomers_for_3d_data, \
               gr.update(choices=dropdown_choices, value=None, visible=True, interactive=True), \
               "Please select an isomer from the list above to view its 3D structure.", \
               gr.update(value="", visible=False)

    except pcp.PubChemHTTPError as e:
        error_msg = f"Error communicating with PubChem: {e}. Please check your internet connection or try again later."
        print(error_msg)
        print(f"FULL TRACEBACK for PubChemHTTPError: {traceback.format_exc()}")
        return empty_gallery, gr.update(value=error_msg, visible=True), \
               empty_isomer_data_state, empty_3d_dropdown_choices, empty_3d_status_text, empty_3d_html_viewer
    except Exception as e:
        error_msg = f"An unexpected server error occurred: {e}"
        print(f"FULL TRACEBACK for general Exception: {traceback.format_exc()}")
        return empty_gallery, gr.update(value=error_msg, visible=True), \
               empty_isomer_data_state, empty_3d_dropdown_choices, empty_3d_status_text, empty_3d_html_viewer

with gr.Blocks() as demo:
    gr.Markdown("# Alkane Structural Isomer Finder and Viewer")
    gr.Markdown(
        "Enter an alkane name (in English) to display its structural isomers (containing only carbon and hydrogen, "
        "without rings, without multiple bonds, without isotopes, and without disconnected fragments) along with their chemical structures and IUPAC (or common) names.\n"
        "Information is retrieved from the PubChem database, and structures are drawn using RDKit and py3Dmol libraries."
    )

    with gr.Row():
        molecule_input = gr.Textbox(
            label="Enter Alkane Name", 
            placeholder="Example: butane, pentane, hexane",
            info="Enter the desired alkane name in English and lowercase. (e.g.: pentane)",
            scale=3
        )
        submit_btn = gr.Button("Search", scale=1)

    isomer_data_state = gr.State(value=[])

    with gr.Tabs() as tabs:
        with gr.TabItem("2D Isomers", id="tab_2d"):
            gallery_output = gr.Gallery(
                label="Found Isomers (2D)", 
                columns=[3], 
                height="auto", 
                object_fit="contain",
                visible=False
            )
            status_textbox = gr.Textbox(label="Status and Messages", visible=False)

        with gr.TabItem("3D Structure", id="tab_3d"):
            gr.Markdown("### Select Isomer for 3D Display")
            three_d_dropdown = gr.Dropdown(
                label="Select Isomer", 
                choices=[], 
                interactive=True, 
                visible=False,
                info="Select the desired isomer from the list to view its 3D structure."
            )
            three_d_status_text = gr.Textbox(label="3D Display Status", value="Please search for an alkane.", visible=True)
            three_d_viewer = gr.HTML(label="3D Viewer", visible=False, value="")
    
    submit_btn.click(
        fn=find_and_display_isomers,
        inputs=molecule_input,
        outputs=[
            gallery_output, 
            status_textbox, 
            isomer_data_state,
            three_d_dropdown,
            three_d_status_text,
            three_d_viewer
        ]
    )

    three_d_dropdown.change(
        fn=on_3d_dropdown_select,
        inputs=[three_d_dropdown, isomer_data_state],
        outputs=[three_d_status_text, three_d_viewer]
    )

    gr.Examples(
        examples=[
            ["butane"], 
            ["pentane"], 
            ["hexane"],
            ["heptane"], 
            ["octane"]
        ],
        inputs=molecule_input
    )

demo.launch()
