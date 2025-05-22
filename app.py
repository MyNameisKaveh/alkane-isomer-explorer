import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage
from rdkit.Chem import AllChem  # FIX: Import AllChem explicitly for 3D conformer generation
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

def render_3d_molecule_content(smiles_string):
    """
    Generates the HTML content for a 3D molecule viewer using py3Dmol.
    Returns (status_string, html_content_string).
    """
    if not smiles_string:
        return "Please select a molecule to display its 3D structure.", "" # Return empty HTML string

    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if not mol:
            print(f"DEBUG: RDKit could not parse SMILES for 3D: {smiles_string}")
            return f"Error: SMILES '{smiles_string}' cannot be processed.", ""
        
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv2())
        AllChem.MMFFOptimizeMolecule(mol)

        sdf_string = Chem.MolToMolBlock(mol)
        
        if not sdf_string or len(sdf_string.strip()) < 100: # A basic check for content
            print(f"DEBUG: Generated SDF string looks short/empty for SMILES {smiles_string}:\n{sdf_string}")
            return f"Error: Failed to generate valid 3D data for {smiles_string}.", ""

        viewer = py3Dmol.view(width=400, height=400)
        viewer.addModel(sdf_string, 'sdf')
        viewer.setStyle({'stick':{}})
        viewer.zoomTo()
        
        html_content = viewer.replicate_html()

        print(f"DEBUG: Generated 3D HTML content length for '{smiles_string}': {len(html_content)}")
        # print(f"DEBUG: First 500 chars of 3D HTML for '{smiles_string}':\n{html_content[:500]}...") # Too verbose for regular use
        
        return "3D structure of the molecule:", html_content
    except Exception as e:
        print(f"Error rendering 3D molecule for SMILES {smiles_string}: {e}")
        return f"Error in 3D display: {e}", "" # Return empty HTML string on error

def on_3d_dropdown_select(selected_name, all_isomers_data):
    """
    Called when an item is selected from the 3D dropdown.
    Finds the SMILES for the selected name and renders the 3D view.
    """
    if not selected_name:
        # If nothing selected, clear the 3D viewer and reset status
        return gr.update(value="Please select a molecule to display its 3D structure.", visible=True), gr.update(value="", visible=False)
    
    selected_smiles = None
    for name, smiles in all_isomers_data:
        if name == selected_name:
            selected_smiles = smiles
            break
    
    if selected_smiles:
        # Call the content generation function
        status_text_str, html_content_str = render_3d_molecule_content(selected_smiles)
        
        # Return gr.update objects for the outputs
        return gr.update(value=status_text_str, visible=True), gr.update(value=html_content_str, visible=True)
    else:
        return gr.update(value="Error: Selected isomer not found.", visible=True), gr.update(value="", visible=False)


def find_and_display_isomers(molecule_name_input):
    """
    Finds and displays structural alkane isomers for a given molecule name.
    Also prepares data for the 3D viewer.
    """
    # Initialize all outputs with appropriate gr.update() calls for initial visibility
    # For gr.State, just return the value as it's not a visible component
    empty_gallery = gr.update(value=[], visible=False)
    empty_status_textbox = gr.update(value="", visible=False)
    empty_isomer_data_state = [] 
    empty_3d_dropdown_choices = gr.update(choices=[], value=None, visible=False)
    initial_3d_status_text = gr.update(value="Please search for an alkane.", visible=True) # Set initial visible status
    empty_3d_html_viewer = gr.update(value="", visible=False)

    if not molecule_name_input or not molecule_name_input.strip():
        return empty_gallery, gr.update(value="Please enter a molecule name.", visible=True), \
               empty_isomer_data_state, empty_3d_dropdown_choices, initial_3d_status_text, empty_3d_html_viewer

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
                   empty_isomer_data_state, empty_3d_dropdown_choices, initial_3d_status_text, empty_3d_html_viewer
        
        print(f"Found {len(compounds)} potential matches for '{molecule_name}'. Checking them for standard alkane properties...")
        
        for i, c in enumerate(compounds):
            cid = c.cid
            common_name = c.synonyms[0] if c.synonyms else "N/A"
            actual_formula = c.molecular_formula if hasattr(c, 'molecular_formula') else None
            
            print(f"  Checking main compound candidate {i+1}: CID {cid}, Name: '{common_name}', Formula: '{actual_formula}'")

            is_standard_alkane_candidate = True 

            if not actual_formula:
                is_standard_alkane_candidate = False
                print(f"    Main candidate CID {cid} has no molecular formula.")
                
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
                                    print(f"    Main candidate CID {cid} has non-standard isotope: {atom.GetSymbol()}{atom.GetIdx()+1}")
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
                   empty_isomer_data_state, empty_3d_dropdown_choices, initial_3d_status_text, empty_3d_html_viewer
        
        print(f"Proceeding with main compound: CID {main_compound_obj.cid}, Formula: {molecular_formula}")
        print(f"Searching for isomers with formula: {molecular_formula} (up to 200 candidates)...")
        isomers_found_raw = pcp.get_compounds(molecular_formula, 'formula', listkey_count=200) 

        if not isomers_found_raw:
            status_message = f"No isomer found for formula {molecular_formula}."
            print(status_message)
            return empty_gallery, gr.update(value=status_message, visible=True), \
                   empty_isomer_data_state, empty_3d_dropdown_choices, initial_3d_status_text, empty_3d_html_viewer

        print(f"Found {len(isomers_found_raw)} potential isomer entries from PubChem. Filtering for true structural alkane isomers...")
        
        valid_structural_alkanes_entries = [] 
        unique_accepted_smiles = set() 
        all_isomers_for_3d_data = [] # To store (name, smiles) for 3D dropdown

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
        
        isomer_outputs_final_2d = [] # For 2D gallery
        
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
                all_isomers_for_3d_data.append((isomer_display_name, smiles_to_draw)) # Store for 3D
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
        all_isomers_for_3d_data.sort(key=lambda x: x[0]) # Sort 3D data by name too
        
        dropdown_choices = [name for name, smiles in all_isomers_for_3d_data]
        
        # Return all outputs, ensuring gr.update() is used where visibility/initial value needs control
        return gr.update(value=isomer_outputs_final_2d, visible=True), \
               gr.update(value=status_message, visible=True), \
               all_isomers_for_3d_data, \
               gr.update(choices=dropdown_choices, value=None, visible=True, interactive=True), \
               gr.update(value="Please select an isomer from the list above to view its 3D structure.", visible=True), \
               gr.update(value="", visible=False) # Clear 3D viewer on new search

    except pcp.PubChemHTTPError as e:
        error_msg = f"Error communicating with PubChem: {e}. Please check your internet connection or try again later."
        print(error_msg)
        print(f"FULL TRACEBACK for PubChemHTTPError: {traceback.format_exc()}")
        return empty_gallery, gr.update(value=error_msg, visible=True), \
               empty_isomer_data_state, empty_3d_dropdown_choices, initial_3d_status_text, empty_3d_html_viewer
    except Exception as e:
        error_msg = f"An unexpected server error occurred: {e}"
        print(f"FULL TRACEBACK for general Exception: {traceback.format_exc()}")
        return empty_gallery, gr.update(value=error_msg, visible=True), \
               empty_isomer_data_state, empty_3d_dropdown_choices, initial_3d_status_text, empty_3d_html_viewer

# --- Gradio Interface ---
with gr.Blocks() as demo: # Theme argument removed as requested
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

    # Hidden state component to store isomer data for 3D viewer.
    # gr.State does not accept 'visible' argument as it's a non-visual component.
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
            # Dropdown for 3D selection, initially hidden and empty
            three_d_dropdown = gr.Dropdown(
                label="Select Isomer", 
                choices=[], 
                interactive=True, 
                visible=False,
                info="Select the desired isomer from the list to view its 3D structure."
            )
            # Textbox for 3D status (e.g., instructions or errors)
            # Initial message for 3D status, visible when app starts
            three_d_status_text = gr.Textbox(label="3D Display Status", value="Please search for an alkane.", visible=True)
            # HTML component for the 3D viewer
            three_d_viewer = gr.HTML(label="3D Viewer", visible=False, value="")
    
    # Define how interactions trigger functions
    submit_btn.click(
        fn=find_and_display_isomers,
        inputs=molecule_input,
        outputs=[
            gallery_output, 
            status_textbox, 
            isomer_data_state, # This output updates the hidden state
            three_d_dropdown, # This output updates the dropdown choices
            three_d_status_text, # Updates 3D tab status text
            three_d_viewer # Clears/hides 3D viewer initially
        ]
    )

    # When an isomer is selected from the 3D dropdown, trigger 3D rendering
    three_d_dropdown.change(
        fn=on_3d_dropdown_select,
        inputs=[three_d_dropdown, isomer_data_state], # Pass selected name and all isomer data
        outputs=[three_d_status_text, three_d_viewer] # Update 3D status and viewer
    )

    # Optional: Example inputs
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
