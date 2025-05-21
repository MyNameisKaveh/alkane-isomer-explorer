import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage
import gradio as gr
import traceback
import py3Dmol # New import

# --- Existing functions (slightly modified returns) ---
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

# --- New function for 3D rendering ---
def render_3d_molecule(smiles_string):
    """
    Renders a 3D molecule viewer using py3Dmol for the given SMILES string.
    """
    if not smiles_string:
        return "لطفا یک مولکول را انتخاب کنید تا ساختار سه‌بعدی آن نمایش داده شود.", gr.update(value="", visible=False)

    try:
        # Generate a basic 3D structure (e.g., using ETKDG for conformer generation)
        mol = Chem.MolFromSmiles(smiles_string)
        if not mol:
            return f"خطا: SMILES '{smiles_string}' قابل پردازش نیست.", gr.update(value="", visible=False)
        
        # Add hydrogens for better 3D representation
        mol = Chem.AddHs(mol)
        
        # Generate 3D conformer
        Chem.AllChem.EmbedMolecule(mol, Chem.AllChem.ETKDGv2())
        Chem.AllChem.MMFFOptimizeMolecule(mol) # Optional: optimize geometry

        sdf_string = Chem.MolToMolBlock(mol)

        # Create py3Dmol viewer
        viewer = py3Dmol.view(width=400, height=400) # Smaller size for clarity
        viewer.addModel(sdf_string, 'sdf')
        viewer.setStyle({'stick':{}}) # Display as sticks
        viewer.zoomTo()
        
        html_content = viewer.replicate_html()
        return "ساختار سه‌بعدی مولکول:", gr.update(value=html_content, visible=True)
    except Exception as e:
        print(f"Error rendering 3D molecule for SMILES {smiles_string}: {e}")
        return f"خطا در نمایش سه‌بعدی: {e}", gr.update(value="", visible=False)

# --- Function to update 3D viewer based on dropdown selection ---
def on_3d_dropdown_select(selected_name, all_isomers_data):
    """
    Called when an item is selected from the 3D dropdown.
    Finds the SMILES for the selected name and renders the 3D view.
    """
    if not selected_name:
        return "", gr.update(value="", visible=False) # Clear 3D viewer if no selection
    
    selected_smiles = None
    for name, smiles in all_isomers_data:
        if name == selected_name:
            selected_smiles = smiles
            break
    
    if selected_smiles:
        status_text, html_content = render_3d_molecule(selected_smiles)
        return status_text, html_content
    else:
        return "خطا: ایزومر انتخاب شده یافت نشد.", gr.update(value="", visible=False)

# --- Modified main function to return isomer data for 3D ---
def find_and_display_isomers(molecule_name_input):
    """
    Finds and displays structural alkane isomers for a given molecule name.
    Now also returns data for 3D viewer.
    """
    # Initialize outputs for gr.update(). We need to return values for all interface outputs.
    # The state component will store (name, smiles) tuples.
    empty_gallery = gr.update(value=[], visible=False)
    initial_status_text = ""
    initial_status_text_visible = False
    empty_isomer_data_state = gr.update(value=[], visible=False) # The actual state component
    empty_3d_dropdown_choices = gr.update(choices=[], value=None, visible=False) # Dropdown for 3D
    empty_3d_status_text = "" # Text for 3D tab
    empty_3d_html_viewer = gr.update(value="", visible=False) # 3D viewer itself

    if not molecule_name_input or not molecule_name_input.strip():
        initial_status_text = "لطفا نام یک مولکول را وارد کنید."
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
            status_message = f"مولکول '{molecule_name}' در PubChem یافت نشد. لطفا املای آن را بررسی کنید."
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
            status_message = f"آلکان ساختاری استاندارد با نام '{molecule_name}' در PubChem یافت نشد. (این ابزار تنها آلکان‌های شامل کربن و هیدروژن، بدون حلقه، بدون پیوند دوگانه/سه‌گانه و بدون ایزوتوپ را جستجو می‌کند.)"
            print(status_message)
            return empty_gallery, gr.update(value=status_message, visible=True), \
                   empty_isomer_data_state, empty_3d_dropdown_choices, empty_3d_status_text, empty_3d_html_viewer
        
        print(f"Proceeding with main compound: CID {main_compound_obj.cid}, Formula: {molecular_formula}")
        print(f"Searching for isomers with formula: {molecular_formula} (up to 200 candidates)...")
        isomers_found_raw = pcp.get_compounds(molecular_formula, 'formula', listkey_count=200) 

        if not isomers_found_raw:
            status_message = f"ایزومری برای فرمول {molecular_formula} یافت نشد."
            print(status_message)
            return empty_gallery, gr.update(value=status_message, visible=True), \
                   empty_isomer_data_state, empty_3d_dropdown_choices, empty_3d_status_text, empty_3d_html_viewer

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
            status_message = "ایزومر ساختاری آلکان استاندارد و قابل رسمی برای مولکول وارد شده پیدا نشد."
            if len(valid_structural_alkanes_entries) > 0:
                status_message += " (برخی ایزومرهای شناسایی شده در مرحله رسم ناموفق بودند یا در فیلترهای نهایی رد شدند.)"
        else:
            status_message = (
                f"{len(isomer_outputs_final_2d)} ایزومر ساختاری آلکان برای '{molecule_name_input}' "
                f"(فرمول: {molecular_formula}) پیدا و نمایش داده شد. "
                f"توجه: این ابزار تنها ایزومرهای شامل کربن و هیدروژن، بدون حلقه، بدون پیوند چندگانه و بدون ایزوتوپ را شناسایی می‌کند. "
                f"(ممکن است ایزومرهای بیشتری نیز در PubChem وجود داشته باشند که در این جستجو دریافت نشده‌اند.)"
            )
        
        isomer_outputs_final_2d.sort(key=lambda x: x[1])
        all_isomers_for_3d_data.sort(key=lambda x: x[0]) # Sort 3D data by name too
        
        # Prepare 3D dropdown choices
        dropdown_choices = [name for name, smiles in all_isomers_for_3d_data]
        
        # Return all outputs, including the state and 3D dropdown updates
        return gr.update(value=isomer_outputs_final_2d, visible=True), \
               gr.update(value=status_message, visible=True), \
               gr.update(value=all_isomers_for_3d_data, visible=True), \
               gr.update(choices=dropdown_choices, value=None, visible=True, interactive=True), \
               "لطفا یک ایزومر از لیست بالا انتخاب کنید تا ساختار سه‌بعدی آن را ببینید.", \
               gr.update(value="", visible=False)


    except pcp.PubChemHTTPError as e:
        error_msg = f"خطا در ارتباط با PubChem: {e}. لطفا اتصال اینترنت خود را بررسی کنید یا بعداً امتحان کنید."
        print(error_msg)
        print(f"FULL TRACEBACK for PubChemHTTPError: {traceback.format_exc()}")
        return empty_gallery, gr.update(value=error_msg, visible=True), \
               empty_isomer_data_state, empty_3d_dropdown_choices, empty_3d_status_text, empty_3d_html_viewer
    except Exception as e:
        error_msg = f"یک خطای غیرمنتظره در سرور رخ داد: {e}"
        print(f"FULL TRACEBACK for general Exception: {traceback.format_exc()}")
        return empty_gallery, gr.update(value=error_msg, visible=True), \
               empty_isomer_data_state, empty_3d_dropdown_choices, empty_3d_status_text, empty_3d_html_viewer


# --- Gradio Interface ---
with gr.Blocks(theme=gr.themes.Soft()) as demo:
    gr.Markdown("# یابنده و نمایشگر ایزومرهای ساختاری آلکان")
    gr.Markdown(
        "نام یک آلکان (به انگلیسی) را وارد کنید تا ایزومرهای ساختاری آن (تنها شامل کربن و هیدروژن، "
        "بدون حلقه، بدون پیوند چندگانه، بدون ایزوتوپ و بدون قطعات جدا شده) به همراه ساختار شیمیایی و نام IUPAC (یا رایج) نمایش داده شوند.\n"
        "اطلاعات از دیتابیس PubChem دریافت شده و ساختارها با استفاده از کتابخانه‌های RDKit و py3Dmol رسم می‌شوند."
    )

    with gr.Row():
        molecule_input = gr.Textbox(
            label="نام آلکان را وارد کنید", 
            placeholder="مثال: butane, pentane, hexane",
            info="نام آلکان مورد نظر خود را به انگلیسی و با حروف کوچک وارد کنید. (مانند: pentane)",
            scale=3
        )
        submit_btn = gr.Button("جستجو", scale=1)

    # Hidden state component to store isomer data for 3D viewer
    isomer_data_state = gr.State(value=[], visible=False) 

    with gr.Tabs() as tabs:
        with gr.TabItem("ایزومرهای دو بعدی", id="tab_2d"):
            gallery_output = gr.Gallery(
                label="ایزومرهای یافت شده (2D)", 
                columns=[3], 
                height="auto", 
                object_fit="contain",
                visible=False
            )
            status_textbox = gr.Textbox(label="وضعیت و پیام‌ها", visible=False)

        with gr.TabItem("ساختار سه‌بعدی", id="tab_3d"):
            gr.Markdown("### انتخاب ایزومر برای نمایش سه‌بعدی")
            # Dropdown for 3D selection, initially hidden and empty
            three_d_dropdown = gr.Dropdown(
                label="انتخاب ایزومر", 
                choices=[], 
                interactive=True, 
                visible=False,
                info="ایزومر مورد نظر را از لیست انتخاب کنید تا ساختار سه‌بعدی آن نمایش داده شود."
            )
            # Textbox for 3D status (e.g., instructions or errors)
            three_d_status_text = gr.Textbox(label="وضعیت نمایش سه‌بعدی", value="لطفا یک آلکان را جستجو کنید.", visible=True)
            # HTML component for the 3D viewer
            three_d_viewer = gr.HTML(label="نمایش سه‌بعدی", visible=False, value="")
    
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
