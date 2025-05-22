import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import gradio as gr
import traceback
import py3Dmol
import io
import base64

def draw_molecule(smiles_string):
    """
    Renders a 2D molecule image from a SMILES string.
    """
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if mol:
            img = Draw.MolToImage(mol, size=(300, 300))
            return img
        else:
            print(f"Could not parse SMILES: {smiles_string}")
            return None
    except Exception as e:
        print(f"Error drawing molecule for SMILES {smiles_string}: {e}")
        return None

def generate_3d_view(smiles_string):
    """
    Generates a 3D visualization of a molecule from a SMILES string using py3Dmol.
    Returns raw HTML for Gradio.
    """
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if not mol:
            return "<div>Error: Invalid SMILES string.</div>"
        mol = Chem.AddHs(mol)  # Add hydrogens for proper 3D structure
        AllChem.EmbedMolecule(mol, randomSeed=42)  # Generate 3D coordinates
        AllChem.MMFFOptimizeMolecule(mol)  # Optimize geometry
        xyz = Chem.MolToXYZBlock(mol)  # Convert to XYZ format

        view = py3Dmol.view(width=400, height=400)
        view.addModel(xyz, "xyz")
        view.setStyle({'stick': {}})
        view.setBackgroundColor('white')
        view.zoomTo()
        
        # Generate raw HTML without requiring IPython
        html = view._repr_html_()
        return f"""
        <div style='width: 400px; height: 400px; position: relative;'>
            {html}
        </div>
        """
    except Exception as e:
        print(f"Error generating 3D view for SMILES {smiles_string}: {e}")
        return f"<div>Error generating 3D view: {str(e)}</div>"

def find_and_display_isomers(molecule_name_input, selected_isomer_index=0, isomer_data_state=None):
    """
    Finds and displays structural alkane isomers with 2D images and a single 3D viewer.
    Returns gallery, status, 3D view HTML, dropdown choices, selected isomer index, and updated isomer_data_state.
    """
    if not molecule_name_input or not molecule_name_input.strip():
        return gr.update(value=[], visible=True), gr.update(value="لطفا نام یک مولکول را وارد کنید.", visible=True), gr.update(value="", visible=True), gr.update(choices=[], visible=True), gr.update(value=0, visible=True), gr.update(value=[])

    molecule_name = molecule_name_input.strip().lower()
    print(f"Processing request for: '{molecule_name}'")
    
    status_message = ""
    isomer_data = []  # Store (name, SMILES) for dropdown

    try:
        print(f"Searching for compound: '{molecule_name}' in PubChem (up to 10 candidates)...")
        compounds = pcp.get_compounds(molecule_name, 'name', listkey_count=10) 
        
        main_compound_obj = None
        molecular_formula = None

        if not compounds:
            status_message = f"مولکول '{molecule_name}' در PubChem یافت نشد. لطفا املای آن را بررسی کنید."
            print(status_message)
            return gr.update(value=[], visible=True), gr.update(value=status_message, visible=True), gr.update(value="", visible=True), gr.update(choices=[], visible=True), gr.update(value=0, visible=True), gr.update(value=[])
        
        print(f"Found {len(compounds)} potential matches for '{molecule_name}'. Checking them for standard alkane properties...")
        
        for i, c in enumerate(compounds):
            cid = c.cid
            common_name = c.synonyms[0] if c.synonyms else "N/A"
            actual_formula = c.molecular_formula if hasattr(c, 'molecular_formula') else None
            
            print(f"  Checking main compound candidate {i+1}: CID {cid}, Name: '{common_name}', Formula: '{actual_formula}'")

            is_standard_alkane_candidate = True 

            if not actual_formula or not c.canonical_smiles:
                is_standard_alkane_candidate = False
                print(f"    Main candidate CID {cid} lacks formula or SMILES.")
            
            if is_standard_alkane_candidate:
                try:
                    mol_obj = Chem.MolFromSmiles(c.canonical_smiles)
                    if mol_obj:
                        if len(Chem.GetMolFrags(mol_obj)) > 1:
                            print(f"    Main candidate CID {cid} is disconnected.")
                            is_standard_alkane_candidate = False
                        atom_symbols_main = set(atom.GetSymbol() for atom in mol_obj.GetAtoms())
                        if not atom_symbols_main.issubset({'C', 'H'}) or any(atom.GetIsotope() != 0 for atom in mol_obj.GetAtoms()):
                            print(f"    Main candidate CID {cid} is not CH only or has isotopes.")
                            is_standard_alkane_candidate = False
                        for bond in mol_obj.GetBonds():
                            if bond.GetBondType() != Chem.BondType.SINGLE:
                                print(f"    Main candidate CID {cid} has non-single bond.")
                                is_standard_alkane_candidate = False
                                break
                        if Chem.GetSymmSSSR(mol_obj):
                            print(f"    Main candidate CID {cid} has rings.")
                            is_standard_alkane_candidate = False
                    else:
                        print(f"    Main candidate CID {cid} SMILES could not be parsed.")
                        is_standard_alkane_candidate = False
                except Exception as rdkit_ex:
                    print(f"    RDKit error for CID {cid}: {rdkit_ex}")
                    is_standard_alkane_candidate = False
            
            if is_standard_alkane_candidate:
                current_compound_name_matches_input = molecule_name in [syn.lower() for syn in c.synonyms]
                if current_compound_name_matches_input:
                    main_compound_obj = c
                    molecular_formula = actual_formula
                    print(f"  SELECTED main compound: CID {main_compound_obj.cid}, Formula: {molecular_formula}")
                    break
                if not main_compound_obj:
                    main_compound_obj = c
                    molecular_formula = actual_formula
                    print(f"  TENTATIVELY selected: CID {main_compound_obj.cid}, Formula: {molecular_formula}")
        
        if not main_compound_obj or not molecular_formula:
            status_message = f"آلکان ساختاری استاندارد با نام '{molecule_name}' یافت نشد."
            print(status_message)
            return gr.update(value=[], visible=True), gr.update(value=status_message, visible=True), gr.update(value="", visible=True), gr.update(choices=[], visible=True), gr.update(value=0, visible=True), gr.update(value=[])
        
        print(f"Searching for isomers with formula: {molecular_formula} (up to 200 candidates)...")
        isomers_found_raw = pcp.get_compounds(molecular_formula, 'formula', listkey_count=200)

        if not isomers_found_raw:
            status_message = f"ایزومری برای فرمول {molecular_formula} یافت نشد."
            print(status_message)
            return gr.update(value=[], visible=True), gr.update(value=status_message, visible=True), gr.update(value="", visible=True), gr.update(choices=[], visible=True), gr.update(value=0, visible=True), gr.update(value=[])

        print(f"Found {len(isomers_found_raw)} potential isomer entries. Filtering for valid alkane isomers...")
        
        valid_structural_alkanes_entries = []
        unique_accepted_smiles = set()

        for isomer_entry in isomers_found_raw:
            smiles = isomer_entry.canonical_smiles
            if not smiles:
                print(f"  Skipping isomer without SMILES: CID {isomer_entry.cid}")
                continue

            try:
                mol_iso = Chem.MolFromSmiles(smiles)
                if not mol_iso:
                    print(f"  FILTERED (Invalid SMILES): CID {isomer_entry.cid}, SMILES: {smiles}")
                    continue

                is_valid_alkane_isomer = True
                if len(Chem.GetMolFrags(mol_iso)) > 1:
                    print(f"  FILTERED (Disconnected): CID {isomer_entry.cid}, SMILES: {smiles}")
                    is_valid_alkane_isomer = False
                atom_symbols = set(atom.GetSymbol() for atom in mol_iso.GetAtoms())
                if not atom_symbols.issubset({'C', 'H'}) or any(atom.GetIsotope() != 0 for atom in mol_iso.GetAtoms()):
                    print(f"  FILTERED (Non-CH or isotopes): CID {isomer_entry.cid}, SMILES: {smiles}")
                    is_valid_alkane_isomer = False
                for bond in mol_iso.GetBonds():
                    if bond.GetBondType() != Chem.BondType.SINGLE:
                        print(f"  FILTERED (Non-single bond): CID {isomer_entry.cid}, SMILES: {smiles}")
                        is_valid_alkane_isomer = False
                        break
                if Chem.GetSymmSSSR(mol_iso):
                    print(f"  FILTERED (Has rings): CID {isomer_entry.cid}, SMILES: {smiles}")
                    is_valid_alkane_isomer = False

                if is_valid_alkane_isomer:
                    canonical_smiles = Chem.MolToSmiles(mol_iso, isomericSmiles=False)
                    if canonical_smiles not in unique_accepted_smiles:
                        print(f"  ACCEPTED: CID {isomer_entry.cid}, SMILES: {smiles}")
                        valid_structural_alkanes_entries.append(isomer_entry)
                        unique_accepted_smiles.add(canonical_smiles)
                    else:
                        print(f"  Skipping (Duplicate SMILES): CID {isomer_entry.cid}, SMILES: {smiles}")
            except Exception as rdkit_iso_ex:
                print(f"  RDKit error for CID {isomer_entry.cid}: {rdkit_iso_ex}, SMILES: {smiles}")
                continue
        
        print(f"Found {len(valid_structural_alkanes_entries)} unique valid isomers.")
        
        isomer_outputs_final = []
        isomer_data = []
        valid_isomers_count_final = 0

        for final_isomer_entry in valid_structural_alkanes_entries:
            smiles_to_draw = final_isomer_entry.canonical_smiles
            isomer_display_name = final_isomer_entry.iupac_name
            if not isomer_display_name and final_isomer_entry.synonyms:
                simple_alkane_names = [
                    s for s in final_isomer_entry.synonyms
                    if s.lower().endswith("ane") and not any(char.isdigit() for char in s) and '-' not in s
                ]
                isomer_display_name = min(simple_alkane_names, key=len) if simple_alkane_names else final_isomer_entry.synonyms[0]
            if not isomer_display_name:
                isomer_display_name = f"Alkane (CID: {final_isomer_entry.cid})"
            isomer_display_name = isomer_display_name.capitalize()

            mol_image = draw_molecule(smiles_to_draw)
            if mol_image:
                isomer_outputs_final.append((mol_image, f"{isomer_display_name}\nSMILES: {smiles_to_draw}"))
                isomer_data.append((isomer_display_name, smiles_to_draw))
                valid_isomers_count_final += 1
            else:
                print(f"  Failed to draw image for CID {final_isomer_entry.cid}, SMILES: {smiles_to_draw}")

        print(f"Displayed {valid_isomers_count_final} isomers in the gallery.")

        if not isomer_outputs_final:
            status_message = "ایزومر ساختاری آلکان استاندارد و قابل رسمی پیدا نشد."
            if len(valid_structural_alkanes_entries) > 0:
                status_message += " (برخی ایزومرها در رسم ناموفق بودند.)"
            return gr.update(value=[], visible=True), gr.update(value=status_message, visible=True), gr.update(value="", visible=True), gr.update(choices=[], visible=True), gr.update(value=0, visible=True), gr.update(value=[])
        else:
            status_message = (
                f"{len(isomer_outputs_final)} ایزومر ساختاری آلکان برای '{molecule_name_input}' "
                f"(فرمول: {molecular_formula}) پیدا و نمایش داده شد. "
                f"از منوی زیر یک ایزومر را برای مشاهده سه‌بعدی انتخاب کنید."
            )

        isomer_outputs_final.sort(key=lambda x: x[1])
        isomer_data.sort(key=lambda x: x[0])
        dropdown_choices = [name for name, _ in isomer_data]
        
        # Generate 3D view for the selected isomer
        selected_index = min(selected_isomer_index, len(isomer_data) - 1) if isomer_data else 0
        three_d_view = generate_3d_view(isomer_data[selected_index][1]) if isomer_data else ""

        return (
            gr.update(value=isomer_outputs_final, visible=True),
            gr.update(value=status_message, visible=True),
            gr.update(value=three_d_view, visible=True),
            gr.update(choices=dropdown_choices, visible=True),
            gr.update(value=selected_index, visible=True),
            gr.update(value=isomer_data)
        )

    except pcp.PubChemHTTPError as e:
        error_msg = f"خطا در ارتباط با PubChem: {e}. لطفا اتصال اینترنت خود را بررسی کنید."
        print(f"FULL TRACEBACK: {traceback.format_exc()}")
        return gr.update(value=[], visible=True), gr.update(value=error_msg, visible=True), gr.update(value="", visible=True), gr.update(choices=[], visible=True), gr.update(value=0, visible=True), gr.update(value=[])
    except Exception as e:
        error_msg = f"خطای غیرمنتظره: {e}"
        print(f"FULL TRACEBACK: {traceback.format_exc()}")
        return gr.update(value=[], visible=True), gr.update(value=error_msg, visible=True), gr.update(value="", visible=True), gr.update(choices=[], visible=True), gr.update(value=0, visible=True), gr.update(value=[])

def get_isomer_index(selected_name, isomer_data_state):
    """
    Maps selected isomer name to its index in isomer_data_state.
    """
    if not isomer_data_state or not selected_name:
        return 0
    for i, (name, _) in enumerate(isomer_data_state):
        if name == selected_name:
            return i
    return 0

# --- Gradio Interface ---
with gr.Blocks(theme=gr.themes.Soft()) as iface:
    gr.Markdown(
        """
        # یابنده و نمایشگر ایزومرهای ساختاری آلکان
        نام یک آلکان (به انگلیسی) را وارد کنید تا ایزومرهای ساختاری آن (تنها شامل کربن و هیدروژن، بدون حلقه، بدون پیوند چندگانه، بدون ایزوتوپ) به همراه ساختار شیمیایی و نام IUPAC (یا رایج) نمایش داده شوند.  
        اطلاعات از دیتابیس PubChem دریافت شده و ساختارها با استفاده از RDKit و py3Dmol رسم می‌شوند.
        """
    )
    
    molecule_input = gr.Textbox(
        label="نام آلکان را وارد کنید",
        placeholder="مثال: butane, pentane, hexane",
        info="نام آلکان را به انگلیسی و با حروف کوچک وارد کنید."
    )
    
    with gr.Row():
        with gr.Column():
            gallery = gr.Gallery(
                label="ایزومرهای یافت شده (2D)",
                columns=[3],
                height="auto",
                object_fit="contain",
                visible=False
            )
            status = gr.Textbox(label="وضعیت و پیام‌ها", visible=False)
        with gr.Column():
            three_d_viewer = gr.HTML(label="نمایش سه‌بعدی ایزومر", visible=False)
            isomer_dropdown = gr.Dropdown(
                label="انتخاب ایزومر برای نمایش سه‌بعدی",
                choices=[],
                value=None,
                visible=False,
                interactive=True
            )
    
    submit_button = gr.Button("جستجوی ایزومرها")
    
    # States to track selected isomer index and isomer data
    selected_isomer_state = gr.State(value=0)
    isomer_data_state = gr.State(value=[])
    
    # Define interactions
    submit_button.click(
        fn=find_and_display_isomers,
        inputs=[molecule_input, selected_isomer_state, isomer_data_state],
        outputs=[gallery, status, three_d_viewer, isomer_dropdown, selected_isomer_state, isomer_data_state]
    )
    
    isomer_dropdown.change(
        fn=get_isomer_index,
        inputs=[isomer_dropdown, isomer_data_state],
        outputs=selected_isomer_state
    ).then(
        fn=find_and_display_isomers,
        inputs=[molecule_input, selected_isomer_state, isomer_data_state],
        outputs=[gallery, status, three_d_viewer, isomer_dropdown, selected_isomer_state, isomer_data_state]
    )

    gr.Examples(
        examples=[["butane"], ["pentane"], ["hexane"], ["heptane"], ["octane"]],
        inputs=[molecule_input]
    )

if __name__ == '__main__':
    iface.launch()
