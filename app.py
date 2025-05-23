import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage
import gradio as gr
import traceback
# import py3Dmol # فعلا استفاده نمی‌شود
import os
# import time # فعلا استفاده نمی‌شود

# تابع رسم مولکول دوبعدی (بدون تغییر)
def draw_molecule(smiles_string):
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

# تابع اصلی پردازش آلکان (باید کد کامل شما باشد)
# این تابع فقط خروجی‌های مربوط به گالری، وضعیت و CIDs را برمی‌گرداند
def find_and_display_isomers_simple(molecule_name_input):
    if not molecule_name_input or not molecule_name_input.strip():
        return [], "لطفا نام یک مولکول را وارد کنید.", [] 

    molecule_name = molecule_name_input.strip().lower()
    print(f"Processing request for: '{molecule_name}' (Simple Mode)")
    status_message = ""

    try:
        print(f"Searching for compound: '{molecule_name}' in PubChem...")
        compounds = pcp.get_compounds(molecule_name, 'name', listkey_count=5)
        main_compound_obj = None
        molecular_formula = None

        if not compounds:
            status_message = f"مولکول '{molecule_name}' در PubChem یافت نشد."
            print(status_message)
            return [], status_message, []
        
        print(f"Found {len(compounds)} potential matches for '{molecule_name}'. Checking them...")
        # ... (کد انتخاب main_compound_obj و molecular_formula مانند قبل) ...
        for i, c in enumerate(compounds):
            cid = c.cid
            common_name = c.synonyms[0] if c.synonyms else "N/A"
            actual_formula = c.molecular_formula if hasattr(c, 'molecular_formula') else None
            if actual_formula:
                is_standard_hydrocarbon = True 
                if c.canonical_smiles:
                    try:
                        mol_obj = Chem.MolFromSmiles(c.canonical_smiles)
                        if mol_obj:
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
            status_message = f"آلکان استاندارد با نام '{molecule_name}' یافت نشد یا با معیارها مطابقت ندارد."
            return [], status_message, []

        print(f"Proceeding with main compound: CID {main_compound_obj.cid}, Formula: {molecular_formula}")
        isomers_found_raw = pcp.get_compounds(molecular_formula, 'formula', listkey_count=50)
        if not isomers_found_raw:
            status_message = f"ایزومری برای فرمول {molecular_formula} یافت نشد."
            return [], status_message, []

        valid_structural_alkanes_entries = []
        unique_accepted_smiles = set()
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
        
        processed_isomers_for_gallery = []
        for final_isomer_entry in valid_structural_alkanes_entries:
            smiles_to_draw = final_isomer_entry.canonical_smiles
            iupac_name = final_isomer_entry.iupac_name or (final_isomer_entry.synonyms[0] if final_isomer_entry.synonyms else f"Alkane (CID: {final_isomer_entry.cid})")
            mol_image = draw_molecule(smiles_to_draw)
            if mol_image:
                caption = f"{iupac_name.capitalize()}\nSMILES: {smiles_to_draw}\nCID: {final_isomer_entry.cid}"
                processed_isomers_for_gallery.append({"image": mol_image, "caption": caption, "cid": final_isomer_entry.cid, "sort_key": iupac_name.lower()})
        
        processed_isomers_for_gallery.sort(key=lambda x: x["sort_key"])
        isomer_outputs_final = [(item["image"], item["caption"]) for item in processed_isomers_for_gallery]
        final_cids_ordered = [item["cid"] for item in processed_isomers_for_gallery]

        if not isomer_outputs_final:
            status_message = "ایزومر آلکان استاندارد و قابل رسمی پیدا نشد."
        else:
            status_message = f"{len(isomer_outputs_final)} ایزومر برای '{molecule_name_input}' پیدا شد."
        
        return isomer_outputs_final, status_message, final_cids_ordered

    except Exception as e:
        error_msg = f"خطای کلی در find_and_display_isomers_simple: {e}"
        print(error_msg)
        print(traceback.format_exc())
        return [], error_msg, []


# تابع کنترل‌کننده ساده شده برای انتخاب آیتم از گالری
def handle_gallery_selection_simplified(current_cids_in_gallery_state, evt: gr.SelectData):
    if evt is None or evt.index is None:
        return "هیچ ایزومری انتخاب نشده است" 

    selected_index = evt.index
    if not current_cids_in_gallery_state or selected_index >= len(current_cids_in_gallery_state):
        print(f"خطا در handle_gallery_selection_simplified: ایندکس نامعتبر.")
        return "خطا: ایزومر نامعتبر انتخاب شده است."

    selected_cid = current_cids_in_gallery_state[selected_index]
    print(f"Gallery item selected (Simplified). Index: {selected_index}, CID: {selected_cid}")
    
    # فقط یک رشته ساده حاوی CID را برگردانید
    return f"ایزومر انتخاب شده با CID: {selected_cid}"


# --- بخش Gradio Interface با استفاده از gr.Blocks ---
with gr.Blocks(theme=gr.themes.Soft()) as iface:
    gr.Markdown(
        "# یابنده و نمایشگر ایزومرهای آلکان (تست ساده‌سازی شده گالری)\n"
        "این نسخه برای تست مشکل 'Too many arguments' با خروجی ساده از گالری است."
    )

    # کامپوننت‌های مخفی فعلاً استفاده نمی‌شوند
    # sdf_data_store = gr.Textbox(label="SDF Data", visible=False, elem_id="sdf_data_store_elem")
    # js_trigger = gr.Textbox(label="JS Trigger", visible=False, elem_id="js_trigger_elem")

    with gr.Row():
        molecule_name_input = gr.Textbox(
            label="نام آلکان را وارد کنید",
            placeholder="مثال: butane, pentane",
            info="نام آلکان مورد نظر خود را به انگلیسی و با حروف کوچک وارد کنید.",
            scale=3
        )
        submit_button = gr.Button("جستجوی ایزومرها", variant="primary", scale=1)

    status_output = gr.Textbox(label="وضعیت و پیام‌ها", lines=2, interactive=False)
    isomer_cids_state = gr.State([]) 

    # تکست باکس جدید برای نمایش خروجی ساده گالری
    test_output_textbox = gr.Textbox(label="خروجی تست انتخاب گالری", interactive=False)

    with gr.Row():
        gallery_output = gr.Gallery(
            label="ایزومرهای یافت شده (2D)",
            columns=[3], height="auto", object_fit="contain"
        )
        # html_3d_output فعلاً کامنت شده یا یک مقدار استاتیک دارد
        html_3d_output_placeholder = gr.HTML(
            label="نمایش سه‌بعدی (غیرفعال در این تست)",
            value="<p style='text-align:center; color:grey;'>نمایشگر سه‌بعدی در این حالت تست غیرفعال است.</p>"
        )

    gr.Examples(
        examples=[["butane"], ["pentane"]],
        inputs=molecule_name_input,
    )
    
    # رویداد بارگذاری اولیه صفحه فعلاً کاری انجام نمی‌دهد
    # iface.load(initial_3d_viewer_html_setup_with_listeners, inputs=None, outputs=[html_3d_output])

    # اتصال رویداد کلیک دکمه به تابع ساده شده
    submit_button.click(
        fn=find_and_display_isomers_simple, # تابع ساده شده
        inputs=molecule_name_input,
        outputs=[gallery_output, status_output, isomer_cids_state] # فقط خروجی‌های مربوط به نمایش دوبعدی
    )
    molecule_name_input.submit(
        fn=find_and_display_isomers_simple, # تابع ساده شده
        inputs=molecule_name_input,
        outputs=[gallery_output, status_output, isomer_cids_state] # فقط خروجی‌های مربوط به نمایش دوبعدی
    )

    # اتصال رویداد select گالری به تابع ساده شده و تکست باکس تست
    gallery_output.select(
        fn=handle_gallery_selection_simplified, # تابع ساده شده
        inputs=[isomer_cids_state], 
        outputs=[test_output_textbox] # فقط یک خروجی به تکست باکس تست
    )

if __name__ == '__main__':
    iface.launch(debug=True)
