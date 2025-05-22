import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage
import gradio as gr
import traceback
import py3Dmol # اضافه شده
import os # اضافه شده

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

# تابع جدید برای ایجاد نمایش سه‌بعدی و برگرداندن HTML
def generate_3d_html_from_cid(cid):
    if cid is None:
        return "<p style='color:orange;'>CID برای نمایش سه‌بعدی ارائه نشده است.</p>"
    
    print(f"در حال ایجاد نمایش سه‌بعدی برای CID: {cid}...")
    temp_sdf_file = f'temp_3d_structure_{cid}.sdf' # نام فایل موقت با CID
    sdf_content = None

    try:
        # دانلود فایل SDF سه‌بعدی
        pcp.download('SDF', temp_sdf_file, str(cid), 'cid', record_type='3d', overwrite=True)

        # خواندن محتوای فایل SDF
        with open(temp_sdf_file, 'r') as f:
            sdf_content = f.read()

        if not sdf_content or sdf_content.strip() == "$$$$\n": # بررسی خالی بودن فایل SDF
            print(f"فایل SDF دانلود شده برای CID {cid} خالی است یا ساختار سه‌بعدی ندارد.")
            return f"<p style='color:red;'>ساختار سه‌بعدی (SDF) برای CID {cid} در PubChem یافت نشد یا خالی است.</p>"

    except pcp.NotFoundError:
        print(f"ساختار سه‌بعدی (SDF) برای CID {cid} در PubChem یافت نشد.")
        return f"<p style='color:red;'>ساختار سه‌بعدی (SDF) برای CID {cid} در PubChem یافت نشد.</p>"
    except Exception as e:
        error_msg = f"خطا در دانلود یا خواندن فایل SDF برای CID {cid}: {e}"
        print(error_msg)
        print(f"FULL TRACEBACK for SDF download/read: {traceback.format_exc()}")
        return f"<p style='color:red;'>{error_msg}</p>"
    finally:
        # پاک کردن فایل موقت پس از استفاده
        if os.path.exists(temp_sdf_file):
            try:
                os.remove(temp_sdf_file)
            except Exception as e_rem:
                print(f"خطا در حذف فایل موقت {temp_sdf_file}: {e_rem}")

    if sdf_content:
        try:
            print("در حال نمایش ساختار سه‌بعدی با py3Dmol...")
            viewer = py3Dmol.view(width=500, height=400)
            viewer.addModel(sdf_content, 'sdf')
            viewer.setStyle({'stick': {}}) # سبک نمایش
            viewer.setBackgroundColor('0xeeeeee')
            viewer.zoomTo()
            # _make_html یک صفحه کامل HTML برمی‌گرداند، Gradio آن را در iframe قرار می‌دهد.
            html_output = viewer._make_html()
            return html_output
        except Exception as e_render:
            error_msg = f"خطا در رندر کردن ساختار سه‌بعدی برای CID {cid} با py3Dmol: {e_render}"
            print(error_msg)
            print(f"FULL TRACEBACK for py3Dmol render: {traceback.format_exc()}")
            return f"<p style='color:red;'>{error_msg}</p>"
    else:
        # این حالت نباید رخ دهد اگر بررسی‌های قبلی انجام شده باشند
        return f"<p style='color:red;'>محتوای SDF برای CID {cid} پس از دانلود در دسترس نبود (وضعیت غیرمنتظره).</p>"

# تابع اصلی پردازش آلکان (قبلاً find_and_display_isomers)
# اکنون لیست CIDها و HTML اولیه برای نمایشگر سه‌بعدی را نیز برمی‌گرداند.
def find_and_display_isomers_and_cids(molecule_name_input):
    if not molecule_name_input or not molecule_name_input.strip():
        return [], "لطفا نام یک مولکول را وارد کنید.", [], "<p style='color:orange;'>نام مولکول وارد نشده است.</p>"

    molecule_name = molecule_name_input.strip().lower()
    print(f"Processing request for: '{molecule_name}'")
    status_message = ""

    try:
        print(f"Searching for compound: '{molecule_name}' in PubChem...")
        compounds = pcp.get_compounds(molecule_name, 'name', listkey_count=5)
        
        main_compound_obj = None
        molecular_formula = None

        if not compounds:
            status_message = f"مولکول '{molecule_name}' در PubChem یافت نشد."
            print(status_message)
            return [], status_message, [], f"<p style='color:red;'>{status_message}</p>"
        
        print(f"Found {len(compounds)} potential matches for '{molecule_name}'. Checking them...")
        # ... (بخش انتخاب main_compound_obj و molecular_formula مانند کد اصلی شما) ...
        # این بخش از کد شما برای انتخاب ترکیب اصلی باید بدون تغییر باقی بماند
        # فقط مطمئن شوید که main_compound_obj و molecular_formula به درستی مقداردهی می‌شوند.
        # کپی کردن بخش مربوط به انتخاب main_compound_obj از کد اصلی شما:
        for i, c in enumerate(compounds):
            cid = c.cid
            common_name = c.synonyms[0] if c.synonyms else "N/A"
            actual_formula = c.molecular_formula if hasattr(c, 'molecular_formula') else None
            
            print(f"  Checking main compound candidate {i+1}: CID {cid}, Name: '{common_name}', Formula: '{actual_formula}'")

            if actual_formula:
                is_standard_hydrocarbon = True 
                if c.canonical_smiles:
                    try:
                        mol_obj = Chem.MolFromSmiles(c.canonical_smiles)
                        if mol_obj:
                            if len(Chem.GetMolFrags(mol_obj)) > 1:
                                print(f"    Main candidate CID {cid} is disconnected.")
                                is_standard_hydrocarbon = False
                            if is_standard_hydrocarbon:
                                atom_symbols_main = set()
                                for atom in mol_obj.GetAtoms():
                                    atom_symbols_main.add(atom.GetSymbol())
                                    if atom.GetIsotope() != 0:
                                        is_standard_hydrocarbon = False; break
                                if not atom_symbols_main.issubset({'C', 'H'}):
                                    is_standard_hydrocarbon = False
                            if is_standard_hydrocarbon:
                                for bond in mol_obj.GetBonds():
                                    if bond.GetBondType() != Chem.BondType.SINGLE:
                                        is_standard_hydrocarbon = False; break
                                if Chem.GetSymmSSSR(mol_obj):
                                    is_standard_hydrocarbon = False
                        else: is_standard_hydrocarbon = False
                    except Exception: is_standard_hydrocarbon = False 
                else: is_standard_hydrocarbon = False 

                if not is_standard_hydrocarbon: continue

                current_compound_name_matches_input = molecule_name in [syn.lower() for syn in c.synonyms]
                if current_compound_name_matches_input: 
                    main_compound_obj = c
                    molecular_formula = actual_formula
                    break 
                if not main_compound_obj: 
                    main_compound_obj = c 
                    molecular_formula = actual_formula 
        
        if not main_compound_obj or not molecular_formula: 
            status_message = f"آلکان استاندارد با نام '{molecule_name}' در PubChem یافت نشد یا با معیارهای آلکان (فقط C و H، پیوندهای یگانه، بدون حلقه، بدون ایزوتوپ غیر استاندارد، متصل) مطابقت ندارد."
            print(status_message)
            return [], status_message, [], f"<p style='color:red;'>{status_message}</p>"
        
        print(f"Proceeding with main compound: CID {main_compound_obj.cid}, Formula: {molecular_formula}")
        print(f"Searching for isomers with formula: {molecular_formula}...")
        isomers_found_raw = pcp.get_compounds(molecular_formula, 'formula', listkey_count=50)

        if not isomers_found_raw:
            status_message = f"ایزومری برای فرمول {molecular_formula} یافت نشد."
            print(status_message)
            return [], status_message, [], f"<p style='color:orange;'>{status_message}</p>"

        print(f"Found {len(isomers_found_raw)} potential isomer entries from PubChem. Filtering for true structural alkane isomers...")
        
        valid_structural_alkanes_entries = [] 
        unique_accepted_smiles = set()

        # ... (بخش فیلتر کردن ایزومرها مانند کد اصلی شما) ...
        # این بخش از کد شما برای فیلتر کردن ایزومرها باید بدون تغییر باقی بماند.
        # کپی کردن بخش فیلتر کردن ایزومرها از کد اصلی شما:
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
                        valid_structural_alkanes_entries.append(isomer_entry)
                        unique_accepted_smiles.add(canonical_smiles_for_uniqueness)
            except Exception: continue
        
        print(f"Found {len(valid_structural_alkanes_entries)} unique, valid structural alkane isomers after filtering.")
        
        processed_isomers_for_gallery = []

        for final_isomer_entry in valid_structural_alkanes_entries:
            smiles_to_draw = final_isomer_entry.canonical_smiles
            iupac_name = final_isomer_entry.iupac_name
            
            if not iupac_name and final_isomer_entry.synonyms:
                chosen_synonym = final_isomer_entry.synonyms[0]
                # (منطق انتخاب نام شما در اینجا)
                simple_names = [s for s in final_isomer_entry.synonyms if s.lower().endswith("ane") and not any(char.isdigit() for char in s.split('-')[0]) and '-' not in s.split(' ')[0]]
                if simple_names:
                    chosen_synonym = min(simple_names, key=len)
                else:
                    non_iupac_synonyms = [s for s in final_isomer_entry.synonyms if s != final_isomer_entry.iupac_name]
                    if non_iupac_synonyms:
                        chosen_synonym = min(non_iupac_synonyms, key=len)
                iupac_name = chosen_synonym.capitalize()
            elif not iupac_name:
                iupac_name = f"Alkane (CID: {final_isomer_entry.cid})"

            mol_image = draw_molecule(smiles_to_draw)
            if mol_image:
                # اضافه کردن CID به کپشن برای اطلاعات بیشتر و استفاده در مرتب‌سازی اگر لازم باشد
                caption = f"{iupac_name}\nSMILES: {smiles_to_draw}\nCID: {final_isomer_entry.cid}"
                processed_isomers_for_gallery.append({
                    "image": mol_image,
                    "caption": caption,
                    "cid": final_isomer_entry.cid,
                    "sort_key": iupac_name.lower() # برای مرتب‌سازی بر اساس نام
                })
            else:
                print(f"  Failed to draw image for accepted isomer: CID {final_isomer_entry.cid}, SMILES: {smiles_to_draw}")

        # مرتب‌سازی ایزومرها بر اساس نام (یا کلید مرتب‌سازی دیگر)
        processed_isomers_for_gallery.sort(key=lambda x: x["sort_key"])

        # آماده‌سازی خروجی‌ها برای Gradio
        isomer_outputs_final = [(item["image"], item["caption"]) for item in processed_isomers_for_gallery]
        final_cids_ordered = [item["cid"] for item in processed_isomers_for_gallery]
        
        valid_isomers_count_final = len(isomer_outputs_final)
        print(f"Displayed {valid_isomers_count_final} isomers in the gallery.")
        print(f"CIDs for gallery (ordered): {final_cids_ordered}")

        initial_3d_view_html = ""
        if not isomer_outputs_final:
            status_message = "ایزومر آلکان استاندارد و قابل رسمی پیدا نشد."
            if len(valid_structural_alkanes_entries) > 0 and not processed_isomers_for_gallery:
                 status_message += " (برخی در مرحله رسم ناموفق بودند یا کاندیدای معتبری نبودند)."
            initial_3d_view_html = f"<p style='color:orange;'>{status_message}</p>"
        else:
            status_message = f"{len(isomer_outputs_final)} ایزومر ساختاری آلکان برای '{molecule_name_input}' (فرمول: {molecular_formula}) پیدا و نمایش داده شد."
            initial_3d_view_html = "<p>یک ایزومر را از گالری بالا برای نمایش سه‌بعدی انتخاب کنید.</p>"
            # اختیاری: نمایش خودکار اولین ایزومر در حالت سه‌بعدی
            # if final_cids_ordered:
            #     initial_3d_view_html = generate_3d_html_from_cid(final_cids_ordered[0])
        
        return isomer_outputs_final, status_message, final_cids_ordered, initial_3d_view_html

    except pcp.PubChemHTTPError as e:
        error_msg = f"خطا در ارتباط با PubChem: {e}."
        print(error_msg)
        print(f"FULL TRACEBACK for PubChemHTTPError: {traceback.format_exc()}")
        return [], error_msg, [], f"<p style='color:red;'>{error_msg}</p>"
    except Exception as e:
        error_msg = f"یک خطای غیرمنتظره در سرور رخ داد: {e}"
        print(f"FULL TRACEBACK for general Exception: {traceback.format_exc()}")
        return [], error_msg, [], f"<p style='color:red;'>{error_msg}</p>"

# تابع کنترل‌کننده برای انتخاب آیتم از گالری
def handle_gallery_selection(current_cids_in_gallery_state, evt: gr.SelectData):
    if evt is None or evt.index is None : # اگر چیزی انتخاب نشده باشد یا رویداد نامعتبر باشد
        return "<p>برای نمایش سه‌بعدی، یک ایزومر را از گالری انتخاب کنید.</p>"
    
    selected_index = evt.index
    if not current_cids_in_gallery_state or selected_index >= len(current_cids_in_gallery_state):
        print(f"خطا: ایندکس {selected_index} خارج از محدوده لیست CIDها است یا لیست CID خالی است.")
        return "<p style='color:red;'>خطا: اطلاعات CID برای ایزومر انتخاب شده یافت نشد.</p>"
        
    selected_cid = current_cids_in_gallery_state[selected_index]
    print(f"Gallery item selected. Index: {selected_index}, CID: {selected_cid}")
    return generate_3d_html_from_cid(selected_cid)

# --- بخش Gradio Interface با استفاده از gr.Blocks ---
with gr.Blocks(theme=gr.themes.Soft()) as iface:
    gr.Markdown(
        "# یابنده و نمایشگر ایزومرهای آلکان (دو بعدی و سه بعدی)\n"
        "نام یک آلکان (به انگلیسی) را وارد کنید تا ایزومرهای آن به همراه ساختار شیمیایی دوبعدی، نام، SMILES و CID نمایش داده شوند.\n"
        "با کلیک بر روی هر ایزومر در گالری، ساختار سه‌بعدی آن (در صورت وجود در PubChem) نمایش داده خواهد شد.\n"
        "اطلاعات از دیتابیس PubChem دریافت شده و ساختارها با استفاده از کتابخانه‌های RDKit (2D) و py3Dmol (3D) رسم می‌شوند."
    )

    with gr.Row():
        molecule_name_input = gr.Textbox(
            label="نام آلکان را وارد کنید",
            placeholder="مثال: butane, pentane, hexane",
            info="نام آلکان مورد نظر خود را به انگلیسی و با حروف کوچک وارد کنید.",
            scale=3 # تخصیص فضای بیشتر به تکست‌باکس
        )
        submit_button = gr.Button("جستجوی ایزومرها", variant="primary", scale=1)

    status_output = gr.Textbox(label="وضعیت و پیام‌ها", lines=2)
    
    # State برای نگهداری CID های ایزومرهای نمایش داده شده در گالری
    # این State با هر جستجوی جدید به‌روز می‌شود.
    isomer_cids_state = gr.State([])

    with gr.Row():
        gallery_output = gr.Gallery(
            label="ایزومرهای یافت شده (2D)",
            columns=[3],
            height="auto",
            object_fit="contain",
            # elem_id="isomer_gallery" # اگر نیاز به ارجاع با JS باشد
            scale=1 # تنظیم مقیاس برای تقسیم فضا
        )
        html_3d_output = gr.HTML(
            label="نمایش سه‌بعدی ایزومر (py3Dmol)",
            scale=1 # تنظیم مقیاس برای تقسیم فضا
        )
    
    gr.Examples(
        examples=[["butane"], ["pentane"], ["hexane"], ["heptane"]],
        inputs=molecule_name_input,
        # outputs=[gallery_output, status_output, isomer_cids_state, html_3d_output], # برای Examples، تابع باید بتواند همه خروجی‌ها را تولید کند
        # fn=find_and_display_isomers_and_cids # تابع examples باید بتواند همه خروجی‌ها را برگرداند
    )
    
    # اتصال رویداد کلیک دکمه یا submit تکست‌باکس به تابع اصلی
    # تابع اصلی باید تمام خروجی‌های متصل را برگرداند
    submit_button.click(
        fn=find_and_display_isomers_and_cids,
        inputs=molecule_name_input,
        outputs=[gallery_output, status_output, isomer_cids_state, html_3d_output]
    )
    molecule_name_input.submit( # همچنین با زدن Enter در تکست‌باکس
        fn=find_and_display_isomers_and_cids,
        inputs=molecule_name_input,
        outputs=[gallery_output, status_output, isomer_cids_state, html_3d_output]
    )

    # اتصال رویداد select گالری به تابع کنترل‌کننده برای نمایش سه‌بعدی
    # evt (gr.SelectData) به صورت خودکار به عنوان آخرین آرگومان به تابع پاس داده می‌شود.
    gallery_output.select(
        fn=handle_gallery_selection,
        inputs=[isomer_cids_state], # ورودی صریح از state
        outputs=html_3d_output
    )

if __name__ == '__main__':
    iface.launch()
