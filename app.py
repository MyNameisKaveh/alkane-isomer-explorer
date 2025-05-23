import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage
import gradio as gr
import traceback
import py3Dmol
import os

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

# تابع برای ایجاد نمایش سه‌بعدی و برگرداندن HTML
def generate_3d_html_from_cid(cid):
    if cid is None:
        return "<p style='color:orange;'>CID برای نمایش سه‌بعدی ارائه نشده است.</p>"
    
    print(f"در حال ایجاد نمایش سه‌بعدی برای CID: {cid}...")
    temp_sdf_file_dir = "/tmp" 
    if not os.path.exists(temp_sdf_file_dir):
        try:
            os.makedirs(temp_sdf_file_dir, exist_ok=True)
        except OSError as e:
            print(f"خطا در ایجاد دایرکتوری موقت {temp_sdf_file_dir}: {e}")
            temp_sdf_file_dir = "."

    temp_sdf_file = os.path.join(temp_sdf_file_dir, f'temp_3d_structure_{cid}.sdf')
    sdf_content = None

    try:
        print(f"در حال دانلود ساختار سه‌بعدی (SDF) برای CID {cid} به {temp_sdf_file}...")
        pcp.download('SDF', temp_sdf_file, str(cid), 'cid', record_type='3d', overwrite=True)

        with open(temp_sdf_file, 'r') as f:
            sdf_content = f.read()

        if not sdf_content or sdf_content.strip() == "$$$$\n" or sdf_content.strip() == "":
            print(f"فایل SDF دانلود شده برای CID {cid} خالی است یا ساختار سه‌بعدی ندارد.")
            return f"<p style='color:red;'>ساختار سه‌بعدی (SDF) برای CID {cid} در PubChem یافت نشد یا خالی است.</p>"

    except pcp.NotFoundError:
        print(f"ساختار سه‌بعدی (SDF) برای CID {cid} در PubChem یافت نشد.")
        return f"<p style='color:red;'>ساختار سه‌بعدی (SDF) برای CID {cid} در PubChem یافت نشد.</p>"
    except FileNotFoundError:
        print(f"فایل موقت SDF برای CID {cid} در مسیر {temp_sdf_file} پس از تلاش برای دانلود یافت نشد.")
        return f"<p style='color:red;'>خطای داخلی: فایل ساختار سه‌بعدی ایجاد نشد.</p>"
    except Exception as e:
        error_msg = f"خطا در دانلود یا خواندن فایل SDF برای CID {cid}: {e}"
        print(error_msg)
        print(f"FULL TRACEBACK for SDF download/read: {traceback.format_exc()}")
        return f"<p style='color:red;'>{error_msg}</p>"
    finally:
        if os.path.exists(temp_sdf_file):
            try:
                os.remove(temp_sdf_file)
            except Exception as e_rem:
                print(f"خطا در حذف فایل موقت {temp_sdf_file}: {e_rem}")

    if sdf_content:
        try:
            print("در حال نمایش ساختار سه‌بعدی با py3Dmol...")
            viewer = py3Dmol.view(width=500, height=400) # می‌توانید ابعاد را تغییر دهید
            viewer.addModel(sdf_content, 'sdf')
            viewer.setStyle({'stick': {}})
            viewer.setBackgroundColor('0xeeeeee')
            viewer.zoomTo()
            html_output = viewer._make_html()
            
            # لاگ کردن HTML خروجی برای بررسی
            print("---- 3DMOL HTML OUTPUT ----")
            # فقط بخشی از HTML را پرینت کنید اگر خیلی طولانی است، یا بخش‌های کلیدی
            # print(html_output[:1000] + "..." if len(html_output) > 1000 else html_output)
            print(html_output) # برای بررسی کامل در لاگ
            print("--------------------------")
            
            return html_output
        except Exception as e_render:
            error_msg = f"خطا در رندر کردن ساختار سه‌بعدی برای CID {cid} با py3Dmol: {e_render}"
            print(error_msg)
            print(f"FULL TRACEBACK for py3Dmol render: {traceback.format_exc()}")
            return f"<p style='color:red;'>{error_msg}</p>"
    else:
        return f"<p style='color:red;'>محتوای SDF برای CID {cid} پس از دانلود در دسترس نبود.</p>"

# تابع اصلی پردازش آلکان (بدون تغییر عمده نسبت به نسخه قبلی که کار می‌کرد)
def find_and_display_isomers_and_cids(molecule_name_input):
    if not molecule_name_input or not molecule_name_input.strip():
        return [], "لطفا نام یک مولکول را وارد کنید.", [], "<p style='color:orange;'>نام مولکول وارد نشده است.</p>"

    molecule_name = molecule_name_input.strip().lower()
    print(f"Processing request for: '{molecule_name}'")
    # ... (بقیه کد این تابع مانند نسخه قبلی که دوبعدی آن به درستی کار می‌کرد) ...
    # ... (بخش جستجوی ترکیب اصلی، فیلتر کردن ایزومرها، آماده‌سازی گالری) ...
    # این بخش برای اختصار حذف شده، اما باید همان کد قبلی شما باشد که دوبعدی آن کار می‌کرد.
    # فرض می‌کنیم این تابع به درستی مقادیر زیر را برمی‌گرداند:
    # isomer_outputs_final, status_message, final_cids_ordered, initial_3d_view_html
    # در اینجا یک پیاده‌سازی خلاصه شده برای تست قرار می‌دهم، شما باید کد کامل خود را جایگزین کنید
    
    # --- آغاز بخش خلاصه شده برای find_and_display_isomers_and_cids ---
    # --- شما باید کد کامل و صحیح خود را در اینجا قرار دهید ---
    try:
        # شبیه‌سازی جستجو و یافتن ایزومرها
        print(f"Searching for compound: '{molecule_name}' in PubChem...")
        # ... (منطق جستجوی main_compound_obj و molecular_formula) ...
        molecular_formula = "C4H10" # مثال برای بوتان
        if molecule_name == "butane":
            main_compound_obj_cid = 7843 # n-butane CID
        elif molecule_name == "isobutane": # مثال
             main_compound_obj_cid = 6360 # isobutane CID
        else: # یک پیش‌فرض
            main_compound_obj_cid = 7843
        
        print(f"Proceeding with main compound: CID {main_compound_obj_cid}, Formula: {molecular_formula}")
        print(f"Searching for isomers with formula: {molecular_formula}...")
        
        # شبیه‌سازی یافتن ایزومرها
        # برای بوتان (C4H10)، دو ایزومر داریم: n-butane (CID: 7843) و isobutane (2-methylpropane, CID: 6360)
        isomer_entries_simulated = []
        if molecular_formula == "C4H10":
            # n-butane
            butane_entry = pcp.Compound.from_cid(7843) # شیء واقعی PubChemPy
            # isobutane
            isobutane_entry = pcp.Compound.from_cid(6360) # شیء واقعی PubChemPy
            isomer_entries_simulated = [butane_entry, isobutane_entry]


        if not isomer_entries_simulated:
            status_message = f"ایزومری برای فرمول {molecular_formula} یافت نشد (شبیه‌سازی شده)."
            return [], status_message, [], f"<p style='color:orange;'>{status_message}</p>"

        processed_isomers_for_gallery = []
        for entry in isomer_entries_simulated:
            img = draw_molecule(entry.canonical_smiles)
            if img:
                name = entry.iupac_name if entry.iupac_name else entry.synonyms[0]
                caption = f"{name.capitalize()}\nSMILES: {entry.canonical_smiles}\nCID: {entry.cid}"
                processed_isomers_for_gallery.append({
                    "image": img, "caption": caption, "cid": entry.cid, "sort_key": name.lower()
                })
        
        processed_isomers_for_gallery.sort(key=lambda x: x["sort_key"])
        isomer_outputs_final = [(item["image"], item["caption"]) for item in processed_isomers_for_gallery]
        final_cids_ordered = [item["cid"] for item in processed_isomers_for_gallery]
        
        if not isomer_outputs_final:
            status_message = "ایزومر آلکان استاندارد و قابل رسمی پیدا نشد (شبیه‌سازی شده)."
            initial_3d_view_html = f"<p style='color:orange;'>{status_message}</p>"
        else:
            status_message = f"{len(isomer_outputs_final)} ایزومر ساختاری آلکان برای '{molecule_name_input}' (فرمول: {molecular_formula}) پیدا و نمایش داده شد (شبیه‌سازی شده)."
            initial_3d_view_html = "<p>یک ایزومر را از گالری بالا برای نمایش سه‌بعدی انتخاب کنید.</p>"
        
        return isomer_outputs_final, status_message, final_cids_ordered, initial_3d_view_html

    except Exception as e:
        error_msg = f"خطا در شبیه‌سازی: {e}"
        print(error_msg)
        return [], error_msg, [], f"<p style='color:red;'>{error_msg}</p>"
    # --- پایان بخش خلاصه شده ---
    # --- شما باید کد کامل و صحیح خود را در اینجا قرار دهید ---


# تابع کنترل‌کننده برای انتخاب آیتم از گالری (با تغییر برای تست HTML ساده)
def handle_gallery_selection(current_cids_in_gallery_state, evt: gr.SelectData):
    if evt is None or evt.index is None: # اگر چیزی انتخاب نشده باشد
        return "<p>برای نمایش سه‌بعدی، یک ایزومر را از گالری انتخاب کنید.</p>"
    
    selected_index = evt.index
    if not current_cids_in_gallery_state or selected_index >= len(current_cids_in_gallery_state):
        print(f"خطا: ایندکس {selected_index} خارج از محدوده لیست CIDها ({len(current_cids_in_gallery_state)}) است یا لیست CID خالی است.")
        return "<p style='color:red;'>خطا: اطلاعات CID برای ایزومر انتخاب شده یافت نشد.</p>"
        
    selected_cid = current_cids_in_gallery_state[selected_index]
    print(f"Gallery item selected. Index: {selected_index}, CID: {selected_cid}")
    
    # ***** شروع بخش تست با HTML ساده *****
    # ابتدا این بخش را برای تست فعال نگه دارید.
    # اگر با این HTML ساده، خطای "Too many arguments" در کنسول مرورگر ظاهر نشد،
    # این بخش را کامنت کرده و خط بعدی (generate_3d_html_from_cid) را از کامنت خارج کنید.
    simple_html_test = f"""
    <div style='padding: 20px; border: 2px solid blue; background-color: lightblue;'>
        <h1>تست نمایش سه‌بعدی</h1>
        <p>ایزومر انتخاب شده با CID: <strong>{selected_cid}</strong></p>
        <p>این یک HTML ساده برای بررسی عملکرد به‌روزرسانی کامپوننت است.</p>
        <p>اگر این پیام را بدون خطای "Too many arguments" در کنسول می‌بینید، مشکل از خود HTML پیچیده 3Dmol است.</p>
    </div>
    """
    print(f"Returning simple HTML for CID {selected_cid}")
    return simple_html_test
    # ***** پایان بخش تست با HTML ساده *****

    # اگر تست بالا موفق بود، این خط را از کامنت خارج کنید و بخش تست ساده بالا را کامنت کنید:
    # print(f"Calling generate_3d_html_from_cid for CID {selected_cid}")
    # return generate_3d_html_from_cid(selected_cid)


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
            scale=3
        )
        submit_button = gr.Button("جستجوی ایزومرها", variant="primary", scale=1)

    status_output = gr.Textbox(label="وضعیت و پیام‌ها", lines=2, interactive=False)
    
    isomer_cids_state = gr.State([])

    with gr.Row():
        gallery_output = gr.Gallery(
            label="ایزومرهای یافت شده (2D)",
            columns=[3],
            height="auto",
            object_fit="contain",
        )
        html_3d_output = gr.HTML(
            label="نمایش سه‌بعدی ایزومر (py3Dmol)"
        )
    
    gr.Examples(
        examples=[["butane"], ["pentane"], ["hexane"], ["heptane"]],
        inputs=molecule_name_input,
        # outputs=[gallery_output, status_output, isomer_cids_state, html_3d_output], # برای فعال‌سازی، تابع fn باید همه اینها را برگرداند
        # fn=find_and_display_isomers_and_cids
    )
    
    # اتصال رویدادها
    submit_event_params = {
        "fn": find_and_display_isomers_and_cids,
        "inputs": molecule_name_input,
        "outputs": [gallery_output, status_output, isomer_cids_state, html_3d_output]
    }
    submit_button.click(**submit_event_params)
    molecule_name_input.submit(**submit_event_params)

    gallery_output.select(
        fn=handle_gallery_selection,
        inputs=[isomer_cids_state], # evt به طور خودکار پاس داده می‌شود
        outputs=html_3d_output
    )

if __name__ == '__main__':
    iface.launch()
