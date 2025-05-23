import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage
import gradio as gr
import traceback
import py3Dmol # py3Dmol دیگر برای تولید HTML کامل در هر بار استفاده نمی‌شود، اما ممکن است برای برخی عملیات لازم باشد
import os
import time # برای ایجاد ماشه منحصر به فرد

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

# تابع جدید برای دریافت فقط محتوای SDF
def get_sdf_content(cid):
    if cid is None:
        print("CID برای دریافت محتوای SDF ارائه نشده است.")
        return None, "CID ارائه نشده است."
    
    print(f"در حال دریافت محتوای SDF برای CID: {cid}...")
    temp_sdf_file_dir = "/tmp" 
    if not os.path.exists(temp_sdf_file_dir):
        try:
            os.makedirs(temp_sdf_file_dir, exist_ok=True)
        except OSError as e_mkdir:
            print(f"خطا در ایجاد دایرکتوری موقت {temp_sdf_file_dir}: {e_mkdir}")
            temp_sdf_file_dir = "."

    temp_sdf_file = os.path.join(temp_sdf_file_dir, f'temp_sdf_for_js_{cid}.sdf')
    sdf_content = None
    error_message = None

    try:
        print(f"در حال دانلود ساختار سه‌بعدی (SDF) برای CID {cid} به {temp_sdf_file}...")
        pcp.download('SDF', temp_sdf_file, str(cid), 'cid', record_type='3d', overwrite=True)

        with open(temp_sdf_file, 'r') as f:
            sdf_content = f.read()

        if not sdf_content or sdf_content.strip() == "$$$$\n" or sdf_content.strip() == "":
            msg = f"فایل SDF دانلود شده برای CID {cid} خالی است یا ساختار سه‌بعدی ندارد."
            print(msg)
            return None, msg

    except pcp.NotFoundError:
        msg = f"ساختار سه‌بعدی (SDF) برای CID {cid} در PubChem یافت نشد."
        print(msg)
        return None, msg
    except FileNotFoundError:
        msg = f"فایل موقت SDF برای CID {cid} در مسیر {temp_sdf_file} پس از تلاش برای دانلود یافت نشد."
        print(msg)
        return None, msg
    except Exception as e:
        error_msg = f"خطا در دانلود یا خواندن فایل SDF برای CID {cid}: {e}"
        print(error_msg)
        print(f"FULL TRACEBACK for SDF download/read: {traceback.format_exc()}")
        return None, error_msg
    finally:
        if os.path.exists(temp_sdf_file):
            try:
                os.remove(temp_sdf_file)
            except Exception as e_rem:
                print(f"خطا در حذف فایل موقت {temp_sdf_file}: {e_rem}")
    
    if sdf_content:
        print(f"محتوای SDF برای CID {cid} با موفقیت دریافت شد.")
        return sdf_content, None # برگرداندن محتوای SDF و بدون خطا
    else:
        # این حالت اگر بررسی‌های قبلی کار کنند، نباید رخ دهد
        return None, error_message if error_message else "محتوای SDF پس از دانلود در دسترس نبود."


# تابع اصلی پردازش آلکان (باید کد کامل شما باشد)
def find_and_display_isomers_and_cids(molecule_name_input):
    # ... (کد کامل این تابع از پاسخ قبلی که بخش 2D آن کار می‌کرد کپی شود)
    # این تابع باید مقادیر زیر را برگرداند:
    # isomer_outputs_final, status_message, final_cids_ordered, initial_js_trigger_val
    # initial_js_trigger_val می‌تواند یک رشته خالی یا null اولیه باشد
    if not molecule_name_input or not molecule_name_input.strip():
        return [], "لطفا نام یک مولکول را وارد کنید.", [], "", "" # گالری، وضعیت، CIDs، داده SDF، ماشه JS

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
            return [], status_message, [], "", ""
        
        print(f"Found {len(compounds)} potential matches for '{molecule_name}'. Checking them...")
        # ... (کد انتخاب main_compound_obj و molecular_formula مانند قبل) ...
        for i, c in enumerate(compounds):
            cid = c.cid
            common_name = c.synonyms[0] if c.synonyms else "N/A"
            actual_formula = c.molecular_formula if hasattr(c, 'molecular_formula') else None
            if actual_formula: # ادامه منطق فیلتر کردن ترکیب اصلی
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
            return [], status_message, [], "", ""

        print(f"Proceeding with main compound: CID {main_compound_obj.cid}, Formula: {molecular_formula}")
        isomers_found_raw = pcp.get_compounds(molecular_formula, 'formula', listkey_count=50)
        if not isomers_found_raw:
            status_message = f"ایزومری برای فرمول {molecular_formula} یافت نشد."
            return [], status_message, [], "", ""

        # ... (کد فیلتر کردن ایزومرها و آماده‌سازی گالری مانند قبل) ...
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
        
        # برای بار اول، SDF و ماشه خالی هستند
        return isomer_outputs_final, status_message, final_cids_ordered, "", ""

    except Exception as e:
        error_msg = f"خطای کلی در find_and_display_isomers: {e}"
        print(error_msg)
        print(traceback.format_exc())
        return [], error_msg, [], "", ""


# تابع کنترل‌کننده جدید برای انتخاب آیتم از گالری
# این تابع دیگر HTML برنمی‌گرداند، بلکه داده SDF و یک ماشه را برمی‌گرداند
def handle_gallery_selection_for_sdf(current_cids_in_gallery_state, evt: gr.SelectData):
    if evt is None or evt.index is None:
        # اگر چیزی انتخاب نشده، داده SDF خالی و ماشه بدون تغییر (یا یک مقدار خاص) برگردان
        return "", gr.Textbox.NO_CHANGE # یا None برای ماشه اگر بهتر است

    selected_index = evt.index
    if not current_cids_in_gallery_state or selected_index >= len(current_cids_in_gallery_state):
        print(f"خطا در handle_gallery_selection_for_sdf: ایندکس نامعتبر.")
        return "خطا: ایزومر نامعتبر.", gr.Textbox.NO_CHANGE # یا None

    selected_cid = current_cids_in_gallery_state[selected_index]
    print(f"Gallery item selected for SDF. Index: {selected_index}, CID: {selected_cid}")
    
    sdf_data, error_msg = get_sdf_content(selected_cid)
    
    if error_msg:
        print(f"خطا در دریافت SDF برای CID {selected_cid}: {error_msg}")
        # در صورت خطا، یک پیام خطا به کاربر (مثلاً از طریق status_output دیگر) و ماشه برای پاک کردن نمایشگر
        # یا می‌توان SDF خالی برگرداند تا JS آن را مدیریت کند
        return f"خطا در بارگذاری سه‌بعدی: {error_msg}", str(time.time()) + "_error"


    if sdf_data:
        print(f"SDF data for CID {selected_cid} prepared for JS update.")
        # یک مقدار جدید و منحصر به فرد برای ماشه ایجاد کنید (مثلاً timestamp)
        # تا listener جاوااسکریپت متوجه تغییر شود.
        trigger_value = str(selected_cid) + "_" + str(time.time())
        return sdf_data, trigger_value
    else:
        print(f"SDF data for CID {selected_cid} is None, though no explicit error was set.")
        return "خطا: داده سه‌بعدی یافت نشد.", str(time.time()) + "_nodata"


# تابع برای بارگذاری HTML اولیه نمایشگر سه‌بعدی
def initial_3d_viewer_html_setup():
    # این HTML یک نمایشگر 3Dmol خالی ایجاد می‌کند و تابع window.update3DModel را تعریف می‌کند
    # برای جلوگیری از خطای CORS یا مشکلات دیگر، می‌توانیم 3Dmol.js را هم به صورت رشته در بیاوریم
    # اما استفاده از CDN ساده‌تر است اگر کار کند.
    html_content = """
    <div id="viewer_3dmol_container" style="width: 500px; height: 400px; border: 1px solid #ccc; position: relative;">
        <p id="loading_3dmol_message" style="text-align: center; padding-top: 50px;">در حال بارگذاری نمایشگر سه‌بعدی...</p>
    </div>

    <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.2/3Dmol-min.js"></script>
    <script>
        var global_3d_viewer = null; // نمایشگر سراسری

        // تابعی برای اطمینان از اجرای کد پس از بارگذاری کامل 3Dmol.js
        function init3DViewerWhenReady() {
            if (typeof $3Dmol !== 'undefined' && $3Dmol.createViewer) {
                try {
                    var container = document.getElementById('viewer_3dmol_container');
                    if (!container) {
                        console.error('Container viewer_3dmol_container not found!');
                        return;
                    }
                    // حذف پیام "در حال بارگذاری"
                    var loadingMsg = document.getElementById('loading_3dmol_message');
                    if (loadingMsg) loadingMsg.style.display = 'none';
                    
                    // ایجاد نمایشگر فقط اگر وجود نداشته باشد
                    if (!global_3d_viewer) {
                       global_3d_viewer = $3Dmol.createViewer(container, { backgroundColor: 'white' }); // یا $('div#viewer_3dmol_container')
                       console.log("3Dmol Viewer initialized.");
                    } else {
                       console.log("3Dmol Viewer already initialized.");
                    }
                    global_3d_viewer.render(); // رندر اولیه (خالی)
                } catch (e) {
                    console.error("Error initializing 3Dmol viewer:", e);
                    container.innerHTML = "<p style='color:red'>خطا در مقداردهی اولیه نمایشگر 3Dmol: " + e.message + "</p>";
                }
            } else {
                console.log("3Dmol.js not ready yet, retrying in 100ms...");
                setTimeout(init3DViewerWhenReady, 100); // تلاش مجدد
            }
        }

        // تابع سراسری برای به‌روزرسانی مدل
        window.update3DModel = function(sdfData, modelId) {
            console.log("window.update3DModel called with modelId:", modelId);
            if (!global_3d_viewer) {
                console.error("Global 3D viewer is not initialized. Attempting to initialize now.");
                init3DViewerWhenReady(); // سعی کن دوباره مقداردهی اولیه کنی
                // پس از یک تاخیر کوتاه دوباره تابع به‌روزرسانی را فراخوانی کن
                setTimeout(function() { window.update3DModel(sdfData, modelId); }, 500);
                return;
            }
            if (sdfData && sdfData.trim() !== "" && !sdfData.startsWith("خطا:")) {
                try {
                    console.log("Clearing previous models from viewer.");
                    global_3d_viewer.removeAllModels(); // یا clear()
                    console.log("Adding new model from SDF data (first 100 chars):", sdfData.substring(0,100));
                    global_3d_viewer.addModel(sdfData, "sdf");
                    global_3d_viewer.setStyle({}, {'stick':{}}); // استایل برای همه مدل‌های جدید
                    global_3d_viewer.zoomTo();
                    console.log("Rendering new model.");
                    global_3d_viewer.render();
                } catch (e) {
                    console.error("Error updating 3D model:", e);
                    // نمایش خطا در خود div
                     var container = document.getElementById('viewer_3dmol_container');
                     if (container) container.innerHTML = "<p style='color:red'>خطا در رندر مدل: " + e.message + "</p>";
                }
            } else if (sdfData && sdfData.startsWith("خطا:")) {
                 var container = document.getElementById('viewer_3dmol_container');
                 if (container) container.innerHTML = "<p style='color:red'>" + sdfData + "</p>";
                 if (global_3d_viewer) global_3d_viewer.removeAllModels(); global_3d_viewer.render();
            } else {
                console.log("No valid SDF data provided to update3DModel. Clearing viewer.");
                if (global_3d_viewer) {
                    global_3d_viewer.removeAllModels();
                    global_3d_viewer.render();
                }
            }
        };

        // شروع فرآیند مقداردهی اولیه نمایشگر
        // منتظر بمانید تا DOM آماده شود یا مستقیماً فراخوانی کنید اگر اسکریپت در انتهای body است
        if (document.readyState === "loading") {
            document.addEventListener("DOMContentLoaded", init3DViewerWhenReady);
        } else {
            init3DViewerWhenReady(); // DOM از قبل آماده است
        }
    </script>
    """
    return html_content

# --- بخش Gradio Interface با استفاده از gr.Blocks ---
with gr.Blocks(theme=gr.themes.Soft()) as iface:
    gr.Markdown(
        "# یابنده و نمایشگر ایزومرهای آلکان (دو بعدی و سه بعدی پیشرفته)\n"
        "این نسخه از روش جدیدی برای نمایش سه‌بعدی با استفاده از JS استفاده می‌کند."
    )

    # کامپوننت‌های مخفی برای داده‌های SDF و ماشه به‌روزرسانی JS
    sdf_data_store = gr.Textbox(label="SDF Data", visible=False, elem_id="sdf_data_store_elem")
    # elem_id برای دسترسی از JS لازم است
    js_trigger = gr.Textbox(label="JS Trigger", visible=False, elem_id="js_trigger_elem")

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

    with gr.Row():
        gallery_output = gr.Gallery(
            label="ایزومرهای یافت شده (2D)",
            columns=[3], height="auto", object_fit="contain"
        )
        # کامپوننت HTML که در ابتدا با تابع initial_3d_viewer_html_setup پر می‌شود
        html_3d_output = gr.HTML(label="نمایش سه‌بعدی ایزومر")

    gr.Examples(
        examples=[["butane"], ["pentane"]],
        inputs=molecule_name_input,
    )
    
    # رویداد بارگذاری اولیه صفحه برای تنظیم نمایشگر سه‌بعدی
    iface.load(initial_3d_viewer_html_setup, inputs=None, outputs=[html_3d_output])

    # اتصال رویداد کلیک دکمه
    submit_button.click(
        fn=find_and_display_isomers_and_cids,
        inputs=molecule_name_input,
        # خروجی‌ها شامل کامپوننت‌های مخفی نمی‌شوند، آن‌ها توسط تابع دیگر به‌روز می‌شوند
        outputs=[gallery_output, status_output, isomer_cids_state, sdf_data_store, js_trigger]
    )
    molecule_name_input.submit(
        fn=find_and_display_isomers_and_cids,
        inputs=molecule_name_input,
        outputs=[gallery_output, status_output, isomer_cids_state, sdf_data_store, js_trigger]
    )

    # اتصال رویداد select گالری به تابع جدید که SDF و ماشه را برمی‌گرداند
    gallery_output.select(
        fn=handle_gallery_selection_for_sdf,
        inputs=[isomer_cids_state], # evt به طور خودکار پاس داده می‌شود
        outputs=[sdf_data_store, js_trigger] # به‌روزرسانی کامپوننت‌های مخفی
    )

    # listener جاوااسکریپت برای تغییرات در js_trigger
    # این بخش مهم است و باید به درستی کار کند.
    # ما از یک تابع ساختگی `_js` استفاده می‌کنیم که Gradio ممکن است مستقیماً آن را ارائه ندهد.
    # روش صحیح‌تر استفاده از `gr.HTML(...).then(...)` است اگر بخواهیم کد JS را در پایتون بنویسیم
    # یا مستقیماً در فایل JS خارجی بنویسیم و به تغییرات گوش دهیم.
    # در اینجا یک راه‌حل با استفاده از تغییر خروجی html_3d_output برای اجرای JS امتحان می‌کنیم
    # (این یک هک است و ممکن است بهترین روش نباشد)

    # یک روش بهتر برای اجرای JS پس از تغییر یک کامپوننت، استفاده از `gr.HTML().then()`
    # یا استفاده از یک دکمه نامرئی است که با تغییر js_trigger کلیک می‌شود و آن دکمه یک تابع JS را اجرا می‌کند.

    # ساده‌سازی: فرض می‌کنیم که تابع window.update3DModel در HTML اولیه تعریف شده است.
    # و ما فقط باید آن را با استفاده از یک ترفند JS فراخوانی کنیم.
    # Gradio اجازه اجرای مستقیم JS دلخواه به عنوان خروجی یک تابع را نمی‌دهد.
    # باید یک کامپوننت HTML را با یک اسکریپت جدید به‌روز کنیم.

    # این بخش نیاز به بازبینی دارد. ساده‌ترین راه این است که JS listener به تغییرات
    # `js_trigger_elem` و `sdf_data_store_elem` در HTML اولیه گوش دهد.

    # با توجه به محدودیت‌های Gradio برای اجرای مستقیم JS پیچیده از پایتون،
    # listener جاوااسکریپت در `initial_3d_viewer_html_setup` باید به گونه‌ای نوشته شود
    # که به تغییرات مقادیر در #sdf_data_store_elem و #js_trigger_elem گوش دهد.
    # Gradio به طور خودکار این مقادیر را در DOM به‌روز می‌کند.

    # اصلاح HTML اولیه برای شامل شدن listener ها:
    def initial_3d_viewer_html_setup_with_listeners():
        html_content = """
        <div id="viewer_3dmol_container" style="width: 500px; height: 400px; border: 1px solid #ccc; position: relative;">
            <p id="loading_3dmol_message" style="text-align: center; padding-top: 50px;">در حال بارگذاری نمایشگر سه‌بعدی...</p>
        </div>
        <!-- اینها باید مقادیرشان توسط Gradio به‌روز شود -->
        <!-- <textarea id="sdf_data_store_elem_input" style="display:none;"></textarea> -->
        <!-- <input type="text" id="js_trigger_elem_input" style="display:none;"> -->


        <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.2/3Dmol-min.js"></script>
        <script>
            var global_3d_viewer = null;

            function init3DViewerWhenReady() {
                // ... (کد init3DViewerWhenReady مانند قبل) ...
                if (typeof $3Dmol !== 'undefined' && $3Dmol.createViewer) {
                    try {
                        var container = document.getElementById('viewer_3dmol_container');
                        if (!container) { console.error('Container viewer_3dmol_container not found!'); return; }
                        var loadingMsg = document.getElementById('loading_3dmol_message');
                        if (loadingMsg) loadingMsg.style.display = 'none';
                        if (!global_3d_viewer) {
                           global_3d_viewer = $3Dmol.createViewer(container, { backgroundColor: 'white' });
                           console.log("3Dmol Viewer initialized.");
                        } else { console.log("3Dmol Viewer already initialized.");}
                        global_3d_viewer.render();
                    } catch (e) { console.error("Error initializing 3Dmol viewer:", e); container.innerHTML = "<p style='color:red'>خطا: " + e.message + "</p>";}
                } else { setTimeout(init3DViewerWhenReady, 100); }
            }

            window.update3DModel = function(sdfData, modelId) {
                // ... (کد window.update3DModel مانند قبل) ...
                console.log("window.update3DModel called with modelId:", modelId);
                if (!global_3d_viewer) {
                    console.error("Global 3D viewer is not initialized. Retrying init.");
                    init3DViewerWhenReady(); 
                    setTimeout(function() { if(global_3d_viewer) window.update3DModel(sdfData, modelId); else console.error("Still no viewer after retry"); }, 500);
                    return;
                }
                if (sdfData && sdfData.trim() !== "" && !sdfData.startsWith("خطا:")) {
                    try {
                        global_3d_viewer.removeAllModels(); 
                        global_3d_viewer.addModel(sdfData, "sdf");
                        global_3d_viewer.setStyle({}, {'stick':{}}); 
                        global_3d_viewer.zoomTo();
                        global_3d_viewer.render();
                        console.log("Model updated for ID:", modelId);
                    } catch (e) { console.error("Error updating 3D model:", e); }
                } else if (sdfData && sdfData.startsWith("خطا:")) {
                     var container = document.getElementById('viewer_3dmol_container');
                     if (container) container.innerHTML = "<p style='color:red'>" + sdfData + "</p>";
                     if (global_3d_viewer) global_3d_viewer.removeAllModels(); global_3d_viewer.render();
                } else {
                    if (global_3d_viewer) { global_3d_viewer.removeAllModels(); global_3d_viewer.render(); }
                }
            };

            // Listener برای تغییرات در کامپوننت‌های مخفی
            // Gradio به طور خودکار مقدار value این input ها را به‌روز می‌کند وقتی از خروجی تابع پایتون می‌آیند.
            // ما به تغییرات آنها گوش می‌دهیم.
            
            // پیدا کردن المان‌های واقعی که Gradio برای Textbox ایجاد می‌کند.
            // این المان‌ها معمولا یک input یا textarea درون یک div با elem_id هستند.
            function setupMutationObserver() {
                const triggerInput = document.querySelector("#js_trigger_elem textarea"); // یا input اگر تک خطی است
                
                if (triggerInput) {
                    console.log("JS Trigger input found. Setting up MutationObserver.");
                    const observer = new MutationObserver(function(mutationsList, observer) {
                        for(let mutation of mutationsList) {
                            if (mutation.type === 'characterData' || (mutation.type === 'attributes' && mutation.attributeName === 'value') || mutation.target.value !== (mutation.oldValue || '')) {
                                console.log('JS Trigger changed:', triggerInput.value);
                                const sdfDataInput = document.querySelector("#sdf_data_store_elem textarea"); // یا input
                                if (sdfDataInput) {
                                    console.log("SDF Data Store value (on trigger):", sdfDataInput.value.substring(0,100) + "...");
                                    window.update3DModel(sdfDataInput.value, triggerInput.value); // triggerInput.value همان modelId است
                                } else {
                                    console.error("SDF Data Store element not found on trigger.");
                                }
                                return; // فقط برای اولین تغییر عمل کن
                            }
                        }
                    });
                    // برای textarea، باید به تغییرات فرزند text node آن گوش دهیم.
                    // یا به تغییرات value اگر Gradio آن را به عنوان attribute تغییر می‌دهد.
                    // یک راه ساده‌تر ممکن است استفاده از یک event listener 'change' یا 'input' باشد.
                    // اما MutationObserver قابل اعتمادتر است برای تغییرات برنامه نویسی شده.
                    // observer.observe(triggerInput, { attributes: true, childList: true, subtree: true, characterData: true });
                    
                    // یک روش ساده‌تر با event listener:
                    triggerInput.addEventListener('input', function() { // یا 'change'
                        console.log('JS Trigger input event:', triggerInput.value);
                        const sdfDataInput = document.querySelector("#sdf_data_store_elem textarea");
                        if (sdfDataInput) {
                            window.update3DModel(sdfDataInput.value, triggerInput.value);
                        }
                    });
                     // همچنین برای بار اول
                    const sdfDataInputInitial = document.querySelector("#sdf_data_store_elem textarea");
                    if(triggerInput.value && sdfDataInputInitial && sdfDataInputInitial.value){
                         window.update3DModel(sdfDataInputInitial.value, triggerInput.value);
                    }


                } else {
                    console.error("JS Trigger element (#js_trigger_elem textarea) not found. Retrying observer setup...");
                    setTimeout(setupMutationObserver, 500); // تلاش مجدد
                }
            }


            if (document.readyState === "loading") {
                document.addEventListener("DOMContentLoaded", function() {
                    init3DViewerWhenReady();
                    setupMutationObserver();
                });
            } else {
                init3DViewerWhenReady();
                setupMutationObserver();
            }
        </script>
        """
        return html_content

    # استفاده از تابع جدید برای بارگذاری HTML اولیه
    iface.load(initial_3d_viewer_html_setup_with_listeners, inputs=None, outputs=[html_3d_output])


if __name__ == '__main__':
    iface.launch(debug=True) # debug=True برای دیدن لاگ‌های بیشتر در کنسول مرورگر
