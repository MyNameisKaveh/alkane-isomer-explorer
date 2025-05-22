# قدم 1: اطمینان از نصب کتابخانه‌های مورد نیاز
# این خطوط را نیازی نیست در app.py بگذارید، بلکه در فایل requirements.txt باید باشند.
# pip install pubchempy rdkit-pypi gradio py3dmol

import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage
import gradio as gr
import traceback
import os 
import tempfile 
import uuid # برای تولید شناسه‌های یکتا برای div های 3D viewer

# --- تابع کمکی برای رسم مولکول 2D ---
def draw_molecule(smiles_string):
    """
    مولکول را از SMILES رسم کرده و یک تصویر PIL Image برمی‌گرداند.
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

# --- تابع برای نمایش سه‌بعدی (اصلاح شده برای Gradio و خطای f-string) ---
def get_3d_viewer_html(cid, style='stick'):
    """
    ساختار سه‌بعدی یک مولکول را بر اساس CID آن از PubChem دریافت کرده و HTML نمایشگر 3Dmol.js را برمی‌گرداند.
    این HTML شامل یک div و یک تگ script است که مدل 3D را در مرورگر کاربر بارگذاری می‌کند.
    """
    if cid is None or cid == "" or cid == "N/A": 
        return "<p style='text-align: center; color: gray;'>برای نمایش ساختار سه‌بعدی، یک ایزومر را از لیست بالا انتخاب کنید.</p>"

    html_output = f"<p style='text-align: center;'>در حال بارگذاری ساختار سه‌بعدی برای CID: {cid}...</p>"
    sdf_content = None
    temp_sdf_path = None # مقداردهی اولیه برای اطمینان از پاک شدن فایل

    try:
        # ایجاد یک فایل موقت برای ذخیره SDF
        with tempfile.NamedTemporaryFile(delete=False, suffix='.sdf') as temp_sdf_file:
            temp_sdf_path = temp_sdf_file.name

        # دانلود فایل SDF از PubChem
        pcp.download('SDF', temp_sdf_path, cid, 'cid', record_type='3d', overwrite=True)

        # خواندن محتوای فایل SDF
        with open(temp_sdf_path, 'r') as f:
            sdf_content = f.read()

        if not sdf_content:
            html_output = f"<p style='color: red; text-align: center;'>فایل 3D SDF برای CID {cid} خالی بود. ممکن است ساختار سه‌بعدی در دسترس نباشد.</p>"
        else:
            # گام 1: Escape کردن کاراکترهای خاص برای JavaScript (newline, backslash, single quote)
            # این کاراکترها اگر در رشته جاوااسکریپت (literal) استفاده شوند، باید escape شوند.
            js_sdf_content_for_js_string = sdf_content.replace('\\', '\\\\').replace('\n', '\\n').replace("'", "\\'")
            
            # گام 2: Escape کردن کاراکترهای '{' و '}' برای Python f-string
            # اگر محتوای SDF شامل '{' یا '}' باشد (که می‌تواند باشد، مثلاً در PROPERTY block),
            # و این محتوا داخل یک f-string پایتون قرار گیرد، مفسر f-string به مشکل می‌خورد.
            # با تبدیل '{' به '{{' و '}' به '}}', این کاراکترها به صورت literal در f-string در نظر گرفته می‌شوند.
            js_sdf_content_final = js_sdf_content_for_js_string.replace('{', '{{').replace('}', '}}')

            # تولید یک ID یکتا برای div نمایشگر 3D. این برای اطمینان از عدم تداخل هنگام به‌روزرسانی مهم است.
            viewer_div_id = f"viewer_{uuid.uuid4().hex}" 

            # ساخت HTML نهایی شامل div و کد JavaScript برای بارگذاری مدل 3D
            html_output = f"""
            <div id="{viewer_div_id}" style="height: 400px; width: 450px; margin: auto; border: 1px solid #ccc; border-radius: 5px;"></div>
            <script type="text/javascript">
                // یک تابع خوداجرا برای جلوگیری از آلودگی global scope
                (function() {{
                    var element = document.getElementById('{viewer_div_id}');
                    if (!element) {{
                        console.error('3D Viewer Error: Element not found for ID: {viewer_div_id}');
                        return; // اگر عنصر div پیدا نشد، ادامه نده
                    }}
                    element.innerHTML = ''; // پاک کردن محتوای قبلی div در صورت re-render

                    // اطمینان از بارگذاری کتابخانه 3Dmol.js
                    if (typeof $3Dmol === 'undefined') {{
                        console.error('3D Viewer Error: $3Dmol library is not loaded. Ensure 3Dmol-min.js is included.');
                        element.innerHTML = "<p style='color: red; text-align: center;'>خطا: کتابخانه 3Dmol بارگذاری نشده است. لطفاً صفحه را رفرش کنید.</p>";
                        return;
                    }}

                    try {{
                        var viewer = $3Dmol.createViewer( element, {{defaultcolors: $3Dmol.elementColors.default}} );
                        // اضافه کردن مدل با محتوای SDF Escape شده (برای جاوااسکریپت)
                        viewer.addModel('{js_sdf_content_final}', 'sdf'); 
                        
                        // تنظیم استایل نمایش (stick, sphere, line)
                        var style_dict = {{ "{style}": {{}} }}; 
                        if ("{style}" === 'cartoon') {{ // تنظیمات خاص برای استایل 'cartoon'
                            style_dict['cartoon']['color'] = 'spectrum';
                        }}
                        viewer.setStyle(style_dict);

                        viewer.setBackgroundColor('0xeeeeee'); // رنگ پس‌زمینه
                        viewer.zoomTo(); // زوم به اندازه مولکول
                        viewer.render(); // رندر کردن مدل
                    }} catch (js_error) {{
                        console.error('3D Viewer JavaScript Error:', js_error);
                        element.innerHTML = "<p style='color: red; text-align: center;'>خطا در رندرینگ 3D: " + js_error.message + "</p>";
                    }}
                }})(); // پایان تابع خوداجرا
            </script>
            """

    except pcp.NotFoundError:
        html_output = f"<p style='color: orange; text-align: center;'>ساختار 3D SDF برای CID {cid} در PubChem یافت نشد.</p>"
    except Exception as e:
        html_output = f"<p style='color: red; text-align: center;'>خطا در نمایش ساختار 3D (پایتون): {e}</p>"
        print(f"FULL TRACEBACK for 3D rendering (Python): {traceback.format_exc()}")
    finally:
        # پاک کردن فایل موقت SDF
        if temp_sdf_path and os.path.exists(temp_sdf_path):
            os.remove(temp_sdf_path)
            
    return html_output

# --- تابع اصلی find_and_display_isomers (اکنون از gr.update() استفاده می‌کند) ---
def find_and_display_isomers(molecule_name_input):
    """
    ایزومرهای آلکان را پیدا کرده، ساختارهای 2D را آماده می‌کند،
    و داده‌ها را برای Dropdown و نمایش 3D بازمی‌گرداند.
    """
    if not molecule_name_input or not molecule_name_input.strip():
        # بازگشت مقادیر پیش‌فرض با gr.update() برای تنظیم مجدد Dropdown
        return [], gr.update(choices=[], value=None), "<p style='text-align: center; color: gray;'>نام یک آلکان را وارد کنید تا ایزومرها نمایش داده شوند.</p>", "لطفا نام یک مولکول را وارد کنید."

    molecule_name = molecule_name_input.strip().lower()
    print(f"Processing request for: '{molecule_name}'")
    
    status_message = ""
    isomer_outputs_final_2d = [] # برای گالری 2D
    isomer_choices_for_3d = [] # برای دراپ‌داون 3D (نام و CID)

    try:
        print(f"Searching for compound: '{molecule_name}' in PubChem...")
        compounds = pcp.get_compounds(molecule_name, 'name', listkey_count=5) 
        
        main_compound_obj = None
        molecular_formula = None

        if not compounds:
            status_message = f"مولکول '{molecule_name}' در PubChem یافت نشد."
            print(status_message)
            return [], gr.update(choices=[], value=None), "<p style='text-align: center; color: red;'>مولکول یافت نشد.</p>", status_message
        
        print(f"Found {len(compounds)} potential matches for '{molecule_name}'. Checking them...")
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
                                is_standard_hydrocarbon = False
                            atom_symbols_main = set()
                            if is_standard_hydrocarbon:
                                for atom in mol_obj.GetAtoms():
                                    atom_symbols_main.add(atom.GetSymbol())
                                    if atom.GetIsotope() != 0:
                                        is_standard_hydrocarbon = False
                                        break
                                if not is_standard_hydrocarbon: continue

                                if not atom_symbols_main.issubset({'C', 'H'}):
                                    is_standard_hydrocarbon = False
                            if is_standard_hydrocarbon:
                                for bond in mol_obj.GetBonds():
                                    if bond.GetBondType() != Chem.BondType.SINGLE:
                                        is_standard_hydrocarbon = False
                                        break
                                if not is_standard_hydrocarbon: continue
                                
                                if Chem.GetSymmSSSR(mol_obj):
                                    is_standard_hydrocarbon = False
                        else: 
                            is_standard_hydrocarbon = False
                    except Exception as rdkit_ex:
                        print(f"    RDKit error processing SMILES for main candidate CID {cid}: {rdkit_ex}")
                        is_standard_hydrocarbon = False 
                else: 
                    is_standard_hydrocarbon = False 

                if not is_standard_hydrocarbon: continue

                current_compound_name_matches_input = molecule_name in [syn.lower() for syn in c.synonyms]
                if current_compound_name_matches_input: 
                    main_compound_obj = c
                    molecular_formula = actual_formula
                    print(f"  SELECTED main compound: CID {main_compound_obj.cid}, Formula: {molecular_formula}")
                    break 
                
                if not main_compound_obj: 
                    main_compound_obj = c 
                    molecular_formula = actual_formula 
                    print(f"  TENTATIVELY selected main compound: CID {main_compound_obj.cid}, Formula: {molecular_formula}")
        
        if not main_compound_obj or not molecular_formula: 
            status_message = f"آلکان استاندارد با نام '{molecule_name}' در PubChem یافت نشد."
            print(status_message)
            return [], gr.update(choices=[], value=None), "<p style='text-align: center; color: red;'>آلکان استاندارد یافت نشد.</p>", status_message
        
        print(f"Proceeding with main compound: CID {main_compound_obj.cid}, Formula: {molecular_formula}")
        print(f"Searching for isomers with formula: {molecular_formula}...")
        isomers_found_raw = pcp.get_compounds(molecular_formula, 'formula', listkey_count=50) 

        if not isomers_found_raw:
            status_message = f"ایزومری برای فرمول {molecular_formula} یافت نشد."
            print(status_message)
            return [], gr.update(choices=[], value=None), "<p style='text-align: center; color: orange;'>ایزومری یافت نشد.</p>", status_message

        print(f"Found {len(isomers_found_raw)} potential isomer entries from PubChem. Filtering for true structural alkane isomers...")
        
        valid_structural_alkanes_entries = [] 
        unique_accepted_smiles = set()

        for isomer_entry in isomers_found_raw:
            smiles = isomer_entry.canonical_smiles
            if not smiles:
                continue

            try:
                mol_iso = Chem.MolFromSmiles(smiles)
                if not mol_iso:
                    continue

                is_valid_candidate = True
                
                if len(Chem.GetMolFrags(mol_iso)) > 1:
                    is_valid_candidate = False
                
                if is_valid_candidate:
                    atom_symbols = set()
                    for atom in mol_iso.GetAtoms():
                        atom_symbols.add(atom.GetSymbol())
                    if not atom_symbols.issubset({'C', 'H'}):
                        is_valid_candidate = False
                
                if is_valid_candidate:
                    for atom in mol_iso.GetAtoms():
                        if atom.GetSymbol() == 'H' and atom.GetDegree() == 0:
                            is_valid_candidate = False
                            break
                        if atom.GetIsotope() != 0:
                            is_valid_candidate = False
                            break 
                
                if is_valid_candidate:
                    for bond in mol_iso.GetBonds():
                        if bond.GetBondType() != Chem.BondType.SINGLE:
                            is_valid_candidate = False
                            break
                    if not is_valid_candidate: continue

                    if Chem.GetSymmSSSR(mol_iso):
                        is_valid_candidate = False

                if is_valid_candidate:
                    canonical_smiles_for_uniqueness = Chem.MolToSmiles(mol_iso, isomericSmiles=False)
                    if canonical_smiles_for_uniqueness not in unique_accepted_smiles:
                        valid_structural_alkanes_entries.append(isomer_entry)
                        unique_accepted_smiles.add(canonical_smiles_for_uniqueness)
                    else:
                        print(f"  Skipping (Duplicate structure based on non-isomeric SMILES): CID {isomer_entry.cid}, SMILES: {smiles}")

            except Exception as rdkit_iso_ex:
                print(f"  RDKit or processing error for isomer SMILES CID {isomer_entry.cid}: {rdkit_iso_ex}, SMILES: {smiles}")
                continue
        
        print(f"Found {len(valid_structural_alkanes_entries)} unique, valid structural alkane isomers after filtering.")
        
        valid_isomers_count_final = 0

        for final_isomer_entry in valid_structural_alkanes_entries:
            smiles_to_draw = final_isomer_entry.canonical_smiles
            iupac_name = final_isomer_entry.iupac_name
            cid = final_isomer_entry.cid
            
            display_name = ""
            if iupac_name:
                display_name = iupac_name
            elif final_isomer_entry.synonyms:
                chosen_synonym = final_isomer_entry.synonyms[0]
                simple_names = [s for s in final_isomer_entry.synonyms if s.lower().endswith("ane") and not any(char.isdigit() for char in s.split('-')[0]) and '-' not in s.split(' ')[0]]
                if simple_names:
                    chosen_synonym = min(simple_names, key=len)
                else:
                    non_iupac_synonyms = [s for s in final_isomer_entry.synonyms if s != final_isomer_entry.iupac_name]
                    if non_iupac_synonyms:
                        chosen_synonym = min(non_iupac_synonyms, key=len)
                display_name = chosen_synonym.capitalize()
            else:
                display_name = f"Alkane (CID: {cid})"

            mol_image = draw_molecule(smiles_to_draw)
            if mol_image:
                isomer_outputs_final_2d.append((mol_image, f"{display_name}\nSMILES: {smiles_to_draw}"))
                isomer_choices_for_3d.append((display_name, str(cid))) # نام و CID برای دراپ‌داون 3D
                valid_isomers_count_final += 1
            else:
                print(f"  Failed to draw image for accepted isomer: CID {cid}, SMILES: {smiles_to_draw}")

        print(f"Displayed {valid_isomers_count_final} isomers in the gallery.")

        if not isomer_outputs_final_2d:
            status_message = "ایزومر آلکان استاندارد و قابل رسمی پیدا نشد."
            if len(valid_structural_alkanes_entries) > 0:
                status_message += " (برخی در مرحله رسم ناموفق بودند یا کاندیدای معتبری نبودند)."
            return [], gr.update(choices=[], value=None), "<p style='text-align: center; color: orange;'>ایزومرها یافت نشدند یا قابل رسم نبودند.</p>", status_message
        else:
            status_message = f"{len(isomer_outputs_final_2d)} ایزومر ساختاری آلکان برای '{molecule_name_input}' (فرمول: {molecular_formula}) پیدا و نمایش داده شد."
        
        isomer_outputs_final_2d.sort(key=lambda x: x[1])
        isomer_choices_for_3d.sort(key=lambda x: x[0])

        # اولین ایزومر را به عنوان پیش‌فرض برای نمایش 3D انتخاب می‌کنیم
        initial_3d_cid = isomer_choices_for_3d[0][1] if isomer_choices_for_3d else None
        
        # آماده کردن gr.update() برای به‌روزرسانی همزمان choices و value
        dropdown_update = gr.update(
            choices=isomer_choices_for_3d, 
            value=initial_3d_cid # مقدار پیش‌فرض را تنظیم می‌کند
        )

        initial_3d_html = get_3d_viewer_html(initial_3d_cid)

        # ترتیب خروجی‌ها باید با outputs در search_button.click و غیره مطابقت داشته باشد.
        return isomer_outputs_final_2d, dropdown_update, initial_3d_html, status_message

    except pcp.PubChemHTTPError as e:
        error_msg = f"خطا در ارتباط با PubChem: {e}."
        print(error_msg)
        print(f"FULL TRACEBACK for PubChemHTTPError: {traceback.format_exc()}")
        return [], gr.update(choices=[], value=None), f"<p style='text-align: center; color: red;'>خطا در PubChem: {e}</p>", error_msg
    except Exception as e:
        error_msg = f"یک خطای غیرمنتظره در سرور رخ داد: {e}"
        print(f"FULL TRACEBACK for general Exception: {traceback.format_exc()}")
        return [], gr.update(choices=[], value=None), f"<p style='text-align: center; color: red;'>خطای غیرمنتظره: {e}</p>", error_msg

# --- بخش Gradio Interface (با استفاده از gr.Blocks) ---

with gr.Blocks(theme=gr.themes.Soft(), title="یابنده و نمایشگر ایزومرهای آلکان") as demo:
    # این تگ script کتابخانه 3Dmol.js را یک بار در ابتدای بارگذاری صفحه، بارگذاری می‌کند
    # این خط بسیار مهم است و کمک می‌کند 3Dmol در مرورگر اجرا شود.
    gr.HTML('<script src="https://3dmol.org/build/3Dmol-min.js"></script>')

    gr.Markdown(
        """
        # 🧪 یابنده و نمایشگر ایزومرهای آلکان ⌬
        نام یک آلکان (به انگلیسی) را وارد کنید تا ایزومرهای آن به همراه ساختار شیمیایی و نام IUPAC (یا رایج) نمایش داده شوند.
        اطلاعات از دیتابیس PubChem دریافت شده و ساختارها با استفاده از کتابخانه‌های RDKit و Py3Dmol رسم می‌شوند.
        """
    )

    with gr.Row():
        molecule_name_input = gr.Textbox(
            label="نام آلکان را وارد کنید", 
            placeholder="مثال: butane, pentane, hexane",
            info="نام آلکان مورد نظر خود را به انگلیسی و با حروف کوچک وارد کنید."
        )
        search_button = gr.Button("جستجو")

    status_message_output = gr.Textbox(label="وضعیت و پیام‌ها", interactive=False)

    with gr.Tabs():
        with gr.TabItem("📊 ساختارهای 2D"):
            gallery_2d_output = gr.Gallery(
                label="ایزومرهای یافت شده", 
                columns=[3], 
                height="auto", 
                object_fit="contain"
            )

        with gr.TabItem("🔬 ساختارهای 3D"):
            with gr.Row():
                isomer_3d_selector = gr.Dropdown(
                    label="ایزومر مورد نظر را برای نمایش 3D انتخاب کنید",
                    choices=[], # در زمان اجرا پر می‌شوند
                    interactive=True
                )
                style_3d_selector = gr.Dropdown(
                    label="استایل نمایش 3D",
                    choices=['stick', 'sphere', 'line'],
                    value='stick', 
                    interactive=True
                )
            viewer_3d_html = gr.HTML(
                value="<p style='text-align: center; color: gray;'>برای نمایش ساختار سه‌بعدی، یک ایزومر را از لیست بالا انتخاب کنید.</p>",
                elem_id="py3dmol-viewer" # یک ID برای CSS یا دسترسی JavaScript (اختیاری)
            )
    
    # --- مدیریت رویدادها ---
    # لیست کامپوننت‌های خروجی که تابع find_and_display_isomers آن‌ها را به‌روز می‌کند.
    # ترتیب آن‌ها باید دقیقاً با ترتیب return در تابع find_and_display_isomers مطابقت داشته باشد.
    output_components = [
        gallery_2d_output,
        isomer_3d_selector, # این کامپوننت حالا یک gr.update() object دریافت می‌کند
        viewer_3d_html,
        status_message_output
    ]

    search_button.click(
        fn=find_and_display_isomers,
        inputs=[molecule_name_input],
        outputs=output_components,
        show_progress=True
    )
    
    # می‌توان رویداد submit را نیز برای ورودی text_box تعریف کرد.
    molecule_name_input.submit(
        fn=find_and_display_isomers,
        inputs=[molecule_name_input],
        outputs=output_components,
        show_progress=True
    )

    # رویداد تغییر انتخاب دراپ‌داون 3D یا تغییر استایل 3D
    # وقتی ایزومر یا استایل تغییر می‌کند، get_3d_viewer_html با مقادیر جدید فراخوانی می‌شود
    # و HTML نمایشگر 3D را به‌روز می‌کند.
    isomer_3d_selector.change(
        fn=get_3d_viewer_html,
        inputs=[isomer_3d_selector, style_3d_selector], # هم CID (از Dropdown) و هم استایل را می‌گیرد
        outputs=[viewer_3d_html],
        show_progress=True
    )
    style_3d_selector.change(
        fn=get_3d_viewer_html,
        inputs=[isomer_3d_selector, style_3d_selector], # هم CID (از Dropdown) و هم استایل را می‌گیرد
        outputs=[viewer_3d_html],
        show_progress=True
    )

    # مثال‌ها برای راحتی کاربر
    gr.Examples(
        examples=[
            ["butane"], 
            ["pentane"], 
            ["hexane"],
            ["heptane"] 
        ],
        inputs=molecule_name_input,
        outputs=output_components, # استفاده از لیست خروجی‌های یکسان
        fn=find_and_display_isomers,
        cache_examples=False,
        run_on_click=True # اجرا شدن مثال‌ها با کلیک، نه به صورت خودکار
    )

if __name__ == '__main__':
    demo.launch()
