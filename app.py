import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage # دیگر برای 2D اصلی استفاده نمی‌شود، اما شاید برای موارد دیگر لازم باشد
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D # برای رسم SVG
import py3Dmol
import os
import traceback
# from PIL import Image # دیگر برای تصاویر 2D لازم نیست

# --- توابع کمکی ---

def draw_molecule_svg(smiles_string, width=250, height=200): # اندازه‌ها را می‌توان تنظیم کرد
    """یک مولکول را از رشته SMILES به صورت SVG با RDKit رسم می‌کند."""
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if not mol:
            print(f"خطا در پارس کردن SMILES برای SVG: {smiles_string}")
            return None

        # تولید مختصات 2D برای رسم بهتر
        rdDepictor.Compute2DCoords(mol)
        
        # ایجاد drawer SVG
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        
        # گزینه‌های رسم (اختیاری، برای سفارشی‌سازی بیشتر)
        # drawer.drawOptions().atomLabels = {} # عدم نمایش لیبل اتم‌ها (فقط ساختار)
        # drawer.drawOptions().addStereoAnnotation = True
        # drawer.drawOptions().bondLineWidth = 1.5
        # drawer.drawOptions().padding = 0.05 # کمی حاشیه

        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg_output = drawer.GetDrawingText()
        return svg_output
    except Exception as e:
        print(f"خطا در ایجاد SVG برای SMILES {smiles_string}: {e}")
        return None

def get_sdf_content(cid):
    # ... (این تابع بدون تغییر باقی می‌ماند) ...
    if cid is None:
        return None, "CID برای دریافت محتوای SDF ارائه نشده است."
    print(f"در حال دریافت محتوای SDF برای CID: {cid}...")
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
            error_message, sdf_data = f"فایل SDF برای CID {cid} خالی است.", None
    except pcp.NotFoundError: error_message = f"ساختار سه‌بعدی (SDF) برای CID {cid} یافت نشد."
    except FileNotFoundError: error_message = f"فایل موقت SDF برای CID {cid} یافت نشد."
    except Exception as e: error_message = f"خطا در دانلود SDF برای CID {cid}: {e}\n{traceback.format_exc()}"
    finally:
        if os.path.exists(temp_sdf_file):
            try: os.remove(temp_sdf_file)
            except Exception: pass
    if error_message and 'st' in globals() and hasattr(st, 'warning'): st.warning(error_message)
    return sdf_data, error_message

def generate_3d_viewer_html(sdf_data, molecule_name, width=500, height=400):
    # ... (این تابع بدون تغییر باقی می‌ماند) ...
    if not sdf_data: return "<p style='color:orange;'>داده SDF برای نمایش سه‌بعدی موجود نیست.</p>"
    try:
        viewer = py3Dmol.view(width=width, height=height)
        viewer.addModel(sdf_data, 'sdf')
        viewer.setStyle({'stick': {}})
        viewer.setBackgroundColor('0xeeeeee')
        viewer.zoomTo()
        return viewer._make_html()
    except Exception as e:
        if 'st' in globals() and hasattr(st, 'error'): st.error(f"خطا در ایجاد نمایشگر سه‌بعدی برای {molecule_name}: {e}")
        return f"<p style='color:red;'>خطا در رندر سه‌بعدی برای {molecule_name}: {e}</p>"

def process_alkane_request(molecule_name_input):
    if not molecule_name_input or not molecule_name_input.strip():
        return [], "لطفا نام یک مولکول را وارد کنید."
    molecule_name = molecule_name_input.strip().lower()
    status_message, isomer_details_list = f"در حال جستجو برای '{molecule_name}'...", []
    try:
        compounds = pcp.get_compounds(molecule_name, 'name', listkey_count=5)
        main_compound_obj, molecular_formula = None, None
        if not compounds: return [], f"مولکول '{molecule_name}' در PubChem یافت نشد."
        for i, c in enumerate(compounds): # ... (منطق انتخاب ترکیب اصلی) ...
            cid = c.cid
            actual_formula = c.molecular_formula if hasattr(c, 'molecular_formula') else None
            if actual_formula:
                is_standard_hydrocarbon = True 
                if c.canonical_smiles:
                    try:
                        mol_obj = Chem.MolFromSmiles(c.canonical_smiles)
                        if mol_obj:
                            rdDepictor.Compute2DCoords(mol_obj)
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
            return [], f"آلکان استاندارد '{molecule_name}' یافت نشد یا با معیارها مطابقت ندارد."
        status_message = f"یافتن ایزومرها برای {main_compound_obj.iupac_name or molecule_name} ({molecular_formula})..."
        isomers_found_raw = pcp.get_compounds(molecular_formula, 'formula', listkey_count=50)
        if not isomers_found_raw: return [], f"ایزومری برای فرمول {molecular_formula} یافت نشد."
        valid_structural_alkanes_entries, unique_accepted_smiles = [], set()
        for isomer_entry in isomers_found_raw: # ... (منطق فیلتر کردن ایزومرها) ...
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
        if not valid_structural_alkanes_entries:
             return [], f"ایزومر آلکان معتبری برای {molecular_formula} پس از فیلتر کردن یافت نشد."
        
        for entry in sorted(valid_structural_alkanes_entries, key=lambda x: (len(x.canonical_smiles), x.cid)):
            # تغییر: استفاده از draw_molecule_svg
            svg_image_str = draw_molecule_svg(entry.canonical_smiles, width=280, height=220) # تنظیم اندازه مناسب برای گالری
            if svg_image_str:
                name = entry.iupac_name
                if not name and entry.synonyms:
                    simple_names = [s for s in entry.synonyms if s.lower().endswith("ane") and not any(char.isdigit() for char in s.split('-')[0]) and '-' not in s.split(' ')[0]]
                    if simple_names: name = min(simple_names, key=len)
                    else: name = entry.synonyms[0]
                elif not name: name = f"Alkane (CID: {entry.cid})"
                isomer_details_list.append({
                    "cid": entry.cid, 
                    "name": name.capitalize(), 
                    "smiles": entry.canonical_smiles, 
                    "image_svg": svg_image_str # ذخیره رشته SVG
                })
        if isomer_details_list: status_message = f"{len(isomer_details_list)} ایزومر برای '{molecule_name}' ({molecular_formula}) پیدا شد."
        else: status_message = f"ایزومر قابل نمایشی برای '{molecule_name}' یافت نشد."
        return isomer_details_list, status_message
    except pcp.PubChemHTTPError as e: return [], f"خطا در ارتباط با PubChem: {e}"
    except Exception as e: return [], f"خطای غیرمنتظره: {e}\n{traceback.format_exc()}"

# --- رابط کاربری Streamlit ---
st.set_page_config(page_title="نمایشگر ایزومرهای آلکان", layout="wide")
st.title("یابنده و نمایشگر ایزومرهای آلکان (با تصاویر SVG)")

# ... (بخش مقداردهی اولیه session state بدون تغییر) ...
if 'isomer_data' not in st.session_state: st.session_state.isomer_data = []
if 'status_message' not in st.session_state: st.session_state.status_message = ""
if 'selected_cid_for_3d' not in st.session_state: st.session_state.selected_cid_for_3d = None
if 'selected_name_for_3d' not in st.session_state: st.session_state.selected_name_for_3d = "" 
if 'molecule_searched' not in st.session_state: st.session_state.molecule_searched = ""
if 'molecule_name_input' not in st.session_state: st.session_state.molecule_name_input = "" 
if 'run_search_after_example' not in st.session_state: st.session_state.run_search_after_example = False

# ... (بخش ورودی و مثال‌ها در سایدبار بدون تغییر) ...
st.sidebar.header("جستجو و مثال‌ها")
current_molecule_input = st.sidebar.text_input(
    label="نام آلکان را وارد کنید:",
    value=st.session_state.molecule_name_input,
    placeholder="مثال: butane, pentane",
    key="sidebar_molecule_input" 
)
if current_molecule_input != st.session_state.molecule_name_input:
    st.session_state.molecule_name_input = current_molecule_input
    st.session_state.run_search_after_example = False 

example_alkanes = ["butane", "pentane", "hexane", "heptane"]
st.sidebar.subheader("یا یکی از مثال‌ها را انتخاب کنید:")
cols_examples = st.sidebar.columns(2) 
for i, example in enumerate(example_alkanes):
    if cols_examples[i % 2].button(example.capitalize(), key=f"example_{example}"):
        st.session_state.molecule_name_input = example 
        st.session_state.selected_cid_for_3d = None 
        st.session_state.selected_name_for_3d = ""
        st.session_state.run_search_after_example = True 
        st.rerun() 

search_button_sidebar = st.sidebar.button("جستجوی ایزومرها", type="primary", key="sidebar_search_button")

should_run_search = False
if search_button_sidebar and st.session_state.molecule_name_input:
    should_run_search = True
elif st.session_state.run_search_after_example and st.session_state.molecule_name_input:
    should_run_search = True
    st.session_state.run_search_after_example = False 

if should_run_search:
    with st.spinner(f"در حال پردازش برای {st.session_state.molecule_name_input}..."):
        isomers, status_msg = process_alkane_request(st.session_state.molecule_name_input)
        st.session_state.isomer_data = isomers
        st.session_state.status_message = status_msg
        st.session_state.selected_cid_for_3d = None 
        st.session_state.selected_name_for_3d = ""
        st.session_state.molecule_searched = st.session_state.molecule_name_input
elif search_button_sidebar and not st.session_state.molecule_name_input: 
    st.session_state.status_message = "لطفا نام یک مولکول را وارد کنید یا یک مثال انتخاب کنید."
    st.session_state.isomer_data = []
    st.session_state.selected_cid_for_3d = None
    st.session_state.selected_name_for_3d = ""
    st.session_state.molecule_searched = ""

if st.session_state.status_message:
    if "خطا" in st.session_state.status_message or "یافت نشد" in st.session_state.status_message :
        st.sidebar.error(st.session_state.status_message)
    else:
        st.sidebar.info(st.session_state.status_message)

# --- تب‌ها برای نمایش نتایج ---
if st.session_state.isomer_data or st.session_state.selected_cid_for_3d :
    tab1_title = f"گالری ایزومرها ({len(st.session_state.isomer_data)} مورد)" if st.session_state.isomer_data else "گالری ایزومرها"
    tab2_title = f"نمایش سه‌بعدی ({st.session_state.selected_name_for_3d})" if st.session_state.selected_cid_for_3d and st.session_state.selected_name_for_3d else "نمایش سه‌بعدی"
    
    tab_gallery, tab_3d_viewer = st.tabs([tab1_title, tab2_title])

    with tab_gallery:
        if st.session_state.isomer_data:
            st.subheader(f"ایزومرهای یافت شده برای: {st.session_state.molecule_searched}")
            num_columns_gallery = 3 
            gallery_cols = st.columns(num_columns_gallery)
            for i, isomer in enumerate(st.session_state.isomer_data):
                with gallery_cols[i % num_columns_gallery]:
                    container = st.container(border=True) 
                    # تغییر: نمایش SVG با markdown
                    if isomer["image_svg"]:
                        container.markdown(isomer["image_svg"], unsafe_allow_html=True)
                    else:
                        container.markdown("<small>تصویر ناموجود</small>", unsafe_allow_html=True)
                    
                    container.markdown(f"**{isomer['name']}**", unsafe_allow_html=True) # نام ایزومر
                    container.markdown(f"<small>SMILES: {isomer['smiles']}<br>CID: {isomer['cid']}</small>", unsafe_allow_html=True)
                    if container.button(f"نمایش سه‌بعدی", key=f"btn_3d_tab_{isomer['cid']}"):
                        st.session_state.selected_cid_for_3d = isomer['cid']
                        st.session_state.selected_name_for_3d = isomer['name'] 
                        st.rerun() 
        elif st.session_state.molecule_searched:
            st.info("ایزومری برای نمایش در گالری وجود ندارد.")
        else:
            st.info("برای مشاهده گالری، ابتدا یک آلکان جستجو کنید.")

    with tab_3d_viewer:
        # ... (بخش نمایش سه‌بعدی بدون تغییر) ...
        if st.session_state.selected_cid_for_3d:
            st.subheader(f"ساختار سه‌بعدی برای: {st.session_state.selected_name_for_3d}")
            with st.spinner(f"در حال بارگذاری ساختار سه‌بعدی برای {st.session_state.selected_name_for_3d} (CID: {st.session_state.selected_cid_for_3d})..."):
                sdf_data, error = get_sdf_content(st.session_state.selected_cid_for_3d)
                if sdf_data:
                    html_3d = generate_3d_viewer_html(sdf_data, st.session_state.selected_name_for_3d, width=600, height=450)
                    st.components.v1.html(html_3d, height=470, width=620, scrolling=False)
            
            if st.button("پاک کردن نمایش سه‌بعدی", key="clear_3d_view_tab"):
                st.session_state.selected_cid_for_3d = None
                st.session_state.selected_name_for_3d = ""
                st.rerun()
        else:
            st.info("برای مشاهده ساختار سه‌بعدی، یک ایزومر را از تب 'گالری ایزومرها' انتخاب کنید.")
else:
    if st.session_state.molecule_searched and not st.session_state.isomer_data:
        pass 
    else:
        st.info("برای شروع، نام یک آلکان را در سایدبار وارد کرده و جستجو کنید یا یک مثال را انتخاب نمایید.")

st.sidebar.markdown("---")
st.sidebar.markdown("ساخته شده با Streamlit, RDKit, PubChemPy, و py3Dmol.")
