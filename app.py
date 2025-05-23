import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage
import py3Dmol
import os
import traceback
from PIL import Image # برای کار با تصاویر در Streamlit

# --- توابع کمکی (بسیاری از آن‌ها از کد قبلی قابل استفاده مجدد هستند) ---

def draw_molecule_pil(smiles_string, size=(250, 250)): # تغییر اندازه برای نمایش بهتر در ستون‌ها
    """یک مولکول را از رشته SMILES با RDKit رسم کرده و یک شیء PIL Image برمی‌گرداند."""
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if mol:
            img = MolToImage(mol, size=size)
            return img
        else:
            st.warning(f"خطا در پارس کردن SMILES: {smiles_string}")
            return None
    except Exception as e:
        st.error(f"خطا در رسم مولکول برای SMILES {smiles_string}: {e}")
        return None

def get_sdf_content(cid):
    """محتوای فایل SDF سه‌بعدی را برای یک CID مشخص از PubChem دریافت می‌کند."""
    if cid is None:
        return None, "CID برای دریافت محتوای SDF ارائه نشده است."
    
    print(f"در حال دریافت محتوای SDF برای CID: {cid}...") # برای لاگ سمت سرور
    # استفاده از یک مسیر موقت که Streamlit به آن دسترسی نوشتن دارد (معمولا نیازی نیست اگر مستقیم محتوا را می‌خوانیم)
    # اما اگر pcp.download نیاز به فایل فیزیکی دارد، این مهم است.
    # در محیط‌هایی مانند Streamlit Cloud، سیستم فایل ممکن است موقتی باشد.
    # برای سادگی، می‌توانیم تلاش کنیم CID را مستقیماً برای گرفتن رکورد استفاده کنیم.
    
    temp_sdf_file_dir = "/tmp" # یک دایرکتوری موقت رایج در لینوکس
    if not os.path.exists(temp_sdf_file_dir):
        try:
            os.makedirs(temp_sdf_file_dir, exist_ok=True)
        except OSError: # اگر ایجاد نشد، در دایرکتوری فعلی تلاش کن (کمتر ایده‌آل)
             temp_sdf_file_dir = "."

    temp_sdf_file = os.path.join(temp_sdf_file_dir, f'temp_3d_structure_{cid}.sdf')
    sdf_data = None
    error_message = None

    try:
        # ابتدا تلاش برای دریافت مستقیم داده‌های SDF (اگر PubChemPy پشتیبانی کند)
        # compounds = pcp.Compound.from_cid(cid)
        # sdf_3d_record = compounds.get_record('3d') # این API ممکن است وجود نداشته باشد یا به این شکل نباشد
        # اگر روش مستقیمی برای گرفتن رشته SDF سه‌بعدی نیست، از دانلود استفاده می‌کنیم.

        pcp.download('SDF', temp_sdf_file, str(cid), 'cid', record_type='3d', overwrite=True)
        with open(temp_sdf_file, 'r') as f:
            sdf_data = f.read()
        
        if not sdf_data or sdf_data.strip() == "$$$$\n" or sdf_data.strip() == "":
            error_message = f"فایل SDF دانلود شده برای CID {cid} خالی است یا ساختار سه‌بعدی ندارد."
            sdf_data = None
        
    except pcp.NotFoundError:
        error_message = f"ساختار سه‌بعدی (SDF) برای CID {cid} در PubChem یافت نشد."
    except FileNotFoundError:
        error_message = f"فایل موقت SDF برای CID {cid} پس از دانلود یافت نشد."
    except Exception as e:
        error_message = f"خطا در دانلود یا خواندن فایل SDF برای CID {cid}: {e}\n{traceback.format_exc()}"
    finally:
        if os.path.exists(temp_sdf_file):
            try:
                os.remove(temp_sdf_file)
            except Exception:
                pass # عدم موفقیت در حذف فایل موقت نباید برنامه را متوقف کند
    
    if error_message:
        st.warning(error_message) # نمایش هشدار در UI
    return sdf_data, error_message


def generate_3d_viewer_html(sdf_data, cid, width=500, height=400):
    """HTML لازم برای نمایشگر py3Dmol را از داده‌های SDF تولید می‌کند."""
    if not sdf_data:
        return "<p style='color:orange;'>داده SDF برای نمایش سه‌بعدی موجود نیست.</p>"
    try:
        viewer = py3Dmol.view(width=width, height=height)
        viewer.addModel(sdf_data, 'sdf')
        viewer.setStyle({'stick': {}})
        viewer.setBackgroundColor('0xeeeeee')
        viewer.zoomTo()
        # viewer.render() # py3Dmol._make_html() خودش این کار را انجام می‌دهد
        html = viewer._make_html()
        return html
    except Exception as e:
        st.error(f"خطا در ایجاد نمایشگر سه‌بعدی برای CID {cid}: {e}")
        return f"<p style='color:red;'>خطا در رندر سه‌بعدی: {e}</p>"

def process_alkane_request(molecule_name_input):
    """
    نام یک آلکان را پردازش می‌کند، ایزومرهای آن را پیدا کرده و لیستی از دیکشنری‌های
    حاوی اطلاعات ایزومرها (شامل تصویر PIL، نام، CID، SMILES) و یک پیام وضعیت برمی‌گرداند.
    """
    if not molecule_name_input or not molecule_name_input.strip():
        return [], "لطفا نام یک مولکول را وارد کنید."

    molecule_name = molecule_name_input.strip().lower()
    status_message = f"در حال جستجو برای '{molecule_name}'..."
    isomer_details_list = []

    try:
        compounds = pcp.get_compounds(molecule_name, 'name', listkey_count=5)
        main_compound_obj = None
        molecular_formula = None

        if not compounds:
            return [], f"مولکول '{molecule_name}' در PubChem یافت نشد."
        
        # ... (بخش انتخاب main_compound_obj و molecular_formula از کد قبلی شما)
        for i, c in enumerate(compounds):
            cid = c.cid
            actual_formula = c.molecular_formula if hasattr(c, 'molecular_formula') else None
            if actual_formula:
                is_standard_hydrocarbon = True 
                if c.canonical_smiles:
                    try:
                        mol_obj = Chem.MolFromSmiles(c.canonical_smiles)
                        if mol_obj: # ادامه فیلترینگ ترکیب اصلی
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
            return [], f"آلکان استاندارد با نام '{molecule_name}' یافت نشد یا با معیارها مطابقت ندارد."

        status_message = f"یافتن ایزومرها برای {main_compound_obj.iupac_name or molecule_name} (فرمول: {molecular_formula})..."
        isomers_found_raw = pcp.get_compounds(molecular_formula, 'formula', listkey_count=50)

        if not isomers_found_raw:
            return [], f"ایزومری برای فرمول {molecular_formula} یافت نشد."

        valid_structural_alkanes_entries = []
        unique_accepted_smiles = set()
        # ... (بخش فیلتر کردن ایزومرها از کد قبلی شما)
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

        if not valid_structural_alkanes_entries:
             return [], f"ایزومر آلکان معتبری برای {molecular_formula} پس از فیلتر کردن یافت نشد."


        for entry in sorted(valid_structural_alkanes_entries, key=lambda x: x.cid): # مرتب‌سازی برای نمایش یکسان
            pil_image = draw_molecule_pil(entry.canonical_smiles)
            if pil_image:
                name = entry.iupac_name
                if not name and entry.synonyms:
                    # منطق ساده‌تر برای انتخاب نام
                    simple_names = [s for s in entry.synonyms if s.lower().endswith("ane") and not any(char.isdigit() for char in s.split('-')[0]) and '-' not in s.split(' ')[0]]
                    if simple_names: name = min(simple_names, key=len)
                    else: name = entry.synonyms[0]
                elif not name:
                    name = f"Alkane (CID: {entry.cid})"
                
                isomer_details_list.append({
                    "cid": entry.cid,
                    "name": name.capitalize(),
                    "smiles": entry.canonical_smiles,
                    "image": pil_image
                })
        
        if isomer_details_list:
            status_message = f"{len(isomer_details_list)} ایزومر ساختاری آلکان برای '{molecule_name}' (فرمول: {molecular_formula}) پیدا شد."
        else:
            status_message = f"ایزومر قابل نمایشی برای '{molecule_name}' یافت نشد."

        return isomer_details_list, status_message

    except pcp.PubChemHTTPError as e:
        return [], f"خطا در ارتباط با PubChem: {e}",
    except Exception as e:
        return [], f"یک خطای غیرمنتظره رخ داد: {e}\n{traceback.format_exc()}"

# --- رابط کاربری Streamlit ---

st.set_page_config(page_title="نمایشگر ایزومرهای آلکان", layout="wide")

st.title("یابنده و نمایشگر ایزومرهای آلکان (Streamlit)")
st.markdown(
    "نام یک آلکان (به انگلیسی) را وارد کنید تا ایزومرهای آن به همراه ساختار شیمیایی دوبعدی، "
    "نام، SMILES و CID نمایش داده شوند. با کلیک بر روی دکمه مربوط به هر ایزومر، "
    "ساختار سه‌بعدی آن (در صورت وجود در PubChem) نمایش داده خواهد شد."
)

# مقداردهی اولیه session state برای نگهداری وضعیت
if 'isomer_data' not in st.session_state:
    st.session_state.isomer_data = []
if 'status_message' not in st.session_state:
    st.session_state.status_message = ""
if 'selected_cid_for_3d' not in st.session_state:
    st.session_state.selected_cid_for_3d = None
if 'molecule_searched' not in st.session_state: # برای جلوگیری از نمایش نتایج قدیمی
    st.session_state.molecule_searched = ""


# فرم ورودی
with st.form(key="alkane_form"):
    molecule_name_input = st.text_input(
        label="نام آلکان را وارد کنید",
        placeholder="مثال: butane, pentane, hexane",
        help="نام آلکان مورد نظر خود را به انگلیسی و با حروف کوچک وارد کنید."
    )
    submit_button = st.form_submit_button(label="جستجوی ایزومرها")

if submit_button and molecule_name_input:
    with st.spinner(f"در حال پردازش برای {molecule_name_input}..."):
        isomers, status_msg = process_alkane_request(molecule_name_input)
        st.session_state.isomer_data = isomers
        st.session_state.status_message = status_msg
        st.session_state.selected_cid_for_3d = None # پاک کردن انتخاب قبلی سه‌بعدی
        st.session_state.molecule_searched = molecule_name_input
elif submit_button and not molecule_name_input:
    st.session_state.status_message = "لطفا نام یک مولکول را وارد کنید."
    st.session_state.isomer_data = []
    st.session_state.selected_cid_for_3d = None
    st.session_state.molecule_searched = ""


# نمایش پیام وضعیت
if st.session_state.status_message:
    # تشخیص نوع پیام بر اساس محتوا (ساده)
    if "خطا" in st.session_state.status_message or "یافت نشد" in st.session_state.status_message and "ایزومر" not in st.session_state.status_message :
        st.error(st.session_state.status_message)
    elif "پیدا شد" in st.session_state.status_message or "ایزومری برای فرمول" in st.session_state.status_message and "یافت نشد" in st.session_state.status_message:
        st.info(st.session_state.status_message)
    else:
        st.info(st.session_state.status_message)


# تقسیم صفحه به دو ستون: گالری در سمت راست، نمایشگر سه‌بعدی در سمت چپ
col_3d, col_gallery = st.columns([2,3]) # ستون سه‌بعدی 40%، گالری 60%

with col_gallery:
    if st.session_state.isomer_data:
        st.subheader(f"ایزومرهای یافت شده برای: {st.session_state.molecule_searched}")
        
        num_columns = 3 # تعداد ستون‌ها برای گالری
        cols = st.columns(num_columns)
        
        for i, isomer in enumerate(st.session_state.isomer_data):
            with cols[i % num_columns]:
                st.image(isomer["image"], caption=f"{isomer['name']} (CID: {isomer['cid']})", use_column_width=True)
                st.markdown(f"<small>SMILES: {isomer['smiles']}</small>", unsafe_allow_html=True)
                if st.button(f"نمایش سه‌بعدی CID: {isomer['cid']}", key=f"btn_3d_{isomer['cid']}"):
                    st.session_state.selected_cid_for_3d = isomer['cid']
                    # با تنظیم selected_cid_for_3d، Streamlit دوباره اجرا می‌شود و بخش col_3d به‌روز خواهد شد.
    elif st.session_state.molecule_searched: # اگر جستجو انجام شده ولی نتیجه‌ای نیست
        if "خطا" not in st.session_state.status_message and "یافت نشد" not in st.session_state.status_message :
             st.info("ایزومری برای نمایش وجود ندارد.")


with col_3d:
    if st.session_state.selected_cid_for_3d:
        st.subheader(f"نمایش سه‌بعدی برای CID: {st.session_state.selected_cid_for_3d}")
        with st.spinner(f"در حال بارگذاری ساختار سه‌بعدی برای CID: {st.session_state.selected_cid_for_3d}..."):
            sdf_data, error = get_sdf_content(st.session_state.selected_cid_for_3d)
            if sdf_data:
                html_3d = generate_3d_viewer_html(sdf_data, st.session_state.selected_cid_for_3d)
                # ارتفاع باید کمی بیشتر از ارتفاع نمایشگر py3Dmol باشد تا اسکرول ایجاد نشود
                st.components.v1.html(html_3d, height=420, width=520, scrolling=False) 
            else:
                # پیام خطا قبلاً توسط get_sdf_content یا generate_3d_viewer_html نمایش داده شده است.
                # st.error(f"ساختار سه‌بعدی برای CID {st.session_state.selected_cid_for_3d} قابل بارگذاری نبود. {error if error else ''}")
                pass
        
        if st.button("پاک کردن نمایش سه‌بعدی", key="clear_3d_view"):
            st.session_state.selected_cid_for_3d = None
            st.experimental_rerun() # برای پاک کردن فوری نمایشگر سه‌بعدی
            
    else:
        st.info("برای مشاهده ساختار سه‌بعدی، یک ایزومر را از گالری انتخاب کنید.")

st.markdown("---")
st.markdown("ساخته شده با Streamlit, RDKit, PubChemPy, و py3Dmol.")

# برای اجرای این فایل: streamlit run app_streamlit.py
