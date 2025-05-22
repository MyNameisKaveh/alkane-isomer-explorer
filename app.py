
import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import traceback
import logging

# تنظیم لاگ برای دیباگ
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def draw_molecule(smiles_string):
    """
    رندر یه تصویر 2D از مولکول با استفاده از SMILES.
    """
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if mol:
            img = Draw.MolToImage(mol, size=(300, 300))
            logger.info(f"تصویر 2D برای SMILES {smiles_string} با موفقیت تولید شد.")
            return img
        else:
            logger.error(f"ناتوانی در تجزیه SMILES: {smiles_string}")
            return None
    except Exception as e:
        logger.error(f"خطا در رندر تصویر 2D برای SMILES {smiles_string}: {e}")
        return None

def generate_3d_view(smiles_string):
    """
    تولید یه تصویر ثابت 3D از مولکول با استفاده از RDKit.
    یه شیء PIL Image برمی‌گردونه برای Streamlit.
    """
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if not mol:
            logger.error(f"SMILES نامعتبر: {smiles_string}")
            return None
        mol = Chem.AddHs(mol)  # اضافه کردن هیدروژن‌ها برای ساختار 3D
        AllChem.EmbedMolecule(mol, randomSeed=42)  # تولید مختصات 3D
        AllChem.MMFFOptimizeMolecule(mol)  # بهینه‌سازی هندسه
        
        # تولید تصویر 3D
        img = Draw.MolToImage(mol, size=(400, 400), kekulize=False, use3D=True)
        logger.info(f"تصویر 3D برای SMILES {smiles_string} با موفقیت تولید شد.")
        return img
    except Exception as e:
        logger.error(f"خطا در تولید تصویر 3D برای SMILES {smiles_string}: {e}")
        return None

def find_and_display_isomers(molecule_name_input):
    """
    پیدا کردن و نمایش ایزومرهای ساختاری آلکان‌ها با تصاویر 2D و آماده‌سازی داده برای نمایش 3D.
    خروجی: isomer_outputs, status_message, isomer_data.
    """
    if not molecule_name_input or not molecule_name_input.strip():
        logger.warning("ورودی نام مولکول خالی است.")
        return [], "لطفا نام یک مولکول را وارد کنید.", []

    molecule_name = molecule_name_input.strip().lower()
    logger.info(f"پردازش درخواست برای: '{molecule_name}'")
    
    status_message = ""
    isomer_data = []  # ذخیره (نام، SMILES) برای منوی کشویی

    try:
        logger.info(f"جستجو برای ترکیب: '{molecule_name}' در PubChem (حداکثر 10 کاندید)...")
        compounds = pcp.get_compounds(molecule_name, 'name', listkey_count=10) 
        
        main_compound_obj = None
        molecular_formula = None

        if not compounds:
            status_message = f"مولکول '{molecule_name}' در PubChem یافت نشد. لطفا املای آن را بررسی کنید."
            logger.error(status_message)
            return [], status_message, []

        logger.info(f"یافتن {len(compounds)} تطابق احتمالی برای '{molecule_name}'. بررسی خواص آلکان استاندارد...")
        
        for i, c in enumerate(compounds):
            cid = c.cid
            common_name = c.synonyms[0] if c.synonyms else "N/A"
            actual_formula = c.molecular_formula if hasattr(c, 'molecular_formula') else None
            
            logger.info(f"بررسی کاندید ترکیب اصلی {i+1}: CID {cid}, نام: '{common_name}', فرمول: '{actual_formula}'")

            is_standard_alkane_candidate = True 

            if not actual_formula or not c.canonical_smiles:
                is_standard_alkane_candidate = False
                logger.warning(f"کاندید اصلی CID {cid} فاقد فرمول یا SMILES است.")
            
            if is_standard_alkane_candidate:
                try:
                    mol_obj = Chem.MolFromSmiles(c.canonical_smiles)
                    if mol_obj:
                        if len(Chem.GetMolFrags(mol_obj)) > 1:
                            logger.warning(f"کاندید اصلی CID {cid} غیرمتصل است.")
                            is_standard_alkane_candidate = False
                        atom_symbols_main = set(atom.GetSymbol() for atom in mol_obj.GetAtoms())
                        if not atom_symbols_main.issubset({'C', 'H'}) or any(atom.GetIsotope() != 0 for atom in mol_obj.GetAtoms()):
                            logger.warning(f"کاندید اصلی CID {cid} فقط شامل CH نیست یا ایزوتوپ دارد.")
                            is_standard_alkane_candidate = False
                        for bond in mol_obj.GetBonds():
                            if bond.GetBondType() != Chem.BondType.SINGLE:
                                logger.warning(f"کاندید اصلی CID {cid} دارای پیوند غیرتکی است.")
                                is_standard_alkane_candidate = False
                                break
                        if Chem.GetSymmSSSR(mol_obj):
                            logger.warning(f"کاندید اصلی CID {cid} دارای حلقه است.")
                            is_standard_alkane_candidate = False
                    else:
                        logger.warning(f"SMILES کاندید اصلی CID {cid} قابل تجزیه نیست.")
                        is_standard_alkane_candidate = False
                except Exception as rdkit_ex:
                    logger.error(f"خطای RDKit برای CID {cid}: {rdkit_ex}")
                    is_standard_alkane_candidate = False
            
            if is_standard_alkane_candidate:
                current_compound_name_matches_input = molecule_name in [syn.lower() for syn in c.synonyms]
                if current_compound_name_matches_input:
                    main_compound_obj = c
                    molecular_formula = actual_formula
                    logger.info(f"انتخاب ترکیب اصلی: CID {main_compound_obj.cid}, فرمول: {molecular_formula}")
                    break
                if not main_compound_obj:
                    main_compound_obj = c
                    molecular_formula = actual_formula
                    logger.info(f"انتخاب موقت: CID {main_compound_obj.cid}, فرمول: {molecular_formula}")
        
        if not main_compound_obj or not molecular_formula:
            status_message = f"آلکان ساختاری استاندارد با نام '{molecule_name}' یافت نشد."
            logger.error(status_message)
            return [], status_message, []
        
        logger.info(f"جستجو برای ایزومرها با فرمول: {molecular_formula} (حداکثر 200 کاندید)...")
        isomers_found_raw = pcp.get_compounds(molecular_formula, 'formula', listkey_count=200)

        if not isomers_found_raw:
            status_message = f"ایزومری برای فرمول {molecular_formula} یافت نشد."
            logger.error(status_message)
            return [], status_message, []

        logger.info(f"یافتن {len(isomers_found_raw)} ورودی ایزومر احتمالی. فیلتر کردن برای ایزومرهای آلکان معتبر...")
        
        valid_structural_alkanes_entries = []
        unique_accepted_smiles = set()

        for isomer_entry in isomers_found_raw:
            smiles = isomer_entry.canonical_smiles
            if not smiles:
                logger.warning(f"رد کردن ایزومر بدون SMILES: CID {isomer_entry.cid}")
                continue

            try:
                mol_iso = Chem.MolFromSmiles(smiles)
                if not mol_iso:
                    logger.warning(f"فیلتر شده (SMILES نامعتبر): CID {isomer_entry.cid}, SMILES: {smiles}")
                    continue

                is_valid_alkane_isomer = True
                if len(Chem.GetMolFrags(mol_iso)) > 1:
                    logger.warning(f"فیلتر شده (غیرمتصل): CID {isomer_entry.cid}, SMILES: {smiles}")
                    is_valid_alkane_isomer = False
                atom_symbols = set(atom.GetSymbol() for atom in mol_iso.GetAtoms())
                if not atom_symbols.issubset({'C', 'H'}) or any(atom.GetIsotope() != 0 for atom in mol_iso.GetAtoms()):
                    logger.warning(f"فیلتر شده (غیر CH یا ایزوتوپ): CID {isomer_entry.cid}, SMILES: {smiles}")
                    is_valid_alkane_isomer = False
                for bond in mol_iso.GetBonds():
                    if bond.GetBondType() != Chem.BondType.SINGLE:
                        logger.warning(f"فیلتر شده (پیوند غیرتکی): CID {isomer_entry.cid}, SMILES: {smiles}")
                        is_valid_alkane_isomer = False
                        break
                if Chem.GetSymmSSSR(mol_iso):
                    logger.warning(f"فیلتر شده (دارای حلقه): CID {isomer_entry.cid}, SMILES: {smiles}")
                    is_valid_alkane_isomer = False

                if is_valid_alkane_isomer:
                    canonical_smiles = Chem.MolToSmiles(mol_iso, isomericSmiles=False)
                    if canonical_smiles not in unique_accepted_smiles:
                        logger.info(f"پذیرفته شده: CID {isomer_entry.cid}, SMILES: {smiles}")
                        valid_structural_alkanes_entries.append(isomer_entry)
                        unique_accepted_smiles.add(canonical_smiles)
                    else:
                        logger.warning(f"رد کردن (SMILES تکراری): CID {isomer_entry.cid}, SMILES: {smiles}")
            except Exception as rdkit_iso_ex:
                logger.error(f"خطای RDKit برای CID {isomer_entry.cid}: {rdkit_iso_ex}, SMILES: {smiles}")
                continue
        
        logger.info(f"یافتن {len(valid_structural_alkanes_entries)} ایزومر معتبر یکتا.")
        
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
                logger.warning(f"ناتوانی در رندر تصویر برای CID {final_isomer_entry.cid}, SMILES: {smiles_to_draw}")

        logger.info(f"نمایش {valid_isomers_count_final} ایزومر در گالری.")

        if not isomer_outputs_final:
            status_message = "ایزومر ساختاری آلکان استاندارد و قابل رسمی پیدا نشد."
            if len(valid_structural_alkanes_entries) > 0:
                status_message += " (برخی ایزومرها در رسم ناموفق بودند.)"
            logger.error(status_message)
        else:
            status_message = (
                f"{len(isomer_outputs_final)} ایزومر ساختاری آلکان برای '{molecule_name_input}' "
                f"(فرمول: {molecular_formula}) پیدا و نمایش داده شد. "
                f"از منوی زیر یک ایزومر را برای مشاهده سه‌بعدی انتخاب کنید."
            )
            logger.info(status_message)

        isomer_outputs_final.sort(key=lambda x: x[1])
        isomer_data.sort(key=lambda x: x[0])
        return isomer_outputs_final, status_message, isomer_data

    except pcp.PubChemHTTPError as e:
        error_msg = f"خطا در ارتباط با PubChem: {e}. لطفا اتصال اینترنت خود را بررسی کنید."
        logger.error(f"خطای PubChem: {traceback.format_exc()}")
        return [], error_msg, []
    except Exception as e:
        error_msg = f"خطای غیرمنتظره: {e}"
        logger.error(f"خطای غیرمنتظره: {traceback.format_exc()}")
        return [], error_msg, []

# --- رابط Streamlit ---
def main():
    logger.info("شروع اپلیکیشن Streamlit...")
    st.set_page_config(page_title="یابنده ایزومرهای آلکان", layout="wide")
    
    st.markdown(
        """
        # یابنده و نمایشگر ایزومرهای ساختاری آلکان
        نام یک آلکان (به انگلیسی) را وارد کنید تا ایزومرهای ساختاری آن (تنها شامل کربن و هیدروژن، بدون حلقه، بدون پیوند چندگانه، بدون ایزوتوپ) به همراه ساختار شیمیایی و نام IUPAC (یا رایج) نمایش داده شوند.  
        اطلاعات از دیتابیس PubChem دریافت شده و ساختارها با استفاده از RDKit رسم می‌شوند.
        """
    )

    # ورودی متن برای نام مولکول
    molecule_input = st.text_input(
        label="نام آلکان را وارد کنید",
        placeholder="مثال: butane, pentane, hexane",
        help="نام آلکان را به انگلیسی و با حروف کوچک وارد کنید."
    )

    # مقداردهی اولیه وضعیت جلسه
    if 'isomer_outputs' not in st.session_state:
        st.session_state.isomer_outputs = []
    if 'status_message' not in st.session_state:
        st.session_state.status_message = ""
    if 'isomer_data' not in st.session_state:
        st.session_state.isomer_data = []
    if 'selected_isomer' not in st.session_state:
        st.session_state.selected_isomer = None

    # دکمه ارسال
    if st.button("جستجوی ایزومرها"):
        logger.info("دکمه جستجو کلیک شد.")
        with st.spinner("در حال جستجوی ایزومرها..."):
            isomer_outputs, status_message, isomer_data = find_and_display_isomers(molecule_input)
            st.session_state.isomer_outputs = isomer_outputs
            st.session_state.status_message = status_message
            st.session_state.isomer_data = isomer_data
            st.session_state.selected_isomer = isomer_data[0][0] if isomer_data else None
            logger.info("جستجوی ایزومرها کامل شد.")

    # نمایش پیام وضعیت
    if st.session_state.status_message:
        st.text_area("وضعیت و پیام‌ها", st.session_state.status_message, height=100)

    # نمایش ایزومرهای 2D در یه شبکه
    if st.session_state.isomer_outputs:
        st.subheader("ایزومرهای یافت شده (2D)")
        cols = st.columns(3)  # 3 ستون برای گالری
        for i, (img, caption) in enumerate(st.session_state.isomer_outputs):
            with cols[i % 3]:
                st.image(img, caption=caption, use_column_width=True)

    # نمایش نمای 3D با منوی کشویی
    if st.session_state.isomer_data:
        st.subheader("نمایش سه‌بعدی ایزومر")
        selected_isomer = st.selectbox(
            "انتخاب ایزومر برای نمایش سه‌بعدی",
            options=[name for name, _ in st.session_state.isomer_data],
            index=0 if st.session_state.selected_isomer is None else [name for name, _ in st.session_state.isomer_data].index(st.session_state.selected_isomer)
        )
        
        # به‌روزرسانی ایزومر انتخاب‌شده در وضعیت جلسه
        st.session_state.selected_isomer = selected_isomer
        
        # پیدا کردن SMILES برای ایزومر انتخاب‌شده
        selected_smiles = next((smiles for name, smiles in st.session_state.isomer_data if name == selected_isomer), None)
        if selected_smiles:
            logger.info(f"تولید نمای 3D برای ایزومر: {selected_isomer}, SMILES: {selected_smiles}")
            three_d_image = generate_3d_view(selected_smiles)
            if three_d_image:
                st.image(three_d_image, caption=f"نمایش سه‌بعدی: {selected_isomer}", use_column_width=False)
            else:
                st.error(f"خطا در تولید نمایش سه‌بعدی برای {selected_isomer}")
                logger.error(f"ناتوانی در تولید نمای 3D برای {selected_isomer}")

    # مثال‌ها
    st.markdown("### مثال‌ها")
    examples = ["butane", "pentane", "hexane", "heptane", "octane"]
    for example in examples:
        if st.button(example.capitalize()):
            logger.info(f"دکمه مثال {example} کلیک شد.")
            molecule_input = example
            with st.spinner("در حال جستجوی ایزومرها..."):
                isomer_outputs, status_message, isomer_data = find_and_display_isomers(molecule_input)
                st.session_state.isomer_outputs = isomer_outputs
                st.session_state.status_message = status_message
                st.session_state.isomer_data = isomer_data
                st.session_state.selected_isomer = isomer_data[0][0] if isomer_data else None
                st.experimental_rerun()

if __name__ == "__main__":
    logger.info("اجرای اسکریپت app.py شروع شد.")
    main()
```
