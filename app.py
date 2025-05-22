```python
import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import traceback
from PIL import Image

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
    Generates a static 3D visualization of a molecule from a SMILES string using RDKit.
    Returns a PIL Image object for Streamlit.
    """
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if not mol:
            return None
        mol = Chem.AddHs(mol)  # Add hydrogens for proper 3D structure
        AllChem.EmbedMolecule(mol, randomSeed=42)  # Generate 3D coordinates
        AllChem.MMFFOptimizeMolecule(mol)  # Optimize geometry
        
        # Generate a 3D image
        img = Draw.MolToImage(mol, size=(400, 400), kekulize=False, use3D=True)
        return img
    except Exception as e:
        print(f"Error generating 3D view for SMILES {smiles_string}: {e}")
        return None

def find_and_display_isomers(molecule_name_input):
    """
    Finds and displays structural alkane isomers with 2D images and prepares data for 3D view.
    Returns isomer_outputs, status_message, isomer_data.
    """
    if not molecule_name_input or not molecule_name_input.strip():
        return [], "لطفا نام یک مولکول را وارد کنید.", []

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
            return [], status_message, []

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
            return [], status_message, []
        
        print(f"Searching for isomers with formula: {molecular_formula} (up to 200 candidates)...")
        isomers_found_raw = pcp.get_compounds(molecular_formula, 'formula', listkey_count=200)

        if not isomers_found_raw:
            status_message = f"ایزومری برای فرمول {molecular_formula} یافت نشد."
            print(status_message)
            return [], status_message, []

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
        else:
            status_message = (
                f"{len(isomer_outputs_final)} ایزومر ساختاری آلکان برای '{molecule_name_input}' "
                f"(فرمول: {molecular_formula}) پیدا و نمایش داده شد. "
                f"از منوی زیر یک ایزومر را برای مشاهده سه‌بعدی انتخاب کنید."
            )

        isomer_outputs_final.sort(key=lambda x: x[1])
        isomer_data.sort(key=lambda x: x[0])
        return isomer_outputs_final, status_message, isomer_data

    except pcp.PubChemHTTPError as e:
        error_msg = f"خطا در ارتباط با PubChem: {e}. لطفا اتصال اینترنت خود را بررسی کنید."
        print(f"FULL TRACEBACK: {traceback.format_exc()}")
        return [], error_msg, []
    except Exception as e:
        error_msg = f"خطای غیرمنتظره: {e}"
        print(f"FULL TRACEBACK: {traceback.format_exc()}")
        return [], error_msg, []

# --- Streamlit Interface ---
def main():
    st.set_page_config(page_title="یابنده ایزومرهای آلکان", layout="wide")
    
    st.markdown(
        """
        # یابنده و نمایشگر ایزومرهای ساختاری آلکان
        نام یک آلکان (به انگلیسی) را وارد کنید تا ایزومرهای ساختاری آن (تنها شامل کربن و هیدروژن، بدون حلقه، بدون پیوند چندگانه، بدون ایزوتوپ) به همراه ساختار شیمیایی و نام IUPAC (یا رایج) نمایش داده شوند.  
        اطلاعات از دیتابیس PubChem دریافت شده و ساختارها با استفاده از RDKit رسم می‌شوند.
        """
    )

    # Text input for molecule name
    molecule_input = st.text_input(
        label="نام آلکان را وارد کنید",
        placeholder="مثال: butane, pentane, hexane",
        help="نام آلکان را به انگلیسی و با حروف کوچک وارد کنید."
    )

    # Initialize session state for isomers and status
    if 'isomer_outputs' not in st.session_state:
        st.session_state.isomer_outputs = []
    if 'status_message' not in st.session_state:
        st.session_state.status_message = ""
    if 'isomer_data' not in st.session_state:
        st.session_state.isomer_data = []
    if 'selected_isomer' not in st.session_state:
        st.session_state.selected_isomer = None

    # Submit button
    if st.button("جستجوی ایزومرها"):
        with st.spinner("در حال جستجوی ایزومرها..."):
            isomer_outputs, status_message, isomer_data = find_and_display_isomers(molecule_input)
            st.session_state.isomer_outputs = isomer_outputs
            st.session_state.status_message = status_message
            st.session_state.isomer_data = isomer_data
            st.session_state.selected_isomer = isomer_data[0][0] if isomer_data else None

    # Display status message
    if st.session_state.status_message:
        st.text_area("وضعیت و پیام‌ها", st.session_state.status_message, height=100)

    # Display 2D isomers in a grid
    if st.session_state.isomer_outputs:
        st.subheader("ایزومرهای یافت شده (2D)")
        cols = st.columns(3)  # 3 columns for the gallery
        for i, (img, caption) in enumerate(st.session_state.isomer_outputs):
            with cols[i % 3]:
                st.image(img, caption=caption, use_column_width=True)

    # Display 3D viewer with dropdown
    if st.session_state.isomer_data:
        st.subheader("نمایش سه‌بعدی ایزومر")
        selected_isomer = st.selectbox(
            "انتخاب ایزومر برای نمایش سه‌بعدی",
            options=[name for name, _ in st.session_state.isomer_data],
            index=0 if st.session_state.selected_isomer is None else [name for name, _ in st.session_state.isomer_data].index(st.session_state.selected_isomer)
        )
        
        # Update selected isomer in session state
        st.session_state.selected_isomer = selected_isomer
        
        # Find the SMILES for the selected isomer
        selected_smiles = next((smiles for name, smiles in st.session_state.isomer_data if name == selected_isomer), None)
        if selected_smiles:
            three_d_image = generate_3d_view(selected_smiles)
            if three_d_image:
                st.image(three_d_image, caption=f"نمایش سه‌بعدی: {selected_isomer}", use_column_width=False)
            else:
                st.error(f"خطا در تولید نمایش سه‌بعدی برای {selected_isomer}")

    # Examples
    st.markdown("### مثال‌ها")
    examples = ["butane", "pentane", "hexane", "heptane", "octane"]
    for example in examples:
        if st.button(example.capitalize()):
            molecule_input = example
            with st.spinner("در حال جستجوی ایزومرها..."):
                isomer_outputs, status_message, isomer_data = find_and_display_isomers(molecule_input)
                st.session_state.isomer_outputs = isomer_outputs
                st.session_state.status_message = status_message
                st.session_state.isomer_data = isomer_data
                st.session_state.selected_isomer = isomer_data[0][0] if isomer_data else None
                st.experimental_rerun()

if __name__ == "__main__":
    main()
```
