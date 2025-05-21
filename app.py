import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage
from PIL import Image
import io

# تابع برای رسم مولکول از SMILES و بازگرداندن آن به عنوان تصویر PIL
def draw_molecule(smiles_string):
    """
    رشته SMILES را گرفته و تصویر مولکول را به عنوان یک شیء PIL Image برمی‌گرداند.
    اگر SMILES نامعتبر باشد، None برمی‌گرداند.
    """
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if mol:
            # اندازه تصویر را می‌توان تنظیم کرد
            img = MolToImage(mol, size=(300, 300))
            return img
        else:
            print(f"Could not parse SMILES: {smiles_string}")
            return None
    except Exception as e:
        print(f"Error drawing molecule for SMILES {smiles_string}: {e}")
        return None

def find_isomers_data_and_images(molecule_name):
    """
    نام یک مولکول را گرفته، ایزومرهای آن را به همراه نام IUPAC، SMILES
    و تصویر ساختاری (PIL Image) برمی‌گرداند.
    """
    results_with_images = []
    try:
        # 1. جستجوی مولکول اصلی برای یافتن فرمول مولکولی
        print(f"Searching for: {molecule_name}")
        compounds = pcp.get_compounds(molecule_name, 'name', record_type='3d') # record_type='3d' گاهی نتایج بهتری میدهد
        if not compounds:
            print(f"Molecule '{molecule_name}' not found in PubChem.")
            return []
        
        main_compound = compounds[0]
        molecular_formula = main_compound.molecular_formula
        print(f"Found {molecule_name} with formula: {molecular_formula}")

        if not molecular_formula:
            print(f"Could not determine molecular formula for {molecule_name}.")
            return []

        # 2. جستجوی ترکیبات با همان فرمول مولکولی (ایزومرها)
        print(f"Searching for isomers with formula: {molecular_formula}...")
        # برای آلکان‌های ساده، تعداد ایزومرها کم است. برای موارد پیچیده‌تر listkey_count را افزایش دهید.
        isomers = pcp.get_compounds(molecular_formula, 'formula', listkey_count=50) 

        if not isomers:
            print(f"No isomers found for formula {molecular_formula}.")
            return []

        print(f"Found {len(isomers)} potential isomers.")

        for isomer in isomers:
            smiles = isomer.canonical_smiles
            iupac_name = isomer.iupac_name
            
            if not iupac_name and isomer.synonyms:
                iupac_name = isomer.synonyms[0]
            elif not iupac_name:
                iupac_name = "N/A (CID: {})".format(isomer.cid) # از CID برای شناسایی استفاده می‌کنیم

            if smiles:
                mol_image = draw_molecule(smiles)
                if mol_image: # فقط اگر تصویر با موفقیت رسم شد اضافه کن
                    results_with_images.append({
                        "name": iupac_name,
                        "smiles": smiles,
                        "cid": isomer.cid,
                        "image": mol_image
                    })
            # else:
            #     print(f"Warning: No SMILES found for CID {isomer.cid} ({iupac_name})")
        
        # مرتب‌سازی بر اساس نام (اختیاری)
        results_with_images.sort(key=lambda x: x['name'])
        return results_with_images

    except Exception as e:
        print(f"An error occurred: {e}")
        return []

# بخش تست (این بخش در نسخه نهایی با Gradio حذف یا کامنت می‌شود)
if __name__ == "__main__":
    alkane_name = "butane" 
    isomers_info = find_isomers_data_and_images(alkane_name)
    
    if isomers_info:
        print(f"\nIsomers of {alkane_name}:")
        for isomer_data in isomers_info:
            print(f"  Name: {isomer_data['name']}, SMILES: {isomer_data['smiles']}")
            # برای نمایش تصویر در محیطی که کتابخانه تصویر نصب است:
            # isomer_data['image'].show() # این خط در هاگینگ فیس کار نخواهد کرد، فقط برای تست محلی
    else:
        print(f"Could not retrieve isomer information for {alkane_name}.")
