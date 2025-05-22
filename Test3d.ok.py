


import pubchempy as pcp
import py3Dmol
import os # برای مدیریت فایل‌های موقت

# تابع اصلی برای یافتن و نمایش ساختار سه‌بعدی
def find_and_show_3d_structure(molecule_name):
    """
    ساختار سه‌بعدی یک مولکول را بر اساس نام آن از PubChem پیدا کرده و نمایش می‌دهد.
    """
    print(f"در حال جستجو برای مولکول: '{molecule_name}'...")

    try:
        # مرحله 1: پیدا کردن CID (PubChem Compound ID) بر اساس نام
        # PubChemPy به صورت پیش‌فرض اولین نتیجه را برمی‌گرداند که معمولا بهترین تطابق است.
        compounds = pcp.get_compounds(molecule_name, 'name')

        if not compounds:
            print(f"مولکول '{molecule_name}' در PubChem یافت نشد.")
            return

        compound = compounds[0]
        cid = compound.cid
        print(f"CID برای '{molecule_name}': {cid}")

        # مرحله 2: دانلود فایل SDF سه‌بعدی از PubChem
        # PubChemPy می‌تواند مستقیماً فایل SDF سه‌بعدی را دانلود کند.
        temp_sdf_file = f'temp_3d_structure_{cid}.sdf'
        sdf_content = None

        print(f"در حال دانلود ساختار سه‌بعدی (SDF) برای CID {cid}...")
        try:
            pcp.download('SDF', temp_sdf_file, cid, 'cid', record_type='3d', overwrite=True)

            # خواندن محتوای فایل SDF
            with open(temp_sdf_file, 'r') as f:
                sdf_content = f.read()

            if not sdf_content:
                print("فایل SDF دانلود شده خالی بود. ممکن است ساختار سه‌بعدی در دسترس نباشد.")
                sdf_content = None

        except pcp.NotFoundError:
            print(f"ساختار سه‌بعدی (SDF) برای CID {cid} در PubChem یافت نشد.")
        except Exception as e:
            print(f"خطا در دانلود یا خواندن فایل SDF: {e}")
        finally:
            # پاک کردن فایل موقت پس از استفاده
            if os.path.exists(temp_sdf_file):
                os.remove(temp_sdf_file)

        # مرحله 3: نمایش ساختار سه‌بعدی با py3Dmol
        if sdf_content:
            print("در حال نمایش ساختار سه‌بعدی...")
            viewer = py3Dmol.view(width=500, height=400) # ابعاد نمایشگر
            viewer.addModel(sdf_content, 'sdf') # افزودن مدل از محتوای SDF

            # تعیین سبک نمایش (مثلاً 'stick', 'sphere', 'line', 'cartoon')
            viewer.setStyle({'stick': {}}) # می‌توانید این را تغییر دهید یا از کاربر بگیرید
            # viewer.setStyle({'sphere': {'scale': 0.3}}) # مثال دیگر: نمایش کره (Ball-and-stick)
            # viewer.setStyle({'line': {}}) # مثال دیگر: خطی

            viewer.setBackgroundColor('0xeeeeee') # رنگ پس‌زمینه (خاکستری روشن)
            viewer.zoomTo() # زوم کردن روی مولکول
            viewer.show() # نمایش ویوور در Google Colab
        else:
            print(f"ساختار سه‌بعدی برای '{molecule_name}' قابل بازیابی نبود.")

    except pcp.PubChemHTTPError as e:
        print(f"خطا از PubChem: {e}")
    except Exception as e:
        print(f"یک خطای غیرمنتظره رخ داد: {e}")

# --- بخش اجرای کد ---
if __name__ == "__main__":
    # از کاربر نام مولکول را دریافت کنید
    mol_name = input("لطفا نام مولکول (به انگلیسی) را وارد کنید (مثلا: butane, water, aspirin): ")

    # فراخوانی تابع برای نمایش ساختار
    find_and_show_3d_structure(mol_name)
