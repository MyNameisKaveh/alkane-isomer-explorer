# app_streamlit.py

import streamlit as st
# توابع ایمپورت شده از فایل‌های دیگر:
from utils_chem import (
    process_alkane_request, 
    process_general_molecule_search,
    get_sdf_content # اگرچه get_sdf_content می‌تواند در utils_ui هم باشد اگر فقط برای UI 3D است
)
from utils_ui import generate_3d_viewer_html

# --- Streamlit UI (بخش UI تقریباً بدون تغییر باقی می‌ماند، فقط فراخوانی توابع تغییر می‌کند) ---
st.set_page_config(page_title="Chemical Compound Explorer", layout="wide", initial_sidebar_state="expanded")
st.title("Chemical Compound Explorer")

# Initialize session state (بدون تغییر)
if 'alkane_isomer_data' not in st.session_state: st.session_state.alkane_isomer_data = []
if 'alkane_main_molecule_props' not in st.session_state: st.session_state.alkane_main_molecule_props = None
if 'selected_isomer_cid_for_3d' not in st.session_state: st.session_state.selected_isomer_cid_for_3d = None
if 'selected_isomer_name_for_3d' not in st.session_state: st.session_state.selected_isomer_name_for_3d = ""
if 'current_isomer_3d_index' not in st.session_state: st.session_state.current_isomer_3d_index = -1
if 'general_molecule_data' not in st.session_state: st.session_state.general_molecule_data = None
if 'general_molecule_name_input' not in st.session_state: st.session_state.general_molecule_name_input = ""
if 'status_message' not in st.session_state: st.session_state.status_message = ""
if 'selected_3d_style' not in st.session_state: st.session_state.selected_3d_style = 'stick' 
if 'last_search_type' not in st.session_state: st.session_state.last_search_type = None
if 'alkane_name_input' not in st.session_state: st.session_state.alkane_name_input = "" 
if 'run_alkane_search_after_example' not in st.session_state: st.session_state.run_alkane_search_after_example = False


# --- Sidebar (بدون تغییر در منطق، فقط فراخوانی توابع متفاوت است اگر آنها را هم جدا کرده بودیم) ---
st.sidebar.header("Search Options")
# Alkane Isomer Search
st.sidebar.subheader("1. Alkane Isomer Search")
current_alkane_input = st.sidebar.text_input(
    label="Enter alkane name (for isomers):", value=st.session_state.alkane_name_input,
    placeholder="e.g., butane, pentane", key="sidebar_alkane_input",
    help="Finds structural isomers for standard alkanes."
)
if current_alkane_input != st.session_state.alkane_name_input:
    st.session_state.alkane_name_input = current_alkane_input
    st.session_state.run_alkane_search_after_example = False 
example_alkanes = ["butane", "pentane", "hexane", "heptane", "octane"]
st.sidebar.caption("Alkane Examples:")
cols_examples_alkane = st.sidebar.columns(len(example_alkanes) if len(example_alkanes) <=3 else 3) 
for i, example in enumerate(example_alkanes):
    if cols_examples_alkane[i % len(cols_examples_alkane)].button(example.capitalize(), key=f"example_alkane_{example}", use_container_width=True):
        st.session_state.alkane_name_input = example 
        st.session_state.run_alkane_search_after_example = True 
        st.rerun() 
if st.sidebar.button("Search Alkane Isomers", type="primary", key="search_alkane_button", use_container_width=True):
    if st.session_state.alkane_name_input:
        st.session_state.last_search_type = "alkane"
        with st.spinner(f"Searching isomers for {st.session_state.alkane_name_input}..."):
            isomers, status_msg, main_props = process_alkane_request(st.session_state.alkane_name_input) # تابع ایمپورت شده
            st.session_state.alkane_isomer_data, st.session_state.status_message, st.session_state.alkane_main_molecule_props = isomers, status_msg, main_props
            st.session_state.selected_isomer_cid_for_3d, st.session_state.selected_isomer_name_for_3d, st.session_state.current_isomer_3d_index = None, "", -1
            st.session_state.alkane_molecule_searched = st.session_state.alkane_name_input
            st.session_state.general_molecule_data = None 
    else: st.session_state.status_message = "Please enter an alkane name for isomer search."
if st.session_state.run_alkane_search_after_example and st.session_state.alkane_name_input:
    st.session_state.last_search_type = "alkane"; st.session_state.run_alkane_search_after_example = False 
    with st.spinner(f"Searching isomers for {st.session_state.alkane_name_input}..."):
        isomers, status_msg, main_props = process_alkane_request(st.session_state.alkane_name_input) # تابع ایمپورت شده
        st.session_state.alkane_isomer_data, st.session_state.status_message, st.session_state.alkane_main_molecule_props = isomers, status_msg, main_props
        st.session_state.selected_isomer_cid_for_3d, st.session_state.selected_isomer_name_for_3d, st.session_state.current_isomer_3d_index = None, "", -1
        st.session_state.alkane_molecule_searched = st.session_state.alkane_name_input
        st.session_state.general_molecule_data = None 
    st.rerun() 
st.sidebar.markdown("---")
# General Molecule Search
st.sidebar.subheader("2. General Molecule Information")
general_molecule_search_term = st.sidebar.text_input(
    label="Enter any chemical name/identifier:", value=st.session_state.general_molecule_name_input,
    placeholder="e.g., aspirin, caffeine, C6H6", key="sidebar_general_input",
    help="Get information, 2D, and 3D structure for any chemical."
)
if general_molecule_search_term != st.session_state.general_molecule_name_input:
    st.session_state.general_molecule_name_input = general_molecule_search_term
if st.sidebar.button("Search Molecule Info", type="primary", key="search_general_button", use_container_width=True):
    if st.session_state.general_molecule_name_input:
        st.session_state.last_search_type = "general"
        with st.spinner(f"Searching for {st.session_state.general_molecule_name_input}..."):
            mol_data, status_msg = process_general_molecule_search(st.session_state.general_molecule_name_input) # تابع ایمپورت شده
            st.session_state.general_molecule_data, st.session_state.status_message = mol_data, status_msg
            st.session_state.alkane_isomer_data, st.session_state.alkane_main_molecule_props = [], None
            st.session_state.selected_isomer_cid_for_3d, st.session_state.selected_isomer_name_for_3d, st.session_state.current_isomer_3d_index = None, "", -1
    else: st.session_state.status_message = "Please enter a chemical name for general search."

if st.session_state.status_message:
    is_error = any(keyword in st.session_state.status_message.lower() for keyword in ["error", "not found", "empty", "invalid"])
    if is_error: st.sidebar.error(st.session_state.status_message)
    else: st.sidebar.info(st.session_state.status_message)

# --- Main Page Tabs (منطق نمایش تب‌ها و محتوای آن‌ها مانند قبل، فقط توابع از ماژول‌ها فراخوانی می‌شوند) ---
gallery_title, props_title, view2d_title, view3d_title = "Isomer Gallery", "Properties", "2D Structure", "3D View"
# ... (بقیه کد UI برای تب‌ها، گالری، نمایش خواص و نمایش سه‌بعدی مانند پاسخ قبلی، با این تفاوت که توابع از ماژول‌ها فراخوانی می‌شوند)
# ... و تابع generate_3d_viewer_html با پارامترهای width و height اصلاح شده فراخوانی می‌شود.

# Content for Alkane Isomer Search
if st.session_state.last_search_type == "alkane":
    # ... (کد کامل نمایش تب‌های آلکان از پاسخ قبلی)
    # مثال برای فراخوانی تابع ایمپورت شده در بخش نمایش سه‌بعدی ایزومر:
    # html_3d = generate_3d_viewer_html(sdf_data, st.session_state.selected_isomer_name_for_3d, 
    #                                   display_style=st.session_state.selected_3d_style, 
    #                                   component_width=620, component_height=470) # پاس دادن ابعاد کامپوننت
    # st.components.v1.html(html_3d, height=470, width=620, scrolling=False)
    if st.session_state.alkane_isomer_data or st.session_state.alkane_main_molecule_props or st.session_state.selected_isomer_cid_for_3d:
        gallery_title = f"Alkane Isomers ({len(st.session_state.alkane_isomer_data)})" if st.session_state.alkane_isomer_data else "Alkane Isomers"
        props_title = f"Alkane Properties: {st.session_state.alkane_main_molecule_props.get('IUPAC Name', '')[:20]}" if st.session_state.alkane_main_molecule_props else "Alkane Properties"
        view3d_title = f"Isomer 3D: {st.session_state.selected_isomer_name_for_3d[:20]}" if st.session_state.selected_isomer_cid_for_3d else "Isomer 3D View"
        
        tab_gallery, tab_main_props, tab_3d_isomer = st.tabs([gallery_title, props_title, view3d_title])

        with tab_gallery: # محتوای تب گالری آلکان
            if st.session_state.alkane_main_molecule_props: # نمایش اول خواص مولکول اصلی آلکان
                props = st.session_state.alkane_main_molecule_props
                main_mol_name = props.get("IUPAC Name", st.session_state.alkane_molecule_searched.capitalize())
                with st.expander(f"Properties for Searched Alkane: {main_mol_name}", expanded=True):
                    prop_cols = st.columns(2)
                    prop_list = list(props.items())
                    for i_prop, (key, value) in enumerate(prop_list):
                        if value != 'N/A' and value is not None:
                            prop_cols[i_prop % 2].markdown(f"**{key}:** {value}")
                    st.markdown("---")
            if st.session_state.alkane_isomer_data:
                st.subheader(f"Isomers found for: {st.session_state.alkane_molecule_searched.capitalize()}")
                num_columns_gallery = 3
                gallery_cols = st.columns(num_columns_gallery)
                for i, isomer in enumerate(st.session_state.alkane_isomer_data):
                    with gallery_cols[i % num_columns_gallery]:
                        container = st.container(border=True) 
                        container.image(isomer["image"], caption=f"{isomer['name']}", use_container_width=True) 
                        container.markdown(f"<small>SMILES: {isomer['smiles']}<br>CID: {isomer['cid']}</small>", unsafe_allow_html=True)
                        if container.button(f"View 3D", key=f"btn_3d_isomer_{isomer['cid']}"):
                            st.session_state.selected_isomer_cid_for_3d = isomer['cid']
                            st.session_state.selected_isomer_name_for_3d = isomer['name'] 
                            st.session_state.current_isomer_3d_index = i 
                            st.rerun()
            elif st.session_state.alkane_molecule_searched: st.info("No isomers to display.")
            else: st.info("Search for an alkane to see isomers.")

        with tab_main_props: # این تب اکنون در بالا در expander نمایش داده می‌شود، اینجا می‌تواند خالی بماند یا برای اطلاعات بیشتر استفاده شود
             if not st.session_state.alkane_main_molecule_props:
                st.info("Search for an alkane to see its properties here or above the gallery.")
             else: # برای جلوگیری از تکرار، می‌توانیم اینجا را خالی بگذاریم چون در expander نشان داده شده
                st.write("") # یا یک پیام دیگر


        with tab_3d_isomer: # محتوای تب سه‌بعدی ایزومر
            if st.session_state.selected_isomer_cid_for_3d:
                st.subheader(f"3D Structure for Isomer: {st.session_state.selected_isomer_name_for_3d}")
                style_options_map = {'Stick': 'stick', 'Line': 'line', 'Ball and Stick': 'ball_and_stick'}
                # ... (کد انتخابگر سبک و دکمه‌های ناوبری مانند قبل) ...
                style_labels = list(style_options_map.keys())
                try:
                    current_style_label = [k for k, v in style_options_map.items() if v == st.session_state.selected_3d_style][0]
                    current_style_index = style_labels.index(current_style_label)
                except IndexError: 
                    current_style_index = 0 
                    if style_labels: st.session_state.selected_3d_style = style_options_map[style_labels[0]]
                    else: st.session_state.selected_3d_style = 'stick'
                selected_style_label = st.radio("Select Display Style:", options=style_labels, key="radio_3d_style_isomer", horizontal=True, index=current_style_index)
                if style_options_map.get(selected_style_label) != st.session_state.selected_3d_style:
                    st.session_state.selected_3d_style = style_options_map[selected_style_label]
                    st.rerun() 
                nav_col_layout = [1, 0.1, 1.5, 0.1, 1]; prev_col, _, clear_col, _, next_col = st.columns(nav_col_layout)
                with prev_col:
                    if st.button("⬅️ Previous", key="prev_3d_isomer", help="View Previous Isomer", use_container_width=True, disabled=(st.session_state.current_isomer_3d_index <= 0)):
                        st.session_state.current_isomer_3d_index -= 1; prev_isomer = st.session_state.alkane_isomer_data[st.session_state.current_isomer_3d_index]
                        st.session_state.selected_isomer_cid_for_3d, st.session_state.selected_isomer_name_for_3d = prev_isomer['cid'], prev_isomer['name']; st.rerun()
                with clear_col:
                    if st.button("Clear Isomer 3D View", key="clear_3d_isomer", use_container_width=True):
                        st.session_state.selected_isomer_cid_for_3d, st.session_state.selected_isomer_name_for_3d, st.session_state.current_isomer_3d_index = None, "", -1; st.rerun()
                with next_col:
                    if st.button("Next ➡️", key="next_3d_isomer", help="View Next Isomer", use_container_width=True, disabled=(st.session_state.current_isomer_3d_index >= len(st.session_state.alkane_isomer_data) - 1)):
                        st.session_state.current_isomer_3d_index += 1; next_isomer = st.session_state.alkane_isomer_data[st.session_state.current_isomer_3d_index]
                        st.session_state.selected_isomer_cid_for_3d, st.session_state.selected_isomer_name_for_3d = next_isomer['cid'], next_isomer['name']; st.rerun()
                st.markdown("---") 
                with st.spinner(f"Loading 3D structure for {st.session_state.selected_isomer_name_for_3d}..."):
                    sdf_data, error = get_sdf_content(st.session_state.selected_isomer_cid_for_3d) # تابع ایمپورت شده
                    if sdf_data:
                        html_3d = generate_3d_viewer_html(sdf_data, st.session_state.selected_isomer_name_for_3d, 
                                                          display_style=st.session_state.selected_3d_style, 
                                                          component_width=620, component_height=470) # تابع ایمپورت شده
                        st.components.v1.html(html_3d, height=470, width=620, scrolling=False)
            else: st.info("Select an isomer from the gallery to view its 3D structure.")


# Content for General Molecule Search
if st.session_state.last_search_type == "general" and st.session_state.general_molecule_data:
    g_data = st.session_state.general_molecule_data
    # تعریف تب‌ها برای جستجوی عمومی
    tabs_general_titles = []
    if g_data["properties"]: tabs_general_titles.append(f"Properties: {g_data['name'][:20]}")
    if g_data.get("image_2d"): tabs_general_titles.append(f"2D: {g_data['name'][:20]}")
    tabs_general_titles.append(f"3D: {g_data['name'][:20]}")
    
    created_tabs_general = st.tabs(tabs_general_titles)
    
    tab_idx = 0
    if g_data["properties"]:
        with created_tabs_general[tab_idx]:
            st.subheader(f"Properties for: {g_data['name']}")
            prop_cols = st.columns(2)
            prop_list = list(g_data["properties"].items())
            for i_prop, (key, value) in enumerate(prop_list):
                if value != 'N/A' and value is not None:
                    prop_cols[i_prop % 2].markdown(f"**{key}:** {value}")
            st.markdown("---")
        tab_idx +=1

    if g_data.get("image_2d"):
        with created_tabs_general[tab_idx]:
            st.subheader(f"2D Structure for: {g_data['name']}")
            st.image(g_data["image_2d"], use_container_width=True)
        tab_idx +=1
            
    with created_tabs_general[tab_idx]: # 3D Tab
        st.subheader(f"3D Structure for: {g_data['name']}")
        # ... (کد انتخابگر سبک برای مولکول عمومی، مشابه بخش ایزومر اما بدون Next/Previous) ...
        style_options_map_g = {'Stick': 'stick', 'Line': 'line', 'Ball and Stick': 'ball_and_stick'}
        style_labels_g = list(style_options_map_g.keys())
        try:
            current_style_label_g = [k for k, v in style_options_map_g.items() if v == st.session_state.selected_3d_style][0]
            current_style_index_g = style_labels_g.index(current_style_label_g)
        except IndexError: 
            current_style_index_g = 0 
            if style_labels_g: st.session_state.selected_3d_style = style_options_map_g[style_labels_g[0]]
            else: st.session_state.selected_3d_style = 'stick'
        selected_style_label_g = st.radio("Select Display Style:", options=style_labels_g, key="radio_3d_style_general", horizontal=True, index=current_style_index_g)
        if style_options_map_g.get(selected_style_label_g) != st.session_state.selected_3d_style:
            st.session_state.selected_3d_style = style_options_map_g[selected_style_label_g]
            st.rerun() 

        with st.spinner(f"Loading 3D structure for {g_data['name']}..."):
            sdf_data, error = get_sdf_content(g_data['cid']) # تابع ایمپورت شده
            if sdf_data:
                html_3d = generate_3d_viewer_html(sdf_data, g_data['name'], 
                                                  display_style=st.session_state.selected_3d_style, 
                                                  component_width=620, component_height=470) # تابع ایمپورت شده
                st.components.v1.html(html_3d, height=470, width=620, scrolling=False)
            else:
                st.info(f"3D structure could not be loaded. {error if error else ''}")
elif st.session_state.last_search_type is None and not st.session_state.alkane_molecule_searched:
     st.info("To begin, use the search options in the sidebar.")


st.sidebar.markdown("---")
st.sidebar.markdown("Built with Streamlit, RDKit, PubChemPy, and py3Dmol.")
