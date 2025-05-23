def generate_3d_viewer_html(sdf_data, molecule_name, display_style='stick', component_width=600, component_height=450):
    """Generates the HTML for the py3Dmol viewer with adjustments for centering."""
    if not sdf_data: return "<p style='color:orange; text-align:center;'>SDF data not available for 3D view.</p>"
    try:
        # عرض داخلی py3Dmol را همان عرض کامپوننت در نظر بگیرید
        internal_viewer_width = component_width
        internal_viewer_height = component_height 

        viewer = py3Dmol.view(width=internal_viewer_width, height=internal_viewer_height)
        viewer.addModel(sdf_data, 'sdf')
        
        if display_style == 'stick': viewer.setStyle({'stick': {}})
        elif display_style == 'sphere': viewer.setStyle({'sphere': {'scale': 0.35}}) 
        elif display_style == 'line': viewer.setStyle({'line': {'linewidth': 2.0, 'colorscheme': 'blackCarbon'}})
        elif display_style == 'ball_and_stick': viewer.setStyle({'stick': {'radius': 0.08}, 'sphere': {'scale': 0.25}})
        else: viewer.setStyle({'stick': {}})
            
        viewer.setBackgroundColor('0xeeeeee')
        viewer.center() 
        viewer.zoomTo()
        viewer.zoom(0.8) # کاهش بیشتر زوم (از 0.9 به 0.8)
        return viewer._make_html()
    except Exception as e:
        error_msg_html = f"<p style='color:red; text-align:center;'>Error rendering 3D view for {molecule_name}: {str(e)}</p>"
        if 'st' in globals() and hasattr(st, 'error'): 
            st.error(f"Error creating 3D viewer for {molecule_name} with style {display_style}: {e}")
        return error_msg_html
