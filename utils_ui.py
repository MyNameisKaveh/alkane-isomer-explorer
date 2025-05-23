# utils_ui.py

import py3Dmol
# import streamlit as st # فقط اگر می‌خواهید مستقیماً st.error/warning را اینجا صدا بزنید، اما بهتر است خطاها به app اصلی برگردانده شوند.

def generate_3d_viewer_html(sdf_data, molecule_name, display_style='stick', component_width=600, component_height=450):
    """
    Generates the HTML for the py3Dmol viewer with adjustments for centering.
    component_width and component_height are the dimensions of the Streamlit HTML component.
    The internal viewer width will be slightly less to provide padding.
    """
    if not sdf_data: 
        return "<p style='color:orange; text-align:center;'>SDF data not available for 3D view.</p>"
    
    try:
        # Make internal viewer slightly smaller than the component to avoid scrollbars or clipping
        internal_viewer_width = component_width - 25 
        internal_viewer_height = component_height - 10

        viewer = py3Dmol.view(width=internal_viewer_width, height=internal_viewer_height)
        viewer.addModel(sdf_data, 'sdf')
        
        if display_style == 'stick': viewer.setStyle({'stick': {}})
        elif display_style == 'sphere': viewer.setStyle({'sphere': {'scale': 0.35}}) # Consider removing if quality is an issue
        elif display_style == 'line': viewer.setStyle({'line': {'linewidth': 2.0, 'colorscheme': 'blackCarbon'}})
        elif display_style == 'ball_and_stick': viewer.setStyle({'stick': {'radius': 0.08}, 'sphere': {'scale': 0.25}})
        else: viewer.setStyle({'stick': {}}) # Default
            
        viewer.setBackgroundColor('0xeeeeee')
        viewer.center() # Attempt to center the molecule
        viewer.zoomTo()
        # viewer.zoom(0.95) # Optional: Slightly zoom out if parts are still clipped
        return viewer._make_html()
    except Exception as e:
        # It's often better to let the main app handle st.error/st.warning
        # so this function remains purely for HTML generation.
        # However, for debugging, you might print here:
        print(f"Error in generate_3d_viewer_html for {molecule_name}: {e}")
        return f"<p style='color:red; text-align:center;'>Error rendering 3D view for {molecule_name}: {str(e)}</p>"
