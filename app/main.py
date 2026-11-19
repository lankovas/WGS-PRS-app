import streamlit as st
import os
import time
import pandas as pd
import numpy as np  
import re 
import subprocess
from pathlib import Path
import plotly.graph_objects as go
import scipy.stats as stats
import gzip
import shutil
import streamlit.components.v1 as components

from pgs_fetcher import fetch_pgs_panel

if 'panel_ready' not in st.session_state: st.session_state.panel_ready = False
if 'panel_info' not in st.session_state: st.session_state.panel_info = {}
if 'af_mode' not in st.session_state: st.session_state.af_mode = "1000G"
if 'variant_actions' not in st.session_state: st.session_state.variant_actions = pd.DataFrame()
if 'edited_panel_df' not in st.session_state: st.session_state.edited_panel_df = pd.DataFrame()

st.set_page_config(page_title="Clinical PRS Pipeline", layout="wide")

st.markdown(
    """
    <style>
    [data-testid="stNumberInputStepUp"] {display: none;}
    [data-testid="stNumberInputStepDown"] {display: none;}
    </style>
    """,
    unsafe_allow_html=True,
)

def render_dashboard(prs_df, sample_name, final_ancestry_label, somalier_df=None, html_content=None, base_name=None):
    st.markdown("---")
    st.subheader(f"📊 Polygenic Risk Score Results ({sample_name})")
    
    hidden_metrics = ['Population Mean', 'Population SD', 'SD of Raw PRS', 'SD of Z-Score', 
                      'Z-Score Lower 95% Uncertainty Interval', 'Z-Score Upper 95% Uncertainty Interval', 
                      'Z-Score Deterministic Minimum', 'Z-Score Deterministic Maximum',
                      'Z2 Predicted Mean', 'Z2 Predicted SD']
    display_df = prs_df[~prs_df['Metric'].isin(hidden_metrics)]
    
    def highlight_alternating(row):
        return ['background-color: rgba(128, 128, 128, 0.1)' if row.name % 2 == 0 else '' for _ in row]
    
    def format_decimals(val):
        try:
            # Return as a float! PyArrow will crash if we force numeric columns into strings
            return round(float(val), 3)
        except (ValueError, TypeError):
            # If it's a string/text, leave it exactly as is
            return val

    # Best practice: Apply the decimal formatting specifically to the 'Value' column
    styled_df = display_df.style.apply(highlight_alternating, axis=1).format({'Value': format_decimals})
    
    st.dataframe(styled_df, width='stretch', hide_index=True)
    
    csv = prs_df.to_csv(index=False).encode('utf-8')
    col_btn, _ = st.columns([0.3, 0.7])
    with col_btn:
        st.download_button("⬇️ Download PRS Results (.csv)", data=csv, file_name=f"{sample_name}_prs_results.csv", mime="text/csv")
    
    def get_metric(df, metric_name):
        try: return float(df.loc[df['Metric'] == metric_name, 'Value'].iloc[0])
        except: return None
            
    z_score = get_metric(prs_df, 'Categorical Z-Score (Z1)')
    z2_score = get_metric(prs_df, 'PC-Adjusted Z-Score (Z2)')
    z_ci_lower = get_metric(prs_df, 'Z-Score Lower 95% Uncertainty Interval')
    z_ci_upper = get_metric(prs_df, 'Z-Score Upper 95% Uncertainty Interval')
    z_min = get_metric(prs_df, 'Z-Score Deterministic Minimum')
    z_max = get_metric(prs_df, 'Z-Score Deterministic Maximum')

    if z_score is not None:
        st.markdown("---")
        st.subheader(f"Contextualized Population Mapping")
        
        col1, col2, col3, col4, col5 = st.columns(5)
        col1.metric("Deterministic Min", f"{z_min:.3f}" if z_min is not None else "N/A")
        col2.metric("95% UI Lower", f"{z_ci_lower:.3f}" if z_ci_lower is not None else "N/A")
        col3.metric("Categorical Z-Score (Z1)", f"{z_score:.3f}")
        col4.metric("95% UI Upper", f"{z_ci_upper:.3f}" if z_ci_upper is not None else "N/A")
        col5.metric("Deterministic Max", f"{z_max:.3f}" if z_max is not None else "N/A")
        
        if z2_score is not None:
            st.markdown("---")
            col_z2_1, col_z2_2 = st.columns([0.3, 0.7])
            with col_z2_1:
                st.metric("PC-Adjusted Z-Score (Z2)", f"{z2_score:.3f}")
            with col_z2_2:
                st.info("💡 **WGS Admixture Adjustment:** This score uses continuous Principal Component (PC) regression to precisely adjust for genetic admixture, capturing complex genetic backgrounds better than categorical labels.")
        
        min_x = min(-4, z_min) if z_min else -4
        max_x = max(4, z_max) if z_max else 4
        
        x_axis = np.linspace(min_x - 0.5, max_x + 0.5, 1000)
        y_axis = stats.norm.pdf(x_axis, 0, 1)
        
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=x_axis, y=y_axis, mode='lines', name='Population Distribution', line=dict(color='rgba(55,126,184,1)', width=2)))
        
        if z_ci_lower is not None and z_ci_upper is not None:
            ci_x = np.linspace(z_ci_lower, z_ci_upper, 500)
            ci_y = stats.norm.pdf(ci_x, 0, 1)
            fig.add_trace(go.Scatter(x=ci_x, y=ci_y, mode='none', fill='tozeroy', fillcolor='rgba(220,20,60,0.3)', name='Patient 95% CI'))
        
        fig.add_trace(go.Scatter(x=[0, 0], y=[0, stats.norm.pdf(0, 0, 1)], mode='lines', name="Population Mean (0)", line=dict(color='black', width=1,)))
        fig.add_trace(go.Scatter(x=[z_score, z_score], y=[0, stats.norm.pdf(z_score, 0, 1)], mode='lines', name='Patient Z1-Score', line=dict(color='crimson', width=3, dash='dash')))
        
        if z2_score is not None:
            fig.add_trace(go.Scatter(x=[z2_score, z2_score], y=[0, stats.norm.pdf(z2_score, 0, 1)], mode='lines', name='Patient Z2-Score (PC Adj)', line=dict(color='orange', width=3, dash='dashdot')))
        
        if z_min is not None:
            fig.add_trace(go.Scatter(x=[z_min, z_min], y=[0, stats.norm.pdf(z_min, 0, 1)], mode='lines', name='Deterministic Min', line=dict(color='gray', width=2, dash='dot')))
        if z_max is not None:
            fig.add_trace(go.Scatter(x=[z_max, z_max], y=[0, stats.norm.pdf(z_max, 0, 1)], mode='lines', name='Deterministic Max', line=dict(color='gray', width=2, dash='dot')))
        
        fig.update_layout(title=f"Normalized Z-Score Distribution for {final_ancestry_label} Ancestry", xaxis_title="Normalized Z-Score", yaxis_title="Probability Density", hovermode="x unified", margin=dict(t=50, b=40, l=40, r=40))
        
        st.plotly_chart(fig, width='stretch')

    if somalier_df is not None:
        st.markdown("---")
        st.subheader("Somalier Ancestry Prediction")
        
        ancestry_probs = {}
        sample_col = next((col for col in ['sample_id', '#sample_id', 'sample', 'Sample'] if col in somalier_df.columns), None)
        
        if base_name and sample_col and somalier_df[sample_col].str.contains(base_name, na=False).any():
            patient_row = somalier_df[somalier_df[sample_col].str.contains(base_name, na=False)].tail(1)
        else:
            patient_row = somalier_df.tail(1)
            
        if not patient_row.empty:
            prob_cols = [col for col in patient_row.columns if col.endswith('_prob')]
            for col in prob_cols:
                ancestry = col.replace('_prob', '')
                probability = patient_row[col].iloc[0]
                ancestry_probs[ancestry] = probability

        if ancestry_probs:
            st.markdown(f"**Final Predicted Ancestry:** `{final_ancestry_label}`")
            prob_df = pd.DataFrame(list(ancestry_probs.items()), columns=['Ancestry', 'Probability'])
            prob_df['Percentage'] = (prob_df['Probability'] * 100).round(2)
            prob_df = prob_df.sort_values(by='Probability', ascending=False)
            
            col_chart, col_table = st.columns([0.6, 0.4])
            with col_chart:
                fig_prob = go.Figure(data=[go.Bar(x=prob_df['Ancestry'], y=prob_df['Percentage'], text=prob_df['Percentage'].astype(str) + '%', textposition='auto')])
                fig_prob.update_layout(title="Prediction Confidence", yaxis_title="Confidence (%)", yaxis_range=[0, 100], margin=dict(t=30, b=0, l=0, r=0))
                
                st.plotly_chart(fig_prob, width='stretch')
                
            with col_table:
                st.write("Raw Probabilities")
                st.dataframe(prob_df[['Ancestry', 'Percentage']].style.format({'Percentage': '{:.2f}%'}), hide_index=True, width='stretch')

    if html_content is not None:
        with st.expander("Interactive Ancestry PCA Map", expanded=True):
            html_content = html_content.replace("showlegend: false", "showlegend: true")
            html_content = html_content.replace('symbol: "square", size:size}', 'symbol: "star", size: 18, color: "black", line: {color: "white", width: 2}}')
            components.html(html_content, height=800, scrolling=True)

INPUTS_DIR = Path("../Pipeline/inputs")

def extract_trait_from_panel(file_path):
    trait = "Unknown Trait"
    try:
        with gzip.open(file_path, "rt", encoding="utf-8", errors="ignore") as f:
            for _ in range(50): 
                line = f.readline()
                if not line or not line.startswith("#"): break
                if line.startswith("#trait_reported=") or line.startswith("##trait_reported="):
                    trait = line.split("=", 1)[1].strip()
                    break
                elif line.startswith("#trait_mapped=") or line.startswith("##trait_mapped="):
                    trait = line.split("=", 1)[1].strip()
    except Exception: pass
    return trait

def panel_has_embedded_af(file_path):
    if not file_path or not os.path.exists(file_path): return False
    try:
        open_func = gzip.open if str(file_path).endswith('.gz') else open
        with open_func(file_path, "rt", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith("#"): continue
                cols = [c.lower().strip() for c in re.split(r'\t|,', line)]
                return "allelefrequency_effect" in cols
    except Exception: pass
    return False

@st.cache_data(ttl=60)
def get_panel_cache_status(inputs_dir=INPUTS_DIR):
    panels = {}
    if not inputs_dir.exists(): return panels
    
    ref_cache_dir = inputs_dir / "reference_cache"

    for pgs_path in inputs_dir.glob("*.txt.gz"):
        match = re.search(r"((?:PGS\d{6,})|(?:CUSTOM_[A-Za-z0-9_]+))", pgs_path.name)
        if match:
            panel_id = match.group(1)
            trait_name = extract_trait_from_panel(pgs_path)
            
            # 1. Check for calibration.done in reference_cache
            ready_modes = []
            panel_ref_dir = ref_cache_dir / panel_id
            if panel_ref_dir.exists() and panel_ref_dir.is_dir():
                for mode_dir in panel_ref_dir.iterdir():
                    if mode_dir.is_dir() and (mode_dir / "calibration.done").exists():
                        ready_modes.append(mode_dir.name)
            
            # 2. Check for AF cache CSV files (any file matching the pattern)
            af_cache_files = list(inputs_dir.glob(f"*{panel_id}*_af_cache_*.csv"))
            
            if ready_modes:
                has_af = True
                status_label = "✅ Ready"
            elif af_cache_files:
                has_af = True
                status_label = "⚠️ Partially Ready"
            else:
                has_af = False
                status_label = "❌ Cache Missing"

            panels[panel_id] = {
                "trait": trait_name if trait_name != "Unknown Trait" else "Custom Panel",
                "pgs_file": str(pgs_path.resolve()),
                "has_af": has_af,
                "status_label": status_label,
                "ready_modes": ready_modes,
                "af_cache_files": [f.name for f in af_cache_files]
            }
    return panels

st.sidebar.title("Application Mode")
app_mode = st.sidebar.radio("Choose Mode:", ["Run New Pipeline", "Visualize Existing Results"])
st.sidebar.markdown("---")

if app_mode == "Run New Pipeline":
    st.sidebar.title("PRS Pipeline Config")

    # --- SECTION 1: PANEL SELECTION ---
    st.sidebar.markdown(
        "## 1. Panel Selection <a href='https://www.pgscatalog.org/' target='_blank' style='font-size: 0.6em; text-decoration: none;' title='Search the PGS Catalog'>ⓘ</a>", 
        unsafe_allow_html=True
    )
    panel_records = get_panel_cache_status()
    
    # Build dropdown options
    dropdown_options = []
    id_map = {}
    for pid, data in sorted(panel_records.items()):
        # Use the status label generated by our updated function
        display_label = f"{pid}\n{data['trait']} ({data['status_label']})"
        dropdown_options.append(display_label)
        id_map[display_label] = pid

    dropdown_options.append("➕ Enter new PGS ID...")
    dropdown_options.append("⬆️ Upload Custom Panel...")

    # Find the index of the currently active panel (the one we just saved)
    current_pgs = st.session_state.panel_info.get('pgs_id', '')
    default_idx = 0 if panel_records else len(dropdown_options) - 2
    for i, label in enumerate(dropdown_options):
        if current_pgs == id_map.get(label, ""):
            default_idx = i
            break

    selected_label = st.sidebar.selectbox("Select or Search Panel:", options=dropdown_options, index=default_idx)

    custom_panel_uploaded = False
    active_pgs_id = ""

    if selected_label == "➕ Enter new PGS ID...":
        active_pgs_id = st.sidebar.text_input("Enter PGS ID (e.g., PGS004687)").strip().upper()
        if st.sidebar.button("🔍 Preview Panel Info"):
            if not active_pgs_id:
                st.sidebar.warning("Please enter a PGS ID.")
            else:
                with st.spinner(f"Fetching metadata for {active_pgs_id} from PGS Catalog..."):
                    success, result = fetch_pgs_panel(active_pgs_id)
                    if success:
                        st.session_state.panel_ready = True
                        st.session_state.panel_info = result['info']
                        st.session_state.pgs_file = result['file_path'] 
                    else:
                        st.sidebar.error("Failed to fetch panel metadata.")
                        st.session_state.panel_ready = False
                        
    elif selected_label == "⬆️ Upload Custom Panel...":
        custom_panel_name = st.sidebar.text_input("Custom Panel Name:", placeholder="e.g., BRCA_Custom").strip()
        custom_panel_file = st.sidebar.file_uploader("Upload Panel Weights", type=["csv", "tsv", "txt"])
        
        if custom_panel_name and custom_panel_file:
            try:
                df_custom = pd.read_csv(custom_panel_file, sep=None, engine='python')
                cols_lower = {c.lower().strip('\ufeff').strip() for c in df_custom.columns}
                
                has_weight = 'effect_weight' in cols_lower or 'or' in cols_lower or 'hr' in cols_lower
                required_cols = {'chr', 'pos', 'ref', 'alt'}
                
                if not required_cols.issubset(cols_lower) or not has_weight:
                    st.sidebar.error(f"❌ Invalid Format. Missing required coordinates or weight (effect_weight, OR, or HR).")
                else:
                    active_pgs_id = f"CUSTOM_{custom_panel_name}"
                    custom_panel_uploaded = True
                    st.session_state.panel_ready = True
                    
                    input_dir = Path("../Pipeline/inputs")
                    input_dir.mkdir(parents=True, exist_ok=True)
                    custom_file_path = input_dir / f"{active_pgs_id}.txt.gz"
                    
                    df_custom.to_csv(custom_file_path, sep='\t', index=False, compression='gzip')
                    st.session_state.panel_info = {'pgs_id': active_pgs_id, 'trait': f"Custom Upload", 'variants': len(df_custom)}
                    st.session_state.pgs_file = str(custom_file_path)
            except Exception as e:
                st.sidebar.error(f"❌ Error reading panel file: {e}")
    else:
        active_pgs_id = id_map.get(selected_label, "")
        if active_pgs_id and not custom_panel_uploaded:
            if active_pgs_id.startswith("CUSTOM_"):
                st.session_state.panel_ready = True
                st.session_state.panel_info = {'pgs_id': active_pgs_id, 'trait': panel_records[active_pgs_id]['trait'], 'variants': "Unknown"}
                st.session_state.pgs_file = panel_records[active_pgs_id]['pgs_file']
            elif st.session_state.panel_info.get('pgs_id') != active_pgs_id:
                with st.spinner(f"Loading metadata for {active_pgs_id}..."):
                    success, result = fetch_pgs_panel(active_pgs_id)
                    if success:
                        st.session_state.panel_ready = True
                        st.session_state.panel_info = result['info']
                        st.session_state.pgs_file = result['file_path']
                    else:
                        st.sidebar.error(f"Failed to fetch metadata.")
                        st.session_state.panel_ready = False

    # --- PANEL INFO DISPLAY ---
    # --- PANEL INFO DISPLAY ---
    if st.session_state.panel_ready:
        info = st.session_state.panel_info
        st.sidebar.success(f"✅ Panel {info.get('pgs_id')} loaded!")
        
        # --- Display cached modes ---
        pid = info.get('pgs_id')
        if pid in panel_records:
            cached_data = panel_records[pid]
            ready_modes = cached_data.get("ready_modes", [])
            af_csvs = cached_data.get("af_cache_files", [])
            
            if ready_modes or af_csvs:
                with st.sidebar.expander(f"📦 View Cached Modes ({len(ready_modes)} ready)"):
                    if ready_modes:
                        st.markdown("**Fully Calibrated Modes:**")
                        for m in ready_modes:
                            st.markdown(f"- `{m}`")
                    if af_csvs:
                        st.markdown("**Available AF Cache Files:**")
                        for c in af_csvs:
                            st.markdown(f"- `{c}`")
        
        # Only show advanced PGS catalog tabs if it's not a custom upload
        if not str(info.get('pgs_id')).startswith("CUSTOM_"):
            st.sidebar.markdown(f"### Panel Details <a href='https://www.pgscatalog.org/score/{info.get('pgs_id')}/' target='_blank' style='font-size: 0.6em; text-decoration: none;'>ⓘ</a>", unsafe_allow_html=True)
            st.sidebar.markdown(f"**Trait:** {info.get('trait')}\n**Variants:** {info.get('variants')}")
            tab1, tab2, tab3 = st.sidebar.tabs(["Discovery", "Training", "Evaluation"])
            with tab1: st.markdown(f"{info.get('gwas_ancestry')}")
            with tab2: st.markdown(f"{info.get('training_ancestry')}")
            with tab3: st.markdown(f"{info.get('evaluation_ancestry')}")

    actions_df = pd.DataFrame() 
    assume_forward = False
    effect_is_alt = False
    skip_beagle = False
    genome_build_mapped = "GRCh37"
    pos_col_map = "chr_position"
    
    if st.session_state.panel_ready:
        st.session_state.panel_info['has_embedded_af'] = panel_has_embedded_af(st.session_state.pgs_file)
        
        st.sidebar.markdown("### 🛠️ Advanced Panel & AF Editor")
        st.sidebar.info("Edit alleles, force actions, or adjust embedded AFs directly.")
        
        assume_forward = st.sidebar.checkbox(
            "Assume panel is forward read (Skip palindromic flips)", 
            value=False,
            help="Check this if you have verified from the original panel publication that all variants are on the forward strand. The harmonizer will take this into account."
        )
        effect_is_alt = st.sidebar.checkbox(
            "Assume Effect Allele is always ALT (Skip harmonization mapping)", 
            value=False,
            help="Check this if you verified in the publication that the effect allele reported by the PGS catalog is ALWAYS the alternate allele. If you are not sure, do NOT select this option as it can ruin the results."
        )
        skip_beagle = st.sidebar.checkbox("Skip Beagle Imputation (Force AF Imputation only)", value=False)

        try:
            panel_df = pd.read_csv(st.session_state.pgs_file, sep='\t', comment='#')
            
            valid_cols = [
                'rsid', 'hm_rsid', 'chr_name', 'hm_chr', 'chr_position', 'hm_pos', 
                'effect_allele', 'other_allele', 'hm_inferotherallele', 
                'effect_weight', 'or', 'hr', 'allelefrequency_effect', 
                'is_dominant', 'is_recessive', 'dosage_0_weight', 'dosage_1_weight', 'dosage_2_weight', 
                'is_interaction', 'is_haplotype'
            ]
            preview_cols = [c for c in panel_df.columns if c.lower() in valid_cols]
            
            edit_df = panel_df[preview_cols].copy()
            edit_df['Action'] = "Keep" 
            
            st.sidebar.info("💡 **Note:** You can edit the panel directly in the displayed table below and save it under a custom name.")
            
            edited_df = st.sidebar.data_editor(
                edit_df,
                column_config={"Action": st.column_config.SelectboxColumn("Action", options=["Keep", "Drop", "Force Impute"], required=True)},
                hide_index=True, width='stretch'
            )
            
            actions_df = edited_df[edited_df['Action'] != "Keep"]
            st.session_state.edited_panel_df = edited_df

            save_name = st.sidebar.text_input("Save edited panel as:", placeholder="Leave blank for auto-name")
            
            if st.sidebar.button("💾 Save Panel Configuration"):
                
                # 1. Handle empty name by generating a default name
                if not save_name.strip():
                    original_id = st.session_state.panel_info.get('pgs_id', 'Unknown')
                    # Prevent double CUSTOM_ prefix if editing an already custom panel
                    if original_id.startswith("CUSTOM_"):
                        final_save_name = f"{original_id}_edited"
                    else:
                        final_save_name = f"CUSTOM_{original_id}_edited"
                else:
                    # 2. STRICTLY enforce CUSTOM_ prefix for all user inputs
                    final_save_name = save_name.strip()
                    if not final_save_name.startswith("CUSTOM_"):
                        final_save_name = f"CUSTOM_{final_save_name}"

                save_path = INPUTS_DIR / f"{final_save_name}.txt.gz"
                
                # 3. Apply edits
                updated_full_panel = panel_df.copy()
                for col in preview_cols:
                    if col in edited_df.columns:
                        updated_full_panel[col] = edited_df[col]
                        
                # 4. Save file
                updated_full_panel.to_csv(save_path, sep='\t', index=False, compression='gzip')
                
                # 5. Update session state to point the pipeline to the new file immediately
                st.session_state.pgs_file = str(save_path)
                st.session_state.panel_info['pgs_id'] = final_save_name
                st.session_state.panel_info['trait'] = "Custom Edited Panel"
                
                # 6. Clear the cache so the dropdown re-scans the /inputs folder
                get_panel_cache_status.clear()
                
                st.sidebar.success(f"Saved as {final_save_name}! Refreshing menu...")
                time.sleep(1) 
                
                # 7. Force the app to reload top-to-bottom so the dropdown populates
                st.rerun()

        except Exception as e:
            st.sidebar.error(f"Could not load panel for editing: {e}")

    st.sidebar.markdown("---")
    
    st.sidebar.header("2. Global Configuration")
    genome_build_ui = st.sidebar.radio("Select Genome Build:", options=["GRCh38 (hg38)", "GRCh37 (hg19)"])
    genome_build_mapped = "GRCh37" if "GRCh37" in genome_build_ui else "GRCh38"
    pos_col_map = "chr_position" if genome_build_mapped == "GRCh37" else "hm_pos"
    
    st.sidebar.markdown("---")

    st.sidebar.header("3. Population Allele Frequencies")
    cache_exists = False
    
    # Updated to rely on the dynamic check rather than hardcoded path suffixes
    if st.session_state.panel_ready and 'pgs_file' in st.session_state:
        pid = st.session_state.panel_info.get('pgs_id')
        if pid in panel_records and panel_records[pid]["has_af"]:
            cache_exists = True

    af_options = [
        "Use 1000 Genomes Local Cache" if cache_exists else "Generate 1000 Genomes Cache via Pipeline",
        "Upload Custom Allele Frequencies File",
        "Skip Population Z-Score (Raw PRS Only)"
    ]
    if st.session_state.panel_ready and st.session_state.panel_info.get('has_embedded_af'):
        af_options.insert(0, "Use Allele Frequencies from PGS Panel")

    af_option = st.sidebar.radio("Select Allele Frequency Source:", options=af_options)
    uploaded_custom_af = None

    if af_option == "Use Allele Frequencies from PGS Panel":
        st.session_state.af_mode = "PGS_FILE"
    elif "1000 Genomes" in af_option:
        st.session_state.af_mode = "1000G"
        if not cache_exists:
            st.sidebar.info("Cache will be generated natively from PLINK reference files during execution.")
    elif af_option == "Upload Custom Allele Frequencies File":
        st.session_state.af_mode = "CUSTOM"
        uploaded_custom_af = st.sidebar.file_uploader("Upload Custom Allele Frequencies", type=["csv", "tsv"])
    else:
        st.session_state.af_mode = "NONE"

    st.sidebar.markdown("---")

    st.sidebar.header("4. Patient Ancestry")
    ancestry_mode = st.sidebar.radio("How should the pipeline handle ancestry?", options=["Determine from data (Run Somalier)", "I know the patient's ancestry (Skip Somalier)"])
    known_ancestry = None
    if ancestry_mode == "I know the patient's ancestry (Skip Somalier)":
        known_ancestry = st.sidebar.selectbox("Select Patient Ancestry:", options=["European (EUR)", "African (AFR)", "Admixed American (AMR)", "East Asian (EAS)", "South Asian (SAS)", "Other / Unknown"])
    
    st.sidebar.markdown("---")

    st.sidebar.header("5. Patient Data Input")
    seq_type = st.sidebar.radio("Sequencing Type:", options=["Whole Genome Sequencing (WGS)", "Targeted Sequencing (Panel)"])
    if seq_type == "Targeted Sequencing (Panel)":
        st.sidebar.warning("⚠️ **Targeted NGS Notice:** If Beagle fails to converge on targeted panels due to low background density, the pipeline will safely fallback to Allele Frequency imputation automatically.")
        
    custom_sample_name = st.sidebar.text_input("Custom Sample Name (Optional):").strip()
        
    data_input_method = st.sidebar.radio("Provide patient data:", options=["Upload File (For smaller VCFs)", "Server File Path (For large files)"])
    uploaded_patient_file = None
    patient_file_path_input = ""

    if data_input_method == "Upload File (For smaller VCFs)":
        uploaded_patient_file = st.sidebar.file_uploader("Upload patient file (g.vcf, bam)", type=["g.vcf", "vcf.gz", "vcf", "bam"])
    else:
        patient_file_path_input = st.sidebar.text_input("Enter absolute path to patient VCF/BAM:")

    st.sidebar.markdown("---")

    st.sidebar.header("6. QC Thresholds & Run Config")
    with st.sidebar.expander("⚙ Advanced QC Variables"):
        min_dp = st.number_input("Minimum Read Depth (DP)", value=15)
        min_gq = st.number_input("Minimum Genotype Quality (GQ)", value=50)
        min_af = st.number_input("Min. Heterozygous Allele Fraction", value=0.20)
        beagle_dr2 = st.number_input("Beagle Imputation Score (DR2)", value=0.90)
        cores_num = st.number_input("Number of Cores to Run on", value=8)
        
        st.markdown("---")
        st.markdown("**Palindromic Discard Thresholds**")
        palindromic_range = st.slider(
            "Discard ambiguous variants if EAF is within:",
            min_value=0.0, max_value=1.0, value=(0.40, 0.60), step=0.01,
            help="Default is 0.4 - 0.6. Variants with frequencies in this range are impossible to safely strand-flip."
        )
        pal_low, pal_high = palindromic_range

    st.sidebar.markdown("---")

    st.sidebar.header("7. Execute")

    if st.sidebar.button("Execute Pipeline"):
        if not st.session_state.panel_ready:
            st.sidebar.error("❌ Please load a panel first!")
            st.stop()
            
        pgs_file = st.session_state.pgs_file
        cache_path = pgs_file.replace(".txt.gz", "_af_cache.csv") if st.session_state.af_mode != "PGS_FILE" else "NONE"
            
        if st.session_state.af_mode == "CUSTOM" and uploaded_custom_af is None:
            st.sidebar.error("❌ Please upload your custom allele frequency file first!")
            st.stop()

        final_patient_path = None
        if data_input_method == "Upload File (For smaller VCFs)":
            if uploaded_patient_file is None:
                st.sidebar.error("❌ Please upload a patient file.")
                st.stop()
            input_dir = Path("../Pipeline/data/inputs")
            input_dir.mkdir(parents=True, exist_ok=True)
            final_patient_path = str(input_dir / uploaded_patient_file.name)
            with open(final_patient_path, "wb") as f:
                f.write(uploaded_patient_file.getbuffer())
        else:
            if not patient_file_path_input or not os.path.exists(patient_file_path_input.replace("\\", "/")):
                st.sidebar.error(f"❌ File not found on server.")
                st.stop()
            final_patient_path = str(Path(patient_file_path_input.replace("\\", "/")).resolve())

        file_name = os.path.basename(final_patient_path)
        base_name = custom_sample_name if custom_sample_name else file_name.split('.')[0]
        base_name = re.sub(r'[^a-zA-Z0-9_\-]', '_', base_name)
        
        tmp_path = "../../tmp-streamlit"
        chosen_ancestry = known_ancestry if ancestry_mode != "Determine from data (Run Somalier)" else "RUN_SOMALIER"
        
        if seq_type == "Targeted Sequencing (Panel)":
            internal_base_name = f"{base_name}_TARGETED_{st.session_state.panel_info.get('pgs_id')}"
        else:
            internal_base_name = f"{base_name}_WGS"

        final_sample_id = f"{internal_base_name}_{st.session_state.panel_info.get('pgs_id')}"

        st.markdown("---")
        st.subheader(f"Processing Pipeline: {final_sample_id}")
        col_status, col_dag = st.columns([0.6, 0.4])
        with col_status: ui_status = st.empty()
        with col_dag: ui_dag = st.empty()
        with st.expander("🔍 Detailed Terminal Log (Expand to view raw output)", expanded=False):
            terminal_container = st.empty()
            
        log_text = ""
        os.makedirs(f"{tmp_path}/configs", exist_ok=True)
        actions_file_path = f"{tmp_path}/configs/{final_sample_id}_actions.tsv"
        
        if not actions_df.empty:
            actions_df.to_csv(actions_file_path, sep='\t', index=False)
            log_text += f"[INFO] Saved Variant Actions config to {actions_file_path}\n"
        else:
            open(actions_file_path, 'w').close()

        ui_status.info("**⏳ Currently processing:** Preparing Pipeline...")
        if st.session_state.af_mode == "CUSTOM":
            uploaded_custom_af.seek(0)
            pd.read_csv(uploaded_custom_af, sep=None, engine='python').to_csv(cache_path, index=False)

        input_type = "bam" if file_name.endswith(".bam") else "vcf"
            
        snakemake_cmd = [
            "snakemake", "-s", "../Pipeline/Snakefile", "--cores", str(cores_num),
            "--use-conda", "--conda-frontend", "conda", "--rerun-incomplete", "--config",          
            f"input_file={final_patient_path}",
            f"input_type={input_type}",
            f"seq_type={seq_type}",
            f"genome_build={genome_build_mapped}",     
            f"pgs_panel={st.session_state.panel_info.get('pgs_id')}",
            f"pgs_file={os.path.abspath(pgs_file)}",
            f"af_cache={cache_path if cache_path == 'NONE' else os.path.abspath(cache_path)}",
            f"af_mode={st.session_state.af_mode}",
            f"known_ancestry={chosen_ancestry}",
            f"assume_forward={assume_forward}",
            f"effect_is_alt={effect_is_alt}",
            f"skip_beagle={skip_beagle}",
            f"min_dp={min_dp}", f"min_gq={min_gq}", f"min_af={min_af}", f"beagle_dr2={beagle_dr2}",
            f"base_sample_name={internal_base_name}", f"final_sample_name={final_sample_id}",    
            f"wgs_cache_key={base_name}",
            f"tmp_dir={tmp_path}", f"variant_actions={os.path.abspath(actions_file_path)}" 
        ]
        
        log_text += f"[INFO] Launching Snakemake for {final_sample_id}...\n[CMD] {' '.join(snakemake_cmd)}\n\n"
        terminal_container.code(log_text, language="bash")
        
        try:
            process = subprocess.Popen(snakemake_cmd, cwd="../Pipeline", stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
            in_job_stats, scheduled_jobs, active_jobs, last_seen_rule = False, [], {}, None
            
            for line in iter(process.stdout.readline, ''):
                log_text += line
                clean_line = line.strip()
                
                if "Job stats:" in clean_line: in_job_stats = True
                elif in_job_stats and (clean_line == "" or clean_line.startswith("[")):
                    in_job_stats = False
                    if scheduled_jobs:
                        ui_dag.info("**Pipeline Steps:**\n" + "".join([f"- `{j}` (x{c})\n" for j, c in scheduled_jobs]))
                elif in_job_stats:
                    parts = clean_line.split()
                    if len(parts) >= 2 and parts[0] not in ["job", "total", "Job", "count"]: scheduled_jobs.append((parts[0], parts[1]))
                
                rule_match = re.search(r'^(?:local)?rule (.*?):', clean_line)
                if rule_match: last_seen_rule = rule_match.group(1)
                
                jobid_match = re.search(r'^\s*jobid:\s*(\d+)', clean_line)
                if jobid_match and last_seen_rule:
                    active_jobs[jobid_match.group(1)] = last_seen_rule
                    last_seen_rule = None
                
                finished_match = re.search(r'Finished job (\d+)', clean_line)
                if finished_match and finished_match.group(1) in active_jobs:
                    del active_jobs[finished_match.group(1)]

                if active_jobs:
                    ui_status.info("**⏳ Currently Executing:**\n\n" + "".join([f"- `{r}`\n" for r in sorted(list(set(active_jobs.values())))]))
                else: ui_status.info("**⏳ Evaluating DAG...**")

                if log_text.count('\n') > 1500: log_text = "\n".join(log_text.split('\n')[-1500:])
                terminal_container.code(log_text, language="bash")
                
            process.stdout.close()
            process.wait()
            
            if process.returncode == 0:
                ui_status.success(f"**Pipeline Completed Successfully!**")
                
                somalier_tsv = Path(f"{tmp_path}/ancestry_raw_full/{internal_base_name}_ancestry.tsv")
                somalier_html = Path(f"{tmp_path}/ancestry_raw_full/{internal_base_name}.somalier-ancestry.html")
                
                final_ancestry_label = "EUR"
                if ancestry_mode == "I know the patient's ancestry (Skip Somalier)":
                    final_ancestry_label = known_ancestry
                elif somalier_tsv.exists():
                    try:
                        som_res = pd.read_csv(somalier_tsv, sep='\t')
                        if not som_res.empty and 'predicted_ancestry' in som_res.columns:
                            final_ancestry_label = str(som_res['predicted_ancestry'].iloc[0]).strip()
                    except: pass
                
                cache_suffix = ""
                if assume_forward: cache_suffix += "_fwd"
                if effect_is_alt: cache_suffix += "_alt"
                if st.session_state.af_mode == "PGS_FILE": cache_suffix += "_pgsaf"
                
                exact_sample_id = f"{final_sample_id}{cache_suffix}"
                prs_results_file = Path(f"{tmp_path}/prs_imputed_full/{exact_sample_id}_prs_results.csv")
                
                if prs_results_file.exists():
                    prs_df = pd.read_csv(prs_results_file)
                    
                    raw_somalier_tsv = Path(f"{tmp_path}/ancestry_raw_full/{internal_base_name}.somalier-ancestry.tsv")
                    som_df = pd.read_csv(raw_somalier_tsv, sep='\t') if chosen_ancestry == "RUN_SOMALIER" and raw_somalier_tsv.exists() else None
                    html_content = open(somalier_html, 'r', encoding='utf-8').read() if somalier_html.exists() else None
                    
                    render_dashboard(prs_df, exact_sample_id, final_ancestry_label, som_df, html_content, internal_base_name)
                else: st.error(f"⚠ PRS Results CSV was not found at expected path: {prs_results_file}")
            else: ui_status.error("**❌ Pipeline Failed!** Expand the terminal output to investigate the error.")
                
        except Exception as e: st.error(f"❌ Failed to launch pipeline: {e}")

elif app_mode == "Visualize Existing Results":
    st.header("Visualize Existing PRS Results")
    uploaded_files = st.file_uploader("Drag files here (.csv, .tsv, .html)", type=["csv", "tsv", "html", "htm"], accept_multiple_files=True)
    col1, col2 = st.columns(2)
    with col1: sample_name_input = st.text_input("Sample Name", value="Patient_Data")
    with col2: ancestry_fallback = st.selectbox("Known Ancestry", ["EUR", "AFR", "AMR", "EAS", "SAS", "Other / Unknown"])
    
    if st.button("Render Dashboard") and uploaded_files:
        prs_df_up, som_df_up, html_up = None, None, None
        for f in uploaded_files:
            if f.name.endswith(".csv"): prs_df_up = pd.read_csv(f)
            elif f.name.endswith(".tsv"): som_df_up = pd.read_csv(f, sep='\t')
            elif f.name.endswith((".html", ".htm")): html_up = f.read().decode('utf-8')
        if prs_df_up is not None: render_dashboard(prs_df_up, sample_name_input, ancestry_fallback, som_df_up, html_up)
        else: st.error("⚠️ Please include the PRS Results (.csv) file.")