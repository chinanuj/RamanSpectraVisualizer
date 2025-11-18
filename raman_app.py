import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from scipy.signal import find_peaks, savgol_filter, medfilt
from scipy.sparse import csc_matrix, diags
from scipy.sparse.linalg import spsolve
from scipy.stats import linregress
import io

st.set_page_config(page_title="Raman Spectra Analyzer", layout="wide")

st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        font-weight: bold;
        color: #1f77b4;
        text-align: center;
        margin-bottom: 2rem;
    }
    .metric-card {
        background-color: #f0f2f6;
        padding: 1rem;
        border-radius: 0.5rem;
        margin: 0.5rem 0;
    }
</style>
""", unsafe_allow_html=True)

st.markdown('<p class="main-header">🔬 Raman Spectra Analyzer</p>', unsafe_allow_html=True)
st.markdown("### Upload your Raman spectroscopy CSV file for automated processing and visualization")

class RamanProcessor:
    def __init__(self):
        pass
    
    def load_csv_data(self, file):
        try:
            df = pd.read_csv(file, skiprows=98, encoding='latin-1')
            df.columns = df.columns.str.strip()
            for col in df.columns:
                try:
                    df[col] = pd.to_numeric(df[col], errors='coerce')
                except:
                    pass
            df = df.dropna(how='all')
            return df
        except Exception as e:
            st.error(f"Error loading file: {e}")
            return None
    
    def get_columns(self, df):
        return 'Raman Shift', 'Dark Subtracted #1'
    
    def remove_cosmic_rays(self, y, threshold=5.0, window=5):
        y_clean = y.copy()
        median = medfilt(y, kernel_size=window)
        deviation = np.abs(y - median)
        median_deviation = np.median(deviation)
        
        if median_deviation > 0:
            modified_z_scores = 0.6745 * deviation / median_deviation
            spikes = modified_z_scores > threshold
            n_spikes = np.sum(spikes)
            y_clean[spikes] = median[spikes]
            return y_clean, n_spikes
        return y_clean, 0
    
    def baseline_als(self, y, lam=1e6, p=0.001, niter=15):
        L = len(y)
        D = csc_matrix(np.diff(np.eye(L), 2))
        w = np.ones(L)
        
        for i in range(niter):
            W = diags(w, 0, shape=(L, L))
            Z = W + lam * D.dot(D.transpose())
            z = spsolve(Z, w * y)
            w = p * (y > z) + (1 - p) * (y < z)
        return z
    
    def detect_baseline_severity(self, x, y):
        mask = (x >= 1900) & (x <= 2200)
        if np.sum(mask) < 10:
            mask = (x >= np.max(x) - 300) & (x <= np.max(x))
        if np.sum(mask) < 10:
            return 'moderate'
        
        x_baseline = x[mask]
        y_baseline = y[mask]
        
        try:
            slope, _, _, _, _ = linregress(x_baseline, y_baseline)
        except:
            return 'moderate'
        
        y_mean = np.mean(y_baseline)
        if y_mean > 0:
            relative_slope = abs(slope * (x_baseline[-1] - x_baseline[0]) / y_mean)
        else:
            relative_slope = 0
        
        y_signal_region = y[(x >= 400) & (x <= 1800)]
        if len(y_signal_region) > 0:
            intensity_ratio = np.max(y_signal_region) / np.mean(y_baseline) if np.mean(y_baseline) > 0 else np.inf
        else:
            intensity_ratio = 1
        
        if relative_slope < 0.1 and intensity_ratio < 50:
            return 'clean'
        elif relative_slope < 0.5 and intensity_ratio < 200:
            return 'moderate'
        else:
            return 'severe'
    
    def get_adaptive_parameters(self, severity):
        params = {
            'clean': {'lam': 5e5, 'p': 0.005, 'niter': 15},
            'moderate': {'lam': 5e6, 'p': 0.0005, 'niter': 20},
            'severe': {'lam': 1e8, 'p': 0.00001, 'niter': 40}
        }
        return params.get(severity, params['moderate'])
    
    def validate_baseline_correction(self, x, y_corrected):
        mask = (x >= 1900) & (x <= 2200)
        if np.sum(mask) < 10:
            mask = (x >= np.max(x) - 300) & (x <= np.max(x))
        if np.sum(mask) < 10:
            return True, 0.0
        
        baseline_region = y_corrected[mask]
        try:
            slope, _, _, _, _ = linregress(x[mask], baseline_region)
            mean_intensity = np.mean(baseline_region)
            if mean_intensity > 0:
                relative_slope = abs(slope * (x[mask][-1] - x[mask][0]) / mean_intensity)
            else:
                relative_slope = 0
            success = relative_slope < 0.15
            return success, relative_slope
        except:
            return True, 0.0
    
    def adaptive_baseline_correction(self, x, y, max_attempts=5):
        severity = self.detect_baseline_severity(x, y)
        params = self.get_adaptive_parameters(severity)
        
        progress_text = st.empty()
        progress_text.text(f"Baseline severity: {severity.upper()}")
        
        best_corrected = None
        best_slope = float('inf')
        
        for attempt in range(max_attempts):
            baseline = self.baseline_als(y, lam=params['lam'], p=params['p'], niter=params['niter'])
            y_corrected = y - baseline
            y_corrected = np.maximum(y_corrected, 0)
            
            success, residual_slope = self.validate_baseline_correction(x, y_corrected)
            
            if residual_slope < best_slope:
                best_corrected = y_corrected.copy()
                best_slope = residual_slope
            
            if success:
                progress_text.text(f"✓ Baseline corrected successfully (slope: {residual_slope:.4f})")
                return best_corrected, severity, params, True, residual_slope
            else:
                params['lam'] *= 10
                params['p'] /= 5
                params['niter'] = min(params['niter'] + 10, 50)
        
        progress_text.text(f"⚠ Best result achieved (slope: {best_slope:.4f})")
        return best_corrected, severity, params, False, best_slope
    
    def smooth_spectrum(self, y, window=15):
        if len(y) > window:
            return savgol_filter(y, window_length=window, polyorder=2)
        return y
    
    def normalize_spectrum(self, y):
        norm = np.linalg.norm(y)
        if norm > 0:
            return y / norm
        return y
    
    def calculate_snr(self, y, x, signal_region=(800, 1600), noise_region=(1900, 2200)):
        signal_mask = (x >= signal_region[0]) & (x <= signal_region[1])
        noise_mask = (x >= noise_region[0]) & (x <= noise_region[1])
        
        signal_power = np.mean(y[signal_mask]**2) if np.any(signal_mask) else 0
        noise_power = np.std(y[noise_mask])**2 if np.any(noise_mask) else 1e-10
        
        snr = 10 * np.log10(signal_power / noise_power) if noise_power > 0 else 0
        return snr
    
    def detect_major_peaks(self, x, y, prominence_factor=0.12):
        y_filtered = y[(x >= 400) & (x <= 1800)]
        if len(y_filtered) > 0:
            height_threshold = np.percentile(y_filtered, 70)
        else:
            height_threshold = 0
        
        y_max = np.max(y)
        prominence = y_max * prominence_factor
        
        peaks, _ = find_peaks(y, prominence=prominence, height=height_threshold, distance=25, width=5)
        mask_region = (x[peaks] >= 400) & (x[peaks] <= 1800)
        return peaks[mask_region]
    
    def process_spectrum(self, df):
        x_col, y_col = self.get_columns(df)
        x = df[x_col].values
        y = df[y_col].values
        
        mask = ~(np.isnan(x) | np.isnan(y)) & (x >= 0) & (x <= 2500)
        x_clean = x[mask]
        y_clean = y[mask]
        
        y_clean, n_spikes = self.remove_cosmic_rays(y_clean, threshold=5.0, window=5)
        y_corrected, severity, params, success, residual_slope = self.adaptive_baseline_correction(x_clean, y_clean)
        y_smooth = self.smooth_spectrum(y_corrected, window=15)
        y_normalized = self.normalize_spectrum(y_smooth)
        snr = self.calculate_snr(y_normalized, x_clean)
        peaks = self.detect_major_peaks(x_clean, y_normalized)
        
        results = {
            'x': x_clean,
            'y_raw': y_clean,
            'y_processed': y_normalized,
            'peaks': peaks,
            'snr': snr,
            'severity': severity,
            'residual_slope': residual_slope,
            'cosmic_rays': n_spikes,
            'success': success
        }
        
        return results

# Sidebar
with st.sidebar:
    st.header("⚙️ Settings")
    sample_type = st.selectbox("Sample Type", ["Pure Protein", "Patient Saliva", "Other"])
    sample_name = st.text_input("Sample Name (optional)", "")
    show_peaks = st.checkbox("Show Peak Labels", value=True)
    show_baseline = st.checkbox("Show Raw Data Comparison", value=False)
    color_scheme = st.selectbox("Color Scheme", ["Blue", "Red", "Green", "Purple", "Orange"])
    
    color_map = {
        "Blue": "#1f77b4",
        "Red": "#d62728",
        "Green": "#2ca02c",
        "Purple": "#9467bd",
        "Orange": "#ff7f0e"
    }

uploaded_file = st.file_uploader("Choose a CSV file", type=['csv'])

if uploaded_file is not None:
    processor = RamanProcessor()
    
    with st.spinner('Loading data...'):
        df = processor.load_csv_data(uploaded_file)
    
    if df is not None:
        st.success("✓ File loaded successfully!")
        
        with st.spinner('Processing spectrum...'):
            results = processor.process_spectrum(df)
        
        st.markdown("### 📊 Processing Results")
        
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.markdown('<div class="metric-card">', unsafe_allow_html=True)
            st.metric("SNR", f"{results['snr']:.1f} dB")
            st.markdown('</div>', unsafe_allow_html=True)
        
        with col2:
            st.markdown('<div class="metric-card">', unsafe_allow_html=True)
            st.metric("Baseline", results['severity'].upper())
            st.markdown('</div>', unsafe_allow_html=True)
        
        with col3:
            st.markdown('<div class="metric-card">', unsafe_allow_html=True)
            st.metric("Residual Slope", f"{results['residual_slope']:.4f}")
            st.markdown('</div>', unsafe_allow_html=True)
        
        with col4:
            st.markdown('<div class="metric-card">', unsafe_allow_html=True)
            st.metric("Cosmic Rays Removed", results['cosmic_rays'])
            st.markdown('</div>', unsafe_allow_html=True)
        
        st.markdown("### 📈 Raman Spectrum")
        
        # PLOTLY INSTEAD OF MATPLOTLIB
        fig = go.Figure()
        
        selected_color = color_map[color_scheme]
        
        # Processed spectrum
        fig.add_trace(go.Scatter(
            x=results['x'],
            y=results['y_processed'],
            mode='lines',
            name='Processed Spectrum',
            line=dict(color=selected_color, width=2.5)
        ))
        
        # Raw data comparison
        if show_baseline:
            y_raw_norm = results['y_raw'] / np.linalg.norm(results['y_raw'])
            fig.add_trace(go.Scatter(
                x=results['x'],
                y=y_raw_norm,
                mode='lines',
                name='Raw Data',
                line=dict(color='gray', width=1.5, dash='dash'),
                opacity=0.5
            ))
        
        # Peak labels
        if show_peaks and len(results['peaks']) > 0:
            peak_x = results['x'][results['peaks']]
            peak_y = results['y_processed'][results['peaks']]
            
            sorted_indices = np.argsort(peak_y)[::-1]
            top_peaks = sorted_indices[:min(10, len(peak_y))]
            
            peak_x_top = peak_x[top_peaks]
            peak_y_top = peak_y[top_peaks]
            
            fig.add_trace(go.Scatter(
                x=peak_x_top,
                y=peak_y_top,
                mode='markers+text',
                name='Peaks',
                marker=dict(color='red', size=10, symbol='star'),
                text=[f'{int(x)}' for x in peak_x_top],
                textposition='top center',
                textfont=dict(size=10, color='black')
            ))
        
        # Layout
        y_max = np.max(results['y_processed'][(results['x'] >= 400) & (results['x'] <= 1800)])
        
        title = f"{sample_name if sample_name else 'Raman Spectrum'} (SNR: {results['snr']:.1f} dB)"
        
        fig.update_layout(
            title=dict(text=title, font=dict(size=16, color='black')),
            xaxis=dict(title='Raman Shift (cm⁻¹)', range=[400, 1800], showgrid=True),
            yaxis=dict(title='Normalized Intensity', range=[0, y_max * 1.15], showgrid=True),
            template='plotly_white',
            height=600,
            showlegend=True,
            hovermode='x unified'
        )
        
        st.plotly_chart(fig, use_container_width=True)
        
        # Peak table
        if len(results['peaks']) > 0:
            st.markdown("### 🔍 Detected Peaks")
            
            peak_x = results['x'][results['peaks']]
            peak_y = results['y_processed'][results['peaks']]
            
            peak_df = pd.DataFrame({
                'Raman Shift (cm⁻¹)': [int(x) for x in peak_x],
                'Intensity': [f"{y:.4f}" for y in peak_y]
            })
            
            peak_df = peak_df.sort_values('Intensity', ascending=False)
            st.dataframe(peak_df, use_container_width=True)
        
        # Download
        st.markdown("### 💾 Download Processed Data")
        
        processed_df = pd.DataFrame({
            'Raman Shift': results['x'],
            'Raw Intensity': results['y_raw'],
            'Processed Intensity': results['y_processed']
        })
        
        csv = processed_df.to_csv(index=False)
        st.download_button(
            label="Download CSV",
            data=csv,
            file_name=f"processed_{uploaded_file.name}",
            mime="text/csv"
        )

else:
    st.info("👆 Please upload a Raman spectroscopy CSV file to begin analysis")
    
    with st.expander("📋 File Format Requirements"):
        st.markdown("""
        Your CSV file should:
        - Have 98 header rows (will be automatically skipped)
        - Contain columns: `Raman Shift` and `Dark Subtracted #1`
        - Have Raman shift values from 0-2500 cm⁻¹
        """)

st.markdown("---")
st.markdown(
    "<p style='text-align: center; color: gray;'>Raman Spectra Analyzer v1.0</p>",
    unsafe_allow_html=True
)
