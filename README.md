
# 🔬 Raman Spectra Visualizer

A web-based application for automated Raman spectroscopy data processing and visualization, developed for saliva-based tuberculosis detection research.

[![Streamlit](https://img.shields.io/badge/Streamlit-FF4B4B?style=for-the-badge&logo=streamlit&logoColor=white)](https://streamlit.io/)
[![Python](https://img.shields.io/badge/Python-3.11-3776AB?style=for-the-badge&logo=python&logoColor=white)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=for-the-badge)](https://opensource.org/licenses/MIT)

## ✨ Features

- 📊 **Automated Baseline Correction** - Adaptive Asymmetric Least Squares (ALS) algorithm with severity detection
- 🧹 **Cosmic Ray Removal** - Median filtering with modified z-score detection
- 📈 **Peak Detection** - Automatic identification and labeling of prominent Raman peaks
- 📉 **SNR Calculation** - Real-time Signal-to-Noise Ratio computation
- 🎨 **Interactive Visualizations** - Powered by Plotly for zoom, pan, and hover capabilities
- 💾 **Data Export** - Download processed spectra as CSV files
- 🎯 **Quality Metrics** - Baseline severity, residual slope, and processing status

## 🚀 Quick Start

### Installation

```
# Clone the repository
git clone https://github.com/chinanuj/RamanSpectraVisualizer.git
cd RamanSpectraVisualizer

# Create virtual environment (recommended)
python -m venv raman
source raman/bin/activate  # On Windows: raman\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

### Usage

```
# Run the application
streamlit run raman_app.py
```

The app will open automatically in your browser at `http://localhost:8501`

## 📁 File Format

Upload CSV files with the following specifications:

- **Header rows**: 98 (automatically skipped)
- **Required columns**: 
  - `Raman Shift` (cm⁻¹)
  - `Dark Subtracted #1` (intensity values)
- **Raman shift range**: 0-2500 cm⁻¹

### Example Data Format

```
... (98 header rows)
Raman Shift,Dark Subtracted #1
400.5,0.0234
401.0,0.0245
401.5,0.0256
...
```

## 🔄 Processing Pipeline

The application employs a sophisticated multi-stage processing pipeline:

```
Raw Data
    ↓
1. Cosmic Ray Removal (Median filtering, threshold=5.0)
    ↓
2. Baseline Detection (Clean/Moderate/Severe classification)
    ↓
3. Adaptive Baseline Correction (ALS with optimized λ, p, iterations)
    ↓
4. Savitzky-Golay Smoothing (window=15, polyorder=2)
    ↓
5. L2 Normalization
    ↓
6. Peak Detection (prominence-based)
    ↓
Processed Spectrum + Metrics
```

## 🎛️ Configurable Settings

- **Sample Type**: Pure Protein, Patient Saliva, or Other
- **Sample Name**: Custom labeling for outputs
- **Peak Labels**: Toggle peak annotations on/off
- **Raw Data Comparison**: Overlay raw vs. processed spectra
- **Color Schemes**: Choose from 5 professional palettes

## 📊 Output Metrics

The application provides comprehensive quality metrics:

| Metric | Description | Acceptable Range |
|--------|-------------|------------------|
| **SNR** | Signal-to-Noise Ratio (dB) | > 15 dB |
| **Baseline Severity** | Clean/Moderate/Severe | - |
| **Residual Slope** | Post-correction baseline flatness | < 0.15 |
| **Cosmic Rays Removed** | Number of spike artifacts | - |

## 🧪 Use Case: TB Detection

This tool was specifically developed for analyzing saliva-based Raman spectra in tuberculosis biomarker research. It identifies characteristic peaks of the MPT64 protein:

**Key MPT64 Peaks (cm⁻¹)**: 445, 541, 838, 910, 1057, 1122, 1257, 1357, 1453

## 📂 Project Structure

```
RamanSpectraVisualizer/
├── raman_app.py           # Main Streamlit application
├── requirements.txt       # Python dependencies
├── README.md             # This file
├── .gitignore           # Git ignore rules
└── LICENSE              # MIT License
```

## 🛠️ Technologies

- **[Streamlit](https://streamlit.io/)** - Web framework
- **[Plotly](https://plotly.com/python/)** - Interactive plotting
- **[SciPy](https://scipy.org/)** - Signal processing & baseline correction
- **[NumPy](https://numpy.org/)** - Numerical computations
- **[Pandas](https://pandas.pydata.org/)** - Data manipulation
- **[scikit-learn](https://scikit-learn.org/)** - Data preprocessing

## 🔬 Scientific Background

### Baseline Correction Algorithm

The Asymmetric Least Squares (ALS) algorithm minimizes:

```
Σ[w_i(y_i - z_i)²] + λΣ[(Δ²z_i)²]
```

Where:
- `y_i`: Raw intensity
- `z_i`: Baseline estimate
- `w_i`: Asymmetric weights
- `λ`: Smoothness parameter (10⁵ - 10⁸)
- `p`: Asymmetry parameter (10⁻⁵ - 10⁻²)

### Peak Detection

Uses `scipy.signal.find_peaks` with:
- **Prominence**: 12% of maximum intensity
- **Height threshold**: 70th percentile of signal region
- **Distance**: Minimum 25 points between peaks

## 🤝 Contributing

Contributions are welcome! Please follow these steps:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📚 Citation

If you use this tool in your research, please cite:

```
@software{raman_visualizer_2025,
  author = {Chincholikar, Anuj},
  title = {Raman Spectra Visualizer: Web-based Tool for Spectroscopy Analysis},
  year = {2025},
  url = {https://github.com/chinanuj/RamanSpectraVisualizer},
  note = {BTech Final Year Project - Saliva-Based TB Detection}
}
```

## 👤 Author

**Anuj Chincholikar**
- GitHub: [@chinanuj](https://github.com/chinanuj)

## 🔮 Future Enhancements

- [ ] Multi-file batch processing
- [ ] Machine learning-based peak classification
- [ ] Database integration for sample tracking
- [ ] Advanced spectral comparison tools
- [ ] Mobile-responsive design
- [ ] Export to publication-ready formats

---

<div align="center">
  
**⭐ If you find this tool useful, please consider giving it a star! ⭐**

</div>
```bash
git add README.md
git commit -m "Add comprehensive README"
git push origin main
```
