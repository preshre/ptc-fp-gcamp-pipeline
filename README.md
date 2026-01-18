# PTC Fiber Photometry Analysis Code

## Overview
This code performs fiber photometry analysis for Pavlovian Threat Conditioning (PTC) experiments.

**Signal Channels:**
- **465nm (GCaMP)**: Signal channel for calcium-dependent fluorescence
- **405nm (Isobestic)**: Control channel for motion artifact correction

**Outputs:**
- **dF/F (Delta F/F)**: Normalized fluorescence change traces
- **PETH (Peri-Event Time Histogram)**: Mean response aligned to CS onset with SEM
- **AUC (Area Under Curve)**: Quantified responses for Pre-CS, CS, and Post-CS periods
- **Ymax**: Peak response values for each CS event
- **Heatmap**: Visualizes all peri-event dF/F snippets for each CS event as a color-coded matrix

---

## 1. System Requirements

### Software Dependencies
| Package | Version | Purpose |
|---------|---------|---------|
| Python | ≥3.8 | Programming language |
| numpy | ≥1.21.0 | Numerical computations |
| pandas | ≥1.3.0 | Data manipulation |
| scipy | ≥1.7.0 | Signal processing, curve fitting |
| matplotlib | ≥3.4.0 | Plotting and visualization |
| tdt | ≥0.5.0 | TDT data file reading |
| openpyxl | ≥3.0.0 | Excel file generation |

### Operating Systems Tested
- macOS 12.0+ (Monterey and later)
- Windows 10/11

### Hardware Requirements
- **Minimum**: 8 GB RAM, 20 GB free disk space (Recommended to use netowrk drive for the dataset)
- **Recommended**: 16 GB RAM for processing multiple recordings
- No non-standard hardware required

---

## 2. Installation Guide

### Step 1: Install Python
Download and install Python 3.8 or later from [python.org](https://www.python.org/downloads/)

### Step 2: Clone or Download the Code
```bash
# Option A: Download and extract the ZIP file
# Option B: Clone from repository
git clone <repository-url>
cd PTC_FiberPhotometry
```

### Step 3: Create Virtual Environment (Recommended)
```bash
python -m venv venv

# Activate on macOS/Linux:
source venv/bin/activate

# Activate on Windows:
venv\Scripts\activate
```

### Step 4: Install Dependencies
```bash
pip install -r requirements.txt
```

### Typical Installation Time
- **Fresh install**: ~5-10 minutes (depending on internet speed)
- **With existing Python environment**: ~2-3 minutes

---

## 3. Demo

### Demo Data
Sample TDT fiber photometry data is provided in the `demo/sample_data/` folder.

### Instructions to Run Demo
```bash
python pl_fp_analysis_465.py --dir demo/sample_data
```

### Expected Output
After running the demo, the following files will be generated in `demo/expected_output/`:

```
demo/expected_output/
├── dFF/
│   └── <mouse_id>/
│       ├── <mouse_id>_train_dFF.png
│       └── <mouse_id>_train_dFF.pdf
├── PETH Histograms/
│   └── <mouse_id>/
│       ├── <mouse_id>_train_PETH.png
│       └── <mouse_id>_train_PETH.pdf
├── AUC Plots/
│   └── <mouse_id>/
│       ├── <mouse_id>_train_AUC.png
│       └── <mouse_id>_train_AUC.pdf
├── Heatmaps/
│   └── <mouse_id>/
│       ├── <mouse_id>_train_heatmap.png
│       └── <mouse_id>_train_heatmap.pdf
├── auc.txt
├── auc.xlsx
├── yMax.txt
└── yMax.xlsx
```

### Expected Run Time
- **Single mouse recording**: ~30-60 seconds
- **Demo dataset**: ~1-2 minutes

---

## 4. Instructions for Use

### Running on Your Own Data

#### Step 1: Prepare Your Data
Organize your TDT fiber photometry data in a directory structure:
```
your_data/
├── Mouse1_Training/
│   ├── *.tev
│   ├── *.tsq
│   └── ...
├── Mouse1_LTM1d/
├── Mouse2_Training/
└── ...
```

#### Step 2: Run the Analysis
```bash
python pl_fp_analysis_465.py
```
When prompted, enter the path to your data directory.

#### Alternative: Command Line Arguments
```bash
python pl_fp_analysis_465.py --dir /path/to/your/data --signal "_465A" --control "_405A"
```

#### Step 3: View Results
Results are saved to `Results/PTC_Results_<timestamp>/`

### Command Line Options
| Argument | Default | Description |
|----------|---------|-------------|
| `--dir` | (prompts) | Path to TDT data directory |
| `--signal` | `_465A` | Signal channel name (465nm GCaMP) |
| `--control` | `_405A` | Control channel name (405nm isobestic) |

---

## 5. Reproduction Instructions

To reproduce the quantitative results in the manuscript:

1. **Install the software** following Section 2
2. **Obtain the data** from [data repository link]
3. **Run the analysis**:
   ```bash
   python pl_fp_analysis_465.py --dir /path/to/full/dataset
   ```
4. **Compare outputs**: The generated `auc.xlsx` and `yMax.xlsx` files contain the values used for statistical analysis in the manuscript

---

## File Structure
```
PTC_FiberPhotometry/
├── README.md                    # This file
├── requirements.txt             # Python dependencies
├── pl_fp_analysis_465.py        # Main analysis script
├── pl_fp_mouse_465.py           # Mouse data class with signal processing
├── utils/
│   ├── __init__.py
│   ├── utils.py                 # Utility functions
│   └── extract_timestamps.py    # Timestamp extraction from Excel
└── demo/
    ├── sample_data/             # Sample TDT recording (copy your data here)
    └── expected_output/         # Expected results from demo
```

---

## Analysis Pipeline

1. **Data Loading**: Load TDT fiber photometry data using the `tdt` library
2. **Low-pass Filtering**: Apply 6Hz 3rd-order Butterworth filter to remove high-frequency noise
3. **Photobleaching Correction**: Fit double exponential curve (y = a·e^(bx) + c·e^(dx)) and subtract
4. **Motion Correction**: Linear regression of 465nm signal against 405nm control, subtract fitted values
5. **dF/F Calculation**: Compute ΔF/F = 100 × (signal - fitted) / (exponential_fit + ε)
6. **Z-score Calculation**: Standardize with Savitzky-Golay smoothing (window=5000, polyorder=3)
7. **PETH Generation**: Extract peri-event snippets aligned to CS onset (-30s to +60s)
8. **AUC Calculation**: Integrate using trapezoidal rule for Pre-CS (-30 to 0s), CS (0 to 30s), Post-CS (30 to 60s)
9. **Ymax Extraction**: Find peak z-score value within each CS response window

---

## Output Description

| Output | Format | Description |
|--------|--------|-------------|
| dF/F plots | PNG, PDF | Full trace of normalized fluorescence over recording session |
| PETH plots | PNG, PDF | Mean ± SEM response aligned to CS onset |
| AUC plots | PNG, PDF | Bar graph of integrated responses by time window |
| Heatmap plots | PNG, PDF | Color-coded matrix of all peri-event dF/F snippets for each CS event |
| auc.xlsx | Excel | Per-CS AUC values for Pre-CS, CS, and Post-CS periods |
| yMax.xlsx | Excel | Peak response values for each CS/US event |

---

## License
MIT License

---

## Contact
For questions or issues, please contact: Prerana.Shrestha@stonybrook.edu
