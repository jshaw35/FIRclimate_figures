# Code for reproducing figures in "Understanding Far-Infrared Signals of 21st Century Climate Change and Variability"

Python scripts to reproduce all main text figures.

## Quick Start

### 1. Clone the Repository

```bash
git clone https://github.com/jshaw35/FIRclimate_figures.git
cd FIRclimate_figures
```

### 2. Set Up the Conda Environment

Create and activate the environment:

```bash
conda env create -f environment.yml
conda activate firclimate
```

The environment file specifies Python 3.12+ with all required scientific computing packages (xarray, zarr, netCDF4, cartopy, etc.).

### 3. Configure Data Paths

Update paths that point to the unzipped and a directory where figures should be save. Two options:

**a. Custom paths:** Override paths using environment variables before running scripts:

```bash
export FIRCLIMATE_DATA_ROOT=/path/to/your/data
export FIRCLIMATE_OUTPUT_ROOT=/path/to/your/figures
```

**b. Edit `data_config.py` directly to change default paths:**

```python
DATA_ROOT = Path('/your/custom/data/path')
OUTPUT_ROOT = Path('/your/custom/output/path')
```

### 4. Download and unzip data from Zenodo

Zipped data is available at: [10.5281/zenodo.22070286](https://doi.org/10.5281/zenodo.22070286)

```
unzip FIRdata.zip
```

### 5. Run Figure Scripts

Execute any figure script from the repository root:

```bash
python scripts/figure01a_profiles_spectrum.py
python scripts/figure02_2Dhistograms_CESM.py
# ... run other figure scripts
```

All generated figures will be saved to your configured output directory.


## Figures Included

| Script | Figure | Description |
|--------|--------|-------------|
| `figure01a_profiles_spectrum.py` | 1(a) | Atmospheric profiles and PREFIRE spectral response functions |
| `figure01b_radianceintuition.py` | 1(b),D1 | Average Brightness Temperature behavior across humidity regimes |
| `figure02_2Dhistograms_CESM.py` | 2 | 2D histograms of "PREFIRE-like" CESM model output |
| `figure03_cloudintuition.py` | 3 | Cloud Radiative Effects |
| `figure04_2Dhistograms_PREFIREobs.py` | 4 | Regime-composite comparisons between PREFIRE and CESM2 |
| `figure05_modelvalidation.py` | 5 | Zonal-mean comparisons between PREFIRE and CESM2 |
| `figure06_futuretrends.py` | 6 | Future spectral flux trends from CESM SSP585 scenario |
| `figure07_FIRfraction_currentfuture.py` | 7 | Far-infrared fraction in current and future climate |
| `figure08_suppf1_2Dhistograms_amplification.py` | 8,G1 | FIR Antarctic Amplification Explain |
| `figure_supp_2Dhistograms_countscomparison.py` | E1 | Supplementary: Comparison of PREFIRE and CESM2 PWV-TS sampling distributions |
| `figure_supp_clearskyvalidation.py` | F1-3 | Supplementary: Clear-sky validation across channels |


## Configuration

### Path Configuration

The `data_config.py` module manages all path configuration centrally. It provides:

- **Environment variable overrides:** Set custom paths without modifying code
- **Path validation:** Ensures output directory exists


## Troubleshooting

### Import Errors

**Problem:** `ModuleNotFoundError: No module named 'xarray'`

**Solution:** Ensure the `firclimate` conda environment is activated:

```bash
conda activate firclimate
python scripts/figure01a_profiles_spectrum.py
```

### File Not Found Errors

**Problem:** `FileNotFoundError: No files matching pattern 'xyz*.zarr' found in /path/to/data`

**Solution:** Verify that `FIRCLIMATE_DATA_ROOT` points to the correct data directory and that all required files exist:

```bash
# Check data directory configuration
python -c "from data_config import DATA_ROOT; print(f'Data root: {DATA_ROOT}'); print(f'Exists: {DATA_ROOT.exists()}')"

# List available files
ls $FIRCLIMATE_DATA_ROOT | head
```

### Cartopy/Projection Errors

**Problem:** `ValueError: Projection xx not recognized` or cartopy import errors

**Solution:** Cartopy requires system-level geographic data. On some systems, additional dependencies may be needed:

```bash
# On Linux/UNIX
conda install -c conda-forge proj geos

# On macOS
conda install -c conda-forge proj geos
```

For detailed cartopy installation, see [Cartopy documentation](https://scitools.org.uk/cartopy/docs/latest/install.html).

## License

MIT

## Citation

If you use these scripts or figures in your research, please acknowledge this repository, the data, and the manuscript.

## Authors

**Jonah Shaw**  
University of Colorado Boulder / Cooperative Institute for Research in Environmental Science (CIRES)

### Limited References
- PREFIRE Mission: [PREFIRE Mission Overview](https://prefire.jpl.nasa.gov/)
- CESM2: Danabasoglu et al., 2020 (https://doi.org/10.1029/2019MS001916)
