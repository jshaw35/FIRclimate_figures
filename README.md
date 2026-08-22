# FIRclimate Figures - Reproducible Figure Generation

Python scripts to reproduce all publication figures from the FIRclimate project. These scripts process observational data from the PREFIRE satellite mission and model output from CESM simulations to generate publication-ready figures for far-infrared climate research.

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

**Default configuration:** Scripts are pre-configured to look for data and save figures to default locations. No setup needed if these paths are available on your system:

- **Data input:** `/glade/derecho/scratch/jonahshaw/FIRclimate_zenodo_data/`
- **Figure output:** `/glade/derecho/scratch/jonahshaw/FIRclimate_figures_output/`

**Custom paths:** Override paths using environment variables before running scripts:

```bash
export FIRCLIMATE_DATA_ROOT=/path/to/your/data
export FIRCLIMATE_OUTPUT_ROOT=/path/to/your/figures
```

Alternatively, edit `data_config.py` directly to change default paths:

```python
DATA_ROOT = Path('/your/custom/data/path')
OUTPUT_ROOT = Path('/your/custom/output/path')
```

### 4. Run Figure Scripts

Execute any figure script from the repository root:

```bash
python scripts/figure01a_profiles_spectrum.py
python scripts/figure02_2Dhistograms_CESM.py
# ... run other figure scripts
```

All generated figures will be saved to your configured output directory.

## Data Requirements

This repository requires preprocessed data in zarr and netCDF formats. The data directory should contain:

### Observational Data
- PREFIRE satellite observations (brightness temperature grids, binned statistics)
- Clear-sky validation datasets
- Spectral response functions for PREFIRE instruments

### Model Data
- CESM climate model output (PREFIRE historical and SSP585 scenarios)
- Binned model output conditioned on cloud water vapor and skin temperature
- Zonal mean and 2D histogram data products

### Reference Data
- Atmospheric profiles (MODTRAN standard atmospheres)
- Spectral radiative transfer calculations
- FIR fraction climatologies

For data availability and detailed documentation, consult the accompanying Zenodo data repository.

## Figures Included

| Script | Figure | Description |
|--------|--------|-------------|
| `figure01a_profiles_spectrum.py` | 1(a) | Atmospheric profiles and PREFIRE spectral response functions |
| `figure01b_radianceintuition.py` | 1(b) | Radiance interpretation across spectral channels |
| `figure02_2Dhistograms_CESM.py` | 2 | 2D histograms of CESM model output |
| `figure03_cloudintuition.py` | 3 | Cloud physics and intuition illustration |
| `figure04_2Dhistograms_PREFIREobs.py` | 4 | 2D histograms comparing PREFIRE observations to CESM |
| `figure05_modelvalidation.py` | 5 | Model validation against satellite observations |
| `figure06_futuretrends.py` | 6 | Future climate trends from CESM SSP585 scenario |
| `figure07_FIRfraction_currentfuture.py` | 7 | Far-infrared fraction in current and future climate |
| `figure08_suppf1_2Dhistograms_amplification.py` | S1 | Supplementary: 2D histograms showing climate change amplification |
| `figure_supp_2Dhistograms_countscomparison.py` | S2 | Supplementary: Comparison of sampling distributions |
| `figure_supp_clearskyvalidation.py` | S3 | Supplementary: Clear-sky validation across channels |

## Python Dependencies

All dependencies are specified in `environment.yml`. Key packages include:

| Package | Version | Purpose |
|---------|---------|---------|
| `python` | ≥3.12 | Programming language |
| `xarray` | ≥2024.6 | N-dimensional data handling |
| `zarr` | ≥2.18,<3 | Zarr format I/O |
| `netcdf4` | ≥1.6 | NetCDF format I/O |
| `cartopy` | ≥0.23 | Map projections and geographic plotting |
| `matplotlib` | ≥3.10 | Figure creation and visualization |
| `numpy` | ≥2.0 | Numerical arrays |
| `scipy` | ≥1.14 | Scientific functions |
| `pandas` | ≥2.2 | Data manipulation |
| `seaborn` | ≥0.13 | Statistical visualization |
| `dask` | ≥2024.6 | Parallel/lazy computation |
| `cftime` | ≥1.6 | Calendrical date handling |

Install the full environment:

```bash
conda env create -f environment.yml
```

## Configuration

### Path Configuration

The `data_config.py` module manages all path configuration centrally. It provides:

- **Environment variable overrides:** Set custom paths without modifying code
- **Path validation:** Ensures output directory exists
- **Helper functions:** `get_data_file()` for glob pattern matching, `get_output_path()` for output naming

**Using environment variables:**

```bash
# Set paths for current session
export FIRCLIMATE_DATA_ROOT=/scratch/data/firclimate
export FIRCLIMATE_OUTPUT_ROOT=/scratch/output/firclimate
python scripts/figure01a_profiles_spectrum.py
```

**Direct configuration:**

Edit `data_config.py` to change default paths:

```python
DATA_ROOT = Path('/your/data/directory')
OUTPUT_ROOT = Path('/your/output/directory')
```

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

### Memory Issues with Large Files

**Problem:** Script crashes with memory errors when loading large zarr files

**Solution:** The environment includes dask for lazy evaluation. Large files are automatically chunked. If still experiencing issues:

1. Reduce computation scope (subset time periods or spatial domains)
2. Check available system memory: `free -h` (Linux) or `vm_stat` (macOS)
3. Use dask's distributed scheduler for larger systems

### Path Issues on Different Systems

**Problem:** Scripts fail with path errors when run on different machines

**Solution:** Use environment variables to configure paths:

```bash
# Before running scripts
export FIRCLIMATE_DATA_ROOT=/local/path/to/data
export FIRCLIMATE_OUTPUT_ROOT=/local/path/to/output
```

This allows the same repository to work across systems without modification.

## Running All Figures

To generate all figures at once, create a simple batch script:

```bash
#!/bin/bash
cd FIRclimate_figures
conda activate firclimate

for script in scripts/figure*.py; do
    echo "Running $(basename $script)..."
    python "$script"
    if [ $? -ne 0 ]; then
        echo "ERROR: $(basename $script) failed"
        exit 1
    fi
done

echo "All figures generated successfully!"
```

Save as `run_all_figures.sh`, make executable, and run:

```bash
chmod +x run_all_figures.sh
./run_all_figures.sh
```

## Contributing

To report issues, suggest improvements, or contribute enhancements:

1. Open an issue on GitHub describing the problem or suggestion
2. For bug reports, include:
   - Error message and full traceback
   - Python version and operating system
   - Data path configuration
   - Steps to reproduce

3. For contributions, please ensure:
   - Code follows existing style
   - Changes are backwards compatible
   - Path references use `data_config` module
   - All scripts remain data-independent (paths configurable)

## License

[Specify appropriate license - GPL, MIT, CC-BY, etc.]

## Citation

If you use these scripts or figures in your research, please cite:

### Recommended citation:
```
Shaw, J. (2026). FIRclimate Figures - Reproducible Figure Generation. 
GitHub: https://github.com/jshaw35/FIRclimate_figures

Shaw, J. et al. (202X). [Manuscript title]. [Journal Name], XX(X), XX-XX.
DOI: [Manuscript DOI]

Shaw, J. (202X). FIRclimate Dataset: Processed PREFIRE and CESM data products. 
Zenodo. DOI: [Data DOI]
```

## Authors

**Jonah Shaw**  
University of Colorado Boulder / National Center for Atmospheric Research (NCAR)

## References

### Scientific References
- PREFIRE Mission: [PREFIRE Mission Overview](https://prefire.jpl.nasa.gov/)
- CESM Model: Danabasoglu et al., 2020 (https://doi.org/10.1029/2019MS001916)
- Spectral Analysis: [cite spectral radiative transfer references]

### Related Documentation
- [PREFIRE Data Archive](https://daac.ornl.gov/PREFIRE/)
- [CESM Documentation](http://www.cesm.ucar.edu/)
- [xarray Documentation](http://xarray.pydata.org/)
- [Zarr Format](https://zarr-specs.readthedocs.io/)

## Additional Notes

### Performance Notes
- First run may be slow as data is indexed/loaded
- Subsequent runs benefit from caching
- Large zarr files use dask for efficient memory handling
- Cartopy projections are computed on first use and cached

### System Requirements
- **Minimum:** 8 GB RAM (for individual scripts)
- **Recommended:** 16+ GB RAM (for parallel figure generation)
- **Storage:** ~100 GB for complete data archive

### Development Workflow
When developing new figures or modifications:

1. Use the configuration system for all paths
2. Test with default paths first
3. Test with custom paths via environment variables
4. Avoid hardcoding filesystem paths
5. Use `from data_config import DATA_ROOT, OUTPUT_ROOT` pattern

---

**Last Updated:** August 22, 2026  
**Repository:** https://github.com/jshaw35/FIRclimate_figures  
**Questions?** Open an issue on GitHub
