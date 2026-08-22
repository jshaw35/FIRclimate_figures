"""
Configuration module for FIRclimate figures repository.

This module centralizes path configuration for data input and figure output.
Paths can be overridden using environment variables:
  - FIRCLIMATE_DATA_ROOT: Path to input data directory
  - FIRCLIMATE_OUTPUT_ROOT: Path to output figures directory
"""

from pathlib import Path
import os


# Base directories - configurable via environment variables
DATA_ROOT = Path(os.getenv('FIRCLIMATE_DATA_ROOT', 
                            '/glade/derecho/scratch/jonahshaw/FIRclimate_zenodo_data'))
OUTPUT_ROOT = Path(os.getenv('FIRCLIMATE_OUTPUT_ROOT', 
                             '/glade/derecho/scratch/jonahshaw/FIRclimate_figures_output'))

# Ensure output directory exists
OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)


def get_data_file(pattern: str):
    """
    Find data file matching pattern in DATA_ROOT.
    
    Parameters
    ----------
    pattern : str
        Glob pattern to search for (e.g., "file*.zarr" or "*.nc")
    
    Returns
    -------
    str or list
        If single match: returns path as string
        If multiple matches: returns list of paths
    
    Raises
    ------
    FileNotFoundError
        If no files match the pattern
    """
    import glob
    matches = glob.glob(str(DATA_ROOT / pattern))
    if not matches:
        raise FileNotFoundError(
            f"No files matching pattern '{pattern}' found in {DATA_ROOT}"
        )
    return matches[0] if len(matches) == 1 else matches


def get_output_path(filename: str) -> Path:
    """
    Get output path for a figure.
    
    Parameters
    ----------
    filename : str
        Name of the output file
    
    Returns
    -------
    Path
        Full path to output directory + filename
    """
    return OUTPUT_ROOT / filename


if __name__ == "__main__":
    print(f"Data root: {DATA_ROOT}")
    print(f"Output root: {OUTPUT_ROOT}")
    print(f"Data directory exists: {DATA_ROOT.exists()}")
    print(f"Output directory exists: {OUTPUT_ROOT.exists()}")
