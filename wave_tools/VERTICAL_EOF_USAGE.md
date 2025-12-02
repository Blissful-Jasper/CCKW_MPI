# Vertical EOF Analysis with NaN Handling

## Overview

The `wave_tools.eof` module now includes specialized functions for performing EOF (Empirical Orthogonal Function) analysis on vertical velocity data with automatic NaN handling. This is particularly useful for ocean-masked atmospheric data.

## Key Features

✅ **Automatic NaN handling** - Removes NaN points before analysis, reconstructs results  
✅ **Sign alignment** - Ensures consistent EOF signs across experiments  
✅ **Easy comparison** - Built-in functions to compare multiple experiments  
✅ **Preserves grid structure** - Results maintain original dimensions with NaN in masked regions  

## Installation Requirements

```bash
pip install xeofs
# or
conda install -c conda-forge xeofs
```

## Quick Start

### Basic Usage

```python
from wave_tools.eof import vertical_eof_with_nan_handling
import xarray as xr

# Load your data (can contain NaN values over land)
wa_data = xr.open_dataset('omega_cntl_layers/wa_all_levels.nc')['wa']

# Apply ocean mask (creates NaN over land)
wa_data = wa_data.where(ocean_mask, drop=False)

# Perform EOF analysis
eofs, pcs, variance = vertical_eof_with_nan_handling(
    wa_data,
    n_modes=2
)

print(f"Mode 1 explains {variance[0].values:.2%} of variance")
print(f"EOFs shape: {eofs.shape}")  # (n_modes, n_levels)
print(f"PCs shape: {pcs.shape}")    # (n_modes, n_time, n_lat, n_lon)
```

### Comparing Multiple Experiments

```python
from wave_tools.eof import (
    vertical_eof_with_nan_handling,
    compare_vertical_eofs
)

# Analyze multiple experiments
results = {
    'CNTL': vertical_eof_with_nan_handling(wa_cntl.wa, n_modes=2),
    'P4K': vertical_eof_with_nan_handling(wa_p4k.wa, n_modes=2),
    '4CO2': vertical_eof_with_nan_handling(wa_4co2.wa, n_modes=2)
}

# Compare and plot (signs automatically aligned to CNTL)
fig = compare_vertical_eofs(
    results,
    reference_key='CNTL',
    save_path='./figures/eof_comparison.png'
)
```

### Manual Sign Alignment

```python
from wave_tools.eof import align_eof_signs

# Analyze experiments
eofs_cntl, pcs_cntl, var_cntl = vertical_eof_with_nan_handling(wa_cntl.wa, n_modes=2)
eofs_p4k, pcs_p4k, var_p4k = vertical_eof_with_nan_handling(wa_p4k.wa, n_modes=2)

# Align signs manually for each mode
eofs_p4k_aligned = eofs_p4k.copy()
for mode_idx in range(2):
    sign = align_eof_signs(
        eofs_cntl.isel(mode=mode_idx),
        eofs_p4k.isel(mode=mode_idx)
    )
    eofs_p4k_aligned[mode_idx, :] *= sign
    print(f"Mode {mode_idx+1}: sign = {'+' if sign > 0 else '-'}")
```

## Function Reference

### `vertical_eof_with_nan_handling()`

Main function for EOF analysis with NaN handling.

**Parameters:**
- `vert_vel` (xr.DataArray): Vertical velocity data, can contain NaN
- `n_modes` (int): Number of EOF modes to compute (default: 2)
- `zg` (xr.DataArray, optional): Geopotential height for coordinate assignment

**Returns:**
- `eofs` (xr.DataArray): EOF spatial patterns, shape (n_modes, n_levels)
- `pcs` (xr.DataArray): Principal components, shape (n_modes, n_time, n_lat, n_lon)
- `explained_variance` (xr.DataArray): Explained variance ratio for each mode

### `align_eof_signs()`

Align EOF signs based on correlation.

**Parameters:**
- `eof_ref` (xr.DataArray): Reference EOF pattern
- `eof_target` (xr.DataArray): Target EOF to align

**Returns:**
- `sign` (int): +1 or -1 multiplier

### `compare_vertical_eofs()`

Compare and visualize EOFs from multiple experiments.

**Parameters:**
- `eofs_dict` (dict): Dictionary with experiment names and (eofs, pcs, variance) tuples
- `reference_key` (str, optional): Reference experiment for sign alignment
- `figsize` (tuple): Figure size (default: (14, 6))
- `colors` (list, optional): Custom colors for each experiment
- `save_path` (str, optional): Path to save figure

**Returns:**
- `fig` (matplotlib.figure.Figure): Figure object

## How It Works

### NaN Handling Pipeline

1. **Stack dimensions**: Convert (time, lat, lon, level) → (sample, level)
2. **Identify valid points**: Mark samples where at least one level has valid data
3. **Remove NaN samples**: Keep only ocean points for EOF computation
4. **Perform EOF**: Standard EOF decomposition on clean data
5. **Reconstruct**: Map results back to original grid, filling NaN where removed

### Why This Matters

Traditional EOF methods fail when data contains NaN values. This implementation:
- ✅ Handles ocean/land masks automatically
- ✅ Preserves original data structure
- ✅ Maintains physical interpretation
- ✅ Enables cross-experiment comparison

## Output Interpretation

### EOFs (Empirical Orthogonal Functions)
- Represent **vertical structure** of dominant modes
- Shape: `(n_modes, n_levels)`
- Units: Normalized amplitude
- Peak values indicate levels of maximum variance

### PCs (Principal Components)  
- Represent **time-space evolution** of each mode
- Shape: `(n_modes, n_time, n_lat, n_lon)`
- Contains NaN over masked regions (e.g., land)
- Units: Match input data units (scaled by EOF amplitude)

### Explained Variance
- Percentage of total variance explained by each mode
- Sum of all modes < 100% (higher modes truncated)
- Mode 1 typically explains 50-70% for vertical velocity

## Example Notebook

See `wa_xeof.ipynb` for complete working examples including:
- Data loading and masking
- EOF analysis for multiple experiments
- Vertical profile comparison
- Physical interpretation

## Troubleshooting

**Issue**: `ImportError: xeofs library is required`  
**Solution**: `pip install xeofs` or `conda install -c conda-forge xeofs`

**Issue**: All values are NaN in PCs  
**Solution**: Check that input data has valid (non-NaN) values over ocean regions

**Issue**: EOF signs differ between runs  
**Solution**: Use `align_eof_signs()` or `compare_vertical_eofs()` with `reference_key`

**Issue**: Shape mismatch errors  
**Solution**: Ensure input data has dimensions (time, lat, lon, level) or can be stacked

## Citation

If you use these functions in your research, please cite:
- xeofs library: https://github.com/nicrie/xeofs
- Original implementation: wave_tools package (2025)

## Contact

For questions or issues, please contact the wave_tools development team or open an issue on GitHub.
