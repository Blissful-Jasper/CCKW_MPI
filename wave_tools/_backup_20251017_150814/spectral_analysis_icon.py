"""
Spectral Analysis Module for Atmospheric Wave Detection
=======================================================

This module performs spectral analysis on atmospheric data to identify
wave patterns in the tropical atmosphere, including Kelvin waves and 
Mixed Rossby-Gravity (MRG) waves.

Author: jianpu
Email: xianpuji@hhu.edu.cn
Date: October 2025

Main Classes:
------------
- Config: Configuration parameters for the analysis
- SpectralAnalysis: Base spectral analysis class
- SpectralAnalysis_v3: Improved version with better symmetric decomposition

Key Features:
------------
1. Symmetric/Antisymmetric decomposition for equatorial waves
2. Wavenumber-frequency spectrum computation
3. Background spectrum smoothing
4. Automatic handling of missing equator grid point

Usage Example:
-------------
```python
from spectral_analysis import SpectralAnalysis_v3, Config
import xarray as xr

# Load your data
data = xr.open_dataset('your_data.nc')['pr']

# Create analysis instance
config = Config()
analysis = SpectralAnalysis_v3(config)

# Run analysis
analysis.load_data(data=data)
analysis.preprocess_data_v3()
analysis.compute_spectra_v3()
analysis.smooth_background()

# Plot results
analysis.plot_spectra()
```
"""

import numpy as np
import xarray as xr
import scipy.signal as signal
from scipy import fft
import matplotlib.pyplot as plt
import cmaps
import math
import time
from typing import Optional, Tuple


# ================ CONFIGURATION PARAMETERS ===================
class Config:
    """Configuration parameters for spectral analysis"""
    
    # Data parameters
    DEFAULT_DATA_PATH = r"I://olr.day.mean.nc"
    DEFAULT_VARIABLE = "olr"
    DEFAULT_LAT_RANGE = [-15, 15]
    DEFAULT_TIME_RANGE = ('1997', '2014')
    
    # Analysis parameters
    WINDOW_SIZE_DAYS = 96
    WINDOW_SKIP_DAYS = 30
    SAMPLES_PER_DAY = 1
    
    # Filtering parameters
    FREQ_CUTOFF = 1.0 / WINDOW_SIZE_DAYS
    
    # Plotting parameters
    CONTOUR_LEVELS = np.array([0.0, 0.4, 0.6, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.7, 2.0, 2.4, 2.8, 4.0])
    COLORMAP = cmaps.NCV_blu_red
    WAVENUMBER_LIMIT = 15


# ================ UTILITY FUNCTIONS ===================
def smooth_121(array: np.ndarray) -> np.ndarray:
    """
    Apply a 1-2-1 smoothing filter to the input array.
    
    Parameters:
    -----------
    array : numpy.ndarray
        Input array to be smoothed
        
    Returns:
    --------
    numpy.ndarray
        Smoothed array
    """
    weight = np.array([1., 2., 1.]) / 4.0
    return np.convolve(np.r_[array[0], array, array[-1]], weight, 'valid')


def decompose_to_symmetric_antisymmetric(data_array: xr.DataArray) -> xr.DataArray:
    """
    Decompose data into symmetric and antisymmetric components with respect to the equator.
    
    DEPRECATED: Use decompose_to_symmetric_antisymmetric_v3 instead.
    
    Parameters:
    -----------
    data_array : xarray.DataArray
        Input data with a 'lat' dimension
        
    Returns:
    --------
    xarray.DataArray
        Data with symmetric and antisymmetric components
    """
    lat_dim = data_array.dims.index('lat')
    nlat = data_array.shape[lat_dim]
    
    # Calculate symmetric and antisymmetric components
    symmetric = 0.5 * (data_array.values - np.flip(data_array.values, axis=lat_dim))
    antisymmetric = 0.5 * (data_array.values + np.flip(data_array.values, axis=lat_dim))
    
    # Convert back to DataArrays
    symmetric = xr.DataArray(symmetric, dims=data_array.dims, coords=data_array.coords)
    antisymmetric = xr.DataArray(antisymmetric, dims=data_array.dims, coords=data_array.coords)
    
    # Combine components
    result = data_array.copy()
    half = nlat // 2
    
    if nlat % 2 == 0:
        # Even number of latitudes
        result.isel(lat=slice(0, half))[:] = symmetric.isel(lat=slice(0, half))
        result.isel(lat=slice(half, None))[:] = antisymmetric.isel(lat=slice(half, None))
    else:
        # Odd number of latitudes (equator is in the middle)
        result.isel(lat=slice(0, half))[:] = symmetric.isel(lat=slice(0, half))
        result.isel(lat=slice(half+1, None))[:] = antisymmetric.isel(lat=slice(half+1, None))
        result.isel(lat=half)[:] = symmetric.isel(lat=half)
        
    return result


def decompose_to_symmetric_antisymmetric_v3(data_array: xr.DataArray) -> Tuple[xr.DataArray, xr.DataArray]:
    """
    Improved symmetric/antisymmetric decomposition (v3).
    
    Automatically handles missing equator (0°) in latitude grid by interpolation.
    
    **Problem**: Original latitude grid [-15°, -13°, ..., -1°, 1°, ..., 15°] missing equator (0°)
    **Solution**: Automatically interpolate to add equator point for correct decomposition
    
    Symmetric component: (data(lat) + data(-lat)) / 2  → Kelvin waves
    Antisymmetric component: (data(lat) - data(-lat)) / 2  → MRG waves
    
    Parameters:
    -----------
    data_array : xarray.DataArray
        Input data with a 'lat' dimension
        
    Returns:
    --------
    tuple of xarray.DataArray
        (symmetric_component, antisymmetric_component)
    """
    print(f"  Input data shape: {data_array.shape}")
    lat_values = data_array.lat.values
    
    # Check if equator is included
    has_equator = any(abs(lat) < 1e-6 for lat in lat_values)
    print(f"  Original latitudes: {lat_values[:5]}...{lat_values[-5:]}")
    print(f"  Has equator: {has_equator}")
    
    if not has_equator:
        print(f"  ⚠️  Latitude grid missing equator, adding 0°...")
        
        # Find position to insert 0
        lat_min = lat_values.min()
        lat_max = lat_values.max()
        lat_step = abs(np.diff(lat_values)[0])
        
        # Create complete symmetric grid including 0
        negative_lats = np.arange(lat_min, 0, lat_step)
        positive_lats = np.arange(lat_step, lat_max + lat_step/2, lat_step)
        new_lats = np.concatenate([negative_lats, [0.0], positive_lats])
        
        print(f"  New latitude grid: {new_lats}")
        
        # Interpolate to new grid
        data_interp = data_array.interp(lat=new_lats, method='linear')
        
    else:
        data_interp = data_array
        new_lats = lat_values
    
    # Flip along latitude dimension
    data_flipped = data_interp.isel(lat=slice(None, None, -1))
    
    # Reassign latitude coordinates to match original data
    data_flipped = data_flipped.assign_coords(lat=data_interp.lat)
    
    # Calculate symmetric and antisymmetric components
    symmetric = 0.5 * (data_interp + data_flipped)
    antisymmetric = 0.5 * (data_interp - data_flipped)
    
    print(f"  ✅ Decomposition complete")
    print(f"     Symmetric range: {symmetric.min().values:.2e} to {symmetric.max().values:.2e}")
    print(f"     Antisymmetric range: {antisymmetric.min().values:.2e} to {antisymmetric.max().values:.2e}")
    print(f"     Antisymmetric std: {antisymmetric.std().values:.2e}")
    
    # Check values at equator
    if any(abs(lat) < 1e-6 for lat in symmetric.lat.values):
        eq_sym = symmetric.sel(lat=0, method='nearest').mean().values
        eq_asym = antisymmetric.sel(lat=0, method='nearest').mean().values
        print(f"     Symmetric mean at equator: {eq_sym:.2e}")
        print(f"     Antisymmetric mean at equator: {eq_asym:.2e} (should be ~0)")
    
    return symmetric, antisymmetric


def remove_annual_cycle(data: xr.DataArray, samples_per_day: float, freq_cutoff: float) -> xr.DataArray:
    """
    Remove the annual cycle (low-frequency component) from the data.
    
    Parameters:
    -----------
    data : xarray.DataArray
        Input data with a 'time' dimension
    samples_per_day : float
        Number of samples per day
    freq_cutoff : float
        Frequency cutoff for the low-pass filter
        
    Returns:
    --------
    xarray.DataArray
        Data with annual cycle removed
    """
    n_time, _, _ = data.shape
    
    # Remove linear trend
    detrended_data = signal.detrend(data, axis=0)
    
    # Apply FFT
    fourier_transform = fft.rfft(detrended_data, axis=0)
    frequencies = fft.rfftfreq(n_time, d=1. / float(samples_per_day))
    
    # Apply low-frequency filter
    cutoff_index = np.argwhere(frequencies <= freq_cutoff).max()
    if cutoff_index > 1:
        fourier_transform[1:cutoff_index + 1, ...] = 0.0
    
    # Inverse FFT
    filtered_data = fft.irfft(fourier_transform, axis=0, n=n_time)
    
    return xr.DataArray(filtered_data, dims=data.dims, coords=data.coords)


# ================ MAIN ANALYSIS CLASS ===================
class SpectralAnalysis:
    """
    Base class for performing spectral analysis on atmospheric data.
    
    Note: This is the base implementation. For improved results,
    use SpectralAnalysis_v3 instead.
    """
    
    def __init__(self, config: Optional[Config] = None):
        """
        Initialize the spectral analysis with configuration.
        
        Parameters:
        -----------
        config : Config, optional
            Configuration object with analysis parameters
        """
        self.config = config or Config()
        self.NA = np.newaxis
        self.pi = math.pi
        
        # Initialize data structures
        self.raw_data = None
        self.processed_data = None
        self.power_symmetric = None
        self.power_antisymmetric = None
        self.background = None
        self.frequency = None
        self.wavenumber = None
        
    def load_data(self, data: Optional[xr.DataArray] = None, data_path: Optional[str] = None, 
                  variable: Optional[str] = None, lat_range: Optional[list] = None, 
                  time_range: Optional[tuple] = None):
        """
        Load data either from a file or directly from a provided DataArray.
        
        Parameters:
        -----------
        data : xr.DataArray, optional
            Input data array. If provided, other parameters are ignored.
        data_path : str, optional
            Path to NetCDF file
        variable : str, optional
            Variable name to load from file
        lat_range : list, optional
            [lat_min, lat_max]
        time_range : tuple, optional
            (start_time, end_time)
            
        Returns:
        --------
        self
        """
        if isinstance(data, xr.DataArray):
            print("Using externally provided DataArray.")
            self.raw_data = data
            return self
        
        # Load from file
        data_path = data_path or self.config.DEFAULT_DATA_PATH
        variable = variable or self.config.DEFAULT_VARIABLE
        lat_range = lat_range or self.config.DEFAULT_LAT_RANGE
        time_range = time_range or self.config.DEFAULT_TIME_RANGE
        
        ds = xr.open_dataset(data_path).sortby('lat')
        self.raw_data = ds[variable].sel(
            time=slice(*time_range),
            lat=slice(*lat_range)
        ).sortby('lat').transpose('time', 'lat', 'lon')
        print(f"Data loaded successfully: {self.raw_data.shape}")
        return self
  
    def preprocess_data(self):
        """
        Preprocess the data: detrend, remove annual cycle, and decompose.
        
        Returns:
        --------
        self
        """
        print("Preprocessing data...")
        start_time = time.time()
        
        # Detrend the data
        mean_value = self.raw_data.mean(dim='time')
        detrended = signal.detrend(self.raw_data, axis=0, type='linear')
        detrended = xr.DataArray(detrended, dims=self.raw_data.dims, 
                                coords=self.raw_data.coords) + mean_value
        
        # Remove annual cycle
        filtered = remove_annual_cycle(
            detrended, 
            self.config.SAMPLES_PER_DAY, 
            self.config.FREQ_CUTOFF
        )
        
        # Decompose into symmetric and antisymmetric components
        self.processed_data = decompose_to_symmetric_antisymmetric(filtered)
        
        print(f"Preprocessing completed in {time.time() - start_time:.1f} seconds.")
        return self
        
    def compute_spectra(self):
        """
        Compute the power spectra using windowed FFT.
        
        Returns:
        --------
        self
        """
        print("Computing power spectra...")
        start_time = time.time()
        
        # Extract dimensions
        ntim, nlat, nlon = self.processed_data.shape
        
        # Calculate window parameters
        spd = self.config.SAMPLES_PER_DAY
        nDayWin = self.config.WINDOW_SIZE_DAYS
        nDaySkip = self.config.WINDOW_SKIP_DAYS
        nDayTot = ntim / spd
        nSampWin = nDayWin * spd
        nSampSkip = nDaySkip * spd
        nWindow = int((nDayTot * spd - nSampWin) / (nSampSkip + nSampWin)) + 1
        
        print(f"Analysis parameters: Window size: {nDayWin} days, Skip: {nDaySkip} days")
        print(f"Total windows: {nWindow}, Window samples: {nSampWin}")
        
        # Initialize power accumulator
        sumpower = np.zeros((nSampWin, nlat, nlon))
        
        # Process each window
        ntStrt, ntLast = 0, nSampWin
        for nw in range(int(nWindow)):
            # Extract window data
            data = self.processed_data[ntStrt:ntLast, :, :]
            data = signal.detrend(data, axis=0)
            
            # Apply taper window
            window = signal.windows.tukey(nSampWin, 0.1, True)
            data *= window[:, self.NA, self.NA]
            
            # Compute FFT
            power = fft.fft2(data, axes=(0, 2)) / (nlon * nSampWin)
            sumpower += np.abs(power) ** 2
            
            # Move to next window
            ntStrt = ntLast + nSampSkip
            ntLast = ntStrt + nSampWin
            
            if nw % 10 == 0:
                print(f"Processed window {nw+1}/{nWindow}")
        
        # Normalize by number of windows
        sumpower /= nWindow
        
        # Setup frequency and wavenumber arrays
        if nlon % 2 == 0:
            self.wavenumber = fft.fftshift(fft.fftfreq(nlon) * nlon)[1:]
            sumpower = fft.fftshift(sumpower, axes=2)[:, :, nlon:0:-1]
        else:
            self.wavenumber = fft.fftshift(fft.fftfreq(nlon) * nlon)
            sumpower = fft.fftshift(sumpower, axes=2)[:, :, ::-1]
        
        self.frequency = fft.fftshift(fft.fftfreq(nSampWin, d=1. / float(spd)))[nSampWin // 2:]
        sumpower = fft.fftshift(sumpower, axes=0)[nSampWin // 2:, :, :]
        
        # Compute symmetric and antisymmetric power
        self.power_symmetric = 2.0 * sumpower[:, nlat // 2:, :].sum(axis=1)
        self.power_antisymmetric = 2.0 * sumpower[:, :nlat // 2, :].sum(axis=1)
        
        # Calculate background spectrum
        self.background = 0.5 * (self.power_symmetric + self.power_antisymmetric)
        
        # Convert to DataArrays
        self.power_symmetric = xr.DataArray(
            self.power_symmetric,
            dims=("frequency", "wavenumber"),
            coords={"wavenumber": self.wavenumber, "frequency": self.frequency}
        )
        self.power_antisymmetric = xr.DataArray(
            self.power_antisymmetric,
            dims=("frequency", "wavenumber"),
            coords={"wavenumber": self.wavenumber, "frequency": self.frequency}
        )
        self.background = xr.DataArray(
            self.background,
            dims=("frequency", "wavenumber"),
            coords={"wavenumber": self.wavenumber, "frequency": self.frequency}
        )
        
        # Mask zero frequency to avoid division by zero
        self.power_symmetric = self.power_symmetric.where(self.power_symmetric.frequency > 0)
        self.power_antisymmetric = self.power_antisymmetric.where(self.power_antisymmetric.frequency > 0)
        self.background = self.background.where(self.background.frequency > 0)
        
        print(f"Spectrum computation completed in {time.time() - start_time:.1f} seconds.")
        return self
    
    def smooth_background(self, preserve_low_freq: bool = True):
        """
        Smooth the background spectrum with improved handling of low frequencies.
        
        Parameters:
        -----------
        preserve_low_freq : bool, optional
            If True, preserve low frequency values (< 0.1 cpd) from over-smoothing.
            Default is True to avoid losing data in 0.0-0.1 frequency range.
        
        Returns:
        --------
        self
        """
        print("Smoothing background spectrum...")
        start_time = time.time()
        
        # Convert to numpy array for smoothing
        bg_values = self.background.values.copy()
        
        # Store original low-frequency values if preserving
        if preserve_low_freq:
            low_freq_mask = self.frequency < 0.1
            original_low_freq = bg_values[low_freq_mask, :].copy()
        
        # Select wavenumbers to smooth
        wave_smooth_indices = np.where(np.abs(self.wavenumber) <= 27)[0]
        
        # Smooth each frequency with adaptive iterations
        for idx, freq in enumerate(self.frequency):
            # Skip if frequency is masked (NaN)
            if np.isnan(freq):
                continue
                
            # Determine smoothing iterations based on frequency
            if freq < 0.1:
                smooth_iterations = 3  # Reduced from 5 to preserve low-freq data
            elif freq >= 0.1 and freq < 0.2:
                smooth_iterations = 10
            elif freq >= 0.2 and freq < 0.3:
                smooth_iterations = 20
            else:  # freq >= 0.3
                smooth_iterations = 40
                
            # Apply smoothing only to non-NaN values
            if not np.all(np.isnan(bg_values[idx, wave_smooth_indices])):
                for _ in range(smooth_iterations):
                    # Only smooth non-NaN values
                    valid_mask = ~np.isnan(bg_values[idx, wave_smooth_indices])
                    if np.sum(valid_mask) > 2:  # Need at least 3 points for 1-2-1 filter
                        bg_values[idx, wave_smooth_indices] = smooth_121(
                            bg_values[idx, wave_smooth_indices]
                        )
        
        # Smooth across frequencies
        for wave_idx in wave_smooth_indices:
            # Skip if all NaN
            if np.all(np.isnan(bg_values[:, wave_idx])):
                continue
                
            valid_mask = ~np.isnan(bg_values[:, wave_idx])
            if np.sum(valid_mask) > 2:
                for _ in range(10):
                    bg_values[:, wave_idx] = smooth_121(bg_values[:, wave_idx])
        
        # Restore original low-frequency values if preserving
        if preserve_low_freq:
            print(f"  ✅ Preserving low-frequency data (< 0.1 cpd)")
            bg_values[low_freq_mask, :] = original_low_freq
        
        # Update background with smoothed values
        self.background.values[:] = bg_values
        
        print(f"Background smoothing completed in {time.time() - start_time:.1f} seconds.")
        return self
    
    def plot_spectra(self, output_path: Optional[str] = None, title: Optional[str] = None):
        """
        Plot the normalized power spectra.
        
        Parameters:
        -----------
        output_path : str, optional
            Path to save the plot. If None, the plot will be displayed.
        title : str, optional
            Custom title for the plot
            
        Returns:
        --------
        self
        """
        print("Plotting results...")
        
        # Create figure
        fig, axes = plt.subplots(1, 2, figsize=(12, 7), dpi=200)
        
        # Normalized spectra
        norm_sym = self.power_symmetric / self.background
        norm_asym = self.power_antisymmetric / self.background
        
        # Plot symmetric component
        im1 = norm_sym.plot.contourf(
            ax=axes[0], cmap=self.config.COLORMAP, add_colorbar=False,
            extend='both', levels=self.config.CONTOUR_LEVELS
        )
        
        # Plot antisymmetric component
        im2 = norm_asym.plot.contourf(
            ax=axes[1], cmap=self.config.COLORMAP, add_colorbar=False,
            extend='both', levels=self.config.CONTOUR_LEVELS
        )
        
        # Add colorbar
        fig.colorbar(im1, ax=[axes[0], axes[1]], orientation='horizontal',
                    shrink=0.45, aspect=30, pad=0.1)
        
        # Set titles and labels
        axes[0].set_title("Symmetric Component Spectrum")
        axes[1].set_title("Antisymmetric Component Spectrum")
        
        for ax in axes:
            ax.axvline(0, linestyle='dashed', color='lightgray')
            ax.set_xlim([-self.config.WAVENUMBER_LIMIT, self.config.WAVENUMBER_LIMIT])
            ax.set_ylim([0, 0.5])
            ax.set_xlabel("Wavenumber")
            ax.set_ylabel("Frequency (cpd)")
        
        if title:
            plt.suptitle(title, fontsize=14, fontweight='bold')
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"  ✅ Plot saved to: {output_path}")
        
        plt.show()
        return self


# ================ IMPROVED ANALYSIS CLASS ===================
class SpectralAnalysis(SpectralAnalysis):
    """
    Improved spectral analysis class.
    
    Improvements:
    1. Uses v3 symmetric decomposition (automatically adds equator grid point)
    2. Correct handling of NaN values
    3. Separate computation of symmetric and antisymmetric component spectra
    """
    
    def preprocess_data_v3(self):
        """
        Improved preprocessing workflow using v3 decomposition.
        
        Returns:
        --------
        self
        """
        print("Preprocessing data...")
        start_time = time.time()
        
        # Fill NaN values
        print("  Handling missing values...")
        raw_filled = self.raw_data.fillna(0)
        
        # Detrend
        mean_value = raw_filled.mean(dim='time')
        detrended = signal.detrend(raw_filled.values, axis=0, type='linear')
        detrended = xr.DataArray(detrended, dims=raw_filled.dims, 
                                coords=raw_filled.coords) + mean_value
        
        # Remove annual cycle
        filtered = remove_annual_cycle(
            detrended, 
            self.config.SAMPLES_PER_DAY, 
            self.config.FREQ_CUTOFF
        )
        
        # Apply v3 symmetric/antisymmetric decomposition
        print("  Applying v3 symmetric/antisymmetric decomposition...")
        self.symmetric_data, self.antisymmetric_data = decompose_to_symmetric_antisymmetric_v3(filtered)
        
        self.processed_data = filtered
        
        print(f"Preprocessing completed in {time.time() - start_time:.1f} seconds")
        return self
    
    def compute_spectra_v3(self):
        """
        Compute power spectra for symmetric and antisymmetric components.
        
        Returns:
        --------
        self
        """
        print("Computing power spectra...")
        start_time = time.time()
        
        # Compute symmetric component spectrum
        print("\n  📊 Computing symmetric component spectrum (Kelvin waves)...")
        self.power_symmetric = self._compute_single_component_spectrum(self.symmetric_data)
        
        # Compute antisymmetric component spectrum
        print("  📊 Computing antisymmetric component spectrum (MRG waves)...")
        self.power_antisymmetric = self._compute_single_component_spectrum(self.antisymmetric_data)
        
        # Compute background spectrum
        self.background = 0.5 * (self.power_symmetric + self.power_antisymmetric)
        
        # Mask zero frequency
        self.power_symmetric = self.power_symmetric.where(self.power_symmetric.frequency > 0)
        self.power_antisymmetric = self.power_antisymmetric.where(self.power_antisymmetric.frequency > 0)
        self.background = self.background.where(self.background.frequency > 0)
        
        print(f"\nSpectrum computation completed in {time.time() - start_time:.1f} seconds")
        return self
    
    def _compute_single_component_spectrum(self, data: xr.DataArray) -> xr.DataArray:
        """
        Compute spectrum for a single component.
        
        Parameters:
        -----------
        data : xr.DataArray
            Input data array
            
        Returns:
        --------
        xr.DataArray
            Power spectrum
        """
        ntim, nlat, nlon = data.shape
        
        # Window parameters
        spd = self.config.SAMPLES_PER_DAY
        nDayWin = self.config.WINDOW_SIZE_DAYS
        nDaySkip = self.config.WINDOW_SKIP_DAYS
        nDayTot = ntim / spd
        nSampWin = nDayWin * spd
        nSampSkip = nDaySkip * spd
        nWindow = int((nDayTot * spd - nSampWin) / (nSampSkip + nSampWin)) + 1
        
        # Initialize power accumulator
        sumpower = np.zeros((nSampWin, nlat, nlon))
        
        # Process each window
        ntStrt, ntLast = 0, nSampWin
        for nw in range(int(nWindow)):
            window_data = data.values[ntStrt:ntLast, :, :]
            window_data = signal.detrend(window_data, axis=0)
            
            # Apply window
            window = signal.windows.tukey(nSampWin, 0.1, True)
            window_data = window_data * window[:, self.NA, self.NA]
            
            # FFT
            power = fft.fft2(window_data, axes=(0, 2)) / (nlon * nSampWin)
            sumpower += np.abs(power) ** 2
            
            ntStrt = ntLast + nSampSkip
            ntLast = ntStrt + nSampWin
        
        # Normalize
        sumpower /= nWindow
        
        # Setup frequency and wavenumber
        if nlon % 2 == 0:
            self.wavenumber = fft.fftshift(fft.fftfreq(nlon) * nlon)[1:]
            sumpower = fft.fftshift(sumpower, axes=2)[:, :, nlon:0:-1]
        else:
            self.wavenumber = fft.fftshift(fft.fftfreq(nlon) * nlon)
            sumpower = fft.fftshift(sumpower, axes=2)[:, :, ::-1]
        
        self.frequency = fft.fftshift(fft.fftfreq(nSampWin, d=1. / float(spd)))[nSampWin // 2:]
        sumpower = fft.fftshift(sumpower, axes=0)[nSampWin // 2:, :, :]
        
        # Sum over all latitudes
        power_spectrum = sumpower.sum(axis=1)
        
        # Convert to DataArray
        power_da = xr.DataArray(
            power_spectrum,
            dims=("frequency", "wavenumber"),
            coords={"wavenumber": self.wavenumber, "frequency": self.frequency}
        )
        
        return power_da
    
    def plot_comparison(self, other_analysis, labels: Tuple[str, str] = ("CNTL", "4K"), 
                       output_path: Optional[str] = None, title: Optional[str] = None):
        """
        Plot comparison of two spectral analyses.
        
        Parameters:
        -----------
        other_analysis : SpectralAnalysis_v3
            Another analysis instance to compare with
        labels : tuple of str
            Labels for the two analyses (self, other)
        output_path : str, optional
            Path to save the plot
        title : str, optional
            Custom title for the plot
            
        Returns:
        --------
        self
        """
        print("Plotting comparison...")
        
        # Compute normalized spectra
        norm_sym_1 = self.power_symmetric / self.background
        norm_asym_1 = self.power_antisymmetric / self.background
        norm_sym_2 = other_analysis.power_symmetric / other_analysis.background
        norm_asym_2 = other_analysis.power_antisymmetric / other_analysis.background
        
        # Create figure
        fig, axes = plt.subplots(2, 2, figsize=(16, 12), dpi=200)
        
        # Setup levels and colormap
        levels = self.config.CONTOUR_LEVELS
        cmap = self.config.COLORMAP
        
        # === Top-left: First analysis antisymmetric (MRG) ===
        im1 = norm_asym_1.plot.contourf(
            ax=axes[0, 0], cmap=cmap, add_colorbar=False,
            extend='both', levels=levels
        )
        axes[0, 0].set_title(f"{labels[0]} - Antisymmetric (MRG)", 
                            fontsize=14, fontweight='bold')
        
        # === Top-right: First analysis symmetric (Kelvin) ===
        im2 = norm_sym_1.plot.contourf(
            ax=axes[0, 1], cmap=cmap, add_colorbar=False,
            extend='both', levels=levels
        )
        axes[0, 1].set_title(f"{labels[0]} - Symmetric (Kelvin)", 
                            fontsize=14, fontweight='bold')
        
        # === Bottom-left: Second analysis antisymmetric (MRG) ===
        im3 = norm_asym_2.plot.contourf(
            ax=axes[1, 0], cmap=cmap, add_colorbar=False,
            extend='both', levels=levels
        )
        axes[1, 0].set_title(f"{labels[1]} - Antisymmetric (MRG)", 
                            fontsize=14, fontweight='bold')
        
        # === Bottom-right: Second analysis symmetric (Kelvin) ===
        im4 = norm_sym_2.plot.contourf(
            ax=axes[1, 1], cmap=cmap, add_colorbar=False,
            extend='both', levels=levels
        )
        axes[1, 1].set_title(f"{labels[1]} - Symmetric (Kelvin)", 
                            fontsize=14, fontweight='bold')
        
        # Format all axes
        for ax in axes.flat:
            ax.axvline(0, linestyle='--', color='black', linewidth=1.5, alpha=0.7)
            ax.axhline(0.15, linestyle=':', color='gray', linewidth=1, alpha=0.5)
            ax.set_xlim([-self.config.WAVENUMBER_LIMIT, self.config.WAVENUMBER_LIMIT])
            ax.set_ylim([0, 0.5])
            ax.set_xlabel("Wavenumber", fontsize=12)
            ax.set_ylabel("Frequency (cpd)", fontsize=12)
        
        # Add labels to distinguish symmetric/antisymmetric
        for i, label_text in enumerate(['anti-sym', 'sym', 'anti-sym', 'sym']):
            row, col = divmod(i, 2)
            axes[row, col].text(0.02, 0.98, label_text, 
                              transform=axes[row, col].transAxes,
                              fontsize=11, va='top', ha='left',
                              bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.6))
        
        # Add shared colorbar
        cbar = fig.colorbar(im1, ax=axes, orientation='horizontal',
                          shrink=0.8, aspect=50, pad=0.08)
        cbar.set_label('Normalized Power (Power/Background)', 
                      fontsize=13, fontweight='bold')
        cbar.ax.tick_params(labelsize=11)
        
        # Add 1.0 reference line on colorbar
        cbar.ax.axvline(1.0, color='black', linewidth=2.5, linestyle='--', alpha=0.8)
        
        # Set title
        if title is None:
            title = f'Wavenumber-Frequency Spectrum Comparison: {labels[0]} vs {labels[1]}'
        plt.suptitle(title, fontsize=16, fontweight='bold', y=0.995)
        
        # Save or show
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"  ✅ Comparison plot saved to: {output_path}")
        
        plt.show()
        return self


# ================ CONVENIENCE FUNCTIONS ===================
def quick_analysis(data: xr.DataArray, config: Optional[Config] = None, 
                  plot: bool = True, output_path: Optional[str] = None) -> SpectralAnalysis:
    """
    Quick spectral analysis with default settings.
    
    Parameters:
    -----------
    data : xr.DataArray
        Input data array with dimensions (time, lat, lon)
    config : Config, optional
        Configuration object. If None, uses default Config.
    plot : bool, optional
        Whether to plot results. Default is True.
    output_path : str, optional
        Path to save the plot
        
    Returns:
    --------
    SpectralAnalysis_v3
        Completed analysis object
    """
    print("="*70)
    print("🚀 Quick Spectral Analysis")
    print("="*70)
    
    config = config or Config()
    analysis = SpectralAnalysis_v3(config)
    
    analysis.load_data(data=data)
    analysis.preprocess_data_v3()
    analysis.compute_spectra_v3()
    analysis.smooth_background()
    
    if plot:
        analysis.plot_spectra(output_path=output_path)
    
    print("\n✅ Analysis complete!")
    return analysis


if __name__ == "__main__":
    print("Spectral Analysis Module")
    print("=" * 50)
    print("This module provides tools for atmospheric wave spectral analysis.")
    print("\nMain classes:")
    print("  - Config: Configuration parameters")
    print("  - SpectralAnalysis: Base analysis class")
    print("  - SpectralAnalysis_v3: Improved analysis class (recommended)")
    print("\nFor usage examples, see module docstring or example notebook.")
