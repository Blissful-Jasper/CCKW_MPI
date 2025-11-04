"""
EOF Analysis Module for Wave Tools
===================================

This module provides Empirical Orthogonal Function (EOF) analysis 
with two implementation methods:
1. SVD-based (numpy) - Custom implementation
2. xeofs-based - Using the xeofs library

Features:
---------
- Flexible input handling (xr.DataArray)
- Automatic dimension detection and transposition
- Land/ocean masking capabilities
- Climatology removal with harmonic filtering
- Vertical profile plotting
- Result saving/loading

Author: xianpuji
Date: October 2025
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from typing import Dict, Tuple, Optional, Union, List
import os
import pickle

# Optional import for land masking
try:
    from global_land_mask import globe
    _LAND_MASK_AVAILABLE = True
except ImportError:
    _LAND_MASK_AVAILABLE = False
    import warnings
    warnings.warn(
        "global_land_mask not available. Land/ocean masking will be disabled. "
        "Install with: pip install global-land-mask",
        ImportWarning
    )


class EOFAnalyzer:
    """
    Comprehensive EOF analysis class with multiple methods and visualization.
    
    This class supports both SVD-based and xeofs-based EOF decomposition,
    with built-in support for masking, climatology removal, and visualization.
    
    Parameters
    ----------
    method : str, default='svd'
        Decomposition method: 'svd' or 'xeofs'
    apply_land_mask : bool, default=False
        Whether to apply land/ocean masking
    ocean_only : bool, default=True
        If True, keep ocean points; if False, keep land points
    mask_resolution : str, default='c'
        Mask resolution: 'c' (coarse ~50km), 'l' (low ~10km), 
        'i' (intermediate ~2km), 'h' (high ~400m), 'f' (full ~100m)
    
    Examples
    --------
    >>> # Basic SVD-based analysis
    >>> analyzer = EOFAnalyzer(method='svd')
    >>> results = analyzer.fit(data, n_modes=4)
    >>> analyzer.plot_vertical_profiles()
    
    >>> # xeofs-based analysis with ocean masking
    >>> analyzer = EOFAnalyzer(method='xeofs', apply_land_mask=True, ocean_only=True)
    >>> results = analyzer.fit(data, n_modes=4)
    >>> fig = analyzer.plot_vertical_profiles(n_modes=4)
    """
    
    def __init__(self, 
                 method: str = 'svd',
                 apply_land_mask: bool = False, 
                 ocean_only: bool = True, 
                 mask_resolution: str = 'c'):
        """Initialize the EOF analyzer."""
        
        if method not in ['svd', 'xeofs']:
            raise ValueError(f"Method must be 'svd' or 'xeofs', got '{method}'")
        
        self.method = method
        self.apply_land_mask = apply_land_mask
        self.ocean_only = ocean_only
        self.mask_resolution = mask_resolution
        self.land_mask = None
        self.eof_results = {}
        
        # Try to import xeofs if needed
        if method == 'xeofs':
            try:
                import xeofs as xe
                self.xe = xe
                print("✅ xeofs library loaded successfully")
            except ImportError:
                raise ImportError(
                    "xeofs library is required for method='xeofs'. "
                    "Install with: pip install xeofs or conda install -c conda-forge xeofs"
                )
    
    def create_land_mask(self, data: xr.DataArray, resolution: str = 'c') -> xr.DataArray:
        """
        Create land/ocean mask using global-land-mask.
        
        Parameters
        ----------
        data : xr.DataArray
            Input data with lat/lon coordinates
        resolution : str
            Mask resolution (not used with global-land-mask, kept for compatibility)
            
        Returns
        -------
        xr.DataArray
            Boolean mask (True=ocean if ocean_only=True, True=land if ocean_only=False)
        """
        if not _LAND_MASK_AVAILABLE:
            raise ImportError(
                "global_land_mask is required for land/ocean masking. "
                "Install with: pip install global-land-mask"
            )
        
        # Detect coordinate names
        lat_name = None
        lon_name = None
        
        for coord in data.coords:
            if coord.lower() in ['lat', 'latitude', 'y']:
                lat_name = coord
            elif coord.lower() in ['lon', 'longitude', 'x']:
                lon_name = coord
                
        if lat_name is None or lon_name is None:
            raise ValueError("Cannot find latitude/longitude coordinates in data")
        
        # Get coordinate values
        lats = data[lat_name].values
        lons = data[lon_name].values
        
        # Convert longitude from 0-360 to -180-180 if needed
        lons_converted = lons.copy()
        if np.any(lons > 180):
            print(f"   Converting longitude from 0-360 to -180-180 range...")
            lons_converted = np.where(lons > 180, lons - 360, lons)
        
        # Create meshgrid
        lon_grid, lat_grid = np.meshgrid(lons_converted, lats)
        
        try:
            # Use global-land-mask (returns True for land)
            land_mask_2d = globe.is_land(lat_grid, lon_grid)
            print("   ✅ Using global-land-mask for land/ocean detection")
        except Exception as e:
            print(f"   ⚠️  global-land-mask failed ({e}), using fallback")
            land_mask_2d = np.zeros_like(lat_grid, dtype=bool)
        
        # Create DataArray
        mask_coords = {lat_name: lats, lon_name: lons}
        land_mask_da = xr.DataArray(
            land_mask_2d, 
            coords=mask_coords, 
            dims=[lat_name, lon_name]
        )
        
        # Apply ocean_only logic
        if self.ocean_only:
            self.land_mask = ~land_mask_da  # True for ocean
        else:
            self.land_mask = land_mask_da   # True for land
            
        return self.land_mask

    def _detect_dims(self, data: xr.DataArray) -> Dict[str, Optional[str]]:
        """
        Detect common dimension names in an xarray DataArray.

        Returns a dict with keys: time, level, lat, lon, ensemble.
        Values are the detected dimension name or None if not present.
        The detection is case-insensitive and uses substring matching
        to be robust to names like 'plev', 'level', 'lev', 'time_counter', etc.
        """
        dims = list(data.dims)
        dims_lc = [d.lower() for d in dims]

        def find_exact_or_contains(candidates: List[str]) -> Optional[str]:
            for cand in candidates:
                # exact match first
                if cand in dims_lc:
                    return dims[dims_lc.index(cand)]
            # substring matches
            for i, d in enumerate(dims_lc):
                for cand in candidates:
                    if cand in d:
                        return dims[i]
            return None

        time_dim = find_exact_or_contains(['time', 'time_counter', 'date', 'datetime'])
        lat_dim = find_exact_or_contains(['lat', 'latitude', 'y'])
        lon_dim = find_exact_or_contains(['lon', 'longitude', 'x'])
        # level candidates: plev, lev, level, pressure, pfull
        level_dim = find_exact_or_contains(['plev', 'lev', 'level', 'pressure', 'pfull'])
        ensemble_dim = find_exact_or_contains(['member', 'ensemble', 'ens'])

        return {
            'time': time_dim,
            'level': level_dim,
            'lat': lat_dim,
            'lon': lon_dim,
            'ensemble': ensemble_dim
        }
    
    def extract_low_harmonics(self, 
                             data: xr.DataArray, 
                             n_harm: int = 3, 
                             dim: str = 'dayofyear') -> xr.DataArray:
        """
        Extract low-order harmonics using FFT.
        
        Parameters
        ----------
        data : xr.DataArray
            Input climatology data
        n_harm : int
            Number of harmonics to retain
        dim : str
            Dimension for FFT application
            
        Returns
        -------
        xr.DataArray
            Smoothed climatology
        """
        # Fourier transform
        z_fft = np.fft.rfft(data, axis=data.get_axis_num(dim))
        
        # Filter harmonics
        z_fft_filtered = z_fft.copy()
        if n_harm < z_fft.shape[data.get_axis_num(dim)]:
            z_fft_filtered[n_harm, :, :] *= 0.5
            z_fft_filtered[(n_harm+1):, :, :] = 0
        
        # Inverse FFT
        smoothed_data = np.fft.irfft(
            z_fft_filtered, 
            n=data.sizes[dim], 
            axis=data.get_axis_num(dim)
        ).real
        
        # Create output DataArray
        coords = {k: v for k, v in data.coords.items()}
        attrs = {
            "smoothing": f"FFT: {n_harm} harmonics retained",
            "units": data.attrs.get("units", ""),
            "long_name": f"Daily Climatology: {n_harm} harmonics retained",
        }
        
        return xr.DataArray(smoothed_data, coords=coords, dims=data.dims, attrs=attrs)
    
    def _compute_eof_svd(self, data: np.ndarray) -> Dict:
        """
        Compute EOF using SVD (numpy-based).
        
        Parameters
        ----------
        data : np.ndarray
            Data matrix (variables × observations)
            
        Returns
        -------
        dict
            EOF results including patterns, PCs, eigenvalues, variance
        """
        # Handle missing values
        valid_mask = ~np.any(np.isnan(data), axis=0)
        data_valid = data[:, valid_mask]
        
        # SVD decomposition
        u, s, v = np.linalg.svd(data_valid, full_matrices=False)
        
        # EOF patterns and PCs
        eof_patterns = u.T
        pc_series = np.dot(eof_patterns, data_valid)
        
        # Eigenvalues and explained variance
        nt = data_valid.shape[1]
        eigenvalues = s**2 / nt
        explained_variance = eigenvalues / np.sum(eigenvalues) * 100
        
        # Estimate degrees of freedom (North test)
        phi_L, phi_0, dof = self._estimate_dof(data_valid, L=1)
        eigenvalue_errors = explained_variance * np.sqrt(2 / dof)
        
        return {
            'eof_patterns': eof_patterns,
            'pc_series': pc_series,
            'eigenvalues': eigenvalues,
            'explained_variance': explained_variance,
            'eigenvalue_errors': eigenvalue_errors,
            'degrees_of_freedom': dof,
            'valid_mask': valid_mask,
            'phi_0': phi_0,
            'phi_L': phi_L
        }
    
    def _compute_eof_xeofs(self, data: xr.DataArray, n_modes: int = 10) -> Dict:
        """
        Compute EOF using xeofs library.
        
        Parameters
        ----------
        data : xr.DataArray
            Input data with dimensions (time, level, lat, lon)
        n_modes : int
            Number of EOF modes to compute
            
        Returns
        -------
        dict
            EOF results including patterns, PCs, explained variance
        """
        print(f"   🔄 Processing with xeofs (n_modes={n_modes})...")
        
        # Detect dimension names (robust to various CMIP/world variations)
        detected = self._detect_dims(data)
        time_dim = detected['time']
        level_dim = detected['level']
        lat_dim = detected['lat']
        lon_dim = detected['lon']

        if not all([time_dim, lat_dim, lon_dim]):
            raise ValueError(f"Cannot identify required time/lat/lon dimensions in {list(data.dims)}")

        print(f"   Detected: time={time_dim}, level={level_dim}, lat={lat_dim}, lon={lon_dim}")

        # Build expected order; include level if present
        expected_order = [time_dim, lat_dim, lon_dim]
        if level_dim is not None:
            # put level after time for xeofs expectation (time, level, lat, lon)
            expected_order = [time_dim, level_dim, lat_dim, lon_dim]

        if list(data.dims) != expected_order:
            print(f"   Transposing from {list(data.dims)} to {expected_order}")
            data = data.transpose(*expected_order)
        
        # Stack spatial and temporal dimensions
        # Stack spatial and temporal dims into a single sample axis
        # If level exists, keep it as 'feature' (vertical features)
        if level_dim is not None:
            stacked = data.stack(sample=(time_dim, lat_dim, lon_dim))
            stacked = stacked.rename({level_dim: 'feature'})
        else:
            # 3D data (time, lat, lon)
            stacked = data.stack(sample=(time_dim, lat_dim, lon_dim))
            # create a dummy 'feature' dimension with size 1 for compatibility
            stacked = stacked.expand_dims({'feature': 1})
        
        print(f"   Stacked shape: {stacked.shape}")
        
        # Fit EOF model
        model = self.xe.single.EOF(n_modes=n_modes, check_nans=False)
        model.fit(stacked, dim="sample")
        
        # Extract results
        components = model.components()  # EOF patterns
        scores = model.scores()  # PCs
        explained_variance = model.explained_variance_ratio()
        
        print(f"   EOF components shape: {components.shape}")
        print(f"   PC scores shape: {scores.shape}")
        
        # Unstack scores
        scores = scores.unstack('sample')
        
        # Rename for consistency
        eofs = components.rename({'feature': level_dim})
        pcs = scores
        
        # Convert to numpy for consistency with SVD method
        eof_patterns = eofs.values
        explained_variance_pct = explained_variance.values * 100
        
        return {
            'eof_patterns': eof_patterns,
            'pc_series': None,  # Not directly comparable to SVD PCs
            'pc_scores': pcs,   # Keep xarray format
            'eofs': eofs,       # Keep xarray format
            'eigenvalues': None,
            'explained_variance': explained_variance_pct,
            'eigenvalue_errors': None,
            'degrees_of_freedom': None,
            'valid_mask': None,
            'phi_0': None,
            'phi_L': None
        }
    
    def _estimate_dof(self, data: np.ndarray, L: int = 1) -> Tuple[float, float, float]:
        """
        Estimate degrees of freedom for North test.
        
        Parameters
        ----------
        data : np.ndarray
            Input data
        L : int
            Lag for autocorrelation
            
        Returns
        -------
        tuple
            (phi_L, phi_0, degrees_of_freedom)
        """
        nt = data.shape[1]
        
        # Calculate lag-L autocorrelation
        B = 0
        for k in range(L, nt - L):
            B += np.sum(data[:, k] * data[:, k + L])
        
        phi_L = B / (nt - 2 * L)
        phi_0 = np.sum(data**2) / nt
        r_L = phi_L / phi_0
        
        # Degrees of freedom
        dof = (1 - r_L**2) / (1 + r_L**2) * nt
        
        return phi_L, phi_0, dof
    
    def fit(self, 
            data: xr.DataArray,
            time_slice: slice = None,
            lat_slice: slice = None,
            level_slice: slice = None,
            n_harmonics: int = 3,
            n_modes: int = 10) -> Dict:
        """
        Complete EOF analysis pipeline.
        
        This method handles data preprocessing (masking, climatology removal),
        performs EOF decomposition, and stores results.
        
        Parameters
        ----------
        data : xr.DataArray
            Input data with dimensions including time, level, lat, lon
        time_slice : slice, optional
            Time period to analyze (e.g., slice('1980', '1990'))
        lat_slice : slice, optional
            Latitude range (e.g., slice(-15, 15))
        level_slice : slice, optional
            Pressure level range (e.g., slice(1000, 100))
        n_harmonics : int, default=3
            Number of harmonics to retain in climatology
        n_modes : int, default=10
            Number of EOF modes to compute
            
        Returns
        -------
        dict
            Complete analysis results including EOF patterns, PCs, 
            explained variance, and metadata
            
        Examples
        --------
        >>> analyzer = EOFAnalyzer(method='svd')
        >>> data = xr.open_dataarray('omega.nc')
        >>> results = analyzer.fit(data, 
        ...                        time_slice=slice('1980', '1990'),
        ...                        lat_slice=slice(-15, 15),
        ...                        n_harmonics=3,
        ...                        n_modes=4)
        >>> print(f"Explained variance: {results['explained_variance'][:4]}")
        """
        print("="*70)
        print(f"🚀 EOF Analysis (method={self.method})")
        print("="*70)
        
        # Apply selections
        print("📦 Preparing data...")

        # Detect dims early to apply slices in a robust manner
        detected = self._detect_dims(data)
        time_dim = detected['time']
        level_dim = detected['level']
        lat_dim = detected['lat']
        lon_dim = detected['lon']

        try:
            if time_slice is not None and time_dim is not None and time_dim in data.dims:
                data = data.sel({time_dim: time_slice})
            if lat_slice is not None and lat_dim is not None and lat_dim in data.dims:
                data = data.sel({lat_dim: lat_slice})
            if level_slice is not None and level_dim is not None and level_dim in data.dims:
                data = data.sel({level_dim: level_slice})
        except Exception as e:
            print(f"   ⚠️  Selection warning: {e}")

        print(f"   Data shape: {data.shape}")
        print(f"   Dimensions: {list(data.dims)}")

        # Re-detect in case selection changed dimension order/names
        detected = self._detect_dims(data)
        time_dim = detected['time']
        level_dim = detected['level']
        lat_dim = detected['lat']
        lon_dim = detected['lon']

        # Ensure correct dimension order: time, (level), lat, lon
        if time_dim is None or lat_dim is None or lon_dim is None:
            raise ValueError(f"Data must contain time/lat/lon dims; found {list(data.dims)}")

        expected_order = [time_dim, lat_dim, lon_dim]
        if level_dim is not None:
            expected_order = [time_dim, level_dim, lat_dim, lon_dim]

        if list(data.dims) != expected_order:
            print(f"   Transposing from {list(data.dims)} to {expected_order}")
            data = data.transpose(*expected_order)
        
        # Apply land mask if requested
        if self.apply_land_mask:
            print("🗺️  Applying land/ocean mask...")
            # Prepare an indexer for a single time & level (if present)
            mask_indexer = {}
            mask_indexer[time_dim] = 0
            if level_dim is not None:
                mask_indexer[level_dim] = 0

            mask = self.create_land_mask(data.isel(mask_indexer))
            mask_expanded = mask.broadcast_like(data)
            data = data.where(mask_expanded)

            mask_type = "Ocean" if self.ocean_only else "Land"
            # Compute retained valid points for the first time & level (or 2D slice)
            indexer = {time_dim: 0}
            if level_dim is not None:
                indexer[level_dim] = 0

            valid_slice = data.isel(indexer)
            valid_points = (~np.isnan(valid_slice)).sum().values
            total_points = valid_slice.size
            print(f"   {mask_type} points retained: {valid_points}/{total_points} "
                  f"({valid_points/total_points*100:.1f}%)")
        
        # Calculate climatology and anomalies
        print("📊 Calculating climatology and anomalies...")
        
        # Check if time coordinate has datetime-like properties
        try:
            # Try to use dayofyear grouping (for real datetime data)
            clim = data.groupby('time.dayofyear').mean(dim='time')
            clim_smoothed = self.extract_low_harmonics(clim, n_harm=n_harmonics)
            anomalies = data.groupby('time.dayofyear') - clim_smoothed
        except (AttributeError, KeyError):
            # Fallback for non-datetime time coordinates (e.g., synthetic data)
            print("   ⚠️  Time coordinate doesn't have dayofyear, using simple mean")
            clim = data.mean(dim='time')
            anomalies = data - clim
            clim_smoothed = clim
        
        print(f"   Anomalies shape: {anomalies.shape}")
        
        # Perform EOF decomposition based on method
        if self.method == 'svd':
            print("🔬 Performing SVD-based EOF decomposition...")
            
            # Reshape for SVD
            detected_anom = self._detect_dims(anomalies)
            atime = detected_anom['time']
            alevel = detected_anom['level']
            alat = detected_anom['lat']
            alon = detected_anom['lon']

            if alevel is not None:
                nt = anomalies.sizes[atime]
                nlev = anomalies.sizes[alevel]
                nlat = anomalies.sizes[alat]
                nlon = anomalies.sizes[alon]
                data_reshaped = anomalies.values.transpose(1, 0, 2, 3).reshape(
                    nlev, nt * nlat * nlon, order='F')
            else:
                # 3D data (time, lat, lon) -> treat level as singleton
                nt = anomalies.sizes[atime]
                nlat = anomalies.sizes[alat]
                nlon = anomalies.sizes[alon]
                data_reshaped = anomalies.values.transpose(0, 1, 2).reshape(
                    1, nt * nlat * nlon, order='F')

            eof_results = self._compute_eof_svd(data_reshaped)
            
        elif self.method == 'xeofs':
            print("🔬 Performing xeofs-based EOF decomposition...")
            eof_results = self._compute_eof_xeofs(anomalies, n_modes=n_modes)
        
        # Store results
        self.eof_results = {
            **eof_results,
            'pressure_levels': data[level_dim].values if level_dim is not None and level_dim in data.dims else None,
            'original_data': data,
            'anomalies': anomalies,
            'climatology': clim_smoothed,
            'n_harmonics': n_harmonics,
            'mask_applied': self.apply_land_mask,
            'ocean_only': self.ocean_only if self.apply_land_mask else None,
            'method': self.method
        }
        
        print("="*70)
        print("✅ EOF analysis completed!")
        print(f"   Method: {self.method}")
        print(f"   Modes computed: {len(eof_results['explained_variance'])}")
        print(f"   Top 4 explained variance: {eof_results['explained_variance'][:4]}")
        print("="*70)
        
        return self.eof_results
    
    def plot_vertical_profiles(self, 
                               n_modes: int = 4, 
                               figsize: Tuple[int, int] = (12, 10),
                               save_path: Optional[str] = None,
                               normalize: bool = True) -> plt.Figure:
        """
        Plot EOF vertical profiles.
        
        Parameters
        ----------
        n_modes : int
            Number of modes to plot
        figsize : tuple
            Figure size (width, height)
        save_path : str, optional
            Path to save the figure
        normalize : bool
            Whether to normalize EOF patterns
            
        Returns
        -------
        matplotlib.figure.Figure
            The created figure
            
        Examples
        --------
        >>> analyzer.fit(data)
        >>> fig = analyzer.plot_vertical_profiles(n_modes=4, figsize=(12, 8))
        >>> plt.show()
        """
        if not self.eof_results:
            raise ValueError("No EOF results available. Run fit() first.")
        
        # Get results
        eof_patterns = self.eof_results['eof_patterns']
        explained_var = self.eof_results['explained_variance']
        pressure_levels = self.eof_results['pressure_levels']

        if pressure_levels is None:
            raise ValueError("No vertical level coordinate found in EOF results. Cannot plot vertical profiles for 2D data.")
        
        # Normalize if requested
        if normalize:
            eof_patterns_plot = eof_patterns.copy()
            for i in range(eof_patterns_plot.shape[0]):
                max_val = np.abs(eof_patterns_plot[i, :]).max()
                if max_val > 0:
                    eof_patterns_plot[i, :] /= max_val
        else:
            eof_patterns_plot = eof_patterns
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        plt.rcParams.update({'font.size': 14})
        
        # Colors
        colors = ['#2E8B57', '#FF8C00', '#4169E1', '#DC143C']
        
        # Plot EOF profiles
        for i in range(min(n_modes, len(explained_var))):
            label = f'EOF{i+1} ({explained_var[i]:.1f}%)'
            ax.plot(eof_patterns_plot[i, :], pressure_levels, 
                   label=label, color=colors[i], 
                   linewidth=3, marker='o', markersize=4)
        
        # Zero line
        ax.axvline(x=0, color='black', linestyle=':', alpha=0.7, linewidth=1)
        
        # Formatting
        ax.set_ylim(pressure_levels.max(), pressure_levels.min())
        ax.set_xlim(-1.0 if normalize else None, 1.0 if normalize else None)
        ax.set_yscale("log")
        
        # Pressure ticks
        major_ticks = np.array([1000, 850, 700, 600, 500, 400, 300, 200, 100])
        valid_ticks = major_ticks[
            (major_ticks >= pressure_levels.min()) & 
            (major_ticks <= pressure_levels.max())
        ]
        ax.set_yticks(valid_ticks)
        ax.get_yaxis().set_major_formatter(plt.ScalarFormatter())
        
        # Labels
        xlabel = 'Normalized EOF Amplitude' if normalize else 'EOF Amplitude'
        ax.set_xlabel(xlabel, fontsize=16, fontweight='bold')
        ax.set_ylabel('Pressure (hPa)', fontsize=16, fontweight='bold')
        ax.legend(loc='best', fontsize=12, framealpha=0.9)
        
        # Grid
        ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
        
        # Ticks
        ax.tick_params(which='both', direction='inout', length=6, width=1.2)
        ax.tick_params(axis='both', which='major', labelsize=12)
        
        # Title
        mask_info = ""
        if self.eof_results['mask_applied']:
            mask_info = f" ({'Ocean Only' if self.eof_results['ocean_only'] else 'Land Only'})"
        
        title = f'EOF Vertical Profiles - {self.method.upper()} Method{mask_info}'
        ax.set_title(title, fontsize=18, fontweight='bold', pad=20)
        
        # Info box
        info_text = (f"Method: {self.method.upper()}\n"
                    f"Harmonics: {self.eof_results['n_harmonics']}\n"
                    f"Total Variance: {np.sum(explained_var[:n_modes]):.1f}%")
        ax.text(0.02, 0.98, info_text, transform=ax.transAxes, fontsize=10,
               verticalalignment='top', 
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"✅ Figure saved to: {save_path}")
            
        return fig
    
    def save_results(self, save_path: str) -> None:
        """
        Save EOF results to file.
        
        Parameters
        ----------
        save_path : str
            Path to save results (.pkl file)
        """
        if not self.eof_results:
            raise ValueError("No EOF results to save. Run fit() first.")
        
        # Prepare data for saving (exclude large DataArrays)
        save_data = {
            'eof_patterns': self.eof_results['eof_patterns'],
            'explained_variance': self.eof_results['explained_variance'],
            'pressure_levels': self.eof_results['pressure_levels'],
            'n_harmonics': self.eof_results['n_harmonics'],
            'mask_applied': self.eof_results['mask_applied'],
            'ocean_only': self.eof_results['ocean_only'],
            'method': self.eof_results['method'],
        }
        
        # Add method-specific data
        if self.method == 'svd':
            save_data.update({
                'pc_series': self.eof_results['pc_series'],
                'eigenvalues': self.eof_results['eigenvalues'],
                'eigenvalue_errors': self.eof_results['eigenvalue_errors'],
                'degrees_of_freedom': self.eof_results['degrees_of_freedom'],
                'valid_mask': self.eof_results['valid_mask'],
            })
        
        # Save to file
        with open(save_path, 'wb') as f:
            pickle.dump(save_data, f, protocol=pickle.HIGHEST_PROTOCOL)
        
        file_size_mb = os.path.getsize(save_path) / 1024 / 1024
        print(f"✅ EOF results saved to: {save_path}")
        print(f"   File size: {file_size_mb:.2f} MB")
    
    def load_results(self, load_path: str) -> Dict:
        """
        Load EOF results from file.
        
        Parameters
        ----------
        load_path : str
            Path to load results from (.pkl file)
            
        Returns
        -------
        dict
            Loaded EOF results
        """
        if not os.path.exists(load_path):
            raise FileNotFoundError(f"Results file not found: {load_path}")
        
        with open(load_path, 'rb') as f:
            loaded_data = pickle.load(f)
        
        self.eof_results = loaded_data
        self.method = loaded_data.get('method', 'svd')
        
        print(f"✅ EOF results loaded from: {load_path}")
        print(f"   Method: {self.method}")
        print(f"   Pressure levels: {len(loaded_data['pressure_levels'])}")
        print(f"   Top 4 explained variance: {loaded_data['explained_variance'][:4]}")
        
        return loaded_data


def quick_eof_analysis(data: xr.DataArray,
                      method: str = 'svd',
                      n_modes: int = 4,
                      time_slice: slice = None,
                      lat_slice: slice = None,
                      level_slice: slice = None,
                      n_harmonics: int = 3,
                      plot: bool = True,
                      save_fig: Optional[str] = None) -> Tuple[EOFAnalyzer, Dict, Optional[plt.Figure]]:
    """
    Quick one-line EOF analysis with optional plotting.
    
    This convenience function performs EOF analysis and optionally creates
    a vertical profile plot in a single call.
    
    Parameters
    ----------
    data : xr.DataArray
        Input data with dimensions (time, level, lat, lon)
    method : str, default='svd'
        Decomposition method: 'svd' or 'xeofs'
    n_modes : int, default=4
        Number of EOF modes to compute and plot
    time_slice : slice, optional
        Time period to analyze
    lat_slice : slice, optional
        Latitude range
    level_slice : slice, optional
        Pressure level range
    n_harmonics : int, default=3
        Number of harmonics for climatology
    plot : bool, default=True
        Whether to create vertical profile plot
    save_fig : str, optional
        Path to save figure
        
    Returns
    -------
    analyzer : EOFAnalyzer
        Fitted analyzer instance
    results : dict
        EOF analysis results
    fig : matplotlib.figure.Figure or None
        Figure object if plot=True, None otherwise
        
    Examples
    --------
    >>> data = xr.open_dataarray('omega.nc')
    >>> analyzer, results, fig = quick_eof_analysis(
    ...     data, 
    ...     method='svd',
    ...     n_modes=4,
    ...     time_slice=slice('1980', '1990'),
    ...     lat_slice=slice(-15, 15),
    ...     plot=True
    ... )
    >>> plt.show()
    """
    # Create analyzer
    analyzer = EOFAnalyzer(method=method)
    
    # Perform analysis
    results = analyzer.fit(
        data=data,
        time_slice=time_slice,
        lat_slice=lat_slice,
        level_slice=level_slice,
        n_harmonics=n_harmonics,
        n_modes=n_modes
    )
    
    # Create plot if requested
    fig = None
    if plot:
        fig = analyzer.plot_vertical_profiles(
            n_modes=n_modes,
            save_path=save_fig
        )
    
    return analyzer, results, fig


# Convenience aliases for backward compatibility
def eof_svd(data: xr.DataArray, **kwargs) -> Tuple[EOFAnalyzer, Dict, Optional[plt.Figure]]:
    """Convenience function for SVD-based EOF analysis."""
    return quick_eof_analysis(data, method='svd', **kwargs)


def eof_xeofs(data: xr.DataArray, **kwargs) -> Tuple[EOFAnalyzer, Dict, Optional[plt.Figure]]:
    """Convenience function for xeofs-based EOF analysis."""
    return quick_eof_analysis(data, method='xeofs', **kwargs)


if __name__ == "__main__":
    print("EOF Analysis Module for Wave Tools")
    print("===================================")
    print("\nAvailable classes:")
    print("  - EOFAnalyzer: Comprehensive EOF analysis")
    print("\nAvailable functions:")
    print("  - quick_eof_analysis: One-line EOF analysis")
    print("  - eof_svd: SVD-based EOF (convenience)")
    print("  - eof_xeofs: xeofs-based EOF (convenience)")
    print("\nExample usage:")
    print("  >>> from wave_tools.eof import EOFAnalyzer")
    print("  >>> analyzer = EOFAnalyzer(method='svd')")
    print("  >>> results = analyzer.fit(data, n_modes=4)")
    print("  >>> analyzer.plot_vertical_profiles()")
