import matplotlib.ticker as ticker
from matplotlib import gridspec
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from wave_tools.utils import filter_series, save_figure
from wave_tools import ma as mp
from wave_tools.utils import get_curve



class Spectrum:
    def __init__(self,  cpd_lines):
        self.max_wn_plot = None
        self.max_freq_plot = None
        self.cpd_lines=cpd_lines
    def set_ax(self, ax, text_size, freq_lines=True,depth=True):
        ax.axvline(x=0, color='k', linestyle='--')

        
        # ax.text(self.max_wn_plot-2*0.25*self.max_wn_plot,-0.03,'EASTWARD',fontweight='bold',fontsize=text_size-2)
        # ax.text(-self.max_wn_plot+0.25*self.max_wn_plot,-0.03,'WESTWARD',fontweight='bold',fontsize=text_size-2)
        ax.set_xlim((-self.max_wn_plot,self.max_wn_plot))
        ax.set_ylim((0.02,self.max_freq_plot))
        
        if self.max_wn_plot is not None and self.max_freq_plot is not None:
            ax.set_xlim((-self.max_wn_plot, self.max_wn_plot))
            ax.set_ylim((0.02, self.max_freq_plot))

        if freq_lines:
            # Assuming self.freq_lines and self.cpd_lines are defined elsewhere
            for d in self.cpd_lines:
                if (1./d) <= self.max_freq_plot:
                    ax.axhline(y=1./d, color='k', linestyle='--',linewidth=0.5)
                    ax.text(-self.max_wn_plot+0.8, (1./d+0.01), str(d)+' days', color='k',
                            size=text_size-6, bbox={'facecolor': 'w', 'alpha': 0.9, 'edgecolor': 'none'})
                    
        if depth:
            kw_x,kw_y = get_curve()
            ax.plot(kw_x[0], kw_y[0], 'green', linewidth=1.2, linestyle='solid',zorder=5)
            # Define the range for filtering
            ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
            left, width = .25, .5
            bottom, height = .25, .5
            right = left + width
            top = bottom + height
            ax.text(right+10, 0.29  * (bottom+top),'kelvin', ha="center",va="center",size=text_size-6,
                    bbox={'facecolor':'w','alpha':0.9,'edgecolor':'none'})
            ax.text(right+5.6, 0.4 * (bottom+top),'h=90', ha="center",va="center",size=text_size-6,
                    bbox={'facecolor':'w','alpha':0.9,'edgecolor':'none'})
            ax.text(right+10.5, 0.38 * (bottom+top),'h=25', ha="center",va="center",size=text_size-6,
                    bbox={'facecolor':'w','alpha':0.9,'edgecolor':'none'})
            ax.text(right+9.9, 0.2 * (bottom+top),'h=8', ha="center",va="center",size=text_size-6,
                    bbox={'facecolor':'w','alpha':0.9,'edgecolor':'none'})
            
    def plot_cmip_future(self, cmip_f,wn,fq, figsize=None, text_size=12, contour_range=[1, 2, 0.1], he=[90, 25, 8],
                      min_contour_range_lines=1.1, matsuno_lines=True, meridional_modes=[1], cmap = 'bwr',
                      max_freq_plot=None, max_wn_plot=None, freq_lines=True, cpd_lines=[3, 6, 30], 
                      labels=False, norm=None, ofil=None):
        
        self.max_wn_plot = max_wn_plot
        self.max_freq_plot = max_freq_plot
        
        contour_range = np.arange(contour_range[0], contour_range[1], contour_range[2])

        contour_range_lines = contour_range[contour_range >= min_contour_range_lines]
        
        fb = [0, .5]  # frequency bounds for plot

        wavefft = wn
        freqfft = fq
        
        X, Y = np.meshgrid(wavefft, freqfft)
        
        
        gs = gridspec.GridSpec(1, 3)

        fig = plt.figure(figsize=figsize,dpi=200)  # 设置Figure的尺寸
        
        modelnames = ['(a) Good','(b) Poor','(c) GPCP']

        for idx,data in enumerate(cmip_f):
            print(idx)
            # 计算子图的行列索引
            row = idx // 3
            col = idx % 3
            print(row,col)
            # 在当前子图位置创建一个Axes对象
            ax = fig.add_subplot(gs[row, col])
        
            if row==0:
                ax.set_xlabel('Zonal Wavenumber', size=16,)
            else:
                ax.set_xticks([])
            if col==0:
                ax.set_ylabel('Frequency (CPD)', size=16,)
            else:
                ax.set_yticks([])
             
           
            ax.set_title(modelnames[idx], size=12,  loc='left')

            sym = cmip_f[idx].transpose().sel(frequency=slice(*fb), 
                                                               wavenumber=slice(-self.max_wn_plot, self.max_wn_plot)) 
            cset_0  =  sym.T.plot.contourf(ax=ax, levels=contour_range, 
                                           extend='neither',
                                           norm=norm,
                                           cmap=cmap,add_colorbar=False,add_labels=False
                                            )
            cset1_0 = sym.T.plot.contour( ax=ax,levels= contour_range_lines, colors='k',add_colorbar=False,add_labels=False)
            self.set_ax(ax, text_size, freq_lines)
            if matsuno_lines:
                matsuno_modes = mp.matsuno_modes_wk(he=he,n=meridional_modes,max_wn=self.max_wn_plot)
                # print(matsuno_modes)
                kelvin_series_90 = filter_series(matsuno_modes[90]['Kelvin(he=90m)'], 2, 5.2)
                kelvin_series_8 = filter_series(matsuno_modes[8]['Kelvin(he=8m)'], 2.6, 14)
                for key in matsuno_modes:
                    # print(key)
                    ax.plot(matsuno_modes[key]['Kelvin(he={}m)'.format(key)],color='k',linestyle='-')
                    ax.plot(matsuno_modes[key]['ER(n=1,he={}m)'.format(key)],color='k',linestyle='-')
                  
                    ax.plot(matsuno_modes[key]['Kelvin(he={}m)'.format(key)],color='k',linestyle='-')
                    ax.plot(matsuno_modes[key]['ER(n=1,he={}m)'.format(key)],color='k',linestyle='-')
                
                    ax.plot(matsuno_modes[key]['Kelvin(he={}m)'.format(key)],color='k',linestyle='-')
                    ax.plot(matsuno_modes[key]['ER(n=1,he={}m)'.format(key)],color='k',linestyle='-')

                    ax.plot(kelvin_series_90, color='g', linestyle='-')
                    ax.plot(kelvin_series_8, color='g', linestyle='-')

                if labels:
                    key = list(matsuno_modes.keys())[len(list(matsuno_modes.keys()))//2]
                    
                    wn = matsuno_modes[key].index.values
                    
                    k = int((len(wn)/2)+0.3*(len(wn)/2))
                    k, = np.where(wn == wn[k])[0]
                    k = int(0.7*(len(wn)/2))
                    k = np.where(wn == wn[k])[0]
                    ax.text(wn[k]+0.4,matsuno_modes[key]['ER(n=1,he={}m)'.format(key)].iloc[k]+0.02,'ER', \
                    bbox={'facecolor':'w','alpha':0.9,'edgecolor':'none'},fontsize=text_size-6)
                        
                        
        cbar_ax = fig.add_axes([0.25, 0.002, 0.55, 0.02])  # left\bottom\length\height              
      
      
        cbar=plt.colorbar(cset_0, cax= cbar_ax, orientation='horizontal',
                          ticks=(1,1.2,1.4,1.6,1.8,2,),
                          extend='neither', shrink=0.5, pad=0.05)
        cbar.ax.tick_params(which='both', direction='in',length=0)
        plt.tight_layout()
        if ofil is not None:
            save_figure(fig,filename=ofil)