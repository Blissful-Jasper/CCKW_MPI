import matplotlib.ticker as ticker
from matplotlib import gridspec
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from wave_tools.utils import filter_series, save_figure
from wave_tools import ma as mp
from wave_tools.utils import get_curve

class SpectrumPlotter:
    def __init__(self,  cpd_lines,norm=None):
        self.max_wn_plot = None
        self.max_freq_plot = None
        self.cpd_lines=cpd_lines
        self.norm = norm


    def set_ax(self, ax, text_size, freq_lines=True,depth=True):
        ax.axvline(x=0, color='k', linestyle='--')
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
            # ax.axhline(y=1./3,  xmin=0.672,xmax=0.964,color='g', linestyle='-')
            # ax.axhline(y=1./20, xmin=0.57,xmax=0.585,color='g', linestyle='-')
            
            # ax.axvline(x=2, ymin=0.066,ymax=0.23, color='g', linestyle='-',zorder=5)
            # ax.axvline(x=14,ymin=0.52,ymax=0.65,color='g', linestyle='-')
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
                      labels=False,ofil=None):
        
        self.max_wn_plot = max_wn_plot
        self.max_freq_plot = max_freq_plot
        
        contour_range = np.arange(contour_range[0], contour_range[1], contour_range[2])

        contour_range_lines = contour_range[contour_range >= min_contour_range_lines]
        
        fb = [0, .5]  # frequency bounds for plot

        wavefft = wn
        freqfft = fq
        
        X, Y = np.meshgrid(wavefft, freqfft)
        
        
        gs = gridspec.GridSpec(1, 1)
        plt.rcParams['font.family'] = 'Arial'

        fig = plt.figure(figsize=figsize)  # 设置Figure的尺寸
        
        modelnames = ['GPCP-obs']
        kw_x, kw_y = get_curve()

        for idx,data in enumerate(cmip_f):
            print(idx)
            # 计算子图的行列索引
            row = 0
            col = 0
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
            fig_width = fig.get_figwidth()
            fig_height = fig.get_figheight()
            ax.set_title('Westward', size=12, loc='left')
            ax.set_title('Eastward', size=12, loc='right')
            ax.set_title('Symmetric/Background', size=12,)
            ax.text(x=0.5,y=1.08,s= 'GPCP-obs',
                    size=15, horizontalalignment='center',
                        verticalalignment='center', transform=ax.transAxes)

            sym = cmip_f[idx].transpose().sel(frequency=slice(*fb), 
                                                               wavenumber=slice(-self.max_wn_plot, self.max_wn_plot)) 
            cset_0  =  sym.T.plot.contourf(ax=ax, levels=contour_range, 
                                           extend='neither',
                                           norm=self.norm,
                                           cmap=cmap,add_colorbar=False,add_labels=False
                                            )
            cset1_0 = sym.T.plot.contour( ax=ax,levels= contour_range_lines, colors='k',add_colorbar=False,add_labels=False)

            self.set_ax(ax, text_size, freq_lines)
            ax.plot(kw_x[0], kw_y[0], 'g', linewidth=1.2, linestyle='solid',zorder=5)
            if matsuno_lines:
                matsuno_modes = mp.matsuno_modes_wk(he=he,n=meridional_modes,max_wn=self.max_wn_plot)
                # print(matsuno_modes)
                kelvin_series_90 = filter_series(matsuno_modes[90]['Kelvin(he=90m)'], 2, 5.2)
                kelvin_series_8 = filter_series(matsuno_modes[8]['Kelvin(he=8m)'], 2.6, 14)
                
                for key in matsuno_modes:
                    print(key)
                    ax.plot(matsuno_modes[key]['Kelvin(he={}m)'.format(key)],color='k',linestyle='-')
                    ax.plot(matsuno_modes[key]['ER(n=1,he={}m)'.format(key)],color='k',linestyle='-')
                  
                    ax.plot(matsuno_modes[key]['Kelvin(he={}m)'.format(key)],color='k',linestyle='-')
                    ax.plot(matsuno_modes[key]['ER(n=1,he={}m)'.format(key)],color='k',linestyle='-')
                
                    ax.plot(matsuno_modes[key]['Kelvin(he={}m)'.format(key)],color='k',linestyle='-')
                    ax.plot(matsuno_modes[key]['ER(n=1,he={}m)'.format(key)],color='k',linestyle='-')
                    
                    
                    ax.plot(matsuno_modes[key]['WIG(n=1,he={}m)'.format(key)],color='k',linestyle='-')
                    ax.plot(matsuno_modes[key]['WIG(n=1,he={}m)'.format(key)],color='k',linestyle='-')
                    ax.plot(matsuno_modes[key]['WIG(n=1,he={}m)'.format(key)],color='k',linestyle='-')
                    
                    ax.plot(matsuno_modes[key]['EIG(n=1,he={}m)'.format(key)],color='k',linestyle='-')
                    ax.plot(matsuno_modes[key]['EIG(n=1,he={}m)'.format(key)],color='k',linestyle='-')
                    
                    ax.plot(kelvin_series_90, color='g', linestyle='-')
                    ax.plot(kelvin_series_8, color='g', linestyle='-')

                if labels:
                    key = list(matsuno_modes.keys())[len(list(matsuno_modes.keys()))//2]
                    
                    wn = matsuno_modes[key].index.values
                    
                    k = int((len(wn)/2)+0.3*(len(wn)/2))
                    k, = np.where(wn == wn[k])[0]
                    k = int(0.7*(len(wn)/2))
                    k = np.where(wn == wn[k])[0]
                    ax.text(wn[k]+0.4,matsuno_modes[key]['ER(n=1,he={}m)'.format(key)].iloc[k]+0.02,'n=1 ER', \
                    bbox={'facecolor':'w','alpha':1,'edgecolor':'none'},fontsize=text_size-6)
                        
                    ax.text(wn[k]+0.4,matsuno_modes[key]['WIG(n=1,he={}m)'.format(key)].iloc[k]+0.02,'n=1 WIG', \
                        bbox={'facecolor':'w','alpha':1,'edgecolor':'none'},fontsize=text_size-6)
                        
                        
        cbar_ax = fig.add_axes([0.25, 0.002, 0.55, 0.02])  # left\bottom\length\height              
      
      
        cbar=plt.colorbar(cset_0, cax= cbar_ax, orientation='horizontal',
                          ticks=(1,1.2,1.4,1.6,1.8,2,),
                          extend='neither', shrink=0.5, pad=0.05)
        cbar.ax.tick_params(which='both', direction='in',length=0,
                            bottom=True, left=True, right=True, top=True)
        
        if ofil is not None:
            save_figure(fig,filename=ofil)