# cckw_tools.py

import numpy as np
import matplotlib.pyplot as plt
import os
from typing import Optional, List, Tuple
from .utils import get_curve

def plot_cckw_envelope(
    he: Optional[List[float]] = None,
    fmax: Optional[List[float]] = None,
    savepath: Optional[str] = None,
    dpi: int = 200
) -> None:
    """
    绘制 CCKW 包络示意图，可选择保存。
    """
    kw_x, kw_y = get_curve(he=he, fmax=fmax)

    if he is None:
        he_all = np.array([8, 25, 90])
    else:
        he_all = np.array(he)

    g = 9.8
    re = 6371e3
    s2d = 86400
    cp = (g * he_all) ** 0.5
    zwnum_goal = np.pi * re / cp / s2d

    title = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)']

    plt.rcParams.update({'font.size': 6.5})
    fig, axs = plt.subplots(2, 3, figsize=(5.8, 3.9), dpi=dpi)
    plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.15, wspace=0.15, hspace=0.22)

    for idx, ax in enumerate(axs.flat):
        i, v = divmod(idx, 3)

        for dd, d in enumerate([3, 6, 20]):
            ax.plot([-20, 20], [1 / d, 1 / d], 'k', linewidth=0.5, linestyle=':')
            ax.text(-14.8, 1 / d + 0.01, ['3d', '6d', '20d'][dd], fontsize=6)

        for hh in range(len(he_all)):
            ax.plot([0, zwnum_goal[hh]], [0, 0.5], 'grey', linewidth=0.5, linestyle='dashed')

        ax.plot([0, 0], [0, 0.5], 'k', linewidth=0.5, linestyle=':')

        ax.plot(kw_x[0], kw_y[0], 'purple', linewidth=1.2, linestyle='solid')

        ax.set_title(title[idx], pad=3, loc='center', fontsize=9)
        if v == 0:
            ax.set_ylabel('Frequency (1/day)', fontsize=7.5)
        if i == 1:
            ax.set_xlabel('Zonal wavenumber', fontsize=7.5)

        ax.set_xlim([-20, 20])
        ax.set_ylim([0, 0.5])
        ax.set_xticks(np.arange(-20, 21, 5))
        ax.set_yticks(np.arange(0, 0.55, 0.05))
        ax.tick_params(labelsize=6, direction='in', top=True, right=True)

        if v != 0:
            ax.tick_params(labelleft=False)

    if savepath:
        folder = os.path.dirname(savepath)
        if folder and not os.path.exists(folder):
            os.makedirs(folder)
            print(f'Created folder: {folder}')
        fig.savefig(f"{savepath}.pdf", dpi=dpi, bbox_inches='tight', format='pdf')
        print(f'Figure saved at: {savepath}.pdf')

    plt.show()
    plt.close()



