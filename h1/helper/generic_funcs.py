import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MaxNLocator
import numpy as np

def setup_publication_image(fig, height_prop = 1./1.618, single_col = True, fig_width = None, replacement_kwargs = None, fig_height = None):
    cm_to_inch=0.393701
    if replacement_kwargs == None: replacement_kwargs = {}
    mpl.rcParams['font.size']=8.0
    mpl.rcParams['axes.titlesize']=8.0#'medium'
    mpl.rcParams['xtick.labelsize']=8.0
    mpl.rcParams['ytick.labelsize']=8.0
    mpl.rcParams['lines.markersize']=5.0
    mpl.rcParams['savefig.dpi']=300
    for i in replacement_kwargs.keys():
        mpl.rcParams[i]=replacement_kwargs[i]
    if fig_width==None:
        if single_col:
            fig_width = 8.48*cm_to_inch
        else:
            fig_width = 8.48*cm_to_inch*2
    else:
        fig_width = fig_width * cm_to_inch
    fig.set_figwidth(fig_width)
    if fig_height == None:
        fig.set_figheight(fig_width * height_prop)
    else:
        fig.set_figheight(fig_height * cm_to_inch)

def setup_axis_publication(ax, n_xticks = None, n_yticks = None):
    if n_yticks!= None: ax.yaxis.set_major_locator(MaxNLocator(n_yticks))
    if n_xticks!= None: ax.xaxis.set_major_locator(MaxNLocator(n_xticks))

def cbar_ticks(cbar_ax, n_ticks = 5):
    cbar_ax.locator = MaxNLocator(nbins = n_ticks)
    cbar_ax.update_ticks()

def create_cbar_ax(original_ax, pad = 3, loc = "right", prop = 5):
    divider = make_axes_locatable(original_ax)
    return divider.append_axes(loc, "{}%".format(prop), pad="{}%".format(pad))

def modify_axis(ax, modifiers):
    for k,v in modifiers.items():
        tmp = getattr(ax, k)
        tmp(v)

def a_b_c(ax, x_locs = 0.1, y_locs = 0.8, kwargs = None):
    if kwargs == None: kwargs = {}
    if ax.__class__ == np.ndarray: ax = ax.tolist()
    if x_locs.__class__!=list:x_locs = [x_locs]*len(ax)
    if y_locs.__class__!=list:y_locs = [y_locs]*len(ax)
    start = ord('a')
    for i, (ax_tmp, x_loc, y_loc) in enumerate(zip(ax, x_locs, y_locs)):
        x_lim = ax_tmp.get_xlim()
        y_lim = ax_tmp.get_ylim()
        ax_tmp.text(x_lim[0]+x_loc*(x_lim[1]-x_lim[0]),
                    y_lim[0]+y_loc*(y_lim[1]-y_lim[0]),
                    '({})'.format(chr(start+i)), **kwargs)
