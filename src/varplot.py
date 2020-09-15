#!/usr/bin/env python3

""" varplot.py
    Oscar Morata 2020

    Functions to configure and plot histograms and xyplots
"""

import astropy
from astropy.visualization import hist
import matplotlib.pyplot as plt
import numpy as np


def set_cnfgvalue(cfg, tag, default):
    """ Read configuration value of tag, or use default value"""

    if tag in cfg:
        val = cfg[tag]
    else :
        val = default

    return val



def read_plot_configuration(cfg):
    """Read the plot configuration"""
    
    plot_cfg = {}

    plot_cfg['type'] = set_cnfgvalue(cfg,'type', 'xyplot')
    plot_cfg['size'] = set_cnfgvalue(cfg,'size', [7,7])
    plot_cfg['xlog'] = set_cnfgvalue(cfg, 'xlog', False)
    plot_cfg['ylog'] = set_cnfgvalue(cfg, 'ylog', False)
    plot_cfg['xlabel'] = set_cnfgvalue(cfg, 'xlabel', '')
    plot_cfg['ylabel'] = set_cnfgvalue(cfg, 'ylabel', '')
    plot_cfg['title'] = set_cnfgvalue(cfg, 'title', '')
    plot_cfg['colors'] = set_cnfgvalue(cfg, 'colors', ['skyblue', 'red',
                                                      'green'])
    plot_cfg['dpi'] = set_cnfgvalue(cfg, 'dpi', 100)

    return plot_cfg



def read_histo_configuration(cfg):
    """ Read the configuration for astropy histo"""

    hcfg = {}
    
    hcfg['meth'] = set_cnfgvalue(cfg, 'method', 'scott')
    hcfg['httype'] = set_cnfgvalue(cfg, 'type', 'stepfilled')
    hcfg['alpha'] = set_cnfgvalue(cfg, 'alpha', 0.8)
    hcfg['dens'] = set_cnfgvalue(cfg, 'density', False)
    hcfg['stack'] = set_cnfgvalue(cfg, 'stacked', True)
    hcfg['cumul'] = set_cnfgvalue(cfg, 'cumulative', False)
    hcfg['ylog'] = set_cnfgvalue(cfg, 'log', False)
    hcfg['outfile'] = set_cnfgvalue(cfg, 'outfile', 'out_histo.pdf')

    return hcfg



def read_xyplot_configuration(cfg):
    """ Read the configuration for xy plot"""

    xycfg = {}

    xycfg['symbols'] = set_cnfgvalue(cfg, 'symbols', ['.', '.', '.', '.'])
    xycfg['outfile'] = set_cnfgvalue(cfg, 'outfile', 'out_xyplot.pdf')

    return xycfg


    
def plot_histo(data, cfg, plcfg, outd):
    """Plot histogram of data

    Using configuration in cfg, and output to outf.
    It returns outbin for further processing.
    """

    histocfg = read_histo_configuration(cfg)

    fig, ax = plt.subplots(1,1,figsize=plcfg['size'])

    if plcfg['xlog'] :
        plt.xscale('log')

    if plcfg['xlabel']:
        plt.xlabel(plcfg['xlabel'])
    if plcfg['ylabel'] :
        plt.ylabel(plcfg['ylabel'])
    if plcfg['title']:
        plt.title(plcfg['title'])

    plmin = np.min(data[0][0])
    plmax = np.max(data[0][0])

    if histocfg['meth'] == 'scott':
        w, bins = astropy.stats.scott_bin_width(data[0][0], return_bins=True)

    elif histocfg['meth'] == 'freedman' :
        w, bins = astropy.stats.freedman_bin_width(data[0][0],
                                                   return_bins=True)
    elif histocfg['meth']== 'knuth' :
        w, bins = astropy.stats.knuth_bin_width(data[0][0].ravel(),
                                                return_bins=True, quiet=False)

    for ds in range(np.shape(data)[0]):
        shdata = data[ds][0]
        outbin = hist(shdata, bins=bins, ax=ax, histtype=histocfg['httype'],
                      alpha=histocfg['alpha'], density=histocfg['dens'],
                      range=(plmin, plmax),
                      stacked=histocfg['stack'], cumulative=histocfg['cumul'],
                      log=histocfg['ylog'],
                      color=plcfg['colors'][ds])
 
    fig.savefig(outd+cfg['outfile'], dpi=plcfg['dpi'])

    return outbin



def plot_xy(data, cfg, plcfg, outd):
    """ Plot an x,y plot"""

    xycfg = read_xyplot_configuration(cfg)
    
    fig = plt.figure(figsize=plcfg['size'])

    if plcfg['xlog'] :
        plt.xscale('log')

    if plcfg['ylog'] :
        plt.yscale('log')
        
    if plcfg['xlabel']:
        plt.xlabel(plcfg['xlabel'])
    if plcfg['ylabel'] :
        plt.ylabel(plcfg['ylabel'])
    if plcfg['title']:
        plt.title(plcfg['title'])

    datasets = np.shape(data)[0]

    for d in range(datasets):
        plt.plot(data[d][1], data[d][0], xycfg['symbols'][d],
                  color=plcfg['colors'][d])

    fig.savefig(outd+cfg['outfile'], dpi=plcfg['dpi'])
