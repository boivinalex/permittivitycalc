# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 14:24:19 2017

@author: alex
"""
# File input
import tkinter as tk
from tkinter.filedialog import askdirectory
import numpy as np
from uncertainties import unumpy as unp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from cycler import cycler
import os
import datetime

plt.style.use("ggplot")

#%% GLOBAL VARIABLES
DATE = str(datetime.date.today())

#%% FUNCTIONS
def _dirprompt():
    root = tk.Tk()
    root.withdraw()
    save_path = askdirectory(title='Select Directory to Save Figure')
    root.update()
    
    return save_path

def make_plot(xval, yval, plot_type='d', legend_label=None, name=None, \
              plot_title=None, ylabel=None, xlabel=None, spacing=None, \
              buffer=None, xlim=None, ylim=None, figure_size=(16,10), \
              publish=False, round_val=None):
    """
    Takes input from S_Param_Script_V5 and plots calculated permittivity. Can \
    handle multiple data sets. Plots uncertainty countour if plotting single \
    dataset with uncertainty present.
    
    Arguments
    ---------
    xval (array or list): x-axis data. Multiple data sets must be in a list.
    
    yval (array or list): y-axis data. Multiple data sets must be in a list.
    
    plot_type (str): Flag for default plot types. Can be set to 'd' for Real \
        Part, 'lf' for Imaginary Part, 'lt' for Loss Tangent, or 'c' for \
        Custom. If 'c' is used, plot_title, ylabel, xlabel, and round_val \
        must be set.
        
    legend_label (list of str): Plot legend label. Each dataset much have it's \
        own label. Stings must be in a list.
        
    name (str): For use with publish=True. Used in file name of saved figure. 
    
    plot_title (str): For use with plot_type='c'. Title of the plot.
    
    ylabel (str): For use with plot_type='c'. Y-axis label.
    
    xlabel (str): For use with plot_type='c'. X-axis label.
    
    spacing (float): For use with plot_type='c'. Sets the spacing between \
        y-axis tick marks. Currently not implemented.
        
    buffer (float): For use with plot_type='c'. Sets buffer space around the \
        min and max y-axis values. Currently not implemented.
        
    xlim (tuple, float): Manually set x-axis limits. Currently not implemented.
    
    ylim (tuples, float): Manually set y-axis limits. Currently not implemented.
    
    figure_size (tuple, int): Set the matplotlib figsize. Default: (12,9).
    
    publish (bool): If True save figure as .eps file. Default: False
    """
    
    # Default settings for plotting permittivity data   
    if plot_type == 'd': # Real part
        plot_title = 'Dielectric Constant'
        ylabel = '$\epsilon^\prime$'
        xlabel = 'Frequency (Hz)'
        rnd = 1 # decimals to round to for axes determination
    elif plot_type == 'lf': # Imaginary part
        plot_title = 'Loss Factor'
        ylabel = '$\epsilon^{\prime\prime}$'
        xlabel = 'Frequency (Hz)'
        rnd = 2
    elif plot_type == 'lt': # Loss tan
        plot_title = 'Loss Tangent'
        ylabel = '$tan\delta$'
        xlabel = 'Frequency (Hz)'
        rnd = 2
    elif plot_type == 'c': # Custom plot
        plot_title = plot_title
        ylabel = ylabel
        xlabel = xlabel
        rnd = round_val
    else:
        raise Exception('Invalid plot type')
    
    # Checks if input data is in a list and determines number of things to plot
    #   NOTE: Multiple data sets must be in a list      
    if isinstance(xval,list):
        number_to_compare = len(xval)
    else:
        number_to_compare = 1
        xval = [xval]
        yval = [yval]
    
    # If no label is specified, make one  
    if not legend_label:    # If legend_label is None
        legend_label = []
        # Label for first dataset is 'Data 1', next is 'Data 2', etc...
        for n in range(0,len(xval)):
            legend_label.append('Data {}'.format(n+1))
    else:   # If legend_label is a list, make sure no list items are None
        for n in range(0,len(xval)):
            # If a list item is None, make it a label
            if not legend_label[n]:
                legend_label[n] = 'Data {}'.format(n+1)
        
    # Remove 300kHz point
    x = []
    y = []
    for n in range(0,number_to_compare):
        x.append(xval[n][xval[n]>300000])
        y.append(yval[n][xval[n]>300000])
    
    # Determined axes limits    
    x_max = 0
    x_min = 9999999999
    y_max = 0
    y_min = 9999999999
    for n in range(0,number_to_compare):
        x_chk = unp.nominal_values(x[n])
        max_tmp = round(max(x_chk[~np.isnan(x_chk)]),rnd)
        min_tmp = round(min(x_chk[~np.isnan(x_chk)]))
        if max_tmp > x_max:
            x_max = max_tmp
        if min_tmp < x_min:
            x_min = min_tmp
        y_chk = unp.nominal_values(y[n])
        max_tmp = round(max(y_chk[~np.isnan(y_chk)]),rnd)
        min_tmp = round(min(y_chk[~np.isnan(y_chk)]),rnd)
        if max_tmp > y_max:
            y_max = max_tmp
        if min_tmp < y_min:
            y_min = min_tmp
            
    # Determine appropriate buffer and spacing depedning on plot type
    thickness = y_max - y_min
    if plot_type == 'd':
        if thickness < 0.1:
            buffer = 0.1
            spacing = 0.02
        else:
            buffer = 0.2
            spacing = round((thickness + 2*buffer)/9,1)
    elif plot_type in ('lf','lt','c'):
        if thickness < 0.01:
            buffer = 0.01
            spacing = 0.002
        else:
            buffer = 0.02
            spacing = round((thickness + 2*buffer)/9,2)

    # Makes sure the lowest point is 0 if y_min is 0
    if y_min == 0:
        y_min+=buffer
    elif y_min-buffer < 0:
        # Make sure buffer does not make ymin negative
        y_min = buffer
    
    # Plot
    f = plt.figure(figsize=figure_size)
    ax = f.add_subplot(111)
    ax.set_title(plot_title, fontsize=40)
    if number_to_compare > 8:
        ax.set_prop_cycle(cycler('color',\
                            sns.cubehelix_palette(number_to_compare)))
    else:
        ax.set_prop_cycle(cycler('color',\
                            sns.color_palette("Dark2",number_to_compare)))    
    ax.spines["top"].set_visible(False)  
    ax.spines["right"].set_visible(False)
    ax.get_xaxis().tick_bottom()  
    ax.get_yaxis().tick_left()
    ax.set_xscale('log')
    ax.set_ylim(y_min-buffer, y_max+buffer)
    ax.set_yticks(np.arange(y_min-buffer, y_max+buffer, spacing))
    ax.set_ylabel(ylabel, fontsize=40)
    ax.set_xlabel(xlabel, fontsize=40)
    ax.tick_params(axis='both', which='major', labelsize=30)
    if number_to_compare == 1: # plot uncertainty only if plotting one dataset
        ax.plot(unp.nominal_values(x[0]), unp.nominal_values(y[0]), lw=2, \
                label=legend_label[0])
        ax.fill_between(unp.nominal_values(x[0]), unp.nominal_values(y[0]) - \
                        unp.std_devs(y[0]), unp.nominal_values(y[0]) + \
                        unp.std_devs(y[0]), color="#3F5D7D",label='Uncertainty')
    else:
        for n in range(0,number_to_compare):
            ax.plot(unp.nominal_values(x[n]), unp.nominal_values(y[n]), lw=2, \
                    label=legend_label[n])
    ax.legend(fontsize=22,loc='best')
    if publish:
        # Check for directory
#        if not os.path.exists(DATAPATH):
#            os.makedirs(DATAPATH)
#        # Make file name    
        datapath = _dirprompt()
        savename = name.replace(' ','-') + '_' + plot_title.replace(' ','-') \
            + '_' + DATE + '.png'
        filepath = os.path.join(datapath,savename)
        # Save figure to .eps file
        plt.savefig(filepath,dpi=300,format='png',pad_inches=0)
    
def make_sparam_plot(freq,s11,s22,s21,s12,label=None):
    """
    Plot raw S-Parameter data from S_Param_Script_V5. Supports multiple \
    datasets for comparisson. Multilple datasets much be stored in a list and \
    muct have same frequency array and number of data points.
    
    Arguments
    ---------
    freq (array): Frequency points.
    
    s11,s22,s21,s12 (array or list of arrays): Mag and Phase S-Parameter data.
    
    label (list): (Optional) List of labels. Default: None
    """
    # Checks if input data is in a list and determines number of things to plot
    #   NOTE: Multiple data sets must be in a list      
    if isinstance(s11,list):
        number_to_compare = len(s11)
    else:
        number_to_compare = 1
        s11 = [s11]
        s22 = [s22]
        s21 = [s21]
        s12 = [s12]
    
    # Plot    
    f,ax = plt.subplots(4, 2, sharex=True, figsize=(18, 15))
    for n in range(0,number_to_compare):
        if label:
            kwargs = {"label":label[n]}
        else:
            kwargs = {}    
        ax[0,0].plot(freq,unp.nominal_values(s11[n][0]),**kwargs) #s11mag
        ax[0,0].set_title('Magnitude of S11')
        ax[0,1].plot(freq,unp.nominal_values(s11[n][1])) #s11phase
        ax[0,1].set_title('Phase of S11')
        ax[1,0].plot(freq,unp.nominal_values(s22[n][0])) #s22mag
        ax[1,0].set_title('Magnitude of S22')
        ax[1,1].plot(freq,unp.nominal_values(s22[n][1])) #s22phase
        ax[1,1].set_title('Phase of S22')
        ax[2,0].plot(freq,unp.nominal_values(s21[n][0])) #s21mag
        ax[2,0].set_title('Magnitude of S21')
        ax[2,1].plot(freq,unp.nominal_values(s21[n][1])) #s21phase
        ax[2,1].set_title('Phase of S21')
        ax[3,0].plot(freq,unp.nominal_values(s12[n][0])) #s12mag
        ax[3,0].set_title('Magnitude of S12')
        ax[3,1].plot(freq,unp.nominal_values(s12[n][1])) #s12phase
        ax[3,1].set_title('Phase of S12')
    # Hide redundant x-axis tick marks
    plt.setp([a.get_xticklabels() for a in ax[0, :]], visible=False)
    if label:
        ax[0,0].legend(loc=2)