# -*- coding: utf-8 -*-
"""
Fucntions to plot permittivity and permeability results and S-parameters.
"""
# File input
import tkinter as tk
from tkinter.filedialog import askdirectory
# Array math
import numpy as np
from uncertainties import unumpy as unp
# System
import os
import datetime
# Plotting
import matplotlib
try:
    from matplotlib import pyplot as plt
except:
    matplotlib.use('TkAgg',warn=False, force=True)
    import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, LogLocator, EngFormatter, NullFormatter, LogFormatter
import seaborn as sns
from cycler import cycler

plt.style.use("ggplot")

#%% GLOBAL VARIABLES
DATE = str(datetime.date.today())

#%% FUNCTIONS
def _dirprompt():
    root = tk.Tk()
    root.withdraw()
    global save_path_for_plots
    save_path_for_plots = askdirectory(title='Select Directory to Save Figure')
    root.update()
    
    return save_path_for_plots

def make_plot(xval, yval, plot_type='d', y_axis_type='normal', \
              legend_label=None, name=None, plot_title=None, ylabel=None, \
              xlabel=None, xticks=None, yticks=None, figure_size=(16,10), \
              freq_cutoff=None, publish=False):
    """
    Takes input from sparam_data and plots calculated permittivity. Can 
    handle multiple data sets. Plots uncertainty countour if plotting single 
    dataset with uncertainty present and plots error bars every 25 points when
    plotting multiple datasets.
    
    Arguments
    ---------
    xval : array or list 
        x-axis data. Multiple data sets must be in a list.
    
    yval : array or list 
        y-axis data. Multiple data sets must be in a list.
    
    plot_type : str  
        Flag for default plot types. Can be set to 'd' for the real 
        part of epsilon, 'lf' for the imaginary part of epsilon, 'lt' for 
        dielectric loss tangent, 'ur' for the real part of mu, 'ui' for the 
        imaginary part of mu, or 'c' for Custom. If 'c' is used, xticks and
        yticks must be provided.
        
    y_axis_type : str, optional
        Flag for type of axis to use for the y-axis. Can be either 'normal'
        (default) or 'log' for a log axis. If set to log, y tick postions must
        be manually provided to yticks.
        
    legend_label : list of str, optional
        Plot legend label. Each dataset much have it's 
        own label. Stings must be in a list.
        
    name : str, optional
        Required when publish=True. Used in file name of saved figure. 
    
    plot_title : str, optional
        For use when plot_type='c'. Title of the plot.
    
    ylabel : str, optional
        For use when plot_type='c'. Y-axis label.
    
    xlabel : str, optional
        For use when plot_type='c'. X-axis label.
    
    xticks : list, optional
        Manually set x-axis tick locations. Required when plot_type='c'. 
    
    yticks : list, optional
        Manually set y-axis tick locations. Required when plot_type='c' and 
        when y_axis_type='log'. 
    
    figure_size : tuple or int, optional
        Set the matplotlib figsize. Default: (16,10).
        
    freq_cutoff : float, optional
        Data points lower than freq_cutoff will not be plotted.
    
    publish : bool, optional
        If True save figure as .eps file. Default: False.
    """
    
    # Default settings for plotting permittivity data   
    if plot_type == 'd': # Real part
        plot_title = 'Real Part of the Permittivity'
        ylabel = r'$\epsilon^{\prime}_{r}$'
        xlabel = 'Frequency'
        rnd = 1 # decimals to round to for axes determination
    elif plot_type == 'lf': # Imaginary part
        plot_title = 'Imaginary Part of the Permittivity'
        ylabel = r'$\epsilon^{\prime\prime}_{r}$'
        xlabel = 'Frequency'
        rnd = 2
    elif plot_type == 'lt': # Loss tan
        plot_title = 'Loss Tangent'
        ylabel = r'$tan\delta$'
        xlabel = 'Frequency'
        rnd = 2
    elif plot_type == 'ur': # Real part of mu
        plot_title = 'Real Part of the Permeability'
        ylabel = r'$\mu^{\prime}_{r}$'
        xlabel = 'Frequency'
        rnd = 2 
    elif plot_type == 'ui': # Imaginary part of mu
        plot_title = 'Imaginary Part of the Permeability'
        ylabel = r'$\mu^{\prime\prime}_{r}$'
        xlabel = 'Frequency'
        rnd = 2 
    elif plot_type == 'c': # Custom plot
        pass
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
        
    # Remove all points lower than freq_cutoff
    x = []
    y = []
    for n in range(0,number_to_compare):
        if freq_cutoff:
            x.append(xval[n][xval[n]>freq_cutoff])
            y.append(yval[n][xval[n]>freq_cutoff])
        else:
            x.append(xval[n])
            y.append(yval[n])
    
    # Determine axes limits
    if plot_type!='c':  # skip for custom plots
        x_max = -np.inf
        x_min = np.inf
        y_max = -np.inf
        y_min = np.inf
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
    if not yticks: # skip if y ticks are manually provided
        thickness = y_max - y_min
        if plot_type in ('d','ur','ui'):
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
    # Plot labels
    ax.set_title(plot_title, fontsize=40)
    ax.set_ylabel(ylabel, fontsize=40)
    ax.set_xlabel(xlabel, fontsize=40)
    # Colours
    if number_to_compare > 6: # cycle through cubehelix palette if plotting more than 6 things
        ax.set_prop_cycle(cycler('color',\
                            sns.cubehelix_palette(number_to_compare)))
    else: # otherise use Dark2 palette (should be colour blind safe)
        ax.set_prop_cycle(cycler('color',\
                            sns.color_palette("Dark2",number_to_compare)))
    # Plot axes
    ax.set_xscale('log')
    # Set y ticks
    if y_axis_type == 'log': # check if y axis should be log
        ax.set_yscale('log')
        ax.set_ylim(min(yticks),max(yticks))
        majorLocator_y = FixedLocator(yticks)
        majorFormatter_y = LogFormatter()
        minorLocator_y = LogLocator(subs='all') # Use all interger multiples of the log base for minor ticks 
        minorFormatter_y =  NullFormatter() # No minor tick labels
        # Apply y ticks
        ax.get_yaxis().set_major_locator(majorLocator_y)
        ax.get_yaxis().set_major_formatter(majorFormatter_y)
        ax.get_yaxis().set_minor_locator(minorLocator_y)
        ax.get_yaxis().set_minor_formatter(minorFormatter_y)
    elif yticks and y_axis_type == 'normal':
        ax.set_ylim(min(yticks),max(yticks))
        ax.set_yticks(yticks)
    elif not yticks: # auto set
        ax.set_ylim(y_min-buffer, y_max+buffer)
        ax.set_yticks(np.arange(y_min-buffer, y_max+buffer, spacing))
    # Set x ticks
    if xticks:
        x_ticklocs = xticks
    elif not xticks: # auto set
        if x_min == 0: # log of min and max x values
            x_logmin = 0 # don't take log of 0
        else:
            x_logmin = np.log10(x_min) 
        if x_max == 0:
            x_logmax = 0
        else:
            x_logmax = np.log10(x_max)
        x_logticks = np.logspace(x_logmin, x_logmax, num=4) # 4 equaly spaced points in log space
        x_ticklocs = []
        for n in range(len(x_logticks)): # round scientific values and make a list
            x_ticklocs.append(np.float(np.format_float_scientific(x_logticks[n],\
                            precision=0)))
        if len(set(x_ticklocs)) < 4: # check that this produced 4 unique values
            x_ticklocs = [] # if not do it again with precision = 1
            for n in range(len(x_logticks)): 
                x_ticklocs.append(np.float(np.format_float_scientific(x_logticks[n]\
                            ,precision=1)))
    majorLocator_x = FixedLocator(x_ticklocs)
    majorFormatter_x = EngFormatter(unit='Hz') # Format major ticks with units
    minorLocator_x = LogLocator(subs='all') # Use all interger multiples of the log base for minor ticks 
    minorFormatter_x =  NullFormatter() # No minor tick labels 
    # Apply x ticks
    ax.get_xaxis().set_major_locator(majorLocator_x)
    ax.get_xaxis().set_major_formatter(majorFormatter_x)
    ax.get_xaxis().set_minor_locator(minorLocator_x)
    ax.get_xaxis().set_minor_formatter(minorFormatter_x)
    # Format the actual tick marks
    ax.tick_params(which='both', width=1, labelsize=30)
    ax.tick_params(which='major', length=7)
    ax.tick_params(which='minor', length=4)
    # Use smaller line width for minor tick grid lines
    ax.grid(b=True, which='major', color='w', linewidth=1.0)
    ax.grid(b=True, which='minor', color='w', linewidth=0.5)
    # Do the actual plotting
    if number_to_compare == 1: # plot uncertainty only if plotting one dataset
        ax.plot(unp.nominal_values(x[0]), unp.nominal_values(y[0]), lw=2, \
                label=legend_label[0])
        ax.fill_between(unp.nominal_values(x[0]), unp.nominal_values(y[0]) - \
                        unp.std_devs(y[0]), unp.nominal_values(y[0]) + \
                        unp.std_devs(y[0]), color="#3F5D7D",label='Uncertainty')
    else:
        for n in range(0,number_to_compare):
            ax.errorbar(unp.nominal_values(x[n]), unp.nominal_values(y[n]), \
                        yerr=unp.std_devs(y[n]), errorevery=25, elinewidth=1, \
                        capthick=1, capsize=2,lw=2,label=legend_label[n])
    ax.legend(fontsize=30,loc='best')
    if publish:
        # Make file name
        #If save_path_for_plots already exits, use it, otherwise promt for path
        if 'save_path_for_plots' in globals():
            datapath = save_path_for_plots
        else:
            datapath = _dirprompt()     # prompt for save dir
        savename = name.replace(' ','-') + '_' + plot_title.replace(' ','-') \
            + '_' + DATE + '.eps'
        filepath = os.path.join(datapath,savename)
        # Save figure to .eps file
        plt.savefig(filepath,dpi=300,format='eps',pad_inches=0)
    plt.show()
    
def make_sparam_plot(freq,s11,s22,s21,s12,label=None,shorted=False,s11_short=None):
    """
    Plot raw S-Parameter data from S_Param_Script_V5. Supports multiple \
    datasets for comparisson. Multilple datasets much be stored in a list and \
    muct have same frequency array and number of data points.
    
    Arguments
    ---------
    freq : array 
        Frequency points.
    
    s11,s22,s21,s12 : array or list of arrays 
        Mag and Phase S-Parameter data.
    
    label : list, optional 
        List of labels. Default: None
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
        if shorted:
            s11_short = [s11_short]
    
    # Plot    
    f,ax = plt.subplots(4, 2, sharex=True, figsize=(18, 15))
    for n in range(0,number_to_compare):
        if label:
            kwargs = {"label":label[n]}
        else:
            kwargs = {}    
        ax[0,0].plot(freq,unp.nominal_values(s11[n][0]),**kwargs) #s11mag
        if shorted:
            ax[0,0].plot(freq,unp.nominal_values(s11_short[n][0]))
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