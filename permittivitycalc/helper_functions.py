#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Contains helper functions for file input"""

# File input
import tkinter as tk
from tkinter.filedialog import askopenfilename
import codecs
# Array math
import numpy as np
# Plotting
import permittivitycalc.permittivity_plot as pplot

def _prompt():
    """Prompt for VNA Tools II text file"""
    print('\n')
    print('Select the VNA Tools II Output Data Table')
    print('\n')
    root = tk.Tk()
    root.withdraw()
    file = askopenfilename(filetypes=[('text files', '*.txt')],\
                            title='Select the VNA Tools II Output Data Table')
    root.update()
    
    return file
    
def _get_file(airline=None,file_path=None):
    """Return the file path and airline name. Use prompts if needed."""
    L_in = None
    # Figure out the file path and the airline name
    if file_path and airline:
        if airline not in ('VAL','PAL','GAL','7','washer','custom','10cm'):
            raise Exception('Wrong airline name. You used %s' % (airline))
        elif airline == 'custom':
            L_in = input('Enter the length of the airline (cm): ')
            file = file_path
        else:
            file = file_path
    elif file_path and not airline:
        # Prompt for airline
        airline = input('Are you using the "VAL", "PAL", "GAL", "7" mm, "10cm"' + \
                        'or "custom" airline?: ')
        if airline not in ('VAL','PAL','GAL','7','custom','10cm'):
            raise Exception('Wrong input')
        elif airline == 'custom':
            L_in = input('Enter the length of the airline (cm): ')
        file = file_path
    elif airline and not file_path:
        # Check airline
        if airline not in ('VAL','PAL','GAL','7','custom','10cm'):
            raise Exception('Wrong input')
        elif airline == 'custom':
            L_in = input('Enter the length of the airline (cm): ')
        # Prompt for file
        file = _prompt()
    else:   # Prompt for both
        # Airline
        airline = input('Are you using the "VAL", "PAL", "GAL", "7" mm, "10cm"' + \
                        'or "custom" airline?: ')
        if airline not in ('VAL','PAL','GAL','7','custom'):
            raise Exception('Wrong input')
        elif airline == 'custom':
            L_in = input('Enter the length of the airline (cm): ')
        # File
        file = _prompt()
        
    return airline, file, L_in
    
def get_METAS_data(airline=None,file_path=None):
    """
    Return data arrays from METAS VNA Tools II text file based on user input 
    Frequencies must be in Hz. Recomeneded precision in VNA Tools II is f9.
    
    Arguments
    ---------
    airline : str 
        Airline used for measurement. Options are: 'VAL', 
        'PAL', 'GAL', '7'. If not provided, will prompt user. Prompt will 
        allow user to input a custom airline length.
    
    file_path : str 
        If a path is not given, will prompt user for file. 
        Default: None
    
    Return
    ------
    L : float 
        Length of airline in cm.
    
    SparmArray : array 
        S-parameters and their uncertainties (if present).
    """
    
    # Get the file path and the airline name
    airline, file, L_in = _get_file(airline,file_path)

    # Open the file and make array    
    open_file = codecs.open(file, encoding="utf-8")
    dataArray = np.loadtxt(open_file,delimiter="\t",skiprows=1)
    open_file.close()
    
    # Coaxial Airline Length
    if airline == '7':
        L = 9.9873     # 7 mm airline - 10 cm
    elif airline == 'VAL' or airline == 'PAL':
        L = 14.989
        #L = 14.9835    # 14mm airline - 15 cm 
    elif airline == 'GAL':
        L = 14.991
    elif airline == '10cm':
        L = 9.993
    else:
        L = float(L_in)
        
    return L, airline, dataArray, file

def perm_compare(classlist,allplots=False,**kwargs):
    """
    Given a list of AirlineData instances, plot their permittivity results 
    together using permittivity_plot. Fot each 
    item in the list, will plot the first available data type out of: 
    normalized data, corrected (de-embeded) data, uncorrected data.
    
    Will use the smallest freq_cutoff value in the list (if avaialbe) for 
    plotting.
        
    Arguments
    ---------
    classlist : list of sparam_data.AirlineData
        List of instances of AirlineData.
    
    allplots : bool 
        Default: False. If True plot all of real and imaginary parts of the 
        permittivity and the loss tangent. If Flase plot only the real part
        and the loss tangent.
    """
    # Check that classlist is a list
    if isinstance(classlist,list):
        pass
    else:
        raise Exception('AirlineData instances must be provided in a list')
    # Create data lists
    freq = []
    dielec = []
    losstan = []
    labels = []
    cutoffs = []
    for item in classlist:
        freq.append(item.freq)
        labels.append(item.name)
        cutoffs.append(item.freq_cutoff)
        if item.normalize_density: # Check for normalize_density
            dielec.append(item.norm_dielec)
            losstan.append(item.norm_losstan)
        elif item.corr:
            dielec.append(item.corr_avg_dielec)
            losstan.append(item.corr_avg_losstan)
        else:
            dielec.append(item.avg_dielec)
            losstan.append(item.avg_losstan)
    # Also get lossfac if allplots
    if allplots:
        lossfac = []
        for item in classlist:
            if item.normalize_density:
                lossfac.append(item.norm_lossfac)
            elif item.corr:
                lossfac.append(item.corr_avg_lossfac)
            else:
                lossfac.append(item.avg_lossfac)
    # Pass arguments to make_plot
    kwargs["legend_label"] = labels
    if None in cutoffs:     #pass None as freq_cutoff if any of them are None
        kwargs['freq_cutoff'] = None
    else:   #else pass the smallest one
        kwargs['freq_cutoff'] = np.min(cutoffs)
    # Make the plots
    if allplots:
        pplot.make_plot(freq,dielec,'d',**kwargs)
        pplot.make_plot(freq,lossfac,'lf',**kwargs)
        pplot.make_plot(freq,losstan,'lt',**kwargs)
    else:
        pplot.make_plot(freq,dielec,'d',**kwargs)
        pplot.make_plot(freq,losstan,'lt',**kwargs)
        
