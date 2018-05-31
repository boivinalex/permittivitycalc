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
        if airline not in ('VAL','PAL','GAL','7','washer','custom'):
            raise Exception('Wrong airline name. You used %s' % (airline))
        elif airline == 'custom':
            L_in = input('Enter the length of the airline (cm): ')
            file = file_path
        else:
            file = file_path
    elif file_path and not airline:
        # Prompt for airline
        airline = input('Are you using the "VAL", "PAL", "GAL", "7" mm, ' + \
                        'or "custom" airline?: ')
        if airline not in ('VAL','PAL','GAL','7','custom'):
            raise Exception('Wrong input')
        elif airline == 'custom':
            L_in = input('Enter the length of the airline (cm): ')
        file = file_path
    elif airline and not file_path:
        # Check airline
        if airline not in ('VAL','PAL','GAL','7','custom'):
            raise Exception('Wrong input')
        elif airline == 'custom':
            L_in = input('Enter the length of the airline (cm): ')
        # Prompt for file
        file = _prompt()
    else:   # Prompt for both
        # Airline
        airline = input('Are you using the "VAL", "PAL", "GAL", "7" mm, ' + \
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
    Return data arrays from METAS VNA Tools II text file based on user input \
    Frequencies must be in Hz. Recomeneded precision in VNA Tools II is f9.
    
    Arguments
    ---------
    airline (str): Airline used for measurement. Options are: 'VAL', \
        'PAL', 'GAL', '7'. If not provided, will prompt user. Prompt will \
        allow user to input a custom airline length.
    
    file_path (str): If a path is not given, will prompt user for file. \
        Default: None
    
    Return
    ------
    L (float): Length of airline in cm.
    
    SparmArray (array): S-parameters and their uncertainties (if present)
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
        L = 14.979
        #L = 14.9835    # 14mm airline - 15 cm 
    elif airline == 'GAL':
        L = 14.991
    else:
        L = float(L_in)
        
    return L, airline, dataArray, file

def perm_compare(classlist,allplots=False,**kwargs):
    """
    Given a list of AirlineData instances, plot their permittivity results \
        together using permittivity_plot_V1.py
        
    Arguments
    ---------
    classlist (list): List of instances of AirlineData
    
    allplots (bool): If True plot all of dielectric constant, loss factor and \
        loss tangent. If Flase plot only the dielectric constant and the loss \
        tangent. Default: False
    """
    freq = []
    dielec = []
    losstan = []
    labels = []
    for item in classlist:
        freq.append(item.freq)
        if item.normalize_density: # Check for normalize_density
            dielec.append(item.norm_dielec)
            losstan.append(item.norm_losstan)
        else:
            dielec.append(item.avg_dielec)
            losstan.append(item.avg_losstan)
        labels.append(item.name)
    kwargs["legend_label"] = labels
    if allplots:
        lossfac = []
        for item in classlist:
            lossfac.append(item.avg_lossfac)
        pplot.make_plot(freq,dielec,'d',**kwargs)
        pplot.make_plot(freq,lossfac,'lf',**kwargs)
        pplot.make_plot(freq,losstan,'lt',**kwargs)
    else:
        pplot.make_plot(freq,dielec,'d',**kwargs)
        pplot.make_plot(freq,losstan,'lt',**kwargs)
        
