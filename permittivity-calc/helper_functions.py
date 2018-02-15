#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Contains helper functions for file input and plotting multiple data sets"""

# File input
import tkinter as tk
from tkinter.filedialog import askopenfilename
import codecs
# Array math
import numpy as np
# Plotting
import permittivity_plot as pplot
# Make relative path
import os

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
    
def _get_file(airline,file_path):
    """Return the file path and airline name. Use prompts if needed."""
    L_in = None
    # Figure out the file path and the airline name
    if file_path and airline:
        if airline not in ('VAL','PAL','GAL','7','washer'):
            raise Exception('Wrong airline name. You used %s' % (airline))
        else:
            airline = airline
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
        # Prompt for file
        airline = airline
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
    elif airline == 'washer':
        L = 0.15
    else:
        L = float(L_in)
        
    return L, airline, dataArray, file

def multiple_meas(file_path=None,airline=None):
    """
    Generate an instance of AirlineData for every file in a directory. Store \
        the intances in a list, and plot them all using perm_compare.
        
    Arguments
    ---------
    file_path (str): Full path of any file in the source directory. \
        (Optional - will produce file dialog box if not provided.)
    
    airlne (str): Name of airline used. Every measurement must have been made \
        in the same airline. (Optional - will prompt if not provided.)
        
    Return
    ------
    class_list (lst): List of generated class instances of AirlineData.
    """
    # Use _get_file to get the filepath and airline name if not provided
    if file_path:   # If file path provided as argument
        file = file_path
    elif not file_path:   # If file path not provided
        print('\n')
        print("Select any data file in the source folder. All .txt "+\
              "files in the source folder must be METAS data tables.")
        # Get the file path and the airline name
        airline, file, L_in = _get_file(airline,file_path)
    elif not airline:   # If file path is given but airline is not
        airline, file, L_in = _get_file(airline,file_path)
        
    # Get directory path    
    directory = os.path.dirname(file)
    # Use a list to maintain order for plotting
    class_list = []
    # Iterate through all .txt files in the directory and run AirlineData
    for file in os.listdir(directory):
        if file.endswith(".txt"):
            filename = os.path.splitext(file)[0]    # Use file name as plot label
            # Append each new instance to class list
            class_list.append(AirlineData(*get_METAS_data(airline,\
                                os.path.join(directory, file)),name=filename))
    
    # Plot all files        
    perm_compare(class_list)
    
    return class_list       

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
        
def run_default(airline_name='VAL',**kwargs):
    """
    Run AirlineData on get_METAS_data with all the prompts and return the \
        instance.
    """
    return AirlineData(*get_METAS_data(airline=airline_name),**kwargs)

def run_example(flag='single'):
    test = AirlineData(*get_METAS_data(airline='GAL',file_path=DATAPATH + \
                        '2.5hrs.txt'),bulk_density=2.0,temperature=None,\
                         name='Alumina Vac 2.5hrs',date='2017/04/07')
    if flag == 'single':
        atm = AirlineData(*get_METAS_data(airline='VAL',\
            file_path=DATAPATH + 'atm.txt'),bulk_density=None,\
            temperature=None,name='Alumina atm',date=None,corr=True,\
            solid_dielec=None,solid_losstan=None,particle_diameter=None,\
            particle_density=None,nrw=False)
        return test, atm
    elif flag == 'multiple':
        test2 = AirlineData(*get_METAS_data(),name='TRM')
        classlist = [test,test2]
        perm_compare(classlist)
        return test, test2, classlist