# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 14:56:12 2017

@author: alex
"""

import os
import fnmatch
import numpy as np
from uncertainties import unumpy as unp
import S_Param_Script as perm
import permittivity_plot as pplot
import pickle

def avg_calc():
    """
    Calculates average S-parameters by loading several METAS produced \
        S-parameter .txt files and computing the average S-parameters from \
        all the files. Must give two directory paths, one for Port 1 \
        measurements and one for Port 2 measurements. Used to calculate \
        S-parameters for teflon washers.
    """
    # Directory paths
    path1 = '/Users/alex/Google Drive/research_data/Modern_dielectric_data/Alex/2017-07-20_washer_flush/LHKM_TRM_LRM_01_out/DUTs/p1/'
    path2 = '/Users/alex/Google Drive/research_data/Modern_dielectric_data/Alex/2017-07-20_washer_flush/LHKM_TRM_LRM_01_out/DUTs/p2/'
    # Airline used
    airline = 'washer' # Special case in get_METAS_data
    
    # Port 1 data
    classlist1 = []
    labels = []
    # Iterate through directory
    for itemName in os.listdir(path1):
        if fnmatch.fnmatch(itemName, '*.txt'): # Check for .txt
            #Loops over each itemName in the path. Joins the path and the 
            #   itemName and assigns the value to itemName. Uses the current 
            #   file name as a label.
            labels.append(os.path.splitext(itemName)[0])
            itemName = os.path.join(path1, itemName)
            print(itemName)
        else: # If not a .txt file, move on
            continue
        # Append the class instance to the list of data
        if os.path.isfile(itemName):
            classlist1.append(perm.AirlineData(*perm.get_METAS_data(airline,itemName),corr=False))
    
    # Port 2 data        
    classlist2 = []
    for itemName in os.listdir(path2):
        if fnmatch.fnmatch(itemName, '*.txt'):
            labels.append(os.path.splitext(itemName)[0])
            itemName = os.path.join(path2, itemName)
            print(itemName)
        else:
            continue
    
        if os.path.isfile(itemName):
            classlist2.append(perm.AirlineData(*perm.get_METAS_data(airline,itemName),corr=False))
    
    # Get sparam data from class instances
    s11 = []
    s22 = []
    s21 = []
    s12 = []
    freq = []
    # Port 1
    for dataItem in classlist1:
        s11.append(dataItem.s11)
        s22.append(dataItem.s22)
        s21.append(dataItem.s21)
        s12.append(dataItem.s12)
        freq.append(dataItem.freq)
    # Port 2    
    for dataItem in classlist2:
        s11.append(dataItem.s22) # Flip so they can be averaged
        s22.append(dataItem.s11)
        s21.append(dataItem.s12)
        s12.append(dataItem.s21)
        freq.append(dataItem.freq)
        
    # Build uarrays
    s11avg = unp.uarray(np.zeros([2,len(freq[0])]),np.zeros([2,len(freq[0])]))
    s22avg = unp.uarray(np.zeros([2,len(freq[0])]),np.zeros([2,len(freq[0])]))
    s21avg = unp.uarray(np.zeros([2,len(freq[0])]),np.zeros([2,len(freq[0])]))
    s12avg = unp.uarray(np.zeros([2,len(freq[0])]),np.zeros([2,len(freq[0])]))
    # Calculate average for each sparam at each frequency point
    for n in range(0,len(freq[0])):
        s11_mag_temp = []
        s11_phase_temp = []
        s22_mag_temp = []
        s22_phase_temp = []
        s21_mag_temp = []
        s21_phase_temp = []
        s12_mag_temp = []
        s12_phase_temp = []
        for m in range(0,len(s11)):
            s11_mag_temp.append(s11[m][0][n])
            s11_phase_temp.append(s11[m][1][n])
            s22_mag_temp.append(s22[m][0][n])
            s22_phase_temp.append(s22[m][1][n])
            s21_mag_temp.append(s21[m][0][n])
            s21_phase_temp.append(s21[m][1][n])
            s12_mag_temp.append(s12[m][0][n])
            s12_phase_temp.append(s12[m][1][n])
        s11avg[0][n] = sum(s11_mag_temp)/len(s11_mag_temp)
        s11avg[1][n] = sum(s11_phase_temp)/len(s11_mag_temp)
        s22avg[0][n] = sum(s22_mag_temp)/len(s11_mag_temp)
        s22avg[1][n] = sum(s22_phase_temp)/len(s11_mag_temp)
        s21avg[0][n] = sum(s21_mag_temp)/len(s11_mag_temp)
        s21avg[1][n] = sum(s21_phase_temp)/len(s11_mag_temp)
        s12avg[0][n] = sum(s12_mag_temp)/len(s11_mag_temp)
        s12avg[1][n] = sum(s12_phase_temp)/len(s11_mag_temp)

    # Pickle results for use in S_Param_Script
    with open('s11avg.p', 'wb') as f:
        pickle.dump(s11avg, f)
    with open('s22avg.p', 'wb') as f:
        pickle.dump(s22avg, f)
    with open('s21avg.p', 'wb') as f:
        pickle.dump(s21avg, f)
    with open('s12avg.p', 'wb') as f:
        pickle.dump(s12avg, f)
    with open('washer_freq.p', 'wb') as f:
        pickle.dump(freq, f)
        
    return s11avg, s22avg, s21avg, s12avg, labels, freq, s11, s22, s21, s12

#%%
def main():
    s11avg, s22avg, s21avg, s12avg, labels, freq, s11, s22, s21, s12 = avg_calc()
    pplot.make_sparam_plot(freq[0],s11,s22,s21,s12,labels)
    pplot.make_sparam_plot(freq[0],s11avg,s22avg,s21avg,s12avg)
if __name__ == '__main__':
    main()