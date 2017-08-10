# -*- coding: utf-8 -*-
"""
Complex Permittivity Determination using New Non-Iterative Method from \
S-Parameters

## To run wtih default settings use "var_name = run_default()" ##

Alexandre Boivin
University of Toronto Department of Earth Sciences, Solar System Exploration \
Group
PhD Research: Complex Permittivity Measurements of Planetary Regolith Analogs

Based on MATLAB script by Dylan Hickson and Alexandre Boivin

Dylan Hickson
York University Lassonde School of Engineering
MSc/PhD Research: Complex Permittivity Measurements of Granular Materials

Created on Wed Jul  5 12:16:27 2017

@author: alex

Version History (all previous versions based in MATLAB):
    Tue Jul 25 2017 - All version tracking now done with Git
    Wed Jul 19 2017 - Python script Ver. 2.0, Permittivity Calc Ver. 5.0
        - Scirpt re-wrote to enclose data in a Class including both raw and \
            corrected S-parameter data and calculated permittivity data.
        - Now uses new permmittivity plotting script permittivity_plot_V1.py \
            to produce all plots including raw S-parameter plots.
        - Supports plotting permittivity data from multiple Class instances.
    Thu Feb 16 2017 - Python script Ver. 1.2, Permittivity Calc Ver. 4.2
        - Fixed issue where NaN values in some files were causing errors
    Tue Jun 20 2017 - Permittivity Calc Ver. 4.1.2
        - Added new airline nomenclature and new length measurements
    Sun Nov 13 2016 - Python script Ver. 1.1.1, Permittivity Calc Ver. 4.1.1
        - Cut out data at 300kHz due to consistantly bad data at the first \
            data point
        - Fixed mistake in harmonic_resonances where n was going from 0 to 1 \
            instead of from 1 to n cauing a divide by 0 error
        - Changed x and y ticks on plots to look better and be more flexible
    Wed Oct 05 2016 - Python script Version 1.1, Permitivity Calc Version 4.1
        - Removed MATLAB-like clear due to errors in Spyder 3.0
        - Plot resonances based on max frequency
        - Disabled ability to choose both S-parameters in get_METAS_data if \
            argument prompt is set to False for use with permittivity_compare\
            script
    Wed Sep 21 2016 - Python script Version 1.0, Permitivity Calc Version 4
"""
# File input
import tkinter as tk
from tkinter import filedialog
import codecs
import pickle
# Array math
import numpy as np
import uncertainties 
#Citation: Uncertainties: a Python package for calculations with uncertainties,
#    Eric O. LEBIGOT, http://pythonhosted.org/uncertainties/
from uncertainties import unumpy as unp
# Plotting
import permittivity_plot as pplot
# Make relative path
import os

#%% GLOBAL VARIABLES
E_0 = 8.854187817620*10**-12 #Permittivity of vacuum (F/m) 
U_0 = 4*np.pi*10**-7 #Permeability of vacuum (V*s/A*m) 
C = (1)*1/np.sqrt(E_0*U_0) #Speed of light (m/s) 
E_R_AIR = 1.00058986 #Dielectric Constant of Air (STP @ 0.9 MHz) 
C_R_AIR = C / E_R_AIR #Speed of light in air
LAM_C = float('inf') #Cut-off wavelength = infinity
DATAPATH = os.path.dirname(__file__) + '/data/'

#%% CLASSES
class AirlineData:
    """
    S-parameter data from METAS text file output
    
    Attributes:
    ----------
    L (float): Length of airline in cm.
    
    airline_name (str): Name of airline used for measurement.
    
    file (str): Input file path.
    
    freq (array): Frequency points.
    
    s11, s21, s12, s22 (array): Mag and Phase S-Parameters.
    
    *_dielec (array): Real part of the permittivity. Can be avg_dielec, \
        forward_dielec, or reverse_dielec for average, forward, or \
        reverse permittivity.
    
    *_lossfac (array): Imaginary part of the permittivity. Same as above.
    
    *_losstan (array): Loss tangent. Same as above.
    
    corr_* (array): De-embeded version of S-parameters or permittivity data. \
        Only average S-parameters are used for permittivity calculations with \
        corrected S-parameters. Examples: corr_s11, corr_avg_losstan. Only \
        created if corr = True.
    
    bulk_density (float): (Optional) Bulk density of material. Nessesary for \
        bulk density normalization.
        
    temperature (str or float): (Optional) Temperature of measurement.
    
    name (str): (Optional) Name of measurement to be used in plots.
    
    date (str): (Optional) Measurement date.
    
    corr (bool): (Optional) If True, also correct S-parameter data and \
        produce corr_* arrays. Default = True.
    """
    def __init__(self,L,airline,dataArray,file,corr=True,bulk_density=None,\
                 temperature=None,name=None,date=None):
        self.L = L
        self.airline_name = airline
        self.file = file
        self.corr = corr
        # Unpack data into arrays
        self.freq, self.s11, self.s21, self.s12, self.s22 = \
            self._unpack(dataArray)
        # Calculate permittivity
        self.avg_dielec, self.avg_lossfac, self.avg_losstan = \
            self._permittivity_calc('a')
        self.forward_dielec, self.forward_lossfac, self.forward_losstan = \
            self._permittivity_calc('f')
        self.reverse_dielec, self.reverse_lossfac, self.reverse_losstan = \
            self._permittivity_calc('r')
        # Calculate corrected permittivity if array length is 601 only
        # Also check for NaNs and don't run if any
        if corr and len(self.freq) == 601:
            if not np.isnan(unp.nominal_values(self.avg_dielec)).any():
                self.corr_s11, self.corr_s21, self.corr_s12, self.corr_s22 = \
                    self._de_embed()
                self.corr_avg_dielec, self.corr_avg_lossfac, \
                    self.corr_avg_losstan = self._permittivity_calc('a',True)
        # Optional attributes
        self.bulk_density = bulk_density
        self.temperature = temperature
        self.name = name
        self.date = date
        
    def __repr__(self):
        rep = 'AirlineData(*get_METAS_data(airline=%r,file_path=%r),' % \
                (self.airline_name,self.file) + \
                'bulk_density=%r,temperature=%r,name=%r,date=%r,corr=%r)' % \
                (self.bulk_density,self.temperature,self.name,self.date,\
                 self.corr)
        return rep
        
    def __str__(self):
        srep = 'measured in ' + self.airline_name 
        if self.name:
            srep = self.name + ' ' + srep
        if self.date:
            srep += ' on ' + str(self.date)
        if self.temperature:
            srep += ' at ' + str(self.temperature) + ' degrees'
        if self.bulk_density:
            srep += ' with a bulk density of ' + str(self.bulk_density) + ' g/cm^3'
        srep += ' from file: \n' + self.file
        return srep
        
    def _unpack(self,dataArray):
        """See if uncertainty in data and unpack to S-parameter arrays"""
        if dataArray.shape[1] == 17: # Has unc so use unumpy
            freq = dataArray[:,0]
            s11 = unp.uarray([dataArray[:,1],dataArray[:,3]],\
                                     [dataArray[:,2],dataArray[:,4]])
            s21 = unp.uarray([dataArray[:,5],dataArray[:,7]],\
                                     [dataArray[:,6],dataArray[:,8]])
            s12 = unp.uarray([dataArray[:,9],dataArray[:,11]],\
                                      [dataArray[:,10],dataArray[:,12]])
            s22 = unp.uarray([dataArray[:,13],dataArray[:,15]],\
                                     [dataArray[:,14],dataArray[:,16]])
        elif dataArray.shape[1] == 9: # No unc
            freq = dataArray[:,0]
            s11 = np.array([dataArray[:,1],dataArray[:,2]])
            s21 = np.array([dataArray[:,3],dataArray[:,4]])
            s12 = np.array([dataArray[:,5],dataArray[:,6]])
            s22 = np.array([dataArray[:,7],dataArray[:,8]])
        else:
            raise Exception('Input file has the wrong number of columns')
        return freq, s11, s21, s12, s22
    
    def _permittivity_calc(self,s_param,corr=False):
        """
        Return the complex permittivity and loss tangent from S-parameters \
            and their uncertainties (if present).
    
        Uses the New Non-iterative method described in Boughriet, A.H., Legrand, \
            C. & Chapoton, A., 1997. Noniterative stable transmission/reflection \
            method for low-loss material complex permittivity determination. \
            IEEE Transactions on Microwave Theory and Techniques, 45(1), pp.52–57.
        
        Assumes the magnetic permeability of the sample = 1 (i.e non-magnetic).
        
        Arguments
        ---------
        s_param (str): Calculates complex permittivity using either \
            the average ('a'), the forward ('f'), or the reverse ('r') \
            S-Parameters.
            
        Return
        ------
        f_0 (array): Array of measured frequency values.
        
        dielec (array): Real part of the complex permittivity \
            (dielectric constant).
        
        lossfac (array): Imaginary part of the complex permittivity \
            (loss factor).
        
        losstan (array): Loss tangent. Defined as the ratio of the imaginary \
            and the real part of the complex permittivity (lossfac/dielec).
        """
        
        # Check if using corrected S-params
        if corr:
            s11 = self.corr_s11
            s21 = self.corr_s21
            s12 = self.corr_s12
            s22 = self.corr_s22
        else:
            s11 = self.s11
            s21 = self.s21
            s12 = self.s12
            s22 = self.s22
        
        # Strip arrays of their uncertainties
        #   This measure is temporary as the uncertainties module does not 
        #   support complex numbers at the time of writing.
        #   Later, uncertainties may be propagated through and compared with
        #   those calculated in Boughriet et al.
        if s_param == 'a':
            s_r = (unp.nominal_values(s11) + \
                   unp.nominal_values(s22))/2
            s_t = (unp.nominal_values(s21) + \
                   unp.nominal_values(s12))/2
        elif s_param == 'f':
            s_r = unp.nominal_values(s11)
            s_t = unp.nominal_values(s21)
        else:
            s_r = unp.nominal_values(s22)
            s_t = unp.nominal_values(s12)
            
        # Initialize arrays
        Nrows = len(self.freq)
        gam = np.zeros(Nrows,dtype=complex)
        sign = np.empty(Nrows,dtype=bool)
        a = np.zeros(Nrows,dtype=complex)
        
        ### NOTE: change function from np to unp to propagate uncertainties ###
        
        # Calculate complex permittivity
        lam_0 = (C/self.freq)*100 #Free-space wavelength (cm)
            
        # Convert phase to radians
        s_r[1] = np.radians(s_r[1])
        s_t[1] = np.radians(s_t[1])
        
        # Convert to rectangular notation (cast to complex)
        s_reflect = 1j*(s_r[0]*np.sin(s_r[1])); \
        s_reflect += s_r[0]*np.cos(s_r[1])
        s_trans = 1j*s_t[0]*np.sin(s_t[1]); \
        s_trans += s_t[0]*np.cos(s_t[1])
        
        # Calculate X (root of reflection coefficient)  
        x = ((s_reflect**2) - (s_trans**2) + 1)/(2*s_reflect)
        
        # Calculate Gamma (reflection coefficient) 
        gam1 = x + np.sqrt(x**2 - 1)
        gam2 = x - np.sqrt(x**2 - 1)
        
        # Determine correct root with condition |GAMMA| < 1 
        sign[np.absolute(gam1<=1)] = True
        sign[np.absolute(gam1>1)] = False
        
        # Check if sign is true or false and assign coresponding gamma value to gam
        gam[sign] = gam1[sign]      #if sign is true set gam to gam1
        gam[np.invert(sign)] = gam2[np.invert(sign)]      #if false set to gam2
        
        # Calculate T (transmission coefficient)
        t = (s_trans + s_reflect - gam) / (1 - ((s_reflect + s_trans) * gam))
        # Unwrap phase ambiguities in phase angle of T
        t_phaser = np.angle(t) #get phase angle from complex
        t_phase_unwrap = np.copy(t_phaser)
        # Ignore NaN values
        t_phase_unwrap[~np.isnan(t_phase_unwrap)] = \
                       np.unwrap(t_phase_unwrap[~np.isnan(t_phase_unwrap)])
                       
        # Calculate ln(1/T) with correct phases
        ln_1_T = 1j*t_phase_unwrap; ln_1_T += np.log(1/abs(t))
        # Also create new unwrapped T vector
        new_t = 1j*t_phase_unwrap; new_t += np.absolute(t)
        
        # Calculate A 
        a_1 = np.sqrt(-(ln_1_T / (2*np.pi*self.L))**2)
        a_2 = -1 * (np.sqrt(-(ln_1_T / (2*np.pi*self.L))**2))
         
        # Determine correct root with condition Re(1/A) > 0
        a[a_1.real > 0] = 1 / a_1[a_1.real > 0]
        a[a_1.real < 0] = 1 / a_2[a_1.real < 0]
          
        # Find wavelength in empty cell
        lam_og = 1 / np.sqrt(1/lam_0**2 - 1/LAM_C**2)
        
        # Calculated effective electromagnetic parameters
        mu_eff = 1      #Assume mu_eff = mu_r = 1
        ep_eff = (lam_og/a)**2
         
        # Calculate e_r (relative permittivity)
        ep_r = (1 - lam_0**2/LAM_C**2)*ep_eff + \
        (lam_0**2/LAM_C*2)/mu_eff
        
        dielec = ep_r.real
        lossfac = ep_r.imag
        losstan = lossfac/dielec
        
        # Calculate uncertainties
        if isinstance(self.s11[0][0], uncertainties.UFloat): # Check for uncertainties
            delta_dielec, delta_lossfac, delta_losstan = \
                self._calc_uncertainties(s_param,Nrows,sign,x,s_reflect,\
                                         s_trans,gam,lam_og,new_t,mu_eff,\
                                         ep_eff,lam_0,dielec,lossfac,losstan,corr)
            dielec = unp.uarray(dielec,delta_dielec)
            lossfac = unp.uarray(lossfac,delta_lossfac)
            losstan = unp.uarray(losstan,delta_losstan)
        
        return dielec,lossfac,losstan
        
    def _calc_uncertainties(self,s_param,Nrows,sign,x,s_reflect,s_trans,gam,\
                            lam_og,new_t,mu_eff,ep_eff,lam_0,dielec,\
                            lossfac,losstan,corr):
        """
        Calculates permittivity uncertainties from Boughriet et al.
        
        Returns
        -------
        delta_dielec (array): Uncertainty in real part.
        
        delta_lossfac (array): Uncertainty in imaginary part.
        
        delta_losstan (array): Uncertainty in loss tangent.
        """
        # Check if using corrected S-params
        if corr:
            s11 = self.corr_s11
            s21 = self.corr_s21
            s12 = self.corr_s12
            s22 = self.corr_s22
            L = self.L - 0.3
        else:
            s11 = self.s11
            s21 = self.s21
            s12 = self.s12
            s22 = self.s22
            L = self.L
        
        # Strip uarrays of their nominal values
        if s_param == 'a':
            err_s_r = (unp.std_devs(s11) + \
                   unp.std_devs(s22))/2
            err_s_t = (unp.std_devs(s21) + \
                   unp.std_devs(s12))/2
        elif s_param == 'f':
            err_s_r = unp.std_devs(s11)
            err_s_t = unp.std_devs(s21)
        else:
            err_s_r = unp.std_devs(s22)
            err_s_t = unp.std_devs(s12)
        
        # Initialize arrays
        delta_length = 0.001 #in cm
        dgam_dS_reflect = np.zeros(Nrows,dtype=complex)
        dgam_dS_trans = np.zeros(Nrows,dtype=complex)
        
        # Calculate partial derivatives
        # Determine +/- ambiguity by condition |GAMMA| < 1 from before            
        dgam_dS_reflect[sign] = (1 + (x[sign]/np.sqrt((x[sign]**2)-1)))*\
                        (((2*s_reflect[sign]**2)-(2*s_trans[sign]**2)+1)\
                         /2*s_reflect[sign]**2)
        dgam_dS_trans[sign] = (1 + (x[sign]/np.sqrt((x[sign]**2)-1)))*\
                        (-1*(s_trans[sign]/s_reflect[sign]))
        
        dgam_dS_reflect[np.invert(sign)] = (1 + \
                        (x[np.invert(sign)]/np.sqrt((x[np.invert(sign)]**2)-1)))\
                        *(((2*s_reflect[np.invert(sign)]**2)\
                           -(2*s_trans[np.invert(sign)]**2)+1)\
                        /2*s_reflect[np.invert(sign)]**2)
        dgam_dS_trans[np.invert(sign)] = (1 + (x[np.invert(sign)]\
                        /np.sqrt((x[np.invert(sign)]**2)-1)))\
                        *(-1*(s_trans[np.invert(sign)]/s_reflect[np.invert(sign)]))
     
        dT_dS_reflect = (1 - gam**2 + dgam_dS_reflect*(((s_reflect + s_trans)**2) \
                        - 1))/((1 - (s_reflect + s_trans)*gam)**2)
        dT_dS_trans = (1 - gam**2 + dgam_dS_trans*(((s_reflect + s_trans)**2) \
                        - 1))/((1 - (s_reflect + s_trans)*gam)**2)
        #deff_dgam = 0  # 0 in Boughriet, 1997
        deff_dT = 2*((1j*lam_og/(2*np.pi*L))*np.log(new_t))\
                        *((1j*lam_og)/(2*np.pi*L*new_t))
        
        deff_dS_reflect_mag = (deff_dT*dT_dS_reflect)*np.exp(1j*np.imag(s_reflect))
        deff_dS_trans_mag = (deff_dT*dT_dS_trans)*np.exp(1j*np.imag(s_trans))
        
        deff_dS_reflect_phase = 1j*np.real(s_reflect)*deff_dS_reflect_mag
        deff_dS_trans_phase = 1j*np.real(s_trans)*deff_dS_trans_mag
        
        dT_dd = (-2*np.pi*1j*np.sqrt(mu_eff*ep_eff)*\
                 np.exp((-2*np.pi*L*1j*np.sqrt(mu_eff*ep_eff))/lam_og))/lam_og
        
        # From Baker-Jarvis, 1990 can get deff_dL
        deff_dL = deff_dT*dT_dd
        
        # Combine partial derivatives for overal measurement uncertainty
        # In Boughriet, 1997 the last term (length uncertainty) isn't squared. 
        #   Mistake???? Answer: Yes Reason: Sqaured in older references
            
        delta_dielec = abs((1 - ((lam_0**2)/(LAM_C**2)))\
            *np.sqrt((((np.real(deff_dS_reflect_mag))*err_s_r[0])**2) + \
            (((np.real(deff_dS_reflect_phase))*err_s_r[1])**2) + \
            (((np.real(deff_dS_reflect_mag))*err_s_t[0])**2) + \
            (((np.real(deff_dS_trans_phase))*err_s_t[1])**2) + \
            (((np.real(deff_dL))*delta_length)**2)))
        
        delta_lossfac = abs((1 - ((lam_0**2)/(LAM_C**2)))\
            *np.sqrt((((np.imag(deff_dS_reflect_mag))*err_s_r[0])**2) + \
            (((np.imag(deff_dS_reflect_phase))*err_s_r[1])**2) + \
            (((np.imag(deff_dS_trans_mag))*err_s_t[0])**2) + \
            (((np.imag(deff_dS_trans_phase))*err_s_t[1])**2) + \
            (((np.imag(deff_dL))*delta_length)**2)))
        
        delta_losstan = abs(losstan*np.sqrt((delta_dielec/dielec)**2 +\
                                        (delta_lossfac/lossfac)**2))
        
        return delta_dielec, delta_lossfac, delta_losstan
    
    def _de_embed(self):
        """
        Perform S-parameter de-embedding to remove influence of washers on \
            measurement.
        """
        # Get washer S-parameters and split into mag and phase
        # Unpickle files produced by avg_sparam.py
        with open(DATAPATH + 's11avg.p', 'rb') as f:
            sw11 = pickle.load(f)
        with open(DATAPATH + 's22avg.p', 'rb') as f:
            sw22 = pickle.load(f)
        with open(DATAPATH + 's21avg.p', 'rb') as f:
            sw21 = pickle.load(f)
        with open(DATAPATH + 's12avg.p', 'rb') as f:
            sw12 = pickle.load(f)
        with open(DATAPATH + 'washer_freq.p', 'rb') as f:
            wfreq = pickle.load(f)
        # Split washer sparams into mag and phase since uncertainties package 
        #   does not support complex numbers
        # TESTING OUT FLUSH WASHERS
        sw22 = sw11
        sw12 = sw21
        sw11_mag = sw11[0]
        sw22_mag = sw22[0]
        sw21_mag = sw21[0]
        sw12_mag = sw12[0]
        sw11_phase = sw11[1]
        sw22_phase = sw22[1]
        sw21_phase = sw21[1]
        sw12_phase = sw12[1]
        # Cast to complex and unwarp phase
        sw11_complex = 1j*(unp.nominal_values(sw11_mag)*\
                           np.sin(np.unwrap(np.radians(unp.nominal_values(sw11_phase))))); \
        sw11_complex += unp.nominal_values(sw11_mag)*\
                            np.cos(np.unwrap(np.radians(unp.nominal_values(sw11_phase))))
        sw22_complex = 1j*(unp.nominal_values(sw22_mag)*\
                           np.sin(np.unwrap(np.radians(unp.nominal_values(sw22_phase))))); \
        sw22_complex += unp.nominal_values(sw22_mag)*\
                            np.cos(np.unwrap(np.radians(unp.nominal_values(sw22_phase))))
        sw21_complex = 1j*(unp.nominal_values(sw21_mag)*\
                           np.sin(np.unwrap(np.radians(unp.nominal_values(sw21_phase))))); \
        sw21_complex += unp.nominal_values(sw21_mag)*\
                            np.cos(np.unwrap(np.radians(unp.nominal_values(sw21_phase))))
        sw12_complex = 1j*(unp.nominal_values(sw12_mag)*\
                           np.sin(np.unwrap(np.radians(unp.nominal_values(sw12_phase))))); \
        sw12_complex += unp.nominal_values(sw12_mag)*\
                            np.cos(np.unwrap(np.radians(unp.nominal_values(sw12_phase))))

        # Split measured S-parameters into mag and phase
        sm11_mag = self.s11[0]
        sm22_mag = self.s22[0]
        sm21_mag = self.s21[0]
        sm12_mag = self.s12[0]
        sm11_phase = self.s11[1]
        sm22_phase = self.s22[1]
        sm21_phase = self.s21[1]
        sm12_phase = self.s12[1]
        # Cast to complex and unwarp phase
        sm11_complex = 1j*(unp.nominal_values(sm11_mag)*\
                           np.sin(np.unwrap(np.radians(unp.nominal_values(sm11_phase))))); \
        sm11_complex += unp.nominal_values(sm11_mag)*\
                            np.cos(np.unwrap(np.radians(unp.nominal_values(sm11_phase))))
        sm22_complex = 1j*(unp.nominal_values(sm22_mag)*\
                           np.sin(np.unwrap(np.radians(unp.nominal_values(sm22_phase))))); \
        sm22_complex += unp.nominal_values(sm22_mag)*\
                            np.cos(np.unwrap(np.radians(unp.nominal_values(sm22_phase))))
        sm21_complex = 1j*(unp.nominal_values(sm21_mag)*\
                           np.sin(np.unwrap(np.radians(unp.nominal_values(sm21_phase))))); \
        sm21_complex += unp.nominal_values(sm21_mag)*\
                            np.cos(np.unwrap(np.radians(unp.nominal_values(sm21_phase))))
        sm12_complex = 1j*(unp.nominal_values(sm12_mag)*\
                           np.sin(np.unwrap(np.radians(unp.nominal_values(sm12_phase))))); \
        sm12_complex += unp.nominal_values(sm12_mag)*\
                            np.cos(np.unwrap(np.radians(unp.nominal_values(sm12_phase))))
                            
        # Check for equal length and same frequency points
        #   In the future, consider interpolation in some cases
        # np.allclose retunrs True if equal within a low tolerance
        if not np.allclose(wfreq,self.freq):
            # Output the washer freq array to compare to self.freq for 
            #   troubleshooting
            global washer_check
            washer_check = wfreq
            raise Exception('Measurement must have same 601 frequency points')
        
        ## De-embed
        # Convert to T-parameters
        # Washers
        tw11_left = -(sw11_complex*sw22_complex-sw12_complex*\
                          sw21_complex)/sw21_complex
        tw12_left = sw11_complex/sw21_complex
        tw21_left = -sw22_complex/sw21_complex
        tw22_left = 1/sw21_complex
        # Measured
        tm11 = -(sm11_complex*sm22_complex-sm12_complex*\
                 sm21_complex)/sm21_complex
        tm12 = sm11_complex/sm21_complex
        tm21 = -sm22_complex/sm21_complex
        tm22 = 1/sm21_complex
        # Make matrices
        left_matrix = np.dstack([tw11_left,tw12_left,tw21_left,\
                                 tw22_left]).reshape(len(sw11_complex),2,2)
        right_matrix = left_matrix  # Washers are symetrical
        meas_matrix = np.dstack([tm11,tm12,tm21,\
                                 tm22]).reshape(len(sm11_complex),2,2)
        # Perform de-embeding
        corr_Tmatrix = np.matmul(np.linalg.inv(left_matrix),meas_matrix)
        corr_Tmatrix = np.matmul(corr_Tmatrix,np.linalg.inv(right_matrix))
        # Re-convert to S-parameters
        corr_s11_complex = corr_Tmatrix[:,0,1]/corr_Tmatrix[:,1,1]
        corr_s12_complex = (corr_Tmatrix[:,0,0]*corr_Tmatrix[:,1,1]-\
                    corr_Tmatrix[:,0,1]*corr_Tmatrix[:,1,0])/corr_Tmatrix[:,1,1]
        corr_s21_complex = 1/corr_Tmatrix[:,1,1]
        corr_s22_complex = -corr_Tmatrix[:,1,0]/corr_Tmatrix[:,1,1]
        # Re-cast to mag and phase
        # Use arctan2 to compute phase since it keeps track of signs (wraps)
        corr_s11_phase = np.arctan2(corr_s11_complex.imag,corr_s11_complex.real)
        corr_s12_phase = np.arctan2(corr_s12_complex.imag,corr_s12_complex.real)
        corr_s21_phase = np.arctan2(corr_s21_complex.imag,corr_s21_complex.real)
        corr_s22_phase = np.arctan2(corr_s22_complex.imag,corr_s22_complex.real)
        # Final corrected arrays
        corr_s11 = np.array([np.sqrt(corr_s11_complex.real**2 + \
                                     corr_s11_complex.imag**2),\
                                     corr_s11_phase*(180/np.pi)])
        corr_s12 = np.array([np.sqrt(corr_s12_complex.real**2 + \
                                     corr_s12_complex.imag**2),\
                                     corr_s12_phase*(180/np.pi)])
        corr_s21 = np.array([np.sqrt(corr_s21_complex.real**2 + \
                                     corr_s21_complex.imag**2),\
                                     corr_s21_phase*(180/np.pi)])
        corr_s22 = np.array([np.sqrt(corr_s22_complex.real**2 + \
                                     corr_s22_complex.imag**2),\
                                     corr_s22_phase*(180/np.pi)])
            
        return corr_s11, corr_s21, corr_s12, corr_s22
        
    def draw_plots(self,default_setting=True,corr=False,publish=False):
        """
        Plots permittivity data using make_plot from permittivity_plot_V1.
        
        Default is to plot average s-param results.
        
        Arguments
        ---------
        default_setting (bool): if True, plots average dielectric constant, \
            loss factor and loss tangent. If false prompts user to determine \
            wether to plot, average, forward, reverse, or both s-parameter \
            results. Default: True
        
        corr (bool): If True, use corrected sparam data. Default: False
        
        publish (bool): If True, save figures.
        """
        # Check if using corrected data
        if corr and default_setting:
            avg_dielec = self.corr_avg_dielec
            avg_lossfac = self.corr_avg_lossfac
            avg_losstan = self.corr_avg_losstan
        elif corr and not default_setting:
            raise Exception('default_setting must be True if using corrected data')
        else:
            avg_dielec = self.avg_dielec
            avg_lossfac = self.avg_lossfac
            avg_losstan = self.avg_losstan
        # Figure out what to plot
        x = self.freq
        kwargs = {}
        if default_setting:
            y1 = avg_dielec
            y2 = avg_lossfac
            y3 = avg_losstan
        else:
            # Prompt user
            s_plot = input('Please designate "f" for Forward, "r" for ' + \
                           'Reverse, "b" for Both, or "a" for average ' + \
                           'S-Parameters: ')
            if s_plot not in ('f','r','b','a'):
                raise Exception('Wrong input')
            if s_plot == 'a':
                y1 = avg_dielec
                y2 = avg_lossfac
                y3 = avg_losstan
            elif s_plot == 'f':
                y1 = self.forward_dielec
                y2 = self.forward_lossfac
                y3 = self.forward_losstan
            elif s_plot == 'r':
                y1 = self.reverse_dielec
                y2 = self.reverse_lossfac
                y3 = self.reverse_losstan
            else:
                x = [self.freq,self.freq]
                y1 = [self.forward_dielec,self.reverse_dielec]
                y2 = [self.forward_lossfac,self.reverse_lossfac]
                y3 = [self.forward_losstan,self.reverse_losstan]
                kwargs = {"legend_label":['Forward','Reverse']}
        if publish:
            kwargs['publish'] = True
            kwargs['name'] = self.name
        pplot.make_plot(x,y1,'d',**kwargs)
        pplot.make_plot(x,y2,'lf',**kwargs)
        pplot.make_plot(x,y3,'lt',**kwargs)
        
    def s_param_plot(self,corr=False):
        """Plot raw S-Parameter data using make_sparam_plot from \
        permittivity_plot_V1"""
        if corr:
            pplot.make_sparam_plot(self.freq,[self.s11,self.corr_s11],\
                                   [self.s22,self.corr_s22],[self.s21,\
                                   self.corr_s21],[self.s12,self.corr_s12],\
                                   label=['Uncorrected','Corrected'])
        else:
            pplot.make_sparam_plot(self.freq,self.s11,self.s22,self.s21,self.s12)

#%% FUCNTIONS
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
        print('\n')
        print('Select the VNA Tools II Output Data Table')
        print('\n')
        root = tk.Tk()
        root.withdraw()
        file = filedialog.askopenfilename(filetypes=[('text files', '*.txt')],\
                                title='Select the VNA Tools II Output Data Table')
    else:   # Prompt for both
        # Airline
        airline = input('Are you using the "VAL", "PAL", "GAL", "7" mm, ' + \
                        'or "custom" airline?: ')
        if airline not in ('VAL','PAL','GAL','7','custom'):
            raise Exception('Wrong input')
        elif airline == 'custom':
            L_in = input('Enter the length of the airline (cm): ')
        # File
        print('\n')
        print('Select the VNA Tools II Output Data Table')
        print('\n')
        root = tk.Tk()
        root.withdraw()
        file = filedialog.askopenfilename(filetypes=[('text files', '*.txt')],\
                                title='Select the VNA Tools II Output Data Table')

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

def perm_compare(classlist,allplots=False):
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
        dielec.append(item.avg_dielec)
        losstan.append(item.avg_losstan)
        labels.append(item.name)
    kwargs = {"legend_label":labels}
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
        
def run_default():
    """
    Run AirlineData on get_METAS_data with all the prompts and return the \
        instance.
    """
    return AirlineData(*get_METAS_data())
                
#%% MAIN
def main():
    ## Single file example:
    #global test
    #test = AirlineData(*get_METAS_data(airline='GAL',file_path=DATAPATH + \
    #                    '2.5hrs.txt'),bulk_density=2.0,temperature=None,\
    #                     name='Alumina Vac 2.5hrs',date='2017/04/07')
    
    ## Multiple file example:
    #test2 = AirlineData(*get_METAS_data(),name='TRM')
    #classlist = [test,test2]
    #perm_compare(classlist)
    pass    # Comment to run example
    
if __name__ == '__main__':
    #main()
    pass    # Comment to run example
