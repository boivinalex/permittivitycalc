#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Contains the AirlineData class object"""

# Array math
import numpy as np
import uncertainties 
#Citation: Uncertainties: a Python package for calculations with uncertainties,
#    Eric O. LEBIGOT, http://pythonhosted.org/uncertainties/
from uncertainties import unumpy as unp
# Plotting
from . import permittivity_plot as pplot
# Make relative path
import os
# Helper functions
from .helper_functions import get_METAS_data, _get_file, perm_compare

# GLOBAL VARIABLES
E_0 = 8.854187817620*10**-12 #Permittivity of vacuum (F/m) 
U_0 = 4*np.pi*10**-7 #Permeability of vacuum (V*s/A*m) 
C = (1)*1/np.sqrt(E_0*U_0) #Speed of light (m/s)
LAM_C = float('inf') #Cut-off wavelength = infinity
PACKAGEPATH = os.path.dirname(os.path.abspath(__file__))
DATAPATH = os.path.join(PACKAGEPATH, 'data')


class AirlineData:
    """
    Use S-parameter data from a METAS VNA Tools II output text file to 
    calculate the complex relative permittivity of a sample in a coaxial 
    transmission line
    
    Parameters:
    -----------
    L : float 
        Length of airline in cm.
        
    airline : {'VAL','PAL','GAL','7','custom'} 
        Name of airline used for measurement.
        
    dataArray : array
        Array containing raw measurement data.
        
    file : str 
        Input file path.
        
    corr : bool, optional
        Default = False. If True, also correct S-parameter data and produce 
        corr_* arrays. 
            
    nrw : bool, optional 
        Default = False. If True, use Nicolson, Rross, Weir (NRW) algorithm to 
        calculate permittivity and magnetic permeability.
            
    shorted : bool, optional
        Default = False. If True, automatically load Shorted S11 data. File 
        name must have the following format and be in the same folder as 
        original file: file_path/file_name_shorted.txt
        
        Example:
            - file_path/air_atm.txt
            - file_path/air_atm_shorted.txt
            
    normalize_density : bool or float or int, optional 
        Default = False. If True, use either Lichtenecker or 
        Landau-Lifshitz-Looyenga equation to normalize the complex 
        permittivity to a constant density of 1.60 g/cm^3. If float or int 
        is given, normalize to the density provided in g/cm^3. Requires 
        that bulk_density be given.
        
    norm_eqn : {'LI','LLL'}, optional 
        Default = 'LI'. For use with normalize_density. Equation to be used for 
        normalization. Options are 'LI' for the Lichtenecker 
        equation and 'LLL' for the Landau-Lifshitz-Looyenga equation. 
        LI used alpha = 1.92 (Olhoeft, 1985) and LLL uses 
        alpha = 0.307 (Hickson et al., 2017, Lunar samples).
            
    name : str, optional 
        Name of measurement to be used in plots.
            
    bulk_density : float, optional
        Bulk density of material.
        
    solid_dielec : float, optional 
        The solid dielectric constant of the material.
        
    solid_losstan : float, optional
        The solid loss tangent of the material.
        
    particle_diameter : float, optional
        The average particle diameter in the airline in cm.
        
    particle_density : float, optional
        The average (solid) particle density of the material in g/cm^3.
        
    temperature : str or float, optional 
        Temperature of measurement.
        
    date : str, optional
        Date measurement was made.
            
    freq_cutoff : float, optional
        Default: 4e8. Starting frequency for forward and reverse difference 
        calculations.
    
    Attributes:
    ----------
    freq : array 
        Frequency points.
    
    s11, s21, s12, s22 : array 
        Mag and Phase S-Parameters.
    
    *_dielec : array 
        Real part of the permittivity. Can be avg_dielec, forward_dielec, or 
        reverse_dielec for average, forward (S11,S21), or reverse (S22,S12) 
        permittivity.
    
    *_lossfac : array 
        Imaginary part of the permittivity. Same as above.
    
    *_losstan : array 
        Loss tangent. Same as above.
        
    *_mu_real : array
        Real part of the magnetic permeability. Same as above.
        
    *_mu_imag : array
        Imaginary part of the magnetic permeability. Same as above.
    
    corr_* : array 
        De-embeded version of S-parameters or electromagnetic data. Only average 
        S-parameters are used for EM calculations with corrected 
        S-parameters. Examples: corr_s11, corr_avg_losstan, corr_avg_mu_real. 
        Only created if corr = True.
            
    norm_* : array
        Bulk density normalized permittivity data. Uses averaged permittivity 
        data. Only created if normalize_density = True.
        
    res_freq : array 
        Resonant frequencies in the sample.
        
    airline_dimensions : dict 
        Dimensions of the airline in cm. D1 is the diameter of the inner 
        conductor and D4 is the diameter of the outer conductor. D2 and D3 
        bound the sample-airline boundary regions if particle_diameter is 
        provided. airline_dimensions is generated automatically for 
        airlines VAL, PAL, and GAL. Empty otherwise.
        
    bcorr : complex array 
        Avg complex permittivity corrected for boundary effects. Computed 
        automatically if solid_dielec, particle_diameter, particle_density, 
        and bulk_density are present. solid_losstan is optional.
        
    real_part_diff_array, imag_part_diff_array : array
        Absolute difference of forward and reverse results for the real and 
        imaginary parts of the permittivity.
            
    max_real_diff, max_imag_diff, min_real_diff, min_imag_diff : float
        Maximum and minimum differences of the forward and reverse results for 
        the real and imaginary parts of the permittivity.
    """
    def __init__(self,L,airline,dataArray,file,name=None,date=None,\
                 freq_cutoff=4e8,nrw=False,shorted=False,corr=False,\
                 normalize_density=False,norm_eqn='LI',bulk_density=None,\
                 solid_dielec=None,solid_losstan=None,particle_diameter=None,\
                 particle_density=None,temperature=None):
        # Required Attributes
        self.L = L
        self.airline_name = airline
        self.file = file
        
        # Optional attributes
        self.name = name
        self.date = date
        self.freq_cutoff = freq_cutoff
        self.nrw = nrw
        self.shorted = shorted
        self.corr = corr
        self.normalize_density = normalize_density
        self.norm_eqn = norm_eqn
        self.bulk_density = bulk_density
        self.solid_dielec = solid_dielec
        self.solid_losstan = solid_losstan
        self.particle_diameter = particle_diameter
        self.particle_density = particle_density
        self.temperature = temperature
        
        # Get the airline dimnetions
        self.airline_dimensions = self._dims()
        
        # Unpack data into arrays
        self.freq, self.s11, self.s21, self.s12, self.s22 = \
            self._unpack(dataArray)
        if self.shorted:
            #Get path of shorted file
            path = os.path.split(self.file)
            folder = path[0]
            file_name_full = path[1]
            file_name = os.path.splitext(file_name_full)[0]
            file_name_shorted = file_name + '_shorted.txt'
            path_shorted = folder + '/' + file_name_shorted
            print(path_shorted)
            # Check if file exists
            if os.path.isfile(path_shorted):
                dataArray2 = get_METAS_data(airline=self.airline_name,file_path=path_shorted)
                self.freq_short, self.s11_short = self._unpack(dataArray2[2])
            else:
                raise Exception('Shorted file does not exists. Sould be in the form: orig-filename_shorted.txt')
        
        # Calculate permittivity
        if nrw:
            self.forward_dielec, self.forward_lossfac, self.forward_losstan, \
                self.forward_mu_real, self.forward_mu_imag \
                = self._permittivity_calc('f')
            self.reverse_dielec, self.reverse_lossfac, self.reverse_losstan, \
                self.reverse_mu_real, self.reverse_mu_imag \
                = self._permittivity_calc('r')
            self.avg_dielec = (self.forward_dielec + self.reverse_dielec)/2
            self.avg_lossfac = (self.forward_lossfac + self.reverse_lossfac)/2
            self.avg_mu_real = (self.forward_mu_real + self.reverse_mu_real)/2
            self.avg_mu_imag = (self.forward_mu_imag + self.reverse_mu_imag)/2
            self.avg_losstan = self.avg_lossfac/self.avg_dielec
        else:
            self.forward_dielec, self.forward_lossfac, self.forward_losstan = \
                self._permittivity_calc('f')
            self.reverse_dielec, self.reverse_lossfac, self.reverse_losstan = \
                self._permittivity_calc('r')
            self.avg_dielec = (self.forward_dielec + self.reverse_dielec)/2
            self.avg_lossfac = (self.forward_lossfac + self.reverse_lossfac)/2
            self.avg_losstan = self.avg_lossfac/self.avg_dielec

        # Try to calculate corrected (de-embedded) permittivity
        # Fail if NaNs in data
        if corr:
            try:
                if not np.isnan(unp.nominal_values(self.avg_dielec)).any():
                    if isinstance(self.corr, (float,int)):
                        self.Lcorr = self.L - 2*(self.corr) # Remove total washer length 
                    else:
                        self.Lcorr = self.L - 0.34 # Remove total washer length 
                    self.corr_s11, self.corr_s21, self.corr_s12, \
                        self.corr_s22 = self._de_embed()
                    if nrw:
                        self.corr_forward_dielec, self.corr_forward_lossfac, \
                            self.corr_forward_losstan, \
                            self.corr_forward_mu_real, \
                            self.corr_forward_mu_imag \
                            = self._permittivity_calc('f',True)
                        self.corr_reverse_dielec, self.corr_reverse_lossfac, \
                            self.corr_reverse_losstan, \
                            self.corr_reverse_mu_real, \
                            self.corr_reverse_mu_imag \
                            = self._permittivity_calc('r',True)
                        self.corr_avg_dielec = (self.corr_forward_dielec + \
                                                 self.corr_reverse_dielec)/2
                        self.corr_avg_lossfac = (self.corr_forward_lossfac + \
                                                 self.corr_reverse_lossfac)/2
                        self.corr_avg_losstan = \
                            self.corr_avg_lossfac/self.corr_avg_dielec
                        self.corr_avg_mu_real = (self.corr_forward_mu_real + \
                                           self.corr_reverse_mu_real)/2
                        self.corr_avg_mu_imag = (self.corr_forward_mu_imag + \
                                           self.corr_reverse_mu_imag)/2
                    else:
                        self.corr_forward_dielec, self.corr_forward_lossfac, \
                            self.corr_forward_losstan = \
                            self._permittivity_calc('f',True)
                        self.corr_reverse_dielec, self.corr_reverse_lossfac, \
                            self.corr_reverse_losstan = \
                            self._permittivity_calc('r',True)
                        self.corr_avg_dielec = (self.corr_forward_dielec + \
                                                 self.corr_reverse_dielec)/2
                        self.corr_avg_lossfac = (self.corr_forward_lossfac + \
                                                 self.corr_reverse_lossfac)/2
                        self.corr_avg_losstan = self.corr_avg_lossfac/self.corr_avg_dielec
            except:
                raise Warning('S-parameter correction failed. Using uncorrected data. Check if NaNs in data.')
            
        # Calculate percentage difference between forward and reverse permittivity
        self.real_part_diff_array, self.imag_part_diff_array, self.tand_part_diff_array,\
            self.max_real_diff, self.max_imag_diff, self.max_tand_diff, \
            self.med_real_diff, self.med_imag_diff, self.med_tand_diff = \
            self._absdiff()
        # Print maximum and median difference as percentages
        diff_results = '\n' + str(self.name) + '\n' + \
            'The maximum precentage difference between forward (S11/S21) and reverse (S22/S12) calculated permittivity is: \n' + \
            'ε′: {:.2%} '.format(self.max_real_diff) + \
            'ε′′: {:.2%} '.format(self.max_imag_diff) + \
            'tanδ: {:.2%} \n'.format(self.max_tand_diff) + \
            'The median precentage difference between forward (S11/S21) and reverse (S22/S12) calculated permittivity is: \n' + \
            'ε′: {:.2%} '.format(self.med_real_diff) + \
            'ε′′: {:.2%} '.format(self.med_imag_diff) + \
            'tanδ: {:.2%} \n'.format(self.med_tand_diff)
        print(diff_results)
        
        # Combine Type A and Type B Uncertainty (if it exists)
        if isinstance(self.s11[0][0], uncertainties.UFloat) and not self.nrw:
            self.avg_dielec, self.avg_lossfac, self.avg_losstan = \
                self._calcTotalUncertainty(self.avg_dielec,self.avg_lossfac,\
                self.avg_losstan)
            if corr:
                self.corr_avg_dielec, self.corr_avg_lossfac,\
                    self.corr_avg_losstan = self._calcTotalUncertainty(\
                    self.corr_avg_dielec,self.corr_avg_lossfac,\
                    self.corr_avg_losstan)
                    
        # Calculare resonant frequencies in sample
        self.res_freq = self._resonant_freq()
        
        # Calculate an average real permittivity and loss tan value from 
        #   midpoint frequency values between resonant frequencies.
        self.freq_avg_dielec, self.freq_avg_losstan, self.freq_avg_dielec_std, \
            self.freq_avg_losstan_std = self._freq_avg()
            
        # If normalize_density is not False and bulk_density exists, normalize
        if self.normalize_density and self.bulk_density:
            if self.corr:   #use corrected data is present
                complex_dielec = 1j*unp.nominal_values(self.corr_avg_lossfac);
                complex_dielec += unp.nominal_values(self.corr_avg_dielec)
            else:
                complex_dielec = 1j*unp.nominal_values(self.avg_lossfac);
                complex_dielec += unp.nominal_values(self.avg_dielec)
            # If normalize_density is True, normalize to to 1.60 g/cm^3
            #   If normalize_density is a value, normalize to that value
            if isinstance(self.normalize_density, bool):
                norm_val = 1.60
            elif isinstance(self.normalize_density, (float,int)):
                norm_val = self.normalize_density
            if self.norm_eqn == 'LI':
                # Lichtenecker equation
                #   alpha from Carrier et al., 1991
                norm_complex_dielec = complex_dielec*((1.92)**\
                    (norm_val-self.bulk_density))
            elif self.norm_eqn == 'LLL':
                # Landau-Lifshitz-Looyenga equation
                #   alpha from Hickson et al., 2018
                norm_complex_dielec = complex_dielec*((norm_val*0.307 + 1)**3 / \
                                            (self.bulk_density*0.307 + 1)**3)
            # Get uncertainty
            if self.corr:
                unc_real = unp.std_devs(self.corr_avg_dielec)
                unc_imag = unp.std_devs(self.corr_avg_lossfac)
            else:
                unc_real = unp.std_devs(self.avg_dielec)
                unc_imag = unp.std_devs(self.avg_lossfac)
            self.norm_dielec = unp.uarray(np.real(norm_complex_dielec),unc_real)
            self.norm_lossfac = unp.uarray(np.imag(norm_complex_dielec),unc_imag)
            self.norm_losstan = self.norm_lossfac/self.norm_dielec
        elif normalize_density:
            raise Exception('Need bulk desnity to normalize to constant density')
                    
        # If appropriate data provided, correct for boundary effects
        if (solid_dielec and particle_diameter and particle_density and \
            bulk_density):
            self.bcorr_dielec, self.bcorr_losstan = self._boundary_correct()
            
    def __repr__(self):
        rep = 'pc.AirlineData(*pc.get_METAS_data(airline=%r,file_path=%r),' % \
                (self.airline_name,self.file) + \
                'bulk_density=%r,temperature=%r,name=%r,date=%r,corr=%r' % \
                (self.bulk_density,self.temperature,self.name,self.date,\
                 self.corr) + ',solid_dielec=%r,solid_losstan=%r' % \
                 (self.solid_dielec,self.solid_losstan) + \
                 ',particle_diameter=%r,particle_density=%r' % \
                 (self.particle_diameter, self.particle_density) + \
                 ',nrw=%r,normalize_density=%r,norm_eqn=%r' % \
                 (self.nrw, self.normalize_density, self.norm_eqn) + \
                 ',shorted=%r,freq_cutoff=%r)' % \
                 (self.shorted, self.freq_cutoff)
        return rep
        
    def __str__(self):
        srep = 'measured in ' + self.airline_name + ' (L = ' + str(self.L) + ')'
        if self.name:
            srep = self.name + ' ' + srep
        if self.date:
            srep += ' on ' + str(self.date)
        if self.temperature:
            srep += ' at ' + str(self.temperature) + ' degrees'
        if self.bulk_density:
            srep += ' with a bulk density of ' + str(self.bulk_density) + ' g/cm^3'
        srep += ' from file: \n' + self.file
        if self.corr:
            srep += '\n' + 'Corrected (de-embeded) data is available.'
        if self.normalize_density:
            if isinstance(self.normalize_density, bool):
                norm_val = str(1.60)
            elif isinstance(self.normalize_density, (float,int)):
                norm_val = str(self.normalize_density)
            srep += '\n' + 'Normalized data at a bulk density of ' + norm_val \
            + ' g/cm^3'
            if self.norm_eqn == 'LI':
                srep += ' using the Lichtenecker equation'
            elif self.norm_eqn == 'LLL':
                srep += ' using the Landau-Lifshitz-Looyenga equation'
            srep += ' is available.'
        return srep
    
    def _unpack(self,dataArray):
        """See if uncertainty in data and unpack to S-parameter arrays"""
        shorted_flag = False
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
        elif self.shorted and dataArray.shape[1] == 5:
            shorted_flag = True
            freq = dataArray[:,0]
            s11 = unp.uarray([dataArray[:,1],dataArray[:,3]],\
                                     [dataArray[:,2],dataArray[:,4]])
        else:
            raise Exception('Input file has the wrong number of columns')
        
        if shorted_flag:
            return freq, s11
        else:
            return freq, s11, s21, s12, s22
    
    def _dims(self):
        """
        Determine the dimensions of the airline used in cm.
        """
        dims = {}
        # Store inner and outer diameters in a dictionary
        if self.airline_name in ('VAL','PAL'):
            dims['D1'] = 6.204
            dims['D4'] = 14.288
        elif self.airline_name == 'GAL':
            dims['D1'] = 6.19
            dims['D4'] = 14.32
        # If particle diameter is known, calculate D2 and D3 for boundary 
        #   effect correction   
        if dims and self.particle_diameter:
            dims['D2'] = dims['D1'] + self.particle_diameter
            dims['D3'] = dims['D4'] - self.particle_diameter
        return dims
    
    def _absdiff(self):
        """
        Calculate the absolute max and median differences in calculated forward 
        (S11,S21) and reverse (S22,S12) permittivity. Use corrected 
        S-parameters if present.
        """
        if self.corr:
            real_part_diff = np.abs(unp.nominal_values(self.corr_forward_dielec) \
                            - unp.nominal_values(self.corr_reverse_dielec))
            imag_part_diff = np.abs(unp.nominal_values(self.corr_forward_lossfac) \
                            - unp.nominal_values(self.corr_reverse_lossfac))
            tand_part_diff = np.abs(unp.nominal_values(self.corr_forward_losstan) \
                            - unp.nominal_values(self.corr_reverse_losstan))
        else:
            real_part_diff = np.abs(unp.nominal_values(self.forward_dielec) \
                            - unp.nominal_values(self.reverse_dielec))
            imag_part_diff = np.abs(unp.nominal_values(self.forward_lossfac) \
                            - unp.nominal_values(self.reverse_lossfac))
            tand_part_diff = np.abs(unp.nominal_values(self.forward_losstan) \
                            - unp.nominal_values(self.reverse_losstan))
            
        if self.corr:
            real_avg = unp.nominal_values(self.corr_avg_dielec)
            imag_avg = unp.nominal_values(self.corr_avg_lossfac)
            tand_avg = unp.nominal_values(self.corr_avg_losstan)
        else:
            real_avg = unp.nominal_values(self.avg_dielec)
            imag_avg = unp.nominal_values(self.avg_lossfac)
            tand_avg = unp.nominal_values(self.avg_losstan)
            
        if self.freq_cutoff:
            real_part_diff_calc = real_part_diff[self.freq >= self.freq_cutoff]\
                /real_avg[self.freq >= self.freq_cutoff]
            imag_part_diff_calc = imag_part_diff[self.freq >= self.freq_cutoff]\
                /imag_avg[self.freq >= self.freq_cutoff]
            tand_part_diff_calc = tand_part_diff[self.freq >= self.freq_cutoff]\
                /tand_avg[self.freq >= self.freq_cutoff]
        else:
            real_part_diff_calc = real_part_diff/real_avg
            imag_part_diff_calc = imag_part_diff/imag_avg
            tand_part_diff_calc = tand_part_diff/tand_avg
        
        max_real_diff = np.max(real_part_diff_calc)
        max_imag_diff = np.max(imag_part_diff_calc)
        max_tand_diff = np.max(tand_part_diff_calc)
        avg_real_diff = np.median(real_part_diff_calc)
        avg_imag_diff = np.median(imag_part_diff_calc)
        avg_tand_diff = np.median(tand_part_diff_calc)
        
        return real_part_diff, imag_part_diff, tand_part_diff, max_real_diff, \
            max_imag_diff, max_tand_diff, avg_real_diff, avg_imag_diff, \
            avg_tand_diff
            
    def _resonant_freq(self):
        """
        Calculate and return array of resonant frequencies from complex 
        permittivity and/or permeability measurement.
    
        Follows David Stillman PhD Thesis

        """
        
        n = np.arange(1,40) # Max number of resonances (overkill)
        if self.corr:
            measured_dielec = unp.nominal_values(self.corr_avg_dielec)
            measured_lossfac = unp.nominal_values(self.corr_avg_lossfac)
            L = self.Lcorr
            if self.nrw:
                global u_r
                u_r = unp.nominal_values(self.corr_avg_mu_real)
                u_i = unp.nominal_values(self.corr_avg_mu_imag)
        else:
            measured_dielec = unp.nominal_values(self.avg_dielec)
            measured_lossfac = unp.nominal_values(self.avg_lossfac)
            L = self.L
            if self.nrw:
                u_r = unp.nominal_values(self.avg_mu_real)
                u_i = unp.nominal_values(self.avg_mu_imag)
        e_r = np.median(measured_dielec[1::]) # Exclude first data point
        e_i = np.median(measured_lossfac[1::])
        if self.nrw:
            u_r = np.median(u_r[1::])
            u_i = np.median(u_i[1::])
        else:
            u_r = 1
            u_i = 0
        
        res_freq = ((2**(1/2))*C*100)/((2*L/n)*((((u_r*e_r - e_i*u_i)**2 + \
                    (u_i*e_r + e_i*u_r)**2)**(1/2)) + e_r*u_r - e_i*u_i)**(1/2))
        # Restrict res_freq to max freq
        res_freq = res_freq[res_freq<=np.max(self.freq)]
        return res_freq
    
    def _freq_avg(self):
        """
        Calculate an average dielectric constant and loss tangent across all 
        measured frequencies from midpoint frequency values between resonant 
        frequencies as described in Hickson et al., 2018
        
        Return
        ------
        freq_avg_dielec : uncertainties.core.AffineScalarFunc 
            Average dielectric constant calculated from midpoint frequency 
            values between resonant frequencies. Skips first two values due to 
            large uncertainty. 
        
        freq_avg_dielec_std : float 
            Standard deviation of freq_avg_dielec from above.
        
        freq_avg_losstan : uncertainties.core.AffineScalarFunc 
            Average loss tangent calculated from midpoint frequency values 
            between resonant frequencies. Skips first two values due to large 
            uncertainty. 
        
        freq_avg_losstan_std : float 
            Standard deviation of freq_avg_losstan from above.
        
        """
        f_r = self.res_freq
        if self.corr:
            dielec = self.corr_avg_dielec
            lossfac = self.corr_avg_lossfac
        else:
            dielec = self.avg_dielec
            lossfac = self.avg_lossfac
        mid_pts  = np.zeros(len(f_r)-1)#, dtype=int)
        f_0_mids = np.zeros(len(f_r)-1)#, dtype=int)
        dielec_mids = []
        loss_tan_mids = []
        f_r = f_r[f_r < np.max(self.freq)]
        for i in range(0, len(f_r)):
            if i != (len(f_r)-1):
                # Find the midpoint frequencies between resonances
                x = (f_r[i+1] - f_r[i])/2
                mid_pts[i] = f_r[i] + x
                # Find the closest corresponding frequency index in f_0
                tmp = np.abs(self.freq - mid_pts[i])
                idx = np.argmin(tmp) # index of closest value
                f_0_mids[i] = self.freq[idx]
                dielec_mids.append(dielec[idx])
                loss_tan_mids.append(lossfac[idx]/dielec[idx])        
        # Calculate averages past first two values 
        freq_avg_dielec  = np.mean(dielec_mids[2:])
        std_dielec = np.std(unp.nominal_values(dielec_mids[2:]))
        freq_avg_losstan = np.mean(loss_tan_mids[2:])
        std_losstan = np.std(unp.nominal_values(loss_tan_mids[2:]))
        
        return freq_avg_dielec, freq_avg_losstan, std_dielec, std_losstan
    
    def _permittivity_calc(self,s_param,corr=False):
        """
        Return the complex permittivity and loss tangent from S-parameters 
        and their uncertainties (if present).
    
        Uses the New Non-iterative method described in Boughriet, A.H., Legrand, 
        C. & Chapoton, A., 1997. Noniterative stable transmission/reflection 
        method for low-loss material complex permittivity determination. 
        IEEE Transactions on Microwave Theory and Techniques, 45(1), pp.52–57.
        
        Assumes the magnetic permeability of the sample = 1 (i.e non-magnetic).
        
        Arguments
        ---------
        s_param : str 
            Calculates complex permittivity using either 
            the forward ('f') (S11,S21), or the reverse ('r') (S22 S12) 
            S-Parameters.
            
        Return
        ------
        f_0 : array 
            Array of measured frequency values.
        
        dielec : array 
            Real part of the complex relative permittivity.
        
        lossfac : array 
            Imaginary part of the complex permittivity.
        
        losstan : array 
            Loss tangent. Defined as the ratio of the imaginary 
            and the real part of the complex permittivity.
        """
        
        # Check if using corrected S-params
        if corr:
            s11 = self.corr_s11
            s21 = self.corr_s21
            s12 = self.corr_s12
            s22 = self.corr_s22
            L = self.Lcorr
        else:
            s11 = self.s11
            s21 = self.s21
            s12 = self.s12
            s22 = self.s22
            L = self.L
        
        # Strip arrays of their uncertainties
        #   This measure is temporary as the uncertainties module does not 
        #   support complex numbers at the time of writing.
        #   Later, uncertainties may be propagated through and compared with
        #   those calculated in Boughriet et al.
        if s_param == 'f':
            s_r = unp.nominal_values(s11)
            s_t = unp.nominal_values(s21)
        else:
            s_r = unp.nominal_values(s22)
            s_t = unp.nominal_values(s12)
            
        # Initialize arrays
        nrows = len(self.freq)
        gam = np.zeros(nrows,dtype=complex)
        sign = np.empty(nrows,dtype=bool)
        a = np.zeros(nrows,dtype=complex)
        
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
        # From Chen, L. F. et al. (2004). Microwave Electronics: Measurement and Materials Characterization
        t_phaser = np.angle(t) #get phase angle from complex
        t_phase_unwrap = np.copy(t_phaser)
        # Ignore NaN values
        t_phase_unwrap[~np.isnan(t_phase_unwrap)] = \
                       np.unwrap(t_phase_unwrap[~np.isnan(t_phase_unwrap)])
                       
        # Calculate ln(1/T) with correct phases
        ln_1_T = 1j*t_phase_unwrap; ln_1_T += np.log(1/abs(t))
        # Also create new unwrapped T vector
        new_t = 1j*t_phase_unwrap; new_t += np.absolute(t)
        
        # Calculate 1/A
        a_1 = np.sqrt(-(ln_1_T / (2*np.pi*L))**2)
         
        # Calculate A
        a = 1 / a_1
          
        # Find wavelength in empty cell
        lam_og = 1 / np.sqrt(1/lam_0**2 - 1/LAM_C**2)
        
        # Calculated effective electromagnetic parameters
        if self.nrw:
            # Calculate mu_r (relative permeability)
            mu_r = (1+gam)/((a*(1-gam))*(np.sqrt((1/lam_0**2)-(1/LAM_C**2))))
            mu_eff = mu_r
            # Calculate e_r (relative permittivity)
            ep_r = (mu_r*(((1-gam)**2)/((1+gam)**2))*(1-(lam_0**2/LAM_C**2))) \
                    + ((lam_0**2/LAM_C**2)*(1/mu_r))
        else: 
            mu_eff = 1      #Assume mu_eff = mu_r = 1
            ep_eff = (lam_og/a)**2
            # Calculate e_r (relative permittivity)
            ep_r = (1 - lam_0**2/LAM_C**2)*ep_eff + \
            (lam_0**2/LAM_C*2)/mu_eff
        
        dielec = ep_r.real
        lossfac = ep_r.imag
        losstan = lossfac/dielec
        if self.nrw:
            complex_mu = mu_eff
            mu_real = complex_mu.real
            mu_imag = complex_mu.imag
        
        # Calculate uncertainties
        # Check for uncertainties and make sure not using NRW
        if isinstance(self.s11[0][0], uncertainties.UFloat) and not self.nrw: 
            delta_dielec, delta_lossfac = \
                self._calc_uncertainties(s_param,nrows,sign,x,s_reflect,\
                                         s_trans,gam,lam_og,new_t,mu_eff,\
                                         ep_eff,lam_0,dielec,lossfac,losstan,\
                                         corr)
            dielec = unp.uarray(dielec,delta_dielec)
            lossfac = unp.uarray(lossfac,delta_lossfac)
            losstan = losstan = lossfac/dielec  
        elif isinstance(self.s11[0][0], uncertainties.UFloat) and self.nrw: 
            delta_dielec, delta_lossfac, delta_mureal, delta_muimag = \
                self._calc_uncertainties_NRW(s_param,nrows,sign,x,s_reflect,\
                                         s_trans,gam,lam_og,new_t,complex_mu,\
                                         ep_r,lam_0,dielec,lossfac,losstan,\
                                         corr)
            dielec = unp.uarray(dielec,delta_dielec)
            lossfac = unp.uarray(lossfac,delta_lossfac)
            losstan = losstan = lossfac/dielec
            mu_real = unp.uarray(mu_real,delta_mureal)
            mu_imag = unp.uarray(mu_imag,delta_muimag)
        
        if self.nrw:
            return dielec,lossfac,losstan,mu_real,mu_imag
        else:
            return dielec,lossfac,losstan
        
    def _calc_uncertainties(self,s_param,nrows,sign,x,s_reflect,s_trans,gam,\
                            lam_og,new_t,mu_eff,ep_eff,lam_0,dielec,\
                            lossfac,losstan,corr):
        """
        Calculates permittivity uncertainties from Boughriet et al.
        
        Returns
        -------
        delta_dielec : array 
            Uncertainty in real part.
        
        delta_lossfac : array 
            Uncertainty in imaginary part.
        """
        # Check if using corrected S-params
        if corr:
            s11 = self.corr_s11
            s21 = self.corr_s21
            s12 = self.corr_s12
            s22 = self.corr_s22
            L = self.Lcorr
        else:
            s11 = self.s11
            s21 = self.s21
            s12 = self.s12
            s22 = self.s22
            L = self.L
        
        # Strip uarrays of their nominal values
        if s_param == 'f':
            err_s_r = unp.std_devs(s11)
            err_s_t = unp.std_devs(s21)
        else:
            err_s_r = unp.std_devs(s22)
            err_s_t = unp.std_devs(s12)
        
        # Initialize arrays
        delta_length = 0.002 #in cm
        dgam_dS_reflect = np.zeros(nrows,dtype=complex)
        dgam_dS_trans = np.zeros(nrows,dtype=complex)
        
        # Calculate partial derivatives
        # Determine +/- ambiguity by condition |GAMMA| < 1 from before            
        dgam_dS_reflect[sign] = (1 + (x[sign]/np.sqrt((x[sign]**2)-1)))*\
                        (((2*s_reflect[sign]**2)-(2*s_trans[sign]**2)+1)\
                         /2*s_reflect[sign]**2)
        dgam_dS_trans[sign] = (1 + (x[sign]/np.sqrt((x[sign]**2)-1)))*\
                        (-1*(s_trans[sign]/s_reflect[sign]))
        
        dgam_dS_reflect[np.invert(sign)] = (1 - \
                        (x[np.invert(sign)]/np.sqrt((x[np.invert(sign)]**2)-1)))\
                        *(((2*s_reflect[np.invert(sign)]**2)\
                           -(2*s_trans[np.invert(sign)]**2)+1)\
                        /2*s_reflect[np.invert(sign)]**2)
        dgam_dS_trans[np.invert(sign)] = (1 - (x[np.invert(sign)]\
                        /np.sqrt((x[np.invert(sign)]**2)-1)))\
                        *(-1*(s_trans[np.invert(sign)]/s_reflect[np.invert(sign)]))
     
        dT_dS_reflect = (1 - gam**2 + dgam_dS_reflect*(((s_reflect + s_trans)**2) \
                        - 1))/((1 - (s_reflect + s_trans)*gam)**2)
        dT_dS_trans = (1 - gam**2 + dgam_dS_trans*(((s_reflect + s_trans)**2) \
                        - 1))/((1 - (s_reflect + s_trans)*gam)**2)
        #deff_dgam = 0  # 0 in Boughriet, 1997
        deff_dT = 2*((1j*lam_og/(2*np.pi*L))*np.log(new_t))\
                        *((1j*lam_og)/(2*np.pi*L*new_t))
        
        deff_dS_reflect_mag = (deff_dT*dT_dS_reflect)*np.exp(1j*np.angle(s_reflect))
        deff_dS_trans_mag = (deff_dT*dT_dS_trans)*np.exp(1j*np.angle(s_trans))
        
        deff_dS_reflect_phase = 1j*np.absolute(s_reflect)*deff_dS_reflect_mag
        deff_dS_trans_phase = 1j*np.absolute(s_trans)*deff_dS_trans_mag
        
        # From Baker-Jarvis tech note 1992 can get dT_dd (or just calculate it)
        dT_dd = (-2*np.pi*1j*np.sqrt(mu_eff*ep_eff)*\
                 np.exp((-2*np.pi*L*1j*np.sqrt(mu_eff*ep_eff))/lam_og))/lam_og
        
        deff_dL = deff_dT*dT_dd
        
        # Combine partial derivatives for overal measurement uncertainty
        # Final Type B Uncertainty
        delta_dielec = (1 - ((lam_0**2)/(LAM_C**2)))\
            *np.sqrt((((np.real(deff_dS_reflect_mag))*err_s_r[0])**2) + \
            (((np.real(deff_dS_reflect_phase))*np.radians(err_s_r[1]))**2) + \
            (((np.real(deff_dS_trans_mag))*err_s_t[0])**2) + \
            (((np.real(deff_dS_trans_phase))*np.radians(err_s_t[1]))**2) + \
            (((np.real(deff_dL))*delta_length)**2))
        
        delta_lossfac = (1 - ((lam_0**2)/(LAM_C**2)))\
            *np.sqrt((((np.imag(deff_dS_reflect_mag))*err_s_r[0])**2) + \
            (((np.imag(deff_dS_reflect_phase))*np.radians(err_s_r[1]))**2) + \
            (((np.imag(deff_dS_trans_mag))*err_s_t[0])**2) + \
            (((np.imag(deff_dS_trans_phase))*np.radians(err_s_t[1]))**2) + \
            (((np.imag(deff_dL))*delta_length)**2))
        
        return delta_dielec, delta_lossfac
    
    def _calc_uncertainties_NRW(self,s_param,nrows,sign,x,s_reflect,\
                                         s_trans,gam,lam_og,new_t,complex_mu,\
                                         ep_r,lam_0,dielec,lossfac,losstan,\
                                         corr):
        """
        Calculates nrw uncertainties from Baker-Jarvis et al. (1993). 
        Transmission/reflection and short-circuit line methods for measuring 
        permittivity and permeability. NIST Technical Note 1355-R. Boulder, 
        CO: National Institute of Standards and Technology. 
        http://doi.org/10.6028/NIST.TN.1355r
        
        Returns
        -------
        delta_dielec : array 
            Uncertainty in real part.
        
        delta_lossfac : array 
            Uncertainty in imaginary part.
            
        delta_mu : complex array
            Uncertainty in complex permeability.
        """
        # Check if using corrected S-params
        if corr:
            s11 = self.corr_s11
            s21 = self.corr_s21
            s12 = self.corr_s12
            s22 = self.corr_s22
            L = self.Lcorr
        else:
            s11 = self.s11
            s21 = self.s21
            s12 = self.s12
            s22 = self.s22
            L = self.L
        
        # Strip uarrays of their nominal values
        if s_param == 'f':
            err_s_r = unp.std_devs(s11)
            err_s_t = unp.std_devs(s21)
        else:
            err_s_r = unp.std_devs(s22)
            err_s_t = unp.std_devs(s12)
        
        delta_length = 0.001 #in cm
        e_0 = E_0*100
        mu_0 = U_0*100
        mu = complex_mu
        eps = ep_r
        t = new_t
        omega = 2*np.pi*self.freq
        small_gam = (1j*2*np.pi/lam_0)*np.sqrt(eps*mu - \
                    (lam_0/LAM_C)**2)
        small_gam_0 = (1j*2*np.pi/lam_0)*np.sqrt(1- (lam_0/LAM_C)**2)
        
        # (2.95)
        dT_dmu = ((L*eps*(omega**2)*e_0*mu_0)/(2*small_gam))*np.exp(-small_gam*L)
        # (2.94)
        dT_deps = ((L*mu*(omega**2)*e_0*mu_0)/(2*small_gam))*np.exp(-small_gam*L)
        # (2.93)
        dT_dL = -small_gam*np.exp(-small_gam*L)
        # (2.91)
        dgam_deps = (small_gam_0*(mu**2)*e_0*mu_0*(omega**2))\
                        /(small_gam*(small_gam + small_gam_0*mu)**2)
        # (2.92)
        dgam_dmu = eps*dgam_deps/mu + 2*small_gam_0*small_gam/(small_gam + \
                        small_gam_0*mu)**2
        # (2.90)
        dS_trans_dgam = 2*t*gam*(t**2 -1)/(1 - (gam**2)*(t**2))**2
        # (2.89)
        dS_trans_dT = (1 - gam**2)*((t**2)*(gam**2) + 1)/(1 - (gam**2)*(t**2))**2
        # (2.88)
        dS_reflect_dT = 2*gam*t*((gam**2) - 1)/(1 - (gam**2)*(t**2))**2
        # (2.87)
        dS_reflect_dgam = -(1 + (gam**2)*(t**2))*((t**1) - 1)/(1 - \
                          (gam**2)*(t**2))**2
        # (2.74)
        a = dS_reflect_dT*dT_deps + dS_reflect_dgam*dgam_deps
        b = dS_reflect_dT*dT_dmu + dS_reflect_dgam*dgam_dmu
        # (2.75)
        c = dS_trans_dT*dT_deps + dS_trans_dgam*dgam_deps
        d = dS_trans_dT*dT_dmu + dS_trans_dgam*dgam_dmu
        # (2.77)
        e = -dS_reflect_dT*dT_dL
        # (2.78)
        f = -dS_trans_dT*dT_dL
        # (2.79)
        deps_dS_reflect_mag = np.exp(1j*np.angle(s_reflect))/(a - b*c/d)
        # (2.80)
        dmu_dS_reflect_mag = -c*deps_dS_reflect_mag/d
        # trans
        deps_dS_trans_mag = np.exp(1j*np.angle(s_trans))/(c - d*a/b)
        dmu_dS_trans_mag = -a*deps_dS_trans_mag/b
        # (2.81)
        deps_dL = (b*f - d*e)/(b*c - a*d)
        # (2.82)
        dmu_dL = (e - a*deps_dL)/b
        # (2.83)
        deps_dS_reflect_phase = 1j*np.absolute(s_reflect)*deps_dS_reflect_mag
        # (2.84)
        dmu_dS_reflect_phase = 1j*np.absolute(s_reflect)*dmu_dS_reflect_mag
        # (2.85)
        deps_dS_trans_phase = 1j*np.absolute(s_trans)*deps_dS_trans_mag
        # (2.86)
        dmu_dS_trans_phase = 1j*np.absolute(s_trans)*dmu_dS_reflect_mag
        # (2.67)
        delta_dielec = np.sqrt((((np.real(deps_dS_reflect_mag))*err_s_r[0])**2) + \
            (((np.real(deps_dS_reflect_phase))*np.radians(err_s_r[1]))**2) + \
            (((np.real(deps_dS_trans_mag))*err_s_t[0])**2) + \
            (((np.real(deps_dS_trans_phase))*np.radians(err_s_t[1]))**2) + \
            (((np.real(deps_dL))*delta_length)**2))
        # (2.68)
        delta_lossfac = np.sqrt((((np.imag(deps_dS_reflect_mag))*err_s_r[0])**2) + \
            (((np.imag(deps_dS_reflect_phase))*np.radians(err_s_r[1]))**2) + \
            (((np.imag(deps_dS_trans_mag))*err_s_t[0])**2) + \
            (((np.imag(deps_dS_trans_phase))*np.radians(err_s_t[1]))**2) + \
            (((np.imag(deps_dL))*delta_length)**2))
        # mu
        delta_mureal = np.sqrt((((np.real(dmu_dS_reflect_mag))*err_s_r[0])**2) + \
            (((np.real(dmu_dS_reflect_phase))*np.radians(err_s_r[1]))**2) + \
            (((np.real(dmu_dS_trans_mag))*err_s_t[0])**2) + \
            (((np.real(dmu_dS_trans_phase))*np.radians(err_s_t[1]))**2) + \
            (((np.real(dmu_dL))*delta_length)**2))
        delta_muimag = np.sqrt((((np.imag(dmu_dS_reflect_mag))*err_s_r[0])**2) + \
            (((np.imag(dmu_dS_reflect_phase))*np.radians(err_s_r[1]))**2) + \
            (((np.imag(dmu_dS_trans_mag))*err_s_t[0])**2) + \
            (((np.imag(dmu_dS_trans_phase))*np.radians(err_s_t[1]))**2) + \
            (((np.imag(dmu_dL))*delta_length)**2))
        
        return delta_dielec, delta_lossfac, delta_mureal, delta_muimag
    
    def _calcTotalUncertainty(self,dielec,lossfac,losstan):
        # Type A Uncertainty. Note 2 regions for lossfac and losstan
        dielec_TypeA = np.where(self.real_part_diff_array > 0.008,\
                                self.real_part_diff_array,0.008)
        lossfac_TypeA_high = np.where(\
                self.imag_part_diff_array[np.where(self.freq<10**8)] > 0.009, \
                self.imag_part_diff_array[np.where(self.freq<10**8)], 0.009)
        lossfac_TypeA_low = np.where(\
                self.imag_part_diff_array[np.where(self.freq>=10**8)] > 0.002, \
                self.imag_part_diff_array[np.where(self.freq>=10**8)], 0.002)
        lossfac_TypeA = np.concatenate((lossfac_TypeA_high,lossfac_TypeA_low))
        losstan_TypeA_high = np.where(\
                self.tand_part_diff_array[np.where(self.freq<10**8)] > 0.009, \
                self.tand_part_diff_array[np.where(self.freq<10**8)], 0.009)
        losstan_TypeA_low = np.where(\
                self.tand_part_diff_array[np.where(self.freq>=10**8)] > 0.002, \
                self.tand_part_diff_array[np.where(self.freq>=10**8)], 0.002)
        losstan_TypeA = np.concatenate((losstan_TypeA_high,losstan_TypeA_low))
        
        # Type B Uncertainty
        delta_dielec = unp.std_devs(dielec)
        delta_lossfac = unp.std_devs(lossfac)
        delta_losstan = unp.std_devs(losstan)
        
        # Combined Uncertainty
        unc_dielec = np.sqrt(delta_dielec**2 + dielec_TypeA**2)
        unc_lossfac = np.sqrt(delta_lossfac**2 + lossfac_TypeA**2)
        unc_losstan = np.sqrt(delta_losstan**2 + losstan_TypeA**2)
        
        # Final uArrays
        dielec = unp.uarray(unp.nominal_values(dielec),unc_dielec)
        lossfac = unp.uarray(unp.nominal_values(lossfac),unc_lossfac)
        losstan = unp.uarray(unp.nominal_values(losstan),unc_losstan)
        
        return dielec, lossfac, losstan
    
    def _de_embed(self):
        """
        Perform S-parameter de-embedding to remove influence of washers on 
        measurement.
        """
        # Get washer length
        if isinstance(self.corr, (list)):
            L_washer1 = self.corr[0]
            L_washer2 = self.corr[1]
        elif isinstance(self.corr, (float)):
            L_washer = self.corr
        else:
            L_washer = 0.17
        # Create synthetic teflon washer s-parameters
        epsilon = -1j*0.0007;
        epsilon += 2.1
        mu = 1 - 1j*0
        lam_0 = (C/self.freq) # (m)
        small_gam = (1j*2*np.pi/lam_0)*np.sqrt(epsilon*mu - \
                        (lam_0/LAM_C)**2)
        small_gam_0 = (1j*2*np.pi/lam_0)*np.sqrt(1- (lam_0/LAM_C)**2)
        big_gam = (small_gam_0*mu - small_gam) / (small_gam_0*mu + \
              small_gam)
        if not isinstance(self.corr, (list)):
            L = L_washer/100 # (m)
            t = np.exp(-small_gam*L)
            sw11_complex_1 = (big_gam*(1-t**2))/(1-(big_gam**2)*(t**2))
            sw21_complex_1 = t*(1-big_gam**2) / (1-(big_gam**2)*(t**2))
            # Symetrical washers
            sw22_complex_1 = sw11_complex_1
            sw12_complex_1 = sw21_complex_1
        else:
            L_1 = L_washer1/100 # (m)
            t_1 = np.exp(-small_gam*L_1)
            sw11_complex_1 = (big_gam*(1-t_1**2))/(1-(big_gam**2)*(t_1**2))
            sw21_complex_1 = t_1*(1-big_gam**2) / (1-(big_gam**2)*(t_1**2))
            sw22_complex_1 = sw11_complex_1
            sw12_complex_1 = sw21_complex_1
            L_2 = L_washer2/100 # (m)
            t_2 = np.exp(-small_gam*L_2)
            sw11_complex_2 = (big_gam*(1-t_2**2))/(1-(big_gam**2)*(t_2**2))
            sw21_complex_2 = t_2*(1-big_gam**2) / (1-(big_gam**2)*(t_2**2))
            sw22_complex_2 = sw11_complex_2
            sw12_complex_2 = sw21_complex_2

        # Split measured sparams into mag and phase since uncertainties package 
        #   does not support complex numbers
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
                            
        ## De-embed
        # Convert to T-parameters
        # Washers
        tw11_left = -(sw11_complex_1*sw22_complex_1-sw12_complex_1*\
                          sw21_complex_1)/sw21_complex_1
        tw12_left = sw11_complex_1/sw21_complex_1
        tw21_left = -sw22_complex_1/sw21_complex_1
        tw22_left = 1/sw21_complex_1
        if isinstance(self.corr, (list)):
            tw11_right = -(sw11_complex_2*sw22_complex_2-sw12_complex_2*\
                          sw21_complex_2)/sw21_complex_2
            tw12_right = sw11_complex_2/sw21_complex_2
            tw21_right = -sw22_complex_2/sw21_complex_2
            tw22_right = 1/sw21_complex_2
        # Measured
        tm11 = -(sm11_complex*sm22_complex-sm12_complex*\
                 sm21_complex)/sm21_complex
        tm12 = sm11_complex/sm21_complex
        tm21 = -sm22_complex/sm21_complex
        tm22 = 1/sm21_complex
        # Make matrices
        left_matrix = np.dstack([tw11_left,tw12_left,tw21_left,\
                                 tw22_left]).reshape(len(sw11_complex_1),2,2)
        if isinstance(self.corr, (list)):
            right_matrix = np.dstack([tw11_right,tw12_right,tw21_right,\
                                 tw22_right]).reshape(len(sw11_complex_1),2,2)
        else:
            # Washers are symetrical
            right_matrix = left_matrix  
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
        ## Build final corrected arrays
        # Get corrected mag and phase
        corr_s11_mag = np.sqrt(corr_s11_complex.real**2 + \
                                     corr_s11_complex.imag**2)
        corr_s11_phase = corr_s11_phase*(180/np.pi)
        corr_s12_mag = np.sqrt(corr_s12_complex.real**2 + \
                                     corr_s12_complex.imag**2)
        corr_s12_phase = corr_s12_phase*(180/np.pi)
        corr_s21_mag = np.sqrt(corr_s21_complex.real**2 + \
                                     corr_s21_complex.imag**2)
        corr_s21_phase = corr_s21_phase*(180/np.pi)
        corr_s22_mag = np.sqrt(corr_s22_complex.real**2 + \
                                     corr_s22_complex.imag**2)
        corr_s22_phase = corr_s22_phase*(180/np.pi)
        # Get original s-param uncertainties
        # NOTE: This is an underestimate of the real uncertainties since they \
        #   are not properly propagated in de-embedding but they are \
        #   probably still much smaller at high freqencies than the stardard \
        #   deviation of the measurement. To be fixed when the uncertainties \
        #   package can handle complex numbers.
        err_s11_mag = unp.std_devs(self.s11[0])
        err_s11_phase = unp.std_devs(self.s11[1])
        err_s12_mag = unp.std_devs(self.s12[0])
        err_s12_phase = unp.std_devs(self.s12[1])
        err_s21_mag = unp.std_devs(self.s21[0])
        err_s21_phase = unp.std_devs(self.s21[1])
        err_s22_mag = unp.std_devs(self.s22[0])
        err_s22_phase = unp.std_devs(self.s22[1])
        # Make final arrays with uncertainties
        corr_s11 = unp.uarray([corr_s11_mag,corr_s11_phase],\
                              [err_s11_mag,err_s11_phase])
        corr_s12 = unp.uarray([corr_s12_mag,corr_s12_phase],\
                              [err_s12_mag,err_s12_phase])
        corr_s21 = unp.uarray([corr_s21_mag,corr_s21_phase],\
                              [err_s21_mag,err_s21_phase])
        corr_s22 = unp.uarray([corr_s22_mag,corr_s22_phase],\
                              [err_s22_mag,err_s22_phase])
            
        return corr_s11, corr_s21, corr_s12, corr_s22
    
    def _boundary_correct(self):
        """
        Correct calculated sprams for boundary effects in the airline after 
        Hickson et al., 2017. Requires the effective solid permittivity 
        of the material, the average particle size in the airline, and 
        the average particle (solid) density to be supplied to the class
        instance. Uses the Looyenga mixing model to calculate the 
        permittivity in the boundary region.
            
        As of now, produces dubious results for the imaginary part.
        """
        beta = 1.835    # Porosity proportinality constant
        # Use washer corrected data if it exists
        if self.corr:
            measured_dielec = unp.nominal_values(self.corr_avg_dielec)
            measured_lossfac = unp.nominal_values(self.corr_avg_lossfac)
            L = self.Lcorr
        else:
            measured_dielec = unp.nominal_values(self.avg_dielec)
            measured_lossfac = unp.nominal_values(self.avg_lossfac)
            L = self.L
        
        # Determine boundary region porosity
        total_volume = np.pi*((self.airline_dimensions['D4']/2)**2)*L - \
            np.pi*((self.airline_dimensions['D1']/2)**2)*L
        sample_volume = np.pi*((self.airline_dimensions['D3']/2)**2)*L - \
            np.pi*((self.airline_dimensions['D2']/2)**2)*L
        boundary_volume = total_volume - sample_volume
        total_porosity = 1 - (self.bulk_density/self.particle_density)
        boundary_porosity = (beta * total_porosity * total_volume) / \
            (sample_volume + beta*boundary_volume)
            
        # Calculate boundary region permittivity using Looyenga mixing model
        if self.solid_losstan:  #Cast to complex if appropriate
            solid_lossfac = self.solid_dielec*self.solid_losstan
            solid_permittivity = 1j*(self.solid_dielec*np.sin(solid_lossfac))
            solid_permittivity += self.solid_dielec*np.cos(solid_lossfac)
        else:
            solid_permittivity = 1j*0
            solid_permittivity += self.solid_dielec
        # Looyenga eqn. for air (ep_air = 1)
        epsilon_thrid = (1 - boundary_porosity)*solid_permittivity**(1/3) \
            + boundary_porosity
        boundary_permittivity = epsilon_thrid**3
        
        # Cast measured average permittivity to complex
        measured_permittivity = 1j*(measured_dielec*np.sin(measured_lossfac))
        measured_permittivity += measured_dielec*np.cos(measured_lossfac)
        
        # Calculate model corrected sample permittivity
        sample_permittivity = (boundary_permittivity * measured_permittivity \
            * np.log(self.airline_dimensions['D3']/\
            self.airline_dimensions['D2'])) / (boundary_permittivity * \
            np.log(self.airline_dimensions['D4']/self.airline_dimensions['D1']) \
            - measured_permittivity * (np.log(self.airline_dimensions['D2']/\
            self.airline_dimensions['D1'])+np.log(self.airline_dimensions['D4']\
            /self.airline_dimensions['D3'])))
            
        # Unpack complex
        sample_dielec = sample_permittivity.real
        samples_losstan = sample_permittivity.imag / sample_dielec
        return sample_dielec, samples_losstan
    
    def _air_gap_correction(self, D2, D3):
        """
        Calculates air gap corrected complex permittivity for solid samples.
        
        Follows Baker-Jarvis et al., 1993 and Rhode & Schwarz 
        Application Note RAC0607-0019_1_4E
        
        Arguments
        ---------
        Ds2 : float 
            The inner diameter (cm) of the solid toroid sample to be specified 
            by user
        
        Ds3 : float T
            The outer diameter (cm) of the solid toroid sample to be specified 
            by user
        
        Return
        ---------
        corr_dielec : array 
            Corrected real part of measured permittivity
        
        corr_lossfac : array 
            Corrected imaginary part of measured permittivity
        
        corr_losstan : array 
            Corrected loss tangent 
        """
        if self.corr:
            measured_dielec = unp.nominal_values(self.corr_avg_dielec)
            measured_lossfac = unp.nominal_values(self.corr_avg_lossfac)
        else:
            measured_dielec = unp.nominal_values(self.avg_dielec)
            measured_lossfac = unp.nominal_values(self.avg_lossfac)
        
        # Calculate L1, L2, and L3 terms
        L1 = np.log(D2/self.airline_dimensions['D1']) + np.log(self.airline_dimensions['D4']/D3)
        L2 = np.log(D3/D2)
        L3 = np.log(self.airline_dimensions['D4']/self.airline_dimensions['D1'])
        
        # Calculate corr_dielec, corr_lossfac 
        corr_dielec  = measured_dielec*(L2/(L3 - measured_dielec*L1))
        corr_lossfac = (corr_dielec*(measured_lossfac/measured_dielec))*(L3/(L3 - L1 \
                       *measured_dielec*(1 + (measured_lossfac/measured_dielec)**2)))
        corr_losstan = corr_lossfac/corr_dielec
        
        return corr_dielec, corr_lossfac, corr_losstan
        
        
    def draw_plots(self,default_settings=True,publish=False,**kwargs):
        """
        Plots permittivity data using make_plot from permittivity_plot.
        
        If freq_cutoff exists, all frequency points lower than freq_cutoff 
        will not be plotted.
        
        Arguments
        ---------
        default_settings : bool 
            Default: True. If True, plots average real part of the permittivity, 
            imaginary part of the permittivity and loss tangent. If corrected
            or normalized data exist, it will be used in the plot. If False 
            prompts user to determine wether to plot, Average, Forward, 
            Reverse, or All results.
        
        publish : bool 
            Default: False. If True, save figures.
            
        Additional Keyword Arguments
        ----------------------------
        These additional arguments can be supplied if default_settings is False.
        
            If default_settings is False and neither of these are provided, 
            uncorrected and unnormalized data will be ploted! 
            
        corr : bool, optional
            Can be used with any of the plot types.
            If True, use corrected sparam data, otheriwse use uncorrected data.
            
        normalized : bool, optional
            Can only be used with Average plots.
            If True, use normalized permittivity data. Supersedes corr if True.
        """
        # If default_settings
        #   Use first of normalized, corrected, or regular data
        if default_settings and self.normalize_density:
            plot_dielec = self.norm_dielec
            plot_lossfac = self.norm_lossfac
            plot_losstan = self.norm_losstan
        elif default_settings and self.corr:
            plot_dielec = self.corr_avg_dielec
            plot_lossfac = self.corr_avg_lossfac
            plot_losstan = self.corr_avg_losstan
            if self.nrw:
                plot_mureal = self.corr_avg_mu_real
                plot_muimag = self.corr_avg_mu_imag
        elif default_settings:
            plot_dielec = self.avg_dielec
            plot_lossfac = self.avg_lossfac
            plot_losstan = self.avg_losstan
            if self.nrw:
                plot_mureal = self.avg_mu_real
                plot_muimag = self.avg_mu_imag
        
        x = self.freq
        # Figure out what to plot
        plot_kwargs = {}
        if default_settings:
            y1 = plot_dielec
            y2 = plot_lossfac
            y3 = plot_losstan
            plot_kwargs = {"legend_label":[self.name]}
            if self.nrw:
                y4 = plot_mureal
                y5 = plot_muimag
        else:
            # Prompt user
            s_plot = input('Please designate "a" for Average, ' + \
                           '"f" for Forward (S11,S21), "r" for ' + \
                           'Reverse (S22,S12), or "all" for All three: ')
            if s_plot not in ('f','r','all','a'):
                raise Exception('Wrong input')
            if s_plot == 'a' and 'normalized' in kwargs:
                y1 = self.norm_dielec
                y2 = self.norm_lossfac
                y3 = self.norm_losstan
            elif s_plot == 'a' and 'corr' in kwargs:
                y1 = self.corr_avg_dielec
                y2 = self.corr_avg_lossfac
                y3 = self.corr_avg_losstan
                if self.nrw:
                    y4 = self.corr_avg_mu_real
                    y5 = self.corr_avg_mu_imag
            elif s_plot == 'a' :
                y1 = self.avg_dielec
                y2 = self.avg_lossfac
                y3 = self.avg_losstan
                if self.nrw:
                    y4 = self.avg_mu_real
                    y5 = self.avg_mu_imag
            elif 'normalized' in kwargs:
                raise Exception('normalized=True can only be used with Average "a" plots')
            elif s_plot == 'f':
                if 'corr' in kwargs:
                    y1 = self.corr_forward_dielec
                    y2 = self.corr_forward_lossfac
                    y3 = self.corr_forward_losstan
                    if self.nrw:
                        y4 = self.corr_forward_mu_real
                        y5 = self.corr_forward_mu_imag
                else:
                    y1 = self.forward_dielec
                    y2 = self.forward_lossfac
                    y3 = self.forward_losstan
                    if self.nrw:
                        y4 = self.forward_mu_real
                        y5 = self.forward_mu_imag
            elif s_plot == 'r':
                if 'corr' in kwargs:
                    y1 = self.corr_reverse_dielec
                    y2 = self.corr_reverse_lossfac
                    y3 = self.corr_reverse_losstan
                    if self.nrw:
                        y4 = self.corr_reverse_mu_real
                        y5 = self.corr_reverse_mu_imag
                else:
                    y1 = self.reverse_dielec
                    y2 = self.reverse_lossfac
                    y3 = self.reverse_losstan
                    if self.nrw:
                        y4 = self.reverse_mu_real
                        y5 = self.reverse_mu_imag
            elif s_plot == 'all' and 'corr' in kwargs:
                x = [x,x,x]
                y1 = [self.corr_forward_dielec,self.corr_reverse_dielec,self.corr_avg_dielec]
                y2 = [self.corr_forward_lossfac,self.corr_reverse_lossfac,self.corr_avg_lossfac]
                y3 = [self.corr_forward_losstan,self.corr_reverse_losstan,self.corr_avg_losstan]
                if self.nrw:
                    y4 = [self.corr_forward_mu_real,self.corr_reverse_mu_real,self.corr_avg_mu_real]
                    y5 = [self.corr_forward_mu_imag,self.corr_reverse_mu_imag,self.corr_avg_mu_imag]
                plot_kwargs = {"legend_label":['Forward [S11,S21]','Reverse [S22,S12]','Average']}
            else:
                x = [x,x,x]
                y1 = [self.forward_dielec,self.reverse_dielec,self.avg_dielec]
                y2 = [self.forward_lossfac,self.reverse_lossfac,self.avg_lossfac]
                y3 = [self.forward_losstan,self.reverse_losstan,self.avg_losstan]
                if self.nrw:
                    y4 = [self.forward_mu_real,self.reverse_mu_real,self.avg_mu_real]
                    y5 = [self.forward_mu_imag,self.reverse_mu_imag,self.avg_mu_imag]
                plot_kwargs = {"legend_label":['Forward [S11,S21]','Reverse [S22,S12]','Average']}
        # Pass publish arguments        
        if publish:
            plot_kwargs['publish'] = True
            plot_kwargs['name'] = self.name
        # Pass freq_cutoff
        if self.freq_cutoff:
            plot_kwargs['freq_cutoff'] = self.freq_cutoff
        # Make plots
        pplot.make_plot(x,y1,'d',**plot_kwargs)
        pplot.make_plot(x,y2,'lf',**plot_kwargs)
        pplot.make_plot(x,y3,'lt',**plot_kwargs)
        if self.nrw:
            pplot.make_plot(x,y4,'ur',**plot_kwargs)
            pplot.make_plot(x,y5,'ui',**plot_kwargs)
        
    def s_param_plot(self,corr=False):
        """
        Plot raw S-Parameter data using make_sparam_plot from
        permittivity_plot
        """
        if corr:
            pplot.make_sparam_plot(self.freq,[self.s11,self.corr_s11],\
                                   [self.s22,self.corr_s22],[self.s21,\
                                   self.corr_s21],[self.s12,self.corr_s12],\
                                   label=['Uncorrected','Corrected'])
        else:
            if self.shorted:
                pplot.make_sparam_plot(self.freq,self.s11,self.s22,self.s21,self.s12,shorted=True,s11_short=self.s11_short)
            else:
                pplot.make_sparam_plot(self.freq,self.s11,self.s22,self.s21,self.s12)
            
    def difference_plot(self):
        """
        Plot the absolute difference between both ε′ and ε′′ calculated from 
        forward (S11,S21) and reverse (S22,S12) S-Paramaters.
        """
        import matplotlib.pyplot as plt
        
        if self.freq_cutoff:
            freq = self.freq[self.freq>=self.freq_cutoff]
            real_diff = self.real_part_diff_array[self.freq>=self.freq_cutoff]
            imag_diff = self.imag_part_diff_array[self.freq>=self.freq_cutoff]
        else:
            freq = self.freq
            real_diff = self.real_part_diff_array
            imag_diff = self.imag_part_diff_array
        
        f,ax = plt.subplots(2, sharex=True, figsize=(18, 15))
        ax[0].loglog(freq,real_diff,'ko-',label='|$\epsilon^\prime[S_{11},S_{21}]$ - $\epsilon^\prime[S_{22},S_{12}]$|')
        ax[0].set_title('Difference in $\epsilon^\prime$', fontsize=40)
        ax[0].legend(fontsize=30,loc=2)
        ax[0].set_ylim(10e-7, 1)
        ax[0].tick_params(axis='both', which='major', width=2, labelsize=30)
        ax[0].tick_params(axis='both', which='minor', width=1.5)
        ax[1].loglog(freq,imag_diff,'ko-',label='|$\epsilon^{\prime\prime}[S_{11},S_{21}]$ - $\epsilon^{\prime\prime}[S_{22},S_{12}]$|')
        ax[1].set_title('Difference in $\epsilon^{\prime\prime}$', fontsize=40)
        ax[1].legend(fontsize=30,loc=2)
        ax[1].set_xlabel('Frequency',fontsize=40)
        ax[1].set_ylim(10e-7, 1)
        ax[1].tick_params(axis='both', which='major', width=2, labelsize=30)
        ax[1].tick_params(axis='both', which='minor', width=1.5)
        if self.freq_cutoff:
            from matplotlib.ticker import FixedLocator, LogLocator, EngFormatter, NullFormatter
            x_logmin = np.log10(np.min(freq)) # log of min and max x values
            x_logmax = np.log10(np.max(freq))
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
            majorLocator = FixedLocator(x_ticklocs)
            majorFormatter = EngFormatter(unit='Hz') # Format major ticks with units
            minorLocator = LogLocator(subs='all') # Use all interger multiples of the log base for minor ticks 
            minorFormatter =  NullFormatter() # No minor tick labels 
            for n in range(0,2):
                # Apply x ticks
                ax[n].get_xaxis().set_major_locator(majorLocator)
                ax[n].get_xaxis().set_major_formatter(majorFormatter)
                ax[n].get_xaxis().set_minor_locator(minorLocator)
                ax[n].get_xaxis().set_minor_formatter(minorFormatter)
                # Format the actual tick marks
                ax[n].tick_params(which='both', width=1, labelsize=30)
                ax[n].tick_params(which='major', length=7)
                ax[n].tick_params(which='minor', length=4)
                # Use smaller line width for minor tick grid lines
                ax[n].grid(b=True, which='major', color='w', linewidth=1.0)
                ax[n].grid(b=True, which='minor', color='w', linewidth=0.5)
        plt.show()
            
#%% Functions
            
def multiple_meas(file_path=None,airline_name=None):
    """
    Generate an instance of AirlineData for every file in a directory. Store 
    the intances in a list, and plot them all using perm_compare.
        
    Arguments
    ---------
    file_path : str, optional 
        Full path of any file in the source directory. Will produce file dialog 
        box if not provided.
    
    airlne : str, optional
        Name of airline used. Every measurement must have been made 
        in the same airline. Will prompt if not provided.
        
    Return
    ------
    class_list : lst 
        List of generated class instances of AirlineData.
    """
    # Use _get_file to get the filepath and airline name if not provided
    if file_path:   # If file path provided as argument
        file = file_path
    elif not file_path:   # If file path not provided
        print('\n')
        print("Select any data file in the source folder. All .txt "+\
              "files in the source folder must be METAS data tables.")
        # Get the file path and the airline name
        airline_name, file, L_in = _get_file(airline_name,file_path)
    elif not airline_name:   # If file path is given but airline is not
        airline_name, file, L_in = _get_file(airline_name,file_path)
        
    # Get directory path    
    directory = os.path.dirname(file)
    # Use a list to maintain order for plotting
    class_list = []
    # Iterate through all .txt files in the directory and run AirlineData
    for file in os.listdir(directory):
        if file.endswith(".txt"):
            filename = os.path.splitext(file)[0]    # Use file name as plot label
            # Append each new instance to class list
            class_list.append(AirlineData(*get_METAS_data(airline_name,\
                                os.path.join(directory, file)),name=filename))
    
    # Plot all files        
    perm_compare(class_list)
    
    return class_list 

def run_default(airline_name='VAL',file_path=None,**kwargs):
    """
    Run AirlineData on get_METAS_data with all the prompts and return the 
    instance.
    """
    return AirlineData(*get_METAS_data(airline_name,file_path),**kwargs)

def run_example():
    rexolite_path = os.path.join(DATAPATH, 'rexolite_PAL.txt')
    serpentine_path = os.path.join(DATAPATH, 'serpentine_dry.txt')
    rexolite_example = AirlineData(*get_METAS_data(airline='PAL',\
        file_path=rexolite_path),name='Rexolite')
    serpentine_example = AirlineData(*get_METAS_data(airline='VAL',\
        file_path=serpentine_path),bulk_density=1.6,\
        name='Serpentine',normalize_density=False,norm_eqn='LI')
    return rexolite_example, serpentine_example

#%% MAIN
def main():
    ## Run Examples
    global rexolite_example
    global serpentine_example
    rexolite_example, serpentine_example = run_example()
    
if __name__ == '__main__':
    main()  # Comment to supress example
#    pass    # Uncomment to supress example
