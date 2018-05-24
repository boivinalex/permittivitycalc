#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Contains the AirlineData class object"""

#File input
import pickle
# Array math
import numpy as np
import uncertainties 
#Citation: Uncertainties: a Python package for calculations with uncertainties,
#    Eric O. LEBIGOT, http://pythonhosted.org/uncertainties/
from uncertainties import unumpy as unp
# Plotting
import permittivitycalc.permittivity_plot as pplot
# Make relative path
import os
# Helper functions
from permittivitycalc.helper_functions import get_METAS_data, _get_file, perm_compare

# GLOBAL VARIABLES
E_0 = 8.854187817620*10**-12 #Permittivity of vacuum (F/m) 
U_0 = 4*np.pi*10**-7 #Permeability of vacuum (V*s/A*m) 
C = (1)*1/np.sqrt(E_0*U_0) #Speed of light (m/s)
LAM_C = float('inf') #Cut-off wavelength = infinity
PACKAGEPATH = os.path.dirname(os.path.abspath(__file__))
DATAPATH = os.path.join(PACKAGEPATH, 'data')


class AirlineData:
    """
    S-parameter data from METAS text file output
    
    Attributes:
    ----------
    L (float): Length of airline in cm.
    
    airline_name (str): Name of airline used for measurement.
    
    file (str): Input file path.
    
    corr (bool): (Optional) If True, also correct S-parameter data and \
        produce corr_* arrays. Default = True.
    
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
        
    res_freq (array): Resonant frequencies in the sample.
        
    name (str): (Optional) Name of measurement to be used in plots.
    
    bulk_density (float): (Optional) Bulk density of material. Nessesary for \
        bulk density normalization.
        
    normalize_density (bool): (Optional) If True, use either Lichtenecker or \
        Landau-Lifshitz-Looyenga equation to normalize the real part of the \
        permittivity to a constant density of 1.60 g/cm^3. Default: False
        
    norm_eqn (str): For use with normalize_density = True. Equation to be \
        used for normalization. Options are 'LI' (default) for the \
        Lichtenecker equation and 'LLL' for the Landau-Lifshitz-Looyenga \
        equation. LI used alpha = 1.92 (Olhoeft, 1985) and LLL uses \
        alpha = 0.307 (Hickson et al., 2017, Lunar samples).
        
    temperature (str or float): (Optional) Temperature of measurement.
    
    date (str): (Optional) Measurement date.
    
    solid_dielec (float): (Optional) The solid dielectric constant \
        of the material.
        
    solid_losstan (float): (Optional) The solid loss tangent of the material.
        
    particle_diameter (float): (Optional) The average particle diameter in \
        the airline in cm.
        
    particle_density (float): (Optional) The average (solid) particle density \
        of the material in g/cm^3.
        
    airline_dimensions (dict): dimensions of the airline in cm. D1 is the \
        diameter of the inner conductor and D4 is the diameter of the outer \
        conductor. D2 and D3 bound the sample-airline boundary regions if \
        particle_diameter is provided. airline_dimensions is generated \
        automatically for airlines VAL, PAL, and GAL. Empty otherwise.
        
    bcorr (complex array): Avg complex permittivity corrected for boundary \
        effects. Computed automatically if solid_dielec, particle_diameter, \
        particle_density, and bulk_density are present. solid_losstan is \
        optional.
    
    nrw (bool): If True, use Nicholson, Rross, Weir (NRW) algorithm to \
        calculate permittivity and magnetic permeability.
        
    shorted (bool): If True, automatically load Shorted S11 data. File name \
        must have the following format and be in the same folder as original \
        file: file_path/file_name_shorted.txt
        
        Example:
            file_path/air_atm.txt
            file_path/air_atm_shorted.txt
    """
    def __init__(self,L,airline,dataArray,file,corr=False,bulk_density=None,\
                 temperature=None,name=None,date=None,solid_dielec=None,\
                 solid_losstan=None,particle_diameter=None,\
                 particle_density=None,nrw=False,shorted=False,\
                 normalize_density=False,norm_eqn='LI'):
        self.L = L
        self.airline_name = airline
        self.file = file
        self.corr = corr
        self.nrw = nrw
        self.shorted = shorted
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
                raise Exception('Shorted file does not exists.')
        # Calculate permittivity
        if nrw:
            self.avg_dielec, self.avg_lossfac, self.avg_losstan, self.mu = \
                self._permittivity_calc('a')
        else:
            self.avg_dielec, self.avg_lossfac, self.avg_losstan = \
                self._permittivity_calc('a')
            self.forward_dielec, self.forward_lossfac, self.forward_losstan = \
                self._permittivity_calc('f')
            self.reverse_dielec, self.reverse_lossfac, self.reverse_losstan = \
                self._permittivity_calc('r')
        # Try to calculate corrected permittivity if array length is 601 only
        # Also check for NaNs and don't run if any
        if corr and len(self.freq) == 601:
            try:
                if not np.isnan(unp.nominal_values(self.avg_dielec)).any():
                    self.Lcorr = self.L - 0.345 # Remove total washer length 
                    self.corr_s11, self.corr_s21, self.corr_s12, \
                        self.corr_s22 = self._de_embed()
                    if nrw:
                        self.corr_avg_dielec, self.corr_avg_lossfac, \
                            self.corr_avg_losstan, self.corr_avg_mu = \
                            self._permittivity_calc('a',True)
                    else:
                        self.corr_avg_dielec, self.corr_avg_lossfac, \
                            self.corr_avg_losstan = \
                            self._permittivity_calc('a',True)
            except:
                pass
        self.res_freq = self._resonant_freq()
        self.freq_avg_dielec, self.freq_avg_losstan, self.freq_avg_dielec_std, \
            self.freq_avg_losstan_std = self._freq_avg()
            
        # Optional attributes
        self.name = name
        self.bulk_density = bulk_density
        self.normalize_density = normalize_density
        self.norm_eqn = norm_eqn
        self.temperature = temperature
        self.date = date
        self.solid_dielec = solid_dielec
        self.solid_losstan = solid_losstan
        self.particle_diameter = particle_diameter
        self.particle_density = particle_density
        self.airline_dimensions = self._dims()
        # If appropriate data provided, correct for boundary effects
        if (solid_dielec and particle_diameter and particle_density and \
            bulk_density):
            self.bcorr_dielec, self.bcorr_losstan = self._boundary_correct()
        # If normalize_density is True and bulk_density exists, do it
        #   Normalize the density to 1.60 g/cm^3
        if normalize_density and bulk_density:
            complex_dielec = 1j*unp.nominal_values(self.avg_lossfac);
            complex_dielec += unp.nominal_values(self.avg_dielec)
            if self.norm_eqn == 'LI':
                # Lichtenecker equation
                #   alpha from Olhoeft, 1985
                norm_complex_dielec = complex_dielec*(1.92)**\
                    (1.60-self.bulk_density)
            elif self.norm_eqn == 'LLL':
                # Landau-Lifshitz-Looyenga equation
                #   alpha from Hickson et al., 2018
                norm_complex_dielec = complex_dielec*((1.60*0.307 + 1)**3 / \
                                            (self.bulk_density*0.307 + 1)**3)
            self.norm_dielec = np.real(norm_complex_dielec)
            self.norm_lossfac = np.imag(norm_complex_dielec)
            self.norm_losstan = self.norm_lossfac/self.norm_dielec
        elif normalize_density:
            raise Exception('Need bulk desnity to normalize to constant density')
            
    def __repr__(self):
        rep = 'AirlineData(*get_METAS_data(airline=%r,file_path=%r),' % \
                (self.airline_name,self.file) + \
                'bulk_density=%r,temperature=%r,name=%r,date=%r,corr=%r' % \
                (self.bulk_density,self.temperature,self.name,self.date,\
                 self.corr) + ',solid_dielec=%r,solid_losstan=%r' % \
                 (self.solid_dielec,self.solid_losstan) + \
                 ',particle_diameter=%r,particle_density=%r' % \
                 (self.particle_diameter, self.particle_density) + \
                 ',nrw=%r,normalize_density=%r,norm_eqn=%r)' % \
                 (self.nrw, self.normalize_density, self.norm_eqn)
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
        if self.normalize_density:
            srep += '\n' + 'Normalized data at a bulk density of 1.60 g/cm^3 using the'
            if self.norm_eqn == 'LI':
                srep += ' Lichtenecker equation'
            elif self.norm_eqn == 'LLL':
                srep += 'Landau-Lifshitz-Looyenga equation'
            srep += ' is available.'
        return srep
    
    def _resonant_freq(self):
        """
        Calculate and return array of resonant frequencies from complex permittivity \ 
        and/or permeability measurement.
    
        Follows David Stillman PhD Thesis / NIST Technical Note 1355 ??

        """
        
        n = np.arange(1,40) # Max number of resonances (overkill)
        if self.corr:
            measured_dielec = unp.nominal_values(self.corr_avg_dielec)
            measured_lossfac = unp.nominal_values(self.corr_avg_lossfac)
            L = self.Lcorr
        else:
            measured_dielec = unp.nominal_values(self.avg_dielec)
            measured_lossfac = unp.nominal_values(self.avg_lossfac)
            L = self.L
        e_r = np.median(measured_dielec[1::]) # Exclude first data point
        e_i = np.median(measured_lossfac[1::])
        if self.nrw:
            u_r = np.real(unp.nominal_values(self.mu))
            u_r = np.median(u_r[1::])
            u_i = np.imag(unp.nominal_values(self.mu))
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
        Calculate an average dielectric constant and loss tangent across all measured frequencies \ 
        from midpoint frequency values between resonant frequencies as described in: \
            Hickson et al., 2018
        
        Return
        ------
        freq_avg_dielec (uncertainties.core.AffineScalarFunc): Average dielectric \
        constant calculated from midpoint frequency values between resonant \
        frequencies. Skips first two values due to large uncertainty. 
        
        freq_avg_dielec_std (float): Standard deviation of freq_avg_dielec from above.
        
        freq_avg_losstan (uncertainties.core.AffineScalarFunc): Average loss \
        tangent calculated from midpoint frequency values between resonant \
        frequencies. Skips first two values due to large uncertainty. 
        
        freq_avg_losstan_std (float): Standard deviation of freq_avg_losstan from above.
        
        """
        f_r = self.res_freq
        if self.corr:
            dielec = self.corr_avg_dielec
            lossfac = self.corr_avg_lossfac
        else:
            dielec = self.avg_dielec
            lossfac = self.avg_lossfac
        mid_pts  = np.zeros(len(f_r)-1, dtype=int)
        f_0_mids = np.zeros(len(f_r)-1, dtype=int)
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
            dims['D1'] = 6.205
            dims['D4'] = 14.285
        elif self.airline_name == 'GAL':
            dims['D1'] = 6.19
            dims['D4'] = 14.32
        # If particle diameter is known, calculate D2 and D3 for boundary 
        #   effect correction   
        if dims and self.particle_diameter:
            dims['D2'] = dims['D1'] + self.particle_diameter
            dims['D3'] = dims['D4'] - self.particle_diameter
        return dims
    
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
        a_1 = np.sqrt(-(ln_1_T / (2*np.pi*L))**2)
        a_2 = -1 * (np.sqrt(-(ln_1_T / (2*np.pi*L))**2))
         
        # Determine correct root with condition Re(1/A) > 0
        a[a_1.real > 0] = 1 / a_1[a_1.real > 0]
        a[a_1.real < 0] = 1 / a_2[a_1.real < 0]
          
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
        if self.nrw:
            lossfac = -ep_r.imag
            complex_mu = mu_eff
        else:
            lossfac = ep_r.imag
        losstan = lossfac/dielec
        
        # Calculate uncertainties
        # Check for uncertainties and make sure not using NRW
        if isinstance(self.s11[0][0], uncertainties.UFloat) and not self.nrw: 
            delta_dielec, delta_lossfac, delta_losstan = \
                self._calc_uncertainties(s_param,nrows,sign,x,s_reflect,\
                                         s_trans,gam,lam_og,new_t,mu_eff,\
                                         ep_eff,lam_0,dielec,lossfac,losstan,corr)
            dielec = unp.uarray(dielec,delta_dielec)
            lossfac = unp.uarray(lossfac,delta_lossfac)
            losstan = unp.uarray(losstan,delta_losstan)
        
        if self.nrw:
            return dielec,lossfac,losstan,complex_mu
        else:
            return dielec,lossfac,losstan
        
    def _calc_uncertainties(self,s_param,nrows,sign,x,s_reflect,s_trans,gam,\
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
            L = self.Lcorr
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
        dgam_dS_reflect = np.zeros(nrows,dtype=complex)
        dgam_dS_trans = np.zeros(nrows,dtype=complex)
        
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
        
        # Add measurement standard deviation error to measurement uncertainty
        #   in quadrature as calculated with multiple rexolite measurements.
        #   Note 2 regions for lossfac and losstan
        delta_dielec = np.sqrt(delta_dielec**2 + 0.008**2)
        
        delta_lossfac[np.where(self.freq<10**8)] = \
            np.sqrt(delta_lossfac[np.where(self.freq<10**8)]**2 + 0.009**2) 
        delta_lossfac[np.where(self.freq>=10**8)] = \
            np.sqrt(delta_lossfac[np.where(self.freq>=10**8)]**2 + 0.002**2)
        
        delta_losstan[np.where(self.freq<10**8)] = \
            np.sqrt(delta_losstan[np.where(self.freq<10**8)]**2 + 0.009**2)
        delta_losstan[np.where(self.freq>=10**8)] = \
            np.sqrt(delta_losstan[np.where(self.freq>=10**8)]**2 + 0.002**2)
        
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
        # Cast to complex and unwrap phase
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
        Correct calculated sprams for boundary effects in the airline after \
            Hickson et al., 2017. Requires the effective solid permittivity \
            of the material, the average particle size in the airline, and \
            the average particle (solid) density to be supplied to the class\
            instance. Uses the Looyenga mixing model to calculate the \
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
        
        Follows Baker-Jarvis et al., 1993 and Rhode & Schwarz Application Note RAC0607-0019_1_4E
        
        Arguments
        ---------
        Ds2 (float): The inner diameter (cm) of the solid toroid sample to be specified by user
        
        Ds3 (float): The outer diameter (cm) of the solid toroid sample to be specified by user
        
        Return
        ---------
        corr_dielec (array): Corrected real part of measured permittivity
        
        corr_lossfac (array): Corrected imaginary part of measured permittivity
        
        corr_losstan (array): Corrected loss tangent 
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
        
        
    def draw_plots(self,default_setting=True,corr=False,normalized=False,\
                   publish=False):
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
        elif normalized:
            avg_dielec = self.norm_dielec
            avg_lossfac = self.norm_lossfac
            avg_losstan = self.norm_losstan
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
            kwargs = {"legend_label":[self.name]}
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
            
#%% Functions
            
def multiple_meas(file_path=None,airline_name=None):
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
    Run AirlineData on get_METAS_data with all the prompts and return the \
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
