# -*- coding: utf-8 -*-
"""
Created on Mon Jun  4 13:28:19 2018

@author: alex
"""
# Array math
import numpy as np
#Citation: Uncertainties: a Python package for calculations with uncertainties,
#    Eric O. LEBIGOT, http://pythonhosted.org/uncertainties/
from uncertainties import unumpy as unp
# Nonlinear fitting
from lmfit import Minimizer, Parameters, report_fit
# Plotting
import permittivitycalc.permittivity_plot as pplot

# GLOBAL VARIABLES
E_0 = 8.854187817620*10**-12 #Permittivity of vacuum (F/m) 
U_0 = 4*np.pi*10**-7 #Permeability of vacuum (V*s/A*m) 
C = (1)*1/np.sqrt(E_0*U_0) #Speed of light (m/s)
LAM_C = float('inf') #Cut-off wavelength = infinity


class AirlineIter():
    """
    Iterative fit to measured S-parameters using a Cole-Cole model. Can fit 
    both 2-port and 2-port + shorted measurements.
        
    To determine an initial guess for the fit, first fit a Cole-Cole
    model to analytical results (either New Non-iterative, NRW, or both).
    Then use the emcee package with lmfit to refine initial guess before
    performing final fit using lmfit. Then plot and store the results.
        
    Attributes
    ----------
    data_instance : permittivitycalc.AirlineData 
        AirlineData class instance  containing raw S-parameters to be 
        iterated over.
        
    trial_run : bool
        If True only fit Cole-Cole model to analytical results and plot the 
        results. Useful for determining the number of poles to be used
        in the Cole-Cole model before perfoming the more time consuming
        final fit (Default: True).
            
    number_of_poles : int
        Number of poles to be used in the Cole-Cole model for epsilon (Default: 2).
        
    number_of_fits : int
        Number of default emcee fits to perform beofre final fit (Default: 2).
        
    start_freq : float or int
        Start frequency for iteration
        
    fit_mu : bool
        If True, fit both permittivity and permeability (Default: False).
        
    number_of_poles_mu : int
        Number of poles to be used in the Cole-Cole model for mu (Default: 1).
        
    Return
    ------
    epsilon_iter : array
        Complex array containing the results of the iterrative fit for epsilon.
        trail_run must be set to False.
        
    mu_iter : array
        Complex array containing the results of the iterative fit for mu.
        trial_run must be set to False and fit_mu must be set to True.
        
    cc_parmas : dict
        Cole-Cole model parameters for epsilon and mu (if fit_mu is True).
    """
    def __init__(self,data_instance,trial_run=True,number_of_poles=2,\
                 number_of_fits=2,start_freq=None,fit_mu=False,\
                 number_of_poles_mu=1):
        self.meas = data_instance
        self.trial = trial_run
        self.poles = number_of_poles
        self.poles_mu = number_of_poles_mu
        self.fits = number_of_fits
        self.start_freq = start_freq
        self.fit_mu = fit_mu
        # Check if shorted
        if self.meas.shorted:
            self.shorted = True
        else:
            self.shorted = False
        
        if self.trial:
            self._permittivity_iterate()
        elif self.fit_mu:
            self.epsilon_iter, self.mu_iter, cc_params = self._permittivity_iterate()
        else:
            self.epsilon_iter, cc_params = self._permittivity_iterate()
            
    def _colecole(self,number_of_poles,freq,v,mu=False):
        """
        Unpack Cole-Cole paramaterd and retuns Cole-Cole model based on 
        number of poles.
        
        Attributes
        ----------
        number_of_poles : int 
            Number of poles in the Cole-Cole model 
            
        freq : numpy array 
            Frequency vector for model
        
        v : dict 
            Cole-Cole parameters to be used in model
            v must be a dictionary with the following values:
                - k_inf
                - sigma
                - k_dc_n
                - tau_n
                - alpha_n
                
                Where n is the pole number in the Cole-Cole model from 1 to inf.
                
        mu : bool
            If True use mu parameters in v:
                - mu_inf
                - mu_dc_n
                - mutau_n
                - mualpha_n
                
                Where n is the pole number in the Cole-Cole model from 1 to inf.
                
        Return
        ------
        k : array
            Cole-Cole model.
        """
        if mu:
            k = (v['mu_inf'])
            # Make k the length of freq if using 0 poles
            if self.poles_mu == 0:
                freq_vector = np.ones(len(freq))
                k = k * freq_vector
        else:
            k = (v['k_inf']) - 1j*v['sigma']/(2*np.pi*freq*E_0)
        
        for n in range(number_of_poles):
            n+=1    # Start names at 1 intead of 0
            if mu and self.poles_mu != 0:
                k += (v['mu_dc_{}'.format(n)] - v['mu_inf'])/(1 + (1j*2*np.pi*freq*v['mutau_{}'.format(n)])**v['mualpha_{}'.format(n)])
            else:
                k += (v['k_dc_{}'.format(n)] - v['k_inf'])/(1 + (1j*2*np.pi*freq*v['tau_{}'.format(n)])**v['alpha_{}'.format(n)])
        
        return k
    
    def _model_sparams(self,freq,L,epsilon,mu):
        # Calculate predicted sparams
        lam_0 = (C/freq)
        
        small_gam = (1j*2*np.pi/lam_0)*np.sqrt(epsilon*mu - \
                    (lam_0/LAM_C)**2)
        
        small_gam_0 = (1j*2*np.pi/lam_0)*np.sqrt(1- (lam_0/LAM_C)**2)
        
        t = np.exp(-small_gam*L)
        
        big_gam = (small_gam_0*mu - small_gam) / (small_gam_0*mu + \
                  small_gam)
        
        # Use shorted S11 data if present
        if self.shorted:
            # Modified S11
            s11_predicted = big_gam - ( ( (1-big_gam**2)*t**2 ) / (1 - (big_gam*t**2) ) )
        else:
            # Baker-Jarvis S11
            s11_predicted = (big_gam*(1-t**2))/(1-(big_gam**2)*(t**2))
        
        # S21
        s21_predicted = t*(1-big_gam**2) / (1-(big_gam**2)*(t**2))
        
        s12_predicted = s21_predicted
        
        return s11_predicted, s21_predicted, s12_predicted
    
    def _iteration_parameters(self,initial_values=None,mu=False):
        """
        Creates Parameter object to be used in _permittivity_iterate
        
        Attributes
        ----------
        number_of_poles : int 
            Number if poles in the Cole-Cole model.
            
        initial_values : dict (optional)
            Initial guess interation parameters for the Cole-Cole model. If 
                none given, will generate default parameters.
            
            initial_values must be a dictionary with the following values:
                - k_inf
                - sigma
                - k_dc_n
                - tau_n
                - alpha_n
                
                Where n is the pole number in the Cole-Cole model from 1 to inf.
                
            If mu = True then initial_values also must contain:
                - mu_inf
                - mu_dc_n
                - mutau_n
                - mualpha_n
                
        mu : bool
            If True, create seperate parameters for mu (Default: False).
                
        Return
        ------
        params : lmfit.Parameter object
            paramaters for iteration
        """
        # Get default initial values if none given
        if not initial_values:
            initial_values = self._default_initial_values(self.poles)
            if mu:
                initial_values_mu = self._default_initial_values_mu(self.poles_mu)
                initial_values = {**initial_values,**initial_values_mu}
        
        # Create parameters
        params = Parameters()
        if mu:
            params.add('mu_inf',value=initial_values['mu_inf'],min=1)
        params.add('k_inf',value=initial_values['k_inf'],min=1)
        params.add('sigma',value=initial_values['sigma'],min=0)
        
        if mu and self.poles_mu != 0:
            for m in range(self.poles_mu):
                m+=1
                params.add('mu_dc_{}'.format(m),value=initial_values['mu_dc_{}'.format(m)],min=1)
                params.add('mutau_{}'.format(m),value=initial_values['mutau_{}'.format(m)],min=0)
                params.add('mualpha_{}'.format(m),value=initial_values['mualpha_{}'.format(m)],min=0,max=1)   
        for n in range(self.poles):
            n+=1 # start variable names at 1 instead of 0
            params.add('k_dc_{}'.format(n),value=initial_values['k_dc_{}'.format(n)],min=1)
            params.add('tau_{}'.format(n),value=initial_values['tau_{}'.format(n)],min=0)
            params.add('alpha_{}'.format(n),value=initial_values['alpha_{}'.format(n)],min=0,max=1)
            
        return params
    
    def _default_initial_values(self,number_of_poles):
        """
        Creates default initial values for iteration parameters
        """
        initial_values = {'k_inf':2.01,'sigma':0.0001}
        
        for n in range(number_of_poles):
            n+=1 # start variable names at 1 instead of 0
            initial_values['k_dc_{}'.format(n)] = 3 + 10**(n-1)
            initial_values['tau_{}'.format(n)] = 1e-9 * 10**-(2*(n-1))
            initial_values['alpha_{}'.format(n)] = 0.5
            
        return initial_values
    
    def _default_initial_values_mu(self,number_of_poles):
        """
        Creates default initial values for iteration parameters
        """
        initial_values = {'mu_inf':1.01}
        
        for n in range(number_of_poles):
            n+=1 # start variable names at 1 instead of 0
            initial_values['mu_dc_{}'.format(n)] = 0.001 + 10**(n-1)
            initial_values['mutau_{}'.format(n)] = 1e-9 * 10**-(2*(n-1))
            initial_values['mualpha_{}'.format(n)] = 0.5
            
        return initial_values
    
    def _colecole_residuals(self,params,number_of_poles,freq,k,mu=False):
        """
        Cole-Cole model objective function
        """
        v = params.valuesdict()
        
        if mu:
            k_predicted = self._colecole(number_of_poles,freq,v,mu=True)
        else:
            k_predicted = self._colecole(number_of_poles,freq,v)
        
        # Residuals
        resid1 = k_predicted.real - k.real
#        resid2 = (np.abs(k_predicted.imag) - np.abs(k.imag))
        resid2 = k_predicted.imag - k.imag
        
        return np.concatenate((resid1,resid2))
    
    def _iterate_objective_function(self,params,L,freq_0,s11c,s21c,s12c,s22c):
        """
        Objective funtion to minimize from modified Baker-Jarvis (NIST) 
            iterative method (Houtz et al. 2016).
        """
        
        freq = self.meas.freq[self.meas.freq>=freq_0]
        L = L/100 #L in m

        # Unpack parameters
        v = params.valuesdict()
        
        # Calculate predicted mu and epsilon
        if self.fit_mu:    #check if fitting mu 
            mu = self._colecole(self.poles_mu,freq,v,mu=True)
        else:   #set mu=1 if not fitting mu
            mu = 1
        epsilon = self._colecole(self.poles,freq,v)
        
        s11_predicted, s21_predicted, s12_predicted = self._model_sparams(freq,L,epsilon,mu)
        
        # Get uncertainty (weights)
        if self.shorted:
            s11m_sigma = unp.std_devs(self.meas.s11_short[0][self.meas.freq>=freq[0]])
            s11p_sigma = unp.std_devs(self.meas.s11_short[1][self.meas.freq>=freq[0]])
        else:   #NOTE: Update to use S22 for non-shorted case
            s11m_sigma = unp.std_devs(self.meas.s11[0][self.meas.freq>=freq[0]])
            s11p_sigma = unp.std_devs(self.meas.s11[1][self.meas.freq>=freq[0]])
        s21m_sigma = unp.std_devs(self.meas.s21[0][self.meas.freq>=freq[0]])
        s21p_sigma = unp.std_devs(self.meas.s21[1][self.meas.freq>=freq[0]])
        s12m_sigma = unp.std_devs(self.meas.s12[0][self.meas.freq>=freq[0]])
        s12p_sigma = unp.std_devs(self.meas.s12[1][self.meas.freq>=freq[0]])
        
        # Create weighted objective functions for magnitute and phase seperately
        obj_func_real = ((np.absolute(s21c) - np.absolute(s21_predicted))/s21m_sigma + \
                         (np.absolute(s12c) - np.absolute(s12_predicted))/s12m_sigma + \
                         (np.absolute(s11c) - np.absolute(s11_predicted))/s11m_sigma)
        obj_func_imag = ((np.unwrap(np.angle(s21c)) - np.unwrap(np.angle(s21_predicted)))**2/s21p_sigma**2 + \
                         (np.unwrap(np.angle(s12c)) - np.unwrap(np.angle(s12_predicted)))**2/s12p_sigma**2 + \
                         (np.unwrap(np.angle(s11c)) - np.unwrap(np.angle(s11_predicted)))**2/s11p_sigma**2)
#        # Un-weighed objective fucntion
#        obj_func_real = ((np.absolute(s21c) - np.absolute(s21_predicted)) + \
#                         (np.absolute(s12c) - np.absolute(s12_predicted)) + \
#                         (np.absolute(s11c) - np.absolute(s11_predicted)))
#        obj_func_imag = ((np.unwrap(np.angle(s21c)) - np.unwrap(np.angle(s21_predicted)))**2 + \
#                         (np.unwrap(np.angle(s12c)) - np.unwrap(np.angle(s12_predicted)))**2 + \
#                         (np.unwrap(np.angle(s11c)) - np.unwrap(np.angle(s11_predicted)))**2)
        
        return np.concatenate((obj_func_real,obj_func_imag))
    
    def _sparam_iterator(self,params,L,freq_0,s11,s21,s12,s22):
        """
        Perform the s-parameter fit using lmfit and emcee and produce the fit 
            report.
        """
        # Fit data
        minner = Minimizer(self._iterate_objective_function,\
                   params,fcn_args=(L,freq_0,s11,s21,s12,s22),\
                   nan_policy='omit')
        
        result = minner.emcee(steps=300,nwalkers=600,is_weighted=True,burn=50)
        
        report_fit(result)
        
        return result
    
    def _permittivity_iterate(self,corr=False):
        """
        Set up iteration and plot results. Corrected data currently NOT FULLY SUPPORTED.
        """
        number_of_poles = self.poles
        number_of_fits = self.fits
        
        # Get electromagnetic properties
        # Note: does not currently check if using corrected data
        if self.start_freq:     #start iteration from self.start_freq
            freq = self.meas.freq[self.meas.freq>=self.start_freq]
        else:   #use full frequency range
            freq = self.meas.freq
        # Get epsilon    
        epsilon = -1j*unp.nominal_values(self.meas.avg_lossfac);
        epsilon += unp.nominal_values(self.meas.avg_dielec)
        epsilon = epsilon[self.meas.freq>=freq[0]]
        # If ierating for mu, get mu
        if self.fit_mu:
            if self.meas.nrw:   #get epsilon and mu
                mu = self.meas.mu[self.meas.freq>=freq[0]]
                mu_real = mu.real
                mu_imag = mu.imag
                mu = -1j*mu_imag;
                mu += mu_real
            else:   #create mu via nrw
                self.meas.nrw = True
                dielec, lossfac, losstan, mu = \
                    self.meas._permittivity_calc('a')
                mu = mu[self.meas.freq>=freq[0]]
                mu_real = mu.real
                mu_imag = mu.imag
                mu = -1j*mu_imag;
                mu += mu_real
                self.meas.nrw = False    # Reset to previous setting
            
        # Create a set of Parameters to the Cole-Cole model
        if self.fit_mu:
            params = self._iteration_parameters(mu=True)
        else:
            params = self._iteration_parameters()
        
        ## First, fit Cole-Cole model(s) to analytical results to get initial guess
        # Iterate to find parameters       
        miner = Minimizer(self._colecole_residuals,params,\
                          fcn_args=(number_of_poles,freq,epsilon))
        result = miner.minimize()
        if self.fit_mu:
            miner_mu = Minimizer(self._colecole_residuals,params,\
                          fcn_args=(self.poles_mu,freq,mu,True))
            result_mu = miner_mu.minimize()
        
        # Write fit report
        report_fit(result)
        if self.fit_mu:
            report_fit(result_mu)
        
        # Get parameter values
        values = result.params
        values = values.valuesdict()
        if self.fit_mu:
            values_mu = result_mu.params
            values_mu = values_mu.valuesdict()
            # Merge results into single object
            values['mu_inf'] = values_mu['mu_inf']
            for m in range(self.poles_mu):
                m+=1
                values['mu_dc_{}'.format(m)] = values_mu['mu_dc_{}'.format(m)]
                values['mutau_{}'.format(m)] = values_mu['mutau_{}'.format(m)]
                values['mualpha_{}'.format(m)] = values_mu['mualpha_{}'.format(m)]
            print(values)
            
        
        # Calculate model EM parameters
        epsilon_iter = self._colecole(number_of_poles,freq,values)
        if self.fit_mu:
            mu_iter = self._colecole(self.poles_mu,freq,values_mu,mu=True)
        
        # Plot                    
        pplot.make_plot([freq,freq],[epsilon.real,epsilon_iter.real],legend_label=['Analytical','Iterative'])
        pplot.make_plot([freq,freq],[-epsilon.imag,-epsilon_iter.imag],plot_type='lf',legend_label=['Analytical','Iterative'])
        if self.fit_mu:
            pplot.make_plot([freq,freq],[mu.real,mu_iter.real],legend_label=['Analytical mu','Iterative mu'])
            pplot.make_plot([freq,freq],[-mu.imag,-mu_iter.imag],plot_type='lf',legend_label=['Analytical mu','Iterative mu'])
        
        # Find values at 8.5 GHz by finding index where freq is closest to 8.5 GHz
        ep_real = epsilon_iter.real[np.where(freq == freq[np.abs(freq - 8.5e9).argmin()])][0]
        ep_imag = epsilon_iter.imag[np.where(freq == freq[np.abs(freq - 8.5e9).argmin()])][0]
        print(ep_real)
        print(ep_imag)
        if self.fit_mu:
            mu_real = mu_iter.real[np.where(freq == freq[np.abs(freq - 8.5e9).argmin()])][0]
            mu_imag = mu_iter.imag[np.where(freq == freq[np.abs(freq - 8.5e9).argmin()])][0]
            print(mu_real)
            print(mu_imag)
        
        # If not in trial mode (no iterative fitting of sparams), perform iteration
        if not self.trial:
            ## Perform Modified Baker-Jarvis iteration
            # Check if using corrected S-params
            if corr:
                s11 = unp.nominal_values(self.meas.corr_s11)
                s21 = unp.nominal_values(self.meas.corr_s21)
                s12 = unp.nominal_values(self.meas.corr_s12)
                s22 = unp.nominal_values(self.meas.corr_s22)
                L = self.meas.Lcorr
            else:
                # Use shorted S11 if available
                if self.shorted:
                    s11 = unp.nominal_values(self.meas.s11_short)
                else:
                    s11 = unp.nominal_values(self.meas.s11)
                s21 = unp.nominal_values(self.meas.s21)
                s12 = unp.nominal_values(self.meas.s12)
                s22 = unp.nominal_values(self.meas.s22)
                L = self.meas.L
            
            # Start arrays at start_freq
            s11 = np.array((s11[0][self.meas.freq>=freq[0]],s11[1][self.meas.freq>=freq[0]]))
            s21 = np.array((s21[0][self.meas.freq>=freq[0]],s21[1][self.meas.freq>=freq[0]]))
            s12 = np.array((s12[0][self.meas.freq>=freq[0]],s12[1][self.meas.freq>=freq[0]]))
            s22 = np.array((s22[0][self.meas.freq>=freq[0]],s22[1][self.meas.freq>=freq[0]]))
            # Cast measured sparams to complex
            s11c = 1j*s11[0]*np.sin(np.radians(s11[1]));
            s11c += s11[0]*np.cos(np.radians(s11[1]))
            s22c = 1j*s22[0]*np.sin(np.radians(s22[1]));
            s22c += s22[0]*np.cos(np.radians(s22[1]))
            s21c = 1j*s21[0]*np.sin(np.radians(s21[1]));
            s21c += s21[0]*np.cos(np.radians(s21[1]))
            s12c = 1j*s12[0]*np.sin(np.radians(s12[1]));
            s12c += s12[0]*np.cos(np.radians(s12[1]))
            
            ## Perform the fits acording to number_of_fits
            values_sp = values # Use Cole-Cole fit for intial values
            for n in range(number_of_fits):
                # Create a set of Parameters
                initial_values = values_sp
                if self.fit_mu:
                    params = self._iteration_parameters(initial_values,mu=True)
                else:
                    params = self._iteration_parameters(initial_values)
                # Fit data
                result_sp = self._sparam_iterator(params,L,freq[0],s11c,s21c,s12c,s22c)
                # Update initial values for next run
                values_sp = result_sp.params
                values_sp = values_sp.valuesdict()
            
            # Get final parameter values
            values_sp = result_sp.params
            values_sp = values_sp.valuesdict()
            
            # Calculate model EM parameters
            epsilon_iter_sp = self._colecole(number_of_poles,freq,values_sp)
            if self.fit_mu:
                mu_iter_sp = self._colecole(self.poles_mu,freq,values_sp,mu=True)
            
            # Plot                    
            pplot.make_plot([freq,freq],[epsilon.real,epsilon_iter_sp.real],legend_label=['Analytical','Iterative'])
            pplot.make_plot([freq,freq],[-epsilon.imag,-epsilon_iter_sp.imag],plot_type='lf',legend_label=['Analytical','Iterative'])
            if self.fit_mu:
                pplot.make_plot([freq,freq],[mu.real,mu_iter_sp.real],legend_label=['Analytical mu','Iterative mu'])
                pplot.make_plot([freq,freq],[-mu.imag,-mu_iter_sp.imag],plot_type='lf',legend_label=['Analytical mu','Iterative mu'])
        
            # Plot s-params
            s11_predicted, s21_predicted, s12_predicted = self._model_sparams(freq,L/100,epsilon_iter,1)
            import matplotlib.pyplot as plt
            plt.figure()
            plt.plot(freq, np.absolute(s11c),label='Measured')
            plt.plot(freq, np.absolute(s11_predicted),label='Predicted')
            plt.title('s11mag')
            plt.legend()
            plt.figure()
            plt.plot(freq, np.angle(s11c),label='Measured')
            plt.plot(freq, np.angle(s11_predicted),label='Predicted')
            plt.title('s11phase')
            plt.figure()
            plt.plot(freq, np.absolute(s21c),label='Measured')
            plt.plot(freq, np.absolute(s21_predicted),label='Predicted')
            plt.title('s21mag')
            plt.legend()
            plt.figure()
            plt.plot(freq, np.angle(s21c),label='Measured')
            plt.plot(freq, np.angle(s21_predicted),label='Predicted')
            plt.title('s21phase')
            
            # Return results
            if self.fit_mu:
                return epsilon_iter_sp, mu_iter_sp, values_sp
            else:
                return epsilon_iter_sp, values_sp
