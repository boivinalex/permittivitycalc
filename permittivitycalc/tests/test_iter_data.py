# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 16:10:37 2019

@author: alex
"""
import os
if os.environ.get('DISPLAY','') == '':
    print('No display found. Using interactive Qt5Agg backend')
    import matplotlib
    matplotlib.use('Qt5Agg')
import sys
import numpy as np
import unittest
from unittest.mock import patch
import permittivitycalc as pc
import matplotlib.pyplot as plt
plt.ion()

class iter_data_TestCase(unittest.TestCase):
    """Tests for 'iter_data.py"""
    
    
    def setUp(self):
        self.data_path = os.path.join(pc.__path__[0], 'data')
        self.file_path = os.path.join(self.data_path, 'serpentine_dry.txt')
        self.dataset1 = pc.run_default(airline_name='VAL',file_path=self.file_path)
        self.dataset2 = pc.run_default(airline_name='VAL',file_path=self.file_path,nrw=True,temperature=35)
        self.dataset3 = pc.run_default(airline_name='VAL',file_path=self.file_path,nrw=True,corr=True)
        
    def test_run(self):
        print(sys._getframe().f_code.co_name)
        test_iter = pc.piter(self.dataset1)
        self.assertIsNotNone(test_iter.meas)
        plt.close('all')
        return test_iter
        
    def test_run_nrw(self):
        print(sys._getframe().f_code.co_name,)
        test_iter = pc.piter(self.dataset2,fit_mu=True,number_of_poles_mu=0)
        self.assertIsNotNone(test_iter.meas)
        plt.close('all')
        return test_iter
        
    def test_run_nrw_corr(self):
        print(sys._getframe().f_code.co_name)
        test_iter = pc.piter(self.dataset3,fit_mu=True,number_of_poles_mu=0)
        self.assertIsNotNone(test_iter.meas)
        plt.close('all')
        return test_iter
        
    def test_run_nrw_corr_freq_cutoff(self):
        print(sys._getframe().f_code.co_name)
        test_iter = pc.piter(self.dataset3,start_freq=1e7,end_freq=3e9)
        self.assertIsNotNone(test_iter.meas)
        plt.close('all')
        return test_iter
        
    def test_mcmc_run(self):
        print(sys._getframe().f_code.co_name)
        test_iter = pc.piter(self.dataset1,trial_run=False,nsteps=2,nwalkers=5,number_of_poles=0,nburn=1,nthin=1)
        self.assertIsNotNone(test_iter.epsilon_iter)
        plt.close('all')
        return test_iter
        
    def test_mcmc_run_1pole(self):
        print(sys._getframe().f_code.co_name)
        test_iter = pc.piter(self.dataset1,trial_run=False,nsteps=2,nwalkers=10,number_of_poles=1,nburn=1,nthin=1)
        self.assertIsNotNone(test_iter.epsilon_iter)
        plt.close('all')
        return test_iter
        
    def test_mcmc_run_1pole_cond(self):
        print(sys._getframe().f_code.co_name)
        test_iter = pc.piter(self.dataset1,trial_run=False,nsteps=2,nwalkers=10,number_of_poles=1,nburn=1,nthin=1)
        self.assertIsNotNone(test_iter.epsilon_iter)
        plt.close('all')
        return test_iter
        
    def test_mcmc_nrw_corr_run(self):
        print(sys._getframe().f_code.co_name)
        try:
            test_iter = pc.piter(self.dataset3,trial_run=False,nsteps=2,nwalkers=15,number_of_poles=0,nburn=0,nthin=1,fit_mu=True,number_of_poles_mu=0)
            return test_iter
        except Exception as e:
            raise
        plt.close('all')
        
    def test_mcmc_nrw_corr_1pole_run(self):
        print(sys._getframe().f_code.co_name)
        try:
            test_iter = pc.piter(self.dataset3,trial_run=False,nsteps=2,nwalkers=20,number_of_poles=1,nburn=1,nthin=1,fit_mu=True,number_of_poles_mu=1,fit_conductivity=True)
            return test_iter
        except Exception as e:
            raise
        plt.close('all')
        
    # def test_waterpole(self):
    #     print(sys._getframe().f_code.co_name)
    #     try:
    #         test_iter = pc.piter(self.dataset2,trial_run=False,nsteps=5,nwalkers=25,number_of_poles=2,nburn=1,nthin=1,fit_conductivity=True,water_pole=True)
    #         return test_iter
    #     except Exception as e:
    #         raise
    #     plt.close('all')
        
    # def test_waterpole_fail(self):
    #     print(sys._getframe().f_code.co_name)
    #     self.assertRaises(Exception,pc.piter,self.dataset1,trial_run=False,nsteps=5,nwalkers=5,number_of_poles=0,nburn=1,nthin=1,water_pole=True,msg='AirlineData class instance must be given a temperature if using a Debye water pole')
    #     plt.close('all')
        
if __name__ == '__main__':
    unittest.main()