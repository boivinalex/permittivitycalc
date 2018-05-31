#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 18:39:15 2018

@author: alex
"""

import os
import unittest
from unittest.mock import patch
import permittivitycalc as pc
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    import matplotlib
    matplotlib.use('Agg')
import matplotlib.pyplot as plt


class sparam_data_TestCase(unittest.TestCase):
    """Tests for 'sparam_data.py"""
    
    def setUp(self):
        self.data_path = os.path.join(pc.__path__[0], 'data')
        self.file_path = os.path.join(self.data_path, 'serpentine_dry.txt')
#        self.file_path2 = os.path.join(self.data_path, 'rexolite_PAL.txt')
        self.fake_path = os.path.join(self.data_path, 'fake.txt')
        self.dataset1 = pc.run_default(airline_name='VAL',file_path=self.file_path)
#        self.dataset2 = pc.run_default(airline_name='VAL',file_path=self.file_path2)
        self.normdataset = pc.AirlineData(*pc.get_METAS_data(airline='VAL',file_path=self.file_path),bulk_density=1.6,name='Serpentine',normalize_density=True)

    def test_repr(self):
        self.assertIsNotNone(self.dataset1.__repr__())

    def test_str(self):
        self.assertIsNotNone(self.dataset1.__str__())

    def test_file_import(self):
        """Test file import"""
        self.assertIsNotNone(pc.AirlineData\
                (*pc.get_METAS_data(airline='VAL',file_path=self.file_path)))
        
    def test_get_data(self):
        """Test that dielectric data is accessible"""
        self.assertIsNotNone(self.dataset1.avg_dielec)
        self.assertIsNotNone(self.dataset1.avg_lossfac)
        self.assertIsNotNone(self.dataset1.avg_losstan)
        
    def test_nrw(self):
        """Test nrw"""
        test = pc.AirlineData\
                (*pc.get_METAS_data(airline='VAL',file_path=self.file_path),nrw=True)
        self.assertIsNotNone(test)
        self.assertIsNotNone(test.mu)
        self.assertIsNotNone(test.avg_dielec)
        self.assertIsNotNone(test.avg_losstan)
        
    def test_normalize_density(self):
        """Test density normalization"""
        test = pc.AirlineData\
                (*pc.get_METAS_data(airline='VAL',file_path=self.file_path)\
                 ,normalize_density=True,bulk_density=3.5)
        self.assertIsNotNone(test.norm_dielec)
        self.assertIsNotNone(test.norm_losstan)
        with self.assertRaises(Exception):
            pc.AirlineData(*pc.get_METAS_data(airline='VAL',\
                            file_path=self.file_path),normalize_density=True)


    def test_shorted(self):
        """Test fail to locate shorted sample"""
        with self.assertRaises(Exception):
            pc.AirlineData(*pc.get_METAS_data(airline='VAL',\
                            file_path=self.file_path),shorted=True)

            
    def test_boundary_correct(self):
        """Test boundary correction"""
        test = pc.AirlineData\
                (*pc.get_METAS_data(airline='VAL',file_path=self.file_path)\
                 ,normalize_density=True,bulk_density=3.5,solid_dielec=9\
                 ,particle_diameter=0.01,particle_density=3)
        self.assertIsNotNone(test.bcorr_dielec)
        self.assertIsNotNone(test.bcorr_losstan)

    def test_boundary_correct_with_losstan(self):
        """Test boundary correction"""
        test = pc.AirlineData\
                (*pc.get_METAS_data(airline='VAL',file_path=self.file_path)\
                 ,normalize_density=True,bulk_density=3.5,solid_dielec=9\
                 ,solid_losstan=0.1,particle_diameter=0.01,particle_density=3)
        self.assertIsNotNone(test.bcorr_dielec)
        self.assertIsNotNone(test.bcorr_losstan)
        
    def test_optional_attributes(self):
        pass #TODO
        
    def test_dims(self):
        """Test that the proper number of airline dimentions are calculated"""
        self.assertIsNotNone(self.dataset1.airline_dimensions)
        assert len(self.dataset1.airline_dimensions) == 2
        test = pc.AirlineData\
                (*pc.get_METAS_data(airline='VAL',file_path=self.file_path)\
                 ,particle_diameter=0.01)
        assert len(test.airline_dimensions) == 4

    def test_dims_gal(self):
        """Test that the proper number of airline dimentions are calculated"""
        self.assertIsNotNone(self.dataset1.airline_dimensions)
        assert len(self.dataset1.airline_dimensions) == 2
        test = pc.AirlineData\
                (*pc.get_METAS_data(airline='GAL',file_path=self.file_path)\
                 ,particle_diameter=0.01)
        assert len(test.airline_dimensions) == 4
        
    def test_air_gap_correction(self):
        test = pc.AirlineData\
                (*pc.get_METAS_data(airline='VAL',file_path=self.file_path)\
                 ,particle_diameter=0.01)
        D2 = test.airline_dimensions['D2']
        D3 = test.airline_dimensions['D3']
        self.assertIsNotNone(test._air_gap_correction(D2,D3))
        
    def test_res_freq(self):
        """Test resonant_freq"""
        self.assertIsNotNone(self.dataset1.res_freq)
        
    def test_freq_avg(self):
        """Test freq_avg"""
        self.assertIsNotNone(self.dataset1._freq_avg())
        
    def test_draw_plots(self):
        """Test draw_plots"""
        try:
            self.dataset1.draw_plots()
        except Exception as e:
            raise
        plt.close('all')
            
    def test_draw_plots_normalized(self):
        """Test draw_plots with normalized data"""
        try:
            self.normdataset.draw_plots(normalized=True)
        except Exception as e:
            raise
        plt.close('all')

    def test_draw_plots_corr_notdefault(self):
        """Test draw_plots fails when corr and not default_setting"""
        with self.assertRaises(Exception):
            self.dataset1.draw_plots(corr=True,default_setting=False)

    @patch('builtins.input',return_value='a')
    def test_draw_plots_notdefault_a(self,mock):
        """Test draw_plots average plot"""
        try:
            self.dataset1.draw_plots(default_setting=False)
        except Exception as e:
            raise
        plt.close('all')

    @patch('builtins.input',return_value='f')
    def test_draw_plots_notdefault_f(self,mock):
        """Test draw_plots forward plot"""
        try:
            self.dataset1.draw_plots(default_setting=False)
        except Exception as e:
            raise
        plt.close('all')

    @patch('builtins.input',return_value='r')
    def test_draw_plots_notdefault_r(self,mock):
        """Test draw_plots backwards plot"""
        try:
            self.dataset1.draw_plots(default_setting=False)
        except Exception as e:
            raise
        plt.close('all')

    @patch('builtins.input',return_value='b')
    def test_draw_plots_notdefault_b(self,mock):
        """Test draw_plots both plots"""
        try:
            self.dataset1.draw_plots(default_setting=False)
        except Exception as e:
            raise
        plt.close('all')

    @patch('builtins.input',return_value='wrong')
    def test_draw_plots_notdefault_wrong(self,mock):
        """Test draw_plots fails for wrong plot type"""
        with self.assertRaises(Exception):
            self.dataset1.draw_plots(default_setting=False)
         
    def test_s_param_plots(self):
        """Test s_param_plots"""
        try:
            self.dataset1.s_param_plot()
        except Exception as e:
            raise
        plt.close('all')

    def test_multiple_meas(self):
        """Test multiple_meas"""
        try:
            dataset_list = pc.multiple_meas(file_path=self.file_path,airline_name='VAL')
            self.assertIsNotNone(dataset_list)
            assert len(dataset_list) == 2
        except Exception as e:
            raise
        plt.close('all')

    def test_multiple_meas_nofile(self): #FIX THIS
        """Test multiple_meas"""
        try:
            with patch('permittivitycalc.helper_functions._get_file',return_value=['VAL',self.file_path,None]):
                dataset_list = pc.multiple_meas(airline_name='VAL')
                self.assertIsNotNone(dataset_list)
                assert len(dataset_list) == 2
        except Exception as e:
            raise
        plt.close('all')
        
if __name__ == '__main__':
    unittest.main()