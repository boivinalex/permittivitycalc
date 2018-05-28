#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 18:39:15 2018

@author: alex
"""

import os
import unittest
import permittivitycalc as pc


class sparam_data_TestCase(unittest.TestCase):
    """Tests for 'sparam_data.py"""
    
    def setUp(self):
        self.data_path = os.path.join(pc.__path__[0], 'data')
        self.file_path = os.path.join(self.data_path, 'serpentine_dry.txt')
#        self.file_path2 = os.path.join(self.data_path, 'rexolite_PAL.txt')
        self.fake_path = os.path.join(self.data_path, 'fake.txt')
        self.dataset1 = pc.run_default(airline_name='VAL',file_path=self.file_path)
#        self.dataset2 = pc.run_default(airline_name='VAL',file_path=self.file_path2)

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
            
    def test_boundary_corrent(self):
        """Test boundary correction"""
        test = pc.AirlineData\
                (*pc.get_METAS_data(airline='VAL',file_path=self.file_path)\
                 ,normalize_density=True,bulk_density=3.5,solid_dielec=9\
                 ,particle_diameter=0.01,particle_density=3)
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
            
    def test_s_param_plots(self):
        """Test draw_plots"""
        try:
            self.dataset1.s_param_plot()
        except Exception as e:
            raise
        
if __name__ == '__main__':
    unittest.main()