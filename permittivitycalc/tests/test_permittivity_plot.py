#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 17:37:35 2018

@author: alex
"""

import numpy as np
import unittest
import permittivitycalc.permittivity_plot as pp

class permittivity_plot_TestCase(unittest.TestCase):
    """Tests for 'permittivity_plot.py"""
        
    def test_make_plot(self):
        """Test make_plot"""
        try:
            x = np.arange(0,6)*300000
            y = np.arange(0,6)
            pp.make_plot(x,y)
        except Exception as e:
            raise

    def test_make_plot_custom(self):
        """Test make_plot with custom plot_type"""
        try:
            x = np.arange(0,6)*300000
            y = np.arange(0,6)
            pp.make_plot(x,y,plot_type='c',plot_title='test',ylabel='test',xlabel='test',round_val=2)
        except Exception as e:
            raise

    def test_make_plot_invalid_plot_type(self):
        """Test makes_plot with wrong plot_type"""
        self.assertRaises(Exception,pp.make_plot,plot_type='wrong')
            
    def test_make_sparam_plot(self):
        """Test make_sprama_plot"""
        freq_1 = np.arange(0,6)
        s11_1 = np.array([np.arange(0,6),np.arange(0,6)])
        s22_1 = np.array([np.arange(0,6),np.arange(0,6)])
        s21_1 = np.array([np.arange(0,6),np.arange(0,6)])
        s12_1 = np.array([np.arange(0,6),np.arange(0,6)])
        try:
            pp.make_sparam_plot(freq_1,s11_1,s22_1,s21_1,s12_1)
        except Exception as e:
            raise

    def test_make_sparam_plot_2plots(self):
        """Test make_sprama_plot with 2 plots"""
        try:
            s11_1 = np.array([np.arange(0,6),np.arange(0,6)])
            s22_1 = np.array([np.arange(0,6),np.arange(0,6)])
            s21_1 = np.array([np.arange(0,6),np.arange(0,6)])
            s12_1 = np.array([np.arange(0,6),np.arange(0,6)])
            s11_2 = np.array([np.arange(0,6),np.arange(0,6)])*2
            s22_2 = np.array([np.arange(0,6),np.arange(0,6)])*2
            s21_2 = np.array([np.arange(0,6),np.arange(0,6)])*2
            s12_2 = np.array([np.arange(0,6),np.arange(0,6)])*2
            freq = np.arange(0,6)
            s11 = [s11_1,s11_2]
            s22 = [s22_1,s22_2]
            s21 = [s21_1,s21_2]
            s12 = [s12_1,s12_2]
            pp.make_sparam_plot(freq,s11,s22,s21,s12)
        except Exception as e:
            raise
        
    def test_make_sparam_plot_2plots_labels(self):
        """Test make_sprama_plot with 2 plots and labels"""
        try:
            s11_1 = np.array([np.arange(0,6),np.arange(0,6)])
            s22_1 = np.array([np.arange(0,6),np.arange(0,6)])
            s21_1 = np.array([np.arange(0,6),np.arange(0,6)])
            s12_1 = np.array([np.arange(0,6),np.arange(0,6)])
            s11_2 = np.array([np.arange(0,6),np.arange(0,6)])*2
            s22_2 = np.array([np.arange(0,6),np.arange(0,6)])*2
            s21_2 = np.array([np.arange(0,6),np.arange(0,6)])*2
            s12_2 = np.array([np.arange(0,6),np.arange(0,6)])*2
            freq = np.arange(0,6)
            s11 = [s11_1,s11_2]
            s22 = [s22_1,s22_2]
            s21 = [s21_1,s21_2]
            s12 = [s12_1,s12_2]
            pp.make_sparam_plot(freq,s11,s22,s21,s12,label=['test1','test2'])
        except Exception as e:
            raise    
            
if __name__ == '__main__':
    unittest.main()