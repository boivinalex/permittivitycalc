#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 17:37:35 2018

@author: alex
"""
import sys
import numpy as np
import unittest
import permittivitycalc.permittivity_plot as pp
import matplotlib.pyplot as plt
plt.ion()

class permittivity_plot_TestCase(unittest.TestCase):
    """Tests for 'permittivity_plot.py"""
        
    def test_make_plot(self):
        """Test make_plot"""
        print(sys._getframe().f_code.co_name)
        try:
            x = np.arange(0,6)*300000
            y = np.arange(0,6)
            pp.make_plot(x,y)
        except Exception as e:
            raise
        plt.close('all')
            
    def test_make_plot_2plots(self):
        """Test make_plot"""
        print(sys._getframe().f_code.co_name)
        try:
            x1 = np.arange(0,6)*300000
            y1 = np.arange(0,6)
            x2 = np.arange(0,6)*300000
            y2 = np.arange(0,6)
            x = [x1,x2]
            y = [y1,y2]
            pp.make_plot(x,y)
        except Exception as e:
            raise
        plt.close('all')
            
    def test_make_plot_9plots_lf(self):
        """Test make_plot"""
        print(sys._getframe().f_code.co_name)
        try:
            x1 = np.arange(0,6)*300000
            y1 = np.arange(0,6)
            x2 = np.arange(0,6)*300000
            y2 = np.arange(0,6)
            x3 = np.arange(0,6)*300000
            y3 = np.arange(0,6)
            x4 = np.arange(0,6)*300000
            y4 = np.arange(0,6)
            x5 = np.arange(0,6)*300000
            y5 = np.arange(0,6)
            x6 = np.arange(0,6)*300000
            y6 = np.arange(0,6)
            x7 = np.arange(0,6)*300000
            y7 = np.arange(0,6)
            x8 = np.arange(0,6)*300000
            y8 = np.arange(0,6)
            x9 = np.arange(0,6)*300000
            y9 = np.arange(0,6)
            x = [x1,x2,x3,x4,x5,x6,x7,x8,x9]
            y = [y1,y2,y3,y4,y5,y6,y7,y8,y9]
            pp.make_plot(x,y)
        except Exception as e:
            raise
        plt.close('all')
            
    def test_make_plot_cutoff(self):
        """Test make_plot"""
        print(sys._getframe().f_code.co_name)
        try:
            x = np.arange(0,6)*300000
            y = np.arange(0,6)
            pp.make_plot(x,y,freq_cutoff=600000)
        except Exception as e:
            raise
        plt.close('all')

    def test_make_plot_custom(self):
        """Test make_plot with custom plot_type"""
        print(sys._getframe().f_code.co_name)
        try:
            x = np.arange(0,6)*300000
            y = np.arange(0,6)
            pp.make_plot(x,y,plot_type='c',plot_title='test',ylabel='test',xlabel='test',round_val=2)
        except Exception as e:
            raise
        plt.close('all')

    def test_make_plot_invalid_plot_type(self):
        """Test makes_plot with wrong plot_type"""
        print(sys._getframe().f_code.co_name)
        self.assertRaises(Exception,pp.make_plot,plot_type='wrong')
            
    def test_make_sparam_plot(self):
        """Test make_sprama_plot"""
        print(sys._getframe().f_code.co_name)
        freq_1 = np.arange(0,6)
        s11_1 = np.array([np.arange(0,6),np.arange(0,6)])
        s22_1 = np.array([np.arange(0,6),np.arange(0,6)])
        s21_1 = np.array([np.arange(0,6),np.arange(0,6)])
        s12_1 = np.array([np.arange(0,6),np.arange(0,6)])
        try:
            pp.make_sparam_plot(freq_1,s11_1,s22_1,s21_1,s12_1)
        except Exception as e:
            raise
        plt.close('all')

    def test_make_sparam_plot_2plots(self):
        """Test make_sprama_plot with 2 plots"""
        print(sys._getframe().f_code.co_name)
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
        plt.close('all')
        
    def test_make_sparam_plot_2plots_labels(self):
        """Test make_sprama_plot with 2 plots and labels"""
        print(sys._getframe().f_code.co_name)
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
        plt.close('all')
            
if __name__ == '__main__':
    unittest.main()