#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 17:37:35 2018

@author: alex
"""
import os
if os.environ.get('DISPLAY','') == '':
    print('No display found. Using non-interactive Agg backend')
    import matplotlib
    matplotlib.use('Agg')
import sys
import numpy as np
import unittest
from unittest.mock import patch
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
            pp.make_plot(x,y,plot_type='c',plot_title='test',ylabel='test',xlabel='test',xticks=[0,2,4,6],yticks=[0,2,4,6])
        except Exception as e:
            raise
        plt.close('all')
        
    def test_make_plot_ur(self):
        """Test make_plot with ur plot_type"""
        print(sys._getframe().f_code.co_name)
        try:
            x = np.arange(0,6)*300000
            y = np.arange(0,6)
            pp.make_plot(x,y,plot_type='ur')
        except Exception as e:
            raise
        plt.close('all')
        
    def test_make_plot_ui(self):
        """Test make_plot with ui plot_type"""
        print(sys._getframe().f_code.co_name)
        try:
            x = np.arange(0,6)*300000
            y = np.arange(0,6)
            pp.make_plot(x,y,plot_type='ui')
        except Exception as e:
            raise
        plt.close('all')
        
    def test_make_plot_log(self):
        """Test make_plot with log plot"""
        print(sys._getframe().f_code.co_name)
        try:
            x = np.arange(0,6)*300000
            y = np.arange(0,6)
            pp.make_plot(x,y,plot_type='d',y_axis_type='log',xticks=[0,1,2,3], yticks=[0,1,2,3])
        except Exception as e:
            raise
        plt.close('all')
        
    def test_make_plot_log_flat(self):
        """Test make_plot with log plot"""
        print(sys._getframe().f_code.co_name)
        try:
            x = np.arange(0,6)*300000
            y = [1,2,3,4,5,6]
            pp.make_plot(x,y,plot_type='d',y_axis_type='log',xticks=[0,1,2,3], yticks=[0,1,2,3])
            y = [-0.00005,-0.00004,-0.00003,-0.00002,-0.00001,0]
            pp.make_plot(x,y,plot_type='lt',y_axis_type='log',xticks=[0,1,2,3], yticks=[0,1,2,3])
        except Exception as e:
            raise
        plt.close('all')

    def test_make_plot_invalid_plot_type(self):
        """Test makes_plot with wrong plot_type"""
        print(sys._getframe().f_code.co_name)
        x = np.arange(0,6)*300000
        y = np.arange(0,6)
        self.assertRaises(Exception,pp.make_plot,x,y,plot_type='wrong',msg='Invalid plot type')
        
#    def test_get_file_nofile_noline(self,mock):
#        print(sys._getframe().f_code.co_name)
#        with patch('permittivitycalc.helper_functions._prompt',return_value=self.file_path):
#            file = hf._get_file()
#            self.assertIsNotNone(file)
        
    def test_save_plot(self):
        """Test publish"""
        print(sys._getframe().f_code.co_name)
        with patch('permittivitycalc.pplot._dirprompt',return_value='./'):
            try:
                x = np.arange(0,6)*300000
                y = np.arange(0,6)
                pp.make_plot(x,y,publish=True,name='test')
            except Exception as e:
                raise
        plt.close('all')
            
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