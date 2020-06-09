#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 23 16:03:51 2018

@author: alex
"""

import os
if os.environ.get('DISPLAY','') == '':
    print('No display found. Using non-interactive Agg backend')
    import matplotlib
    matplotlib.use('Qt5Agg')
import sys
import unittest
from unittest.mock import patch
import permittivitycalc.helper_functions as hf
import matplotlib.pyplot as plt
plt.ion()


class helper_functions_TestCase(unittest.TestCase):
    """Tests for 'helpper_funtions.py"""
    
    def setUp(self):
        import permittivitycalc as pc
        self.data_path = os.path.join(pc.__path__[0], 'data')
        self.file_path = os.path.join(self.data_path, 'serpentine_dry.txt')
        self.file_path2 = os.path.join(self.data_path, 'rexolite_PAL.txt')
        self.fake_path = os.path.join(self.data_path, 'fake.txt')
        self.dataset1 = pc.run_default(airline_name='VAL',file_path=self.file_path)
        self.dataset2 = pc.run_default(airline_name='VAL',file_path=self.file_path2)
        self.nonedataset = pc.AirlineData(*pc.get_METAS_data(airline='VAL',file_path=self.file_path2),bulk_density=1.6,name='Serpentine',normalize_density=True,freq_cutoff=None)
        self.cutdataset = pc.AirlineData(*pc.get_METAS_data(airline='VAL',file_path=self.file_path2),bulk_density=1.6,name='Serpentine',normalize_density=True,freq_cutoff=2e8)
        self.normdataset = pc.AirlineData(*pc.get_METAS_data(airline='VAL',file_path=self.file_path2),bulk_density=1.6,name='Serpentine',normalize_density=True)
        self.corrdataset = pc.AirlineData(*pc.get_METAS_data(airline='VAL',file_path=self.file_path2),bulk_density=1.6,name='Serpentine',corr=True)
    
    # #patch tkinter modules
    # @patch('tkinter.Tk') 
    # @patch('tkinter.Tk.withdraw')
    # @patch('tkinter.Tk.update')
    # def test_prompt(self,mock1,mock2,mock3):
    #     """Test _prompt"""
    #     with patch('tkinter.filedialog.askopenfilename',autospec=True) as mock:   # patch askopenfilename
    #         with patch('os.fspath',return_value=self.file_path,autospec=True): # patch os.fspath used by tkinter to split path wirh os.path.split
    #             #set return value for askopenfilename
    #             mock().return_value = self.file_path
    #             #run function
    #             result = hf._prompt()
    #             print(self.file_path)
    #             self.assertEqual(result,self.file_path)

    def test_get_file(self):
        """Test file import"""
        print(sys._getframe().f_code.co_name)
        real_file = hf._get_file('VAL',self.file_path)
        fake_file = hf._get_file('VAL',self.fake_path)
        self.assertIsNotNone(real_file)
        self.assertIsNotNone(fake_file)
        assert os.path.isfile(real_file[1])
        assert not os.path.isfile(fake_file[1])
        
    def test_get_file_wrong_airline(self):
        """Test wrong airline"""
        print(sys._getframe().f_code.co_name)
        kwargs = {"airline":'wrong',"file_path":self.file_path}
        self.assertRaises(Exception,hf._get_file,**kwargs)
    
    @patch('builtins.input',return_value=5)    
    def test_get_file_custom_airline(self,mock):
        """Test custom airline"""
        print(sys._getframe().f_code.co_name)
        custom_line = hf._get_file('custom',self.file_path)
        self.assertIsNotNone(custom_line)

    @patch('builtins.input',return_value='VAL')
    def test_get_file_no_airline(self,mock):
        print(sys._getframe().f_code.co_name)
        noline = hf._get_file(file_path=self.file_path)
        self.assertIsNotNone(noline)
        
    @patch('builtins.input',return_value='wrong')
    def test_get_file_no_airline_wrong_input(self,mock):
        print(sys._getframe().f_code.co_name)
        self.assertRaises(Exception,hf._get_file,file_path=self.file_path)
        
    @patch('builtins.input',side_effect=['custom',5])    
    def test_get_file_no_airline_custom(self,mock):
        print(sys._getframe().f_code.co_name)
        noline = hf._get_file(file_path=self.file_path)
        self.assertIsNotNone(noline)

    def test_get_file_no_file(self):
        print(sys._getframe().f_code.co_name)
        with patch('permittivitycalc.helper_functions._prompt',return_value=self.file_path):
            nofile = hf._get_file(airline='VAL')
            self.assertIsNotNone(nofile)

    def test_get_file_no_file_wrong_airline(self):
        print(sys._getframe().f_code.co_name)
        self.assertRaises(Exception,hf._get_file,airline='wrong')

    @patch('builtins.input',return_value=5)
    def test_get_file_no_file_custom_airline(self,mock):
        print(sys._getframe().f_code.co_name)
        with patch('permittivitycalc.helper_functions._prompt',return_value=self.file_path):
            nofile = hf._get_file(airline='custom')
            self.assertIsNotNone(nofile)

    @patch('builtins.input',return_value='VAL')
    def test_get_file_nofile_noline(self,mock):
        print(sys._getframe().f_code.co_name)
        with patch('permittivitycalc.helper_functions._prompt',return_value=self.file_path):
            file = hf._get_file()
            self.assertIsNotNone(file)

    @patch('builtins.input',side_effect=['custom',5])
    def test_get_file_nofile_noline_custom(self,mock):
        print(sys._getframe().f_code.co_name)
        with patch('permittivitycalc.helper_functions._prompt',return_value=self.file_path):
            file = hf._get_file()
            self.assertIsNotNone(file)

    @patch('builtins.input',return_value='wrong')
    def test_get_file_nofile_noline_wrong(self,mock):
        print(sys._getframe().f_code.co_name)
        self.assertRaises(Exception,hf._get_file)
        
    def test_get_metas_data(self):
        """Test data import"""
        print(sys._getframe().f_code.co_name)
        real_file = hf.get_METAS_data('VAL',self.file_path)
        self.assertIsNotNone(real_file)
        assert len(real_file) == 4
        kwargs = {"airline":'VAL',"file_path":self.fake_path}
        self.assertRaises(FileNotFoundError,hf.get_METAS_data,**kwargs)

    @patch('builtins.input',return_value='7')
    def test_get_metas_data_prompt_for_7line(self,mock):
        print(sys._getframe().f_code.co_name)
        file = hf.get_METAS_data(file_path=self.file_path)
        self.assertIsNotNone(file)

    @patch('builtins.input',return_value='GAL')
    def test_get_metas_data_prompt_for_galline(self,mock):
        print(sys._getframe().f_code.co_name)
        file = hf.get_METAS_data(file_path=self.file_path)
        self.assertIsNotNone(file)

    @patch('builtins.input',side_effect=['custom',5])
    def test_get_metas_data_prompt_for_customline(self,mock):
        print(sys._getframe().f_code.co_name)
        file = hf.get_METAS_data(file_path=self.file_path)
        self.assertIsNotNone(file)

    def test_perm_compare(self):
        print(sys._getframe().f_code.co_name)
        datasets = [self.dataset1,self.dataset2]
        try:
            hf.perm_compare(datasets)
            plt.close('all')
        except Exception as e:
            raise

    def test_perm_compare_allplots(self):
        print(sys._getframe().f_code.co_name)
        datasets = [self.dataset1,self.dataset2]
        try:
            hf.perm_compare(datasets,allplots=True)
        except Exception as e:
            raise
        plt.close('all')

    def test_perm_compare_norm(self):
        print(sys._getframe().f_code.co_name)
        datasets = [self.dataset1,self.normdataset]
        try:
            hf.perm_compare(datasets)
        except Exception as e:
            raise
        plt.close('all')
            
    def test_perm_compare_corr(self):
        print(sys._getframe().f_code.co_name)
        datasets = [self.dataset1,self.corrdataset]
        try:
            hf.perm_compare(datasets)
        except Exception as e:
            raise
        plt.close('all')
            
    def test_perm_compare_none(self):
        print(sys._getframe().f_code.co_name)
        datasets = [self.dataset1,self.nonedataset]
        try:
            hf.perm_compare(datasets)
        except Exception as e:
            raise
        plt.close('all')
            
    def test_perm_compare_cut(self):
        print(sys._getframe().f_code.co_name)
        datasets = [self.dataset1,self.cutdataset]
        try:
            hf.perm_compare(datasets)
        except Exception as e:
            raise
        plt.close('all')

    def test_perm_compare_fail(self):
        print(sys._getframe().f_code.co_name)
        fake = 5
        datasets = [self.dataset1,fake]
        self.assertRaises(Exception,hf.perm_compare,datasets)
        
    def test_perm_compare_fail_nolist(self):
        print(sys._getframe().f_code.co_name)
        self.assertRaises(Exception,hf.perm_compare,self.dataset1)

if __name__ == '__main__':
    unittest.main()
        
    
