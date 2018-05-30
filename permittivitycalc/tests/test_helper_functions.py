#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 23 16:03:51 2018

@author: alex
"""

import os
import unittest
from unittest.mock import patch
from unittest.mock import Mock
import permittivitycalc.helper_functions as hf


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
        
    def test_get_file(self):
        """Test file import"""
        real_file = hf._get_file('VAL',self.file_path)
        fake_file = hf._get_file('VAL',self.fake_path)
        self.assertIsNotNone(real_file)
        self.assertIsNotNone(fake_file)
        assert os.path.isfile(real_file[1])
        assert not os.path.isfile(fake_file[1])
        
    def test_get_file_wrong_airline(self):
        """Test wrong airline"""
        kwargs = {"airline":'wrong',"file_path":self.file_path}
        self.assertRaises(Exception,hf._get_file,**kwargs)
    
    @patch('builtins.input',return_value=5)    
    def test_get_file_custom_airline(self,mock):
        """Test custom airline"""
        custom_line = hf._get_file('custom',self.file_path)
        self.assertIsNotNone(custom_line)

    @patch('builtins.input',return_value='VAL')
    def test_get_file_no_airline(self,mock):
        noline = hf._get_file(file_path=self.file_path)
        self.assertIsNotNone(noline)
        
    @patch('builtins.input',return_value='wrong')
    def test_get_file_no_airline_wrong_input(self,mock):
        self.assertRaises(Exception,hf._get_file,file_path=self.file_path)
        
    @patch('builtins.input',side_effect=['custom',5])    
    def test_get_file_no_airline_custom(self,mock):
        noline = hf._get_file(file_path=self.file_path)
        self.assertIsNotNone(noline)
        
    def test_get_metas_data(self):
        """Test data import"""
        try:
            real_file = hf.get_METAS_data('VAL',self.file_path)
            self.assertIsNotNone(real_file)
            assert len(real_file) == 4
        except Exception as e:
            raise
        kwargs = {"airline":'VAL',"file_path":self.fake_path}
        self.assertRaises(FileNotFoundError,hf.get_METAS_data,**kwargs)
        
    def test_perm_compare(self):
        datasets = [self.dataset1,self.dataset2]
        try:
            hf.perm_compare(datasets)
        except Exception as e:
            raise

if __name__ == '__main__':
    unittest.main()
        
    
