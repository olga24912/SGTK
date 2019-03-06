#!/usr/bin/env python

import unittest
import os
import time
from selenium import webdriver

def compileSGTK():
    if os.path.exists("~/tmp/testSGTK"):
        shutil.rmtree("~/tmp/testSGTK")
    
    os.system("PREFIX=~/tmp/testSGTK  ../../compile.sh")
    data_path = "~/CAB/SGTK_/SGTK/resources/data/E.coli/"
    gfa_path = data_path + "assembly_graph_with_scaffolds.gfa"
    first_path = data_path + "s_6_1.fastq.00.0_0.cor.fastq.gz" 
    second_path = data_path + "s_6_2.fastq.00.0_0.cor.fastq.gz"
    ref_path = data_path + "MG1655-K12.fasta"
    long_path = data_path + "filtered_subreads.fastq"
    os.system("python3 ~/tmp/testSGTK/bin/sgtk.py --gfa "+ gfa_path +
              " --fr " + first_path + " " + second_path +
              " --long " + long_path + 
              " --ref " + ref_path +
              " -o ~/tmp/testSGTK/out")

        
class BaseTests(object):
    def test_from_local(self):
        self.driver.get('file:///home/olga/tmp/testSGTK/out/main.html')
        time.sleep(5)
        self.assertEqual('SGTK', self.driver.title)

    def tearDown(self):
        self.driver.quit()


class FirefoxTests(unittest.TestCase, BaseTests):
    def setUp(self):
        self.driver = webdriver.Firefox()


class ChromeTests(unittest.TestCase, BaseTests):
    def setUp(self):
        self.driver = webdriver.Chrome()


class OperaTests(unittest.TestCase, BaseTests):
    def setUp(self):
        self.driver = webdriver.Opera()


#class SafariTests(unittest.TestCase, BaseTests):
#    def setUp(self):
#        self.driver = webdriver.Safari()


if __name__ == '__main__':
    #compileSGTK()
    unittest.main()
