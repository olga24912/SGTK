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
    def test_load(self):
        self.driver.get('file:///home/olga/tmp/testSGTK/out/main.html')
        time.sleep(5)
        self.assertEqual('SGTK', self.driver.title)

        info = self.driver.find_element_by_xpath("//div[@id='extra_info']/p[1]")
        self.assertEqual('Nodes: 850\nEdges: 19279\nChromosomes: 2\nSources: 4', info.text)

        trs_libs = self.driver.find_elements_by_xpath("//table[@id='lib_table']/tr")
        self.assertEqual(len(trs_libs), 4)

        blocks = self.driver.find_elements_by_xpath("//div[@id='show_block']/div[@class='block']")
        page_num = self.driver.find_element_by_id("choose_page").get_attribute('max')
        self.assertGreater(len(blocks), 0)
        self.assertEqual(len(blocks), int(page_num))

        graph = self.driver.find_elements_by_class_name("__________cytoscape_container")
        self.assertEqual(len(graph), 1)
        
        


class FirefoxTests(unittest.TestCase, BaseTests):
    def setUp(self):
        self.driver = webdriver.Firefox()

    def tearDown(self):
        self.driver.quit()


class ChromeTests(unittest.TestCase, BaseTests):
    def setUp(self):
        self.driver = webdriver.Chrome()

    def tearDown(self):
        self.driver.quit()


class OperaTests(unittest.TestCase, BaseTests):
    def setUp(self):
        self.driver = webdriver.Opera()

    def tearDown(self):
        self.driver.quit()


#class SafariTests(unittest.TestCase, BaseTests):
#    def setUp(self):
#        self.driver = webdriver.Safari()
#
#    def tearDown(self):
#        self.driver.quit()


if __name__ == '__main__':
    #compileSGTK()
    unittest.main()
