#!/usr/bin/env python

import unittest
from selenium import webdriver

class Selenium2OnLocal(unittest.TestCase):
    def setUp(self):
        self.driver = webdriver.Firefox()

    def test_from_local(self):
        self.driver.get('file:///home/olga/CAB/SGTK_/SGTK/UI/mainPage.html')
        self.assertEqual('SGTK', self.driver.title)

    def tearDown(self):
        self.driver.quit()


if __name__ == '__main__':
    unittest.main()
