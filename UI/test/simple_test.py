#!/usr/bin/env python

from selenium import webdriver

driver = webdriver.Firefox()
driver.get('file:///home/olga/CAB/SGTK_/SGTK/UI/mainPage.html')
driver.quit()
