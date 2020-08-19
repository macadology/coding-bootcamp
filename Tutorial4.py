# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 19:02:01 2020
Practise session
@author: raarthi
"""

import doctest


def print_grades(score):
    ''' this function can be used to print the grades
  (int) --> NoneType
  >>> print_grades(70)
  You must work harder!
  >>> print_grades(90)
  S
  >>> print_grades(85)
  A'''
    #myscore = '79'
    #print(myscore)
    if 85 < score <= 100:  # this statement tests if the scoee
        print('S')
    elif 75 < score <= 85:
        print('A')
    else:
        print('You must work harder!')


doctest.testmod(verbose=True)

# String operations
'MNKMDLVADVAEKTDLSKAKATEVIDAVFA'.count('DL')
'MNKMDLVADVAEKTDLSKAKATEVIDAVFA'.lower()
'MNKMDLVADVAEKTDLSKAKATEVIDAVFA'.upper()
'MNKMDLVADVAEKTDLSKAKATEVIDAVFA'.replace('DL', 'AD')
'MNKMDLVADVAEKTDLSKAKATEVIDAVFA'.split('DL')
'MNKMDLVADVAEKTDLSKAKATEVIDAVFA'.strip('FA')
'MNKMDLVADVAEKTDLSKAKATEVIDAVFA\nMNKMDLVADVAEKTDLSKAKATEVIDAVFA'.splitlines()
'MNKMDLVADVAEKTDLSKAKATEVIDAVFA'.index('DL')
'MNKMDLVADVAEKTDLSKAKATEVIDAVFA'.find('DL')
'MNKMDLVADVAEKTDLSKAKATEVIDAVFA'.rfind('DL')


# How to find all the indices of substring within a string?
def find_all(a_string, sub):
    ''' this function finds all the indices of substring within a given string'''
    result = []
    kcount = 0
    while kcount < len(a_string):
        kcount = a_string.find(sub, kcount)
        if kcount == -1:
            #print(kcount)
            return result
        else:
            result.append(kcount)
            #print(kcount)
            kcount += 1  #change to kcount += len(sub) to not search overlapping results
    return result
