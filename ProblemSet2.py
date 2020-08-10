# %% [md]
# # Problem Set 2
#
# Questions are based on tutorials 2 and 3.
#
# # Futher reading
# 1. [List comprehension](https://www.pythonforbeginners.com/basics/list-comprehensions-in-python#:~:text=List%20comprehensions%20provide%20a%20concise,be%20anything%2C%20meaning%20you%20can)
# 1. [Breaking a loop](https://www.digitalocean.com/community/tutorials/how-to-use-break-continue-and-pass-statements-when-working-with-loops-in-python-3)
# 1. Check out this [guide](https://www.geeksforgeeks.org/check-multiple-conditions-in-if-statement-python/) to compare multiple conditions.
# 1. [Alpha diversity](https://en.wikipedia.org/wiki/Diversity_index)
# ## Set up
# %%
import pandas as pd
import numpy as np
import scipy
from numpy import random

species = ['E.coli', 'C.acnes', 'K.pneumoniae', 'H.sapiens']
n = 100
a = np.random.uniform(0,1,n) * 0.2
b = a * 3 + np.random.uniform(0,1,n) * 0.5
c = np.random.uniform(0,1,n) * 0.2
d = c * 0.5 + np.random.uniform(0,1,n) * 0.1
absolute_abundance = np.array([a,b,c,d])
relative_abundance = absolute_abundance / absolute_abundance.sum(axis=0)
# %% [markdown]
# # Questions 
# 1. Loop through samples to print the relative abundance of *K.pneumoniae*
#
# 1. If sample 1's *E.coli's* relative abundance is greater than 0.5, print('E.coli's rel abd is greater than 0.5')
#
# 1. If sample 2's *E.coli's* relative abundance is lower than 0.5 AND *H.sapiens'* rel abd is lower than 0.3, print('Both conditions are satisfied')
#
# 1. Find all samples where *H.sapiens* relative abundance is greater than 0.2
#
# 1. Find the ratio of *E.coli* / *C.acnes* for the first 10 samples.
#
# 1. Print rel abd of *C.acnes* starting from sample 0 until *C.acnes* > 0.3
#
# 1. Find geometric mean of all the samples 
#
# 1. Remove *H.sapiens* and rebalance the rest of the abundances so that they sum to 1
#
# 1. Find the alpha diversity of each sample using the Shannon index
#
# 1. Find the alpha diversity of each sample using the Simpson index
#
# 1. Use list comprehension to obtain a list of the squares of the first 10 positive integers. 
#
# 1. Write a Python function that that prints out the first n rows of Pascal's triangle.
#
# 1. Write a Python function that that prints out the first n rows of Pascal's triangle.
#     
#     Expected output for 6 rows:
#     ```
#     [1]                                                                                                           
#     [1, 1]                                                                                                        
#     [1, 2, 1]                                                                                                     
#     [1, 3, 3, 1]                                                                                                  
#     [1, 4, 6, 4, 1]                                                                                               
#     [1, 5, 10, 10, 5, 1]
#     ```
#
# 1. Given a nucleotide sequence: 
#     ATGCAAGCGGATATGCATGGAAAACTTCACGCTGCCTTAGAAGATGGTTTCTTCCTCTTTTTTTGAACAGCAACAACAACCTAACATTTATTATGACACAACCACCGATCAAGAAGAC
#
#     Write a function count_substring(dna, substring), to count how many times a certain substring appears in the string of nucleotides.
#     For example, the function returns 6 when called with the substring 'CAA'and the nucleotide string above. 
#
# 1. Write a Python function that returns a list of boolean where the Shannon index > Simpson index
#     
#     Expected output format:
#     ```
#     [True, False, True, True, True, False...]
#     ```

# %%
