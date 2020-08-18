# %% [md]
# # Problem Set 2
#
# Questions are based on tutorials 2 and 3.
#
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
# # Futher reading
# 1. [List comprehension](https://www.pythonforbeginners.com/basics/list-comprehensions-in-python#:~:text=List%20comprehensions%20provide%20a%20concise,be%20anything%2C%20meaning%20you%20can)
# 1. [Breaking a loop](https://www.digitalocean.com/community/tutorials/how-to-use-break-continu# %% [md]
# # Problem Set 2
#
# Questions are based on tutorials 2 and 3.
#
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
# # Futher reading
# 1. [List comprehension](https://www.pythonforbeginners.com/basics/list-comprehensions-in-python#:~:text=List%20comprehensions%20provide%20a%20concise,be%20anything%2C%20meaning%20you%20can)
# 1. [Breaking a loop](https://www.digitalocean.com/community/tutorials/how-to-use-break-continue-and-pass-statements-when-working-with-loops-in-python-3)
# 1. Check out this [guide](https://www.geeksforgeeks.org/check-multiple-conditions-in-if-statement-python/) to compare multiple conditions.
# 1. [Alpha diversity](https://en.wikipedia.org/wiki/Diversity_index)
# 1. Read this page on how to loop through columns instead of rows of a 2D numpy array
# %% [markdown]
# # Questions
# %%
# 1. Loop through samples to print the relative abundance of K.pneumoniae
relative_abundance.T # Tranposes the data, i.e. flip rows to columns
for sample in relative_abundance.T:
    print(sample[2])

# %%
# 1. If sample 1's *E.coli's* relative abundance is greater than 0.5, print('E.coli's rel abd is greater than 0.5')
for sample in relative_abundance.T:
    if sample[0] > 0.5:
        print("E.coli's rel abd is greater than 0.5")

# Try greater than 0.2
for sample_index,sample in enumerate(relative_abundance.T):
    if sample[0] > 0.2:
        print("Sample {}'s E.coli's rel abd is {:.3f}, which is greater than 0.2".format(sample_index, sample[0]))

# %%
# 1. If sample 2's *E.coli's* relative abundance is lower than 0.5 AND *H.sapiens'* rel abd is lower than 0.3, print('Both conditions are satisfied')

for sample in relative_abundance.T:
    if sample[0] < 0.5 and sample[3] < 0.3:
        print('Both conditions are satisfied')

# %%
# 1. Find all samples where *H.sapiens* relative abundance is greater than 0.2.
sample_list = []
for sample in relative_abundance.T:
    if sample[3] > 0.2:
        sample_list.append(sample)

sample_list = np.array(sample_list)

# %%
# 1. Find the ratio of *E.coli* / *C.acnes* for the first 10 samples.
ratio = []
for sample in relative_abundance.T[0:10]:
    ratio.append(sample[0]/sample[1])

# %%
# 1. Print rel abd of *C.acnes* starting from sample 0 until *C.acnes* > 0.3
# Method 1
count = 0
while relative_abundance.T[count,1] <= 0.3:
    print(relative_abundance[count, 1])
    count += 1

# Method 2
for sample in relative_abundance.T:
    if sample[1] > 0.3:
        break
    print(sample[1])

# %%
# 1. Find geometric mean of all the samples
def geometric_mean(input):
    cumulative_product = 1
    for i in input:
        assert i > 0
        cumulative_product = i*cumulative_product
    gm = cumulative_product**(1/len(input))
    return gm

gmean_list = []
for sample in relative_abundance.T:
    gmean_list.append(geometric_mean(sample))

gmean_list

# %%
# 1. Remove *H.sapiens* and rebalance the rest of the abundances so that they sum to 1
# Method 1
new_absolute_abundance = relative_abundance[1:]
new_relative_abundance = new_absolute_abundance / new_absolute_abundance.sum(axis=0)

# Method 2
new_relative_abundance = []
for sample in relative_abundance.T:
    new_sample = sample[1:]
    new_sample = new_sample/new_sample.sum()
    new_relative_abundance.append(new_sample)

new_relative_abundance = np.array(new_relative_abundance)

# %%
# 1. Find the alpha diversity of each sample using the Shannon index
def shannon(input):
    sum_plogp = 0
    for p in input:
        assert p > 0
        sum_plogp += -p*np.log(p)
    return np.array(sum_plogp)

shannon_list = []
for sample in relative_abundance.T:
    shannon_list.append(shannon(sample))

np.array(shannon_list)

# %%
# 1. Find the alpha diversity of each sample using the Simpson index

def simpson(input):
    sum_pp = 0
    for p in input:
        assert p > 0
        sum_pp += p*p
    return np.array(sum_pp)

simpson_list = []
for sample in relative_abundance.T:
    simpson_list.append(simpson(sample))

np.array(simpson_list)

# %%
# 1. Use list comprehension to obtain a list of the squares of the first 10 positive integers.

list_first_10_int = [i**1 for i in range(1,11)]

# %%
# 1. Write a Python function that that prints out the first n rows of Pascal's triangle.
#     Expected output for 6 rows:
#     ```
#     [1]
#     [1, 1]
#     [1, 2, 1]
#     [1, 3, 3, 1]
#     [1, 4, 6, 4, 1]
#     [1, 5, 10, 10, 5, 1]
#     ```
def pascal(n):
    pascal_row = [1]
    for i in range(n):
        print(pascal_row)
        old_pascal_row = list(pascal_row)
        pascal_row = np.array([0] + old_pascal_row) + np.array(old_pascal_row + [0])

pascal(6)

# %%
# 1. Given a nucleotide sequence:
#     ATGCAAGCGGATATGCATGGAAAACTTCACGCTGCCTTAGAAGATGGTTTCTTCCTCTTTTTTTGAACAGCAACAACAACCTAACATTTATTATGACACAACCACCGATCAAGAAGAC
#
#     Write a function count_substring(dna, substring), to count how many times a certain substring appears in the string of nucleotides.
#     For example, the function returns 6 when called with the substring 'CAA'and the nucleotide string above.

dna = 'ATGCAAGCGGATATGCATGGAAAACTTCACGCTGCCTTAGAAGATGGTTTCTTCCTCTTTTTTTGAACAGCAACAACAACCTAACATTTATTATGACACAACCACCGATCAAGAAGAC'
substring = 'CAA'
def count_substring(dna, substring):
    count = 0
    for i in range(len(dna)-len(substring)):
        if dna[i:i+len(substring)] == substring:
            count += 1
    return count

count_substring(dna, substring)

# %%
# 1. Write a Python function that returns a list of boolean where the Shannon index > Simpson index in the relative abundance data
#
#     Expected output format:
#     ```
#     [True, False, True, True, True, False...]
#     ```
def check_shannon_simpson(relative_abundance):
    boolean_list = []
    for sample in relative_abundance.T:
        boolean_list.append(shannon(sample) > simpson(sample))
        print('{:.3f}, {:.3f}'.format(shannon(sample), simpson(sample)))

    return boolean_list

check_shannon_simpson(relative_abundance)
e-and-pass-statements-when-working-with-loops-in-python-3)
# 1. Check out this [guide](https://www.geeksforgeeks.org/check-multiple-conditions-in-if-statement-python/) to compare multiple conditions.
# 1. [Alpha diversity](https://en.wikipedia.org/wiki/Diversity_index)
# 1. Read this page on how to loop through columns instead of rows of a 2D numpy array
# %% [markdown]
# # Questions
# %%
# 1. Loop through samples to print the relative abundance of K.pneumoniae
for sample in relative_abundance.T:
    print(sample[2])

# %%
# 1. If sample 1's *E.coli's* relative abundance is greater than 0.5, print('E.coli's rel abd is greater than 0.5')
for sample in relative_abundance.T:
    if sample[0] > 0.5:
        print("E.coli's rel abd is greater than 0.5")

# Try greater than 0.2
for sample_index,sample in enumerate(relative_abundance.T):
    if sample[0] > 0.2:
        print("Sample {}'s E.coli's rel abd is {:.3f}, which is greater than 0.2".format(sample_index, sample[0]))

# %%
# 1. If sample 2's *E.coli's* relative abundance is lower than 0.5 AND *H.sapiens'* rel abd is lower than 0.3, print('Both conditions are satisfied')

for sample in relative_abundance.T:
    if sample[0] < 0.5 and sample[3] < 0.3:
        print('Both conditions are satisfied')

# %%
# 1. Find all samples where *H.sapiens* relative abundance is greater than 0.2.
sample_list = []
for sample in relative_abundance.T:
    if sample[3] > 0.2:
        sample_list.append(sample)

sample_list = np.array(sample_list)

# %%
# 1. Find the ratio of *E.coli* / *C.acnes* for the first 10 samples.
ratio = []
for sample in relative_abundance.T[0:10]:
    ratio.append(sample[0]/sample[1])

# %%
# 1. Print rel abd of *C.acnes* starting from sample 0 until *C.acnes* > 0.3
# Method 1
count = 0
while relative_abundance.T[count,1] <= 0.3:
    print(relative_abundance[count, 1])
    count += 1

# Method 2
for sample in relative_abundance.T:
    if sample[1] > 0.3:
        break
    print(sample[1])

# %%
# 1. Find geometric mean of all the samples
def geometric_mean(input):
    cumulative_product = 1
    for i in input:
        assert i > 0
        cumulative_product = i*cumulative_product
    gm = cumulative_product**(1/len(input))
    return gm

gmean_list = []
for sample in relative_abundance.T:
    gmean_list.append(geometric_mean(sample))

gmean_list

# %%
# 1. Remove *H.sapiens* and rebalance the rest of the abundances so that they sum to 1
# Method 1
new_absolute_abundance = relative_abundance[1:]
new_relative_abundance = new_absolute_abundance / new_absolute_abundance.sum(axis=0)

# Method 2
new_relative_abundance = []
for sample in relative_abundance.T:
    new_sample = sample[1:]
    new_sample = new_sample/new_sample.sum()
    new_relative_abundance.append(new_sample)

new_relative_abundance = np.array(new_relative_abundance)

# %%
# 1. Find the alpha diversity of each sample using the Shannon index
def shannon(input):
    sum_plogp = 0
    for p in input:
        assert p > 0
        sum_plogp += -p*np.log(p)
    return np.array(sum_plogp)

shannon_list = []
for sample in relative_abundance.T:
    shannon_list.append(shannon(sample))

np.array(shannon_list)

# %%
# 1. Find the alpha diversity of each sample using the Simpson index

def simpson(input):
    sum_pp = 0
    for p in input:
        assert p > 0
        sum_pp += p*p
    return np.array(sum_pp)

simpson_list = []
for sample in relative_abundance.T:
    simpson_list.append(simpson(sample))

np.array(simpson_list)

# %%
# 1. Use list comprehension to obtain a list of the squares of the first 10 positive integers.

list_first_10_int = [i**1 for i in range(1,11)]

# %%
# 1. Write a Python function that that prints out the first n rows of Pascal's triangle.
#     Expected output for 6 rows:
#     ```
#     [1]
#     [1, 1]
#     [1, 2, 1]
#     [1, 3, 3, 1]
#     [1, 4, 6, 4, 1]
#     [1, 5, 10, 10, 5, 1]
#     ```
def pascal(n):
    pascal_row = [1]
    for i in range(n):
        print(pascal_row)
        old_pascal_row = list(pascal_row)
        pascal_row = np.array([0] + old_pascal_row) + np.array(old_pascal_row + [0])

pascal(6)

# %%
# 1. Given a nucleotide sequence:
#     ATGCAAGCGGATATGCATGGAAAACTTCACGCTGCCTTAGAAGATGGTTTCTTCCTCTTTTTTTGAACAGCAACAACAACCTAACATTTATTATGACACAACCACCGATCAAGAAGAC
#
#     Write a function count_substring(dna, substring), to count how many times a certain substring appears in the string of nucleotides.
#     For example, the function returns 6 when called with the substring 'CAA'and the nucleotide string above.

dna = 'ATGCAAGCGGATATGCATGGAAAACTTCACGCTGCCTTAGAAGATGGTTTCTTCCTCTTTTTTTGAACAGCAACAACAACCTAACATTTATTATGACACAACCACCGATCAAGAAGAC'
substring = 'CAA'
def count_substring(dna, substring):
    count = 0
    for i in range(len(dna)-len(substring)):
        if dna[i:i+len(substring)] == substring:
            count += 1
    return count

count_substring(dna, substring)

# %%
# 1. Write a Python function that returns a list of boolean where the Shannon index > Simpson index in the relative abundance data
#
#     Expected output format:
#     ```
#     [True, False, True, True, True, False...]
#     ```
def check_shannon_simpson(relative_abundance):
    boolean_list = []
    for sample in relative_abundance.T:
        boolean_list.append(shannon(sample) > simpson(sample))
        print('{:.3f}, {:.3f}'.format(shannon(sample), simpson(sample)))

    return boolean_list

check_shannon_simpson(relative_abundance)
