# %% [markdown]
# # Tutorial 2 - Learning loops and conditionals?
# Suppose we have a gene. How can we find its reverse complement

# %%
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
import numpy as np
import random

# %%
# Using the biopython package
my_seq = Seq("AGTACGCAACTGGTA")
my_seq_revcomp = my_seq.reverse_complement()

# Coding it from scratch
gene1 = "AGTACGCAACTGGTA"
gene1_mutate_pos1 = list(gene1)
gene1_mutate_pos1[0] = 'T'
gene1_mutate_pos1[1] = 'C'
# so on and so forth

# %% [markdown]
# What if I want python to change all the A to T automatically?

# # For Loops
# %%
# Looping a list
for i in [1,2,3,4,5]:
    print('Loop number ith')
    print(i)

# Looping a range
for j in range(5):
    print('Loop number jth')
    print(i)

# Looping a cell
for ind in (1,2,3,4,5):
    print('Loop number jth')
    print(i)

# Looping a dictionary (more complicated)
for key, value in {'a':1,'b':2,'c':3}.items():
    print('Loop key {}'.format(key))
    print(value)

# Using for loops to add
a = 0
for i in range(10):
    #a = a+1
    a += 1

print(a)

# Fibonacci series
a = [1, 1]
for i in range(10):
    a.append(a[-1]+a[-2])

print(a)

# Adding 5 telomeres
gene1_telomere = gene1
telomere_rep = 'TTAGGG'
for i in range(5):
    gene1_telomere = gene1_telomere + telomere_rep

print(gene1_telomere)

# %% [markdown]
# # While Loops
# Sometimes, I want the loop to run until a certain condition is met

# %%
a = 0
while a < 10: # a < 10 is a conditional
    a += 1
    print('Current count: {}'.format(count))

# Fibonacci series
a = [1, 1]
count = 0
while count < 10:
    a.append(a[-1]+a[-2])
    count += 1

print(a)

# Adding 5 telomeres
count = 0
gene1_telomere = gene1
telomere_rep = 'TTAGGG'
while count < 5:
    gene1_telomere = gene1_telomere + telomere_rep
    count += 1

print(gene1_telomere)

# %% [markdown]
# # Conditionals

# %%
a = 1
b = 2

a == 1
b == 2
a == b

c = a + b
c == 3
c < 4
c >=3
c > 3

# %%
a = random.randint(0,10)
if (a % 2) == 0:
    print('{} is even'.format(a))
else:
    print('{} is odd'.format(a))

# %%
# Reverse complement
gene1 = "AGTACGCAACTGGTA"
gene1_mut = []
for i in gene1:
    if i == 'A':
        gene1_mut.append('T')
    elif i == 'T':
        gene1_mut.append('A')
    elif i == 'C':
        gene1_mut.append('G')
    elif i == 'G':
        gene1_mut.append('C')

gene1_mut = gene1_mut[::-1] #Reverse
gene1_reversecomplement = ''.join(gene1_mut)
print(gene1_reversecomplement)
