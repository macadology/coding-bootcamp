# %% [markdown]
# # Tutorial 1 - Understanding variables

# ## Variables
# %%
# Assigning variables
a = 1 #integer
b = [1,2,3,4,5] #list
c = range(1,5) #rangeobject
d = '1,2,3,4,5' #string
e = 1.0 #float
f = (d, e) #cell
g = {'gene1':'atcg','gene2':'atcggcc'} #dictionary
h = dict([('gene1',1),('gene2',2)]) #dictionary

# %%
# Variables in variables
gene1 = 'atcgtcgacagtagctagtcagct'
gene2 = 'tgtgcccgcgca'
gene_list = [gene1, gene2]
gene_cell = (gene1, gene2)
gene_dict = {'gene1':gene1,'gene2':gene2}

# %%
# To view values of variables
a
print(a)
print(a,b,c,d,e)
print('Hello!')
print('a is {}'.format(a))

# %%
# To check what type of variable, 
type(a)
type(b)

# %% [markdown]
# ## Accessing and manipulating variables
# %%
# Indexing variables
base_pos0 = gene1[0]   
base_pos5 = # Fill in the blank
base_pos10_to_12 = # Fill in the blank

# %%


# %%
# Change variables
gene1_mutat1 = gene1
gene1_mutat1[0] = 't'

# %%
# Combining variables
gene3 = 'gctgctgtga'
gene1n2_str = gene1 + gene2 + gene3
gene1n2_list = [gene1, gene2] + [gene3]

# %%
# Count variables
num_a = gene1.count('a')
num_t = gene1.count('t')
# How to count number of 'cg' 

# %% [markdown]
# ## Arithmetric
# %%
a = 4
b = 2
a + b
a - b
a / b
a * b

# %%
# How to add 1 to all elements
a = [1,2,3,4,5]
a + 1

import numpy as np
a_array = np.array(a)
a_array + 1
