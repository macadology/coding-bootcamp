# %% [markdown]
# # Tutorial 1 - Understanding variables
#
# ## Variables
#
# From [guru99.com](https://www.guru99.com/variables-in-python.html), a Python variable is a reserved memory location to store values. In other words, a variable in a python program gives data to the computer for processing. Every value in Python has a datatype. Different data types in Pythogn are Numbers, List, Tuple, Strings, Dictionary, etc. Variables can be declared by any name or even alphabets like a, aa, abc, etc.
#
# Below is a list of commonly-used python variables that is included in the base python installation.
# 1. Integer - An integer datatype
# 1. Float - A floating number datatype
# 1. String - A text datatype
# 1. List - A python object/container to hold different types of data
# 1. Cell - Another python object/container to hold different types of data
# 1. Dictionary - Another python object/container to hold different types of data. Each data is also linked to a key.
#
# There are other data types included in add-on packages. Given the object-oriented nature of python, you can also create a new variable type!
# %%
# Assigning variables
a = 1 #integer
b = [1,2,3,4,5] #list
d = '1,2,3,4,5' #string
e = 1.0 #float
f = (d, e) #cell
g = {'gene1':'atcg','gene2':'atcggcc'} #dictionary
h = dict([('gene1',1),('gene2',2)]) #dictionary
print('Hello World!')

# %%
# Nested variables
gene1 = 'atcgtcgacagtagctagtcagct'
gene2 = 'tgtgcccgcgca'
gene_list = [gene1, gene2]
gene_cell = (gene1, gene2)
gene_dict = {'gene1':gene1,'gene2':gene2}

# %% [markdown]
# One important skill when writing and troubleshooting code is to find out the **value, datatype and other properties** of a variable. Only when we know the variable type can we perform meaningful calculations or functions, or troubleshoot an error. For example, we can't add a string to an integer unless we first convert the string to an integer. Below are some common *functions* to parse variables:
#
# 1. `print` is the easiest way to \`print\` out the value of a variable on a console.
# 1. `type` is a way to get the data type.
# 1. `len` is a way to get the length of a list or array.
# 1. `help` retrieves the documentation of a particular variable.

# %%
# To view value of variables
a
print(a)
print(a,b,c,d,e)
print('Hello!')
print('a is {}'.format(a))

# %%
# To check what type of variable,
type(a)
type(b)

# %%
# To check what length of an array
len(b)

# %% [markdown]
# Terminal Shortcuts
# [up] Access recent command
# [down] Access recent command in reverse
# Ctrl + A Move cursor to start
# Ctrl + E Move cursor to end
# Ctrl + K Delete everything
# Ctrl + L Clear screen

# %% [markdown]
# ## Accessing and manipulating variables
# For list arrays, cells and strings, we can access their elements using indices. For dictionary, we can access its values using an appropriate key. The syntax for referencing an index is different depending on the variable type.
# %%
# Accessing variables using indices
base_pos0 = gene1[0]
base_pos5 = # Fill in the blank
base_pos10_to_12 = # Fill in the blank
gene_dict['gene1']
gene_dict['gene2']

# %% [markdown]
# Note the difference between square brackets and curved brackets. To 'call' a function, use curved brackets. To access an element, use square brackets.

# %%
# Convert variables
print(gene1)

gene1_list = list(gene1) # Convert from string to list
print(gene1_list)

gene1_str = ''.join(gene_list) # Convert from list to string
print(gene1_str)

# %%
# Change variables at position 0
gene1_mutat1 = gene1
gene1_mutat1[0] = 't' # This will not work. Look for a different solution

gene1_mutat1 = list(gene1_mutat1)
gene1_mutat1[0]
gene1_mutat1 = ''.join(gene1_mutat1)

# %%
# Combining variables
gene3 = 'gctgctgtga'
gene1n2_str = gene1 + gene2 + gene3 # Combine strings
gene1n2_list = [gene1, gene2] + [gene3] # Combine arrays

# %%
# Count variables
# %% [markdown]

num_a = gene1.count('a') # Count the number of a
num_t = gene1.count('t') # Count the number of t

# %% [markdown]
# Notice that instead of using `count(gene1,'a')`, we use `gene1.count('a')`. The former is how we call a function, the latter is how we call a method. We will discuss the distinction in a future session, but you are encourage to read about the difference online.
#
# ## Arithmetic
# %%
a = 4
b = 2
a + b
a - b
a / b
a * b
a**b
a % b
#floor
#ceil
int

# %%
# How to add 1 to all elements
a = [1,2,3,4,5]
a + 1

# %% [markdown]
# Notice that when we run `a+1`, we run into an error. A `list` cannot be used to perform arithmetic i.e. it is not designed for that. Instead, we need a vector-type variable. Python has a widely used package called `numpy` that adds that functionality.
#
# We first import the library numpy using the code `import numpy as np`. Note that there are multiple ways to import a library. For the purposes of this tutorial, we will stick to the `import ... as ...` syntax.

# %%
import numpy as np
a_array = np.array(a) # Creates a numpy.array variable, which behaves like a vector.
print(a_array)
print(a_array + 1)
