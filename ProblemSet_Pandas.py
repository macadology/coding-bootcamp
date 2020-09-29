# %% [markdown]
# # Problem Set - Pandas
# To learn pandas, it is easiest to try out some problems on a dataset. Instead of the usual relative abundance dataset, we will use some sample datasets from seaborn.

# For additional plotting exercises, check out https://www.w3resource.com/python-exercises/pandas/index.php

# %%
# Set up
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 10) # Set max number of rows displayed
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram
from pathlib import Path

# %% Importing data from excel into pandas dataframe
relative_abundance_data = pd.read_excel('./rel_abd.xlsx') # Import relative abundance data from excel

# %%
# 1. Convert each of the following data to a pandas Series
## name the index a,b,c,d,e
## name the series value_A and value_B

value_A = {'a': 100, 'b': 200, 'c': 300, 'd': 400, 'e': 800}
value_B = np.array([10, 20, 30, 40, 50])
A_series =
B_series =

# 2. Convert the following into a dataframe
## Name the columns A,B and C
np.random.seed(1234)
n = 20
relabdA = abs(0.33+np.random.uniform(size=n))
relabdB = abs(0.33+np.random.uniform(size=n))
relabdC = 1 - relabdA - relabdB
relabd_df =

# 3. Using relabd_df, create a 4th (mean) and 5th (std) column
## Name the columns M and S.

# 4a. Select the rows where the mean is greater than 0.4
# 4b. Count the number of rows where the mean is greater than 0.4

# 5. Import the iris dataframe. When testing a computer model, it is useful to test with a subset of the data. Write a function that generates a random subset of the iris dataset containing 100 samples.
iris = sns.load_dataset("iris")
def gen_subset(dataframe, n=100):
    """
    Generate random subset of n samples
    """
    # Write function
    return dataframe_subset

iris_subset = gen_subset(iris)

# 6. Import the titanic dataframe. Explore the dataset. Try info() and describe().
# sibsp (number of siblings and children)
# parch (number of parents and children)
titanic = sns.load_dataset("titanic")
titanic

# 6a. How many males and females are there?
## (Hint) Use .groupby(), conditionals and/or .count()

male_count =
female_count =

# 6b. How many different classes and embark_town are there?
## (Hint) Use the method .unique()

number_of_classes
number_of_embarktowns

# 6c. Sort the data by fare paid from both smallest to largest and largest to smallest.
## (Hint) https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.sort_values.html

titanic_s2l =
titanic_l2s =

# 6d. Who is the youngest and oldest among those who survived and those who didn't?
survived_oldest
survived_youngest
dead_oldest
dead_youngest

# 6e. Check if the alive column and survived column are one and the same. If not, identify the index of the discrepancy.
## (Hint) Convert survived column to ones and zeros. use a==b conditional to check if the values are the same.


# 6f. Plot a histogram of the fares paid for fares lower than $100. Use 30 bins
## (Hint) https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.plot.hist.html


# 6e. Find the correlation between columns and plot a heatmap of the correlation matrix. Are there any variables that are correlated to survival?
## https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.corr.html


# 6f. Find the rows containing NaN. In the deck column, fill all NaN cells with letter 'Z'
