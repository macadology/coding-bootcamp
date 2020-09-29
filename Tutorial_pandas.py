# %% [markdown]
# Tutorial - Data manipulation with pandas

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

# %% Data Setup
# Import iris dataset
sns.set_style("white")
iris = sns.load_dataset("iris")

sepal_length = iris['sepal_length']
sepal_width = iris['sepal_width']
petal_length = iris['petal_length']
petal_width = iris['petal_width']
species = iris['species']

# %% Importing exporting data
# Setting up path using pathlib
## pathlib is a library specifically for manipulating strings of directories etc. Strings work most of the time, but pathlib offers greater flexibility and functionality.
csv_filename = Path('./iris.csv')
excel_filename = Path('./iris.xlsx')

# Writing .csv files
## .csv files are basic filetypes for storing data. They are the simplest to interpret and can be read by many programs. However, csv files can be slow to read and write huge datasets. For most use cases though, csv files are more than sufficient.
## Note that to write a file, use the .to_csv() method
## to read a file, use the pd.read_csv() function
iris.to_csv(csv_filename)
iris_csv = pd.read_csv(csv_filename)
iris_csv2 = pd.read_csv(csv_filename,index_col=0) # Ignores the index column

# Writing .xlsx files
iris.to_excel(excel_filename)
iris_excel = pd.read_excel(excel_filename)
iris_excel2 = pd.read_excel(excel_filename, index_col=0) # Ignores the index column

# %% Creating dataframe
## Dataframe is to pandas what np.array is to numpy.
## Dataframe has several methods useful for data exploration.

# Create Series
stem_length = abs(np.round(3 + 0.5 * np.random.normal(size = iris.shape[0]),1))
stem_length_series = pd.Series(stem_length, name='stem_length')

# From numpy array
iris
iris_df1 = pd.DataFrame(np.array([sepal_length, sepal_width, petal_length, petal_width, species]).T, columns=['sepal_length', 'sepal_width', 'petal_length', 'petal_width', 'species'])

# From a dictionary
iris_dict = {'sepal_length': sepal_length,
'sepal_width': sepal_width,
'petal_length': petal_length,
'petal_width': petal_width,
'species': species} # Create dictionary
iris_df2 = pd.DataFrame(iris_dict)

# %% Access data
# Get column
iris
iris['species']

# Get column names
iris.columns

# Get row
iris.iloc[3] # Based on row location
iris.loc[3] # Based on name of index
type(iris.iloc[3])

# Get row names
iris.index
list(iris.index)

# Get values of a particular row and column
iris.iloc[3,0]
iris.loc[3,'sepal_length']

# Get values based on conditionals
## Check the output data type
iris['sepal_length'] < 6
boolean_array = iris['sepal_length'] < 6
iris[boolean_array]
iris[boolean_array]['sepal_length']

# Iterate rows
for ind, row in iris.iterrows():
    print(row[4])

# %% Add data
# Add columns
## Note that the input is in the form of a list
iris_newcol = pd.concat([stem_length_series, iris], axis=1)

# Add rows
## Note the difference between .concat and .append
new_sample = pd.Series([5.2,3.5,1.6,0.3,'setosa'], index=['sepal_length', 'sepal_width', 'petal_length', 'petal_width', 'species'])
iris_newrow = pd.concat([iris, new_sample], axis=0) # Does not work
iris_newrow = iris.append(new_sample,ignore_index=True)

## You can combine multiple dataframes as well

# %% Manipulate data
##
# Creating a copy
## A copy lets you manipulate the data without affecting the original dataset.
iris_copy1 = iris.copy()
iris_copy2 = iris.copy()

# Change value
iris_copy1['species'][149]='custom1'
iris_copy1

iris_copy2.loc[149,'species']='custom2'
iris_copy2

# Change NaN data
nandf = pd.Series([np.nan for i in range(iris.shape[0])], name='nan_col')
iris_w_nan = pd.concat([iris,nandf], axis=1)
iris_w_nan['nan_col'].isna() #Get boolean of nan

iris_w_nan['nan_col'] = iris_w_nan['nan_col'].fillna(0.1)

# Round data
iris.round({'sepal_length': 1, 'petal_length': 0})

# Rename columns or rows
iris.rename(columns = {'sepal_length': 'SL', 'petal_length': 'PL'})
iris.rename(index={1:'One'})

# Wide to long form
iris_long = pd.melt(iris, id_vars=['species'], value_vars=['sepal_length','sepal_width','petal_length','petal_width'], var_name='plant_properties', value_name='perperty_values')
iris_long

# Sort by column
iris.sort_values('sepal_length')
iris.sort_values(['sepal_length','petal_length'])

# %% Explore data
# Statistics
iris.mean() # Calculate mean of all relevant columns
iris.std() # Calculate standard deviation of all relevant columns
iris['sepal_length'].var() # Calculate variance of a single column
gmean = lambda x: np.exp(np.mean(np.log(x))) # Geometric mean
gmean(iris['sepal_length']) # Calculate gmean of a particular column
iris_without_species = iris.iloc[:,0:-1]
iris_without_species.apply(gmean, axis=0) # Calculate gmean of columns
iris_without_species.apply(gmean, axis=1) # Calculate gmean of rows


# Fast exploration
## Pandas is commonly used for data exploration. .info() and .describe() are two commonly used methods to explore the data quickly.
iris.info()
iris.describe()

# Groupby
## Groupby creates groups based on categorical values of a particular column.
iris_gb = iris.groupby('species')
iris_gb.groups
iris_gb.get_group('setosa')
iris_gb['sepal_length'].mean()
iris_gb.describe()
