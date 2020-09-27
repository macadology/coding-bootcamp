# %% [markdown]
# Tutorial - Data manipulation with pandas

# %%
# Set up
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
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
csv_filename = Path('./iris.csv')
excel_filename = Path('./iris.xlsx')

# Writing .csv files
iris.to_csv(csv_filename)
iris_csv = pd.read_csv(csv_filename)
iris_csv2 = pd.read_csv(csv_filename,index_col=0) # Ignores the index column

# Writing .xlsx files
iris.to_excel(excel_filename)
iris_excel = pd.read_excel(excel_filename)
iris_excel2 = pd.read_excel(excel_filename, index_col=0) # Ignores the index column

# %% Creating dataframe
# Create Series
stem_length = abs(np.round(3 + 0.5 * np.random.normal(size = iris.shape[0]),1))
stem_length_series = pd.Series(stem_length, name='stem_length')

# From numpy array
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
iris['sepal_width']

# Get column names
iris.columns

# Get row
iris.iloc[3] # Based on row location
iris.loc[3] # Based on name of index
type(iris.iloc[3])

# Get row names
iris.index
list(iris.index)

# Get row and column
iris.iloc[3,0]
iris.loc[3,'sepal_length']

# Get values based on conditionals
## Check the output data type
iris['sepal_length'] < 6
boolean_array = iris['sepal_length'] < 6
iris[boolean_array]
iris[boolean_array]['sepal_length']

# %% Add data
# Add columns
## Note that the input is a form of a list
iris_newcol = pd.concat([iris,stem_length_series], axis=1)

# Add rows
## Note the difference between .concat and .append
new_sample = pd.Series([5.2,3.5,1.6,0.3,'setosa'], index=['sepal_length', 'sepal_width', 'petal_length', 'petal_width', 'species'])
iris_newrow = iris.append(new_sample,ignore_index=True)

## You can combine multiple dataframes as well

# %% Manipulate data
# Creating a copy
iris_copy1 = iris.copy()
iris_copy2 = iris.copy()

# Change value
iris_copy1['species'][149]='custom1'
iris_copy1

iris_copy.loc[149,'species']='custom2'
iris_copy

# Change NaN data
nandf = pd.Series([np.nan for i in range(iris.shape[0])], name='nan_col')
iris_w_nan = pd.concat([iris,nandf], axis=1)
iris_w_nan['nan_col'].isna() #Get boolean of nan
iris_w_nan['nan_col'] = iris_w_nan['nan_col'].fillna(0.1)

# Round data
iris.round({'sepal_length': 0, 'petal_length': 0})

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
iris.mean()
iris.std()
iris['sepal_length'].var()
gmean = lambda x: np.exp(np.mean(np.log(x)))
iris['sepal_length'].apply(gmean)

# Fast exploration
iris.info()
iris.describe()

# Groupby
iris_gb = iris.groupby('species')
iris_gb.groups
iris_gb.get_group('setosa')
iris_gb['sepal_length'].mean()
iris_gb.describe()
