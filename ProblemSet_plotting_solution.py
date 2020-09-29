# %% [markdown]
# # Problem Set - Plotting

# %% Setup
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

# %% [markdown]
# We first import the excel file into python as a pandas object. The pandas.DataFrame object contains 5 columns, the first 4 columns corresponding to the relative abundance, the last column corresponding to the site the sample was collected from.
# If you encounter an error when reading excel files, download and install the xlrd package.

# %% Importing data from excel into pandas dataframe
relative_abundance_data = pd.read_excel('./rel_abd.xlsx') # Import relative abundance data from excel
ecoli_rel_abd = relative_abundance_data['e_coli'] # Get the relative abundance of ecoli

# Extract the relative abundances of c_acnes, k_pneumoniae and h_sapiens


# %% [markdown]
# # Figure 1. Line/Scatter plot
# We are interested in visualizing the relationship between h.sapiens and the 3 different microbes. Plot a scatter plot of h.sapiens (xaxis) vs microbes (yaxis) with three different colors for each respective microbe. I have provided links to several examples depending on the method. Feel free to use any/all of the three methods listed below


# %% Method 1 - Run plt.scatter or plt.plot using the extracted relative abundances as inputs for x and y.
#==============================================================================
# Refer to this page for an example on how to use plt.scatter. [](https://matplotlib.org/3.3.1/api/_as_gen/matplotlib.pyplot.scatter.html)
# Refer to this page for an example on how to use plt.plot. [](https://matplotlib.org/3.3.1/api/_as_gen/matplotlib.pyplot.plot.html#matplotlib.pyplot.plot)
# Here is another example using plt.scatter. [](https://matplotlib.org/3.3.1/gallery/lines_bars_and_markers/scatter_symbol.html#sphx-glr-gallery-lines-bars-and-markers-scatter-symbol-py)
# If using plt.plot, remember to hide the curve and show the marker. [](https://jakevdp.github.io/PythonDataScienceHandbook/04.02-simple-scatter-plots.html#Scatter-Plots-with-plt.plot)
#==============================================================================
ecoli_rel_abd = relative_abundance_data['e_coli']
cacnes_rel_abd = relative_abundance_data['c_acnes']
kpneumoniae_rel_abd = relative_abundance_data['k_pneumoniae']
hsapiens_rel_abd = relative_abundance_data['h_sapiens']
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.scatter(hsapiens_rel_abd, ecoli_rel_abd, label="E.coli")
ax1.scatter(hsapiens_rel_abd, cacnes_rel_abd, label="C.acnes")
ax1.scatter(hsapiens_rel_abd, kpneumoniae_rel_abd, label="K.pneum")
ax1.legend()
ax1.set_xlabel('H Sapien')
ax1.set_ylabel('Relative Abundance')



# %% Method 2 - Run plt.scatter three times using keyword strings
#==============================================================================
# Refer to this guide on how to plot using keyword strings. [](https://matplotlib.org/3.3.1/tutorials/introductory/pyplot.html#plotting-with-keyword-strings)
# The difference between method 1 and 2 is that in method 2, you don't specify the inputs to x and y explicitly with a variable. Instead, you use a keyword that corresponds to the column name of the pandas object.
#==============================================================================
ecoli_rel_abd = relative_abundance_data['e_coli']
cacnes_rel_abd = relative_abundance_data['c_acnes']
kpneumoniae_rel_abd = relative_abundance_data['k_pneumoniae']
hsapiens_rel_abd = relative_abundance_data['h_sapiens']
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.scatter('h_sapiens','e_coli', data=relative_abundance_data, label="E.coli")
ax1.scatter('h_sapiens','c_acnes', data=relative_abundance_data, label="C.acnes")
ax1.scatter('h_sapiens','k_pneumoniae', data=relative_abundance_data, label="K.pneum")
ax1.legend()
ax1.set_xlabel('H Sapiens')
ax1.set_ylabel('Relative Abundance')


# %% Method 3 - Run sns.scatterplot three times using keyword strings
#==============================================================================
# Create a plot using seaborn's version of scatterplot.
# Refer to the following link for several examples on how to plot a scatter plot [](https://seaborn.pydata.org/generated/seaborn.scatterplot.html)
# I have rearranged the data using pd.melt to convert the data from wide form to long form. Read this post to learn about the differences between long form and wide form data. [](https://towardsdatascience.com/long-wide-data-and-how-to-efficiently-plot-them-7a96887e309d)
# Tip: Use the rearranged data below to plot a scatter plot. Call the variable `rearranged_data` to explore what the data looks like after rearrangement.
#==============================================================================
rearranged_data = pd.melt(relative_abundance_data, id_vars=['sites','h_sapiens'], value_vars=['e_coli','c_acnes','k_pneumoniae'], var_name='microbes', value_name='relabd (y)')
rearranged_data = rearranged_data.rename(columns={'h_sapiens':'h_sapiens_relabd (x)'})
rearranged_data
sns.scatterplot(data=rearranged_data, x="h_sapiens_relabd (x)", y="relabd (y)", hue="microbes")



# %% [markdown]
# # Figure 2. Box plot
# Plot a boxplot to show the distribution of relative abundances of all 4 species. Each species should contain 2 boxes for the two sites.

# %% Run sns.boxplot
#==============================================================================
# Use the examples here as a guide on how to make box plots. [](https://seaborn.pydata.org/generated/seaborn.boxplot.html)
# Tip: use the long form data created below
#==============================================================================
rearranged_data = pd.melt(relative_abundance_data, id_vars=['sites'], value_vars=['e_coli','c_acnes','k_pneumoniae','h_sapiens'], var_name='microbes', value_name='relabd')
rearranged_data
sns.boxplot('microbes','relabd',hue='sites',data=rearranged_data)



# %% [markdown]
# # Figure 3. Stacked bar chart
# Plot a stacked bar graph for each sample showing the relative abundance of each species.
# Method 1 is the fastest, but methods 2 and 3 demonstrate how to create a stacked bar plot from scratch.
# Methods 2 and 3 are optional. You may attempt them if keen.

# %% Method 1 - Use df.plot.bar
#==============================================================================
# Refer to this link for an example [](https://pandas.pydata.org/docs/user_guide/visualization.html#bar-plots)
# Note 1: pandas is using matplotlib as the backend figure generator
# Note 2: the df in `df.plot.bar` refers to the variable name of the target data. In our case, replace df the name of the right variable.
#==============================================================================
ax3 = relative_abundance_data.plot.bar(stacked=True)
ax3.set_xlabel('Sample number')
ax3.set_ylabel('Relative abundance')
ax3.get_figure()

# %% Method 2 - Run sns.barplot
#==============================================================================
# Refer to this example ()[https://randyzwitch.com/creating-stacked-bar-chart-seaborn/]
#==============================================================================
sns.set_style("white")
#sns.set_context({"figure.figsize": (24, 10)})
sample_ind = relative_abundance_data.index
total_eckh = relative_abundance_data.sum(axis=1)
total_eck = relative_abundance_data[['e_coli','c_acnes','k_pneumoniae']].sum(axis=1)
total_ec = relative_abundance_data[['e_coli','c_acnes']].sum(axis=1)
total_e = relative_abundance_data['e_coli']

ax = sns.barplot(x = sample_ind, y = total_eckh, color = "red", label = 'h_sapiens')
ax = sns.barplot(x = sample_ind, y = total_eck, color = "green", label = 'k_pneu')
ax = sns.barplot(x = sample_ind, y = total_ec, color = "orange", label = 'c_acnes')
ax = sns.barplot(x = sample_ind, y = total_e, color = "blue", label = 'e_coli')
ax.legend()
ax.set_xlabel('Samples')
ax.set_ylabel('Relative Abundance')
ax

# %% Method 3 - Run plt.bar
#==============================================================================
# Refer to this example ()[https://matplotlib.org/3.1.1/gallery/lines_bars_and_markers/bar_stacked.html]
#==============================================================================
sns.set_style("white")
sns.set_context({"figure.figsize": (12, 6)})
sample_ind = relative_abundance_data.index
total_eckh = relative_abundance_data.sum(axis=1)
total_eck = relative_abundance_data[['e_coli','c_acnes','k_pneumoniae']].sum(axis=1)
total_ec = relative_abundance_data[['e_coli','c_acnes']].sum(axis=1)
total_e = relative_abundance_data['e_coli']
ind = np.arange(len(relative_abundance_data.index))    # the x locations for the groups
width = 0.9

fig = plt.figure()
b1 = plt.bar(ind, total_e, width, label = 'e_coli', figure = fig)
b2 = plt.bar(ind, relative_abundance_data['c_acnes'], width, bottom=total_e, label = 'c_acnes', figure = fig)
b3 = plt.bar(ind, relative_abundance_data['k_pneumoniae'], width, bottom = total_ec,label = 'k_pneumoniae', figure = fig)
b4 = plt.bar(ind, relative_abundance_data['h_sapiens'], width, bottom = total_eck,label = 'h_sapiens', figure = fig)
plt.legend()
plt.xlabel('Samples')
plt.ylabel('Relative Abundance')



# %% Generating Data
# The following code was used to create the excel file used above. I have included the code here for reference sake, but it is not necessary for the completion of the problem set.
#
# species = ['E.coli', 'C.acnes', 'K.pneumoniae', 'H.sapiens']
# n = 100
# site = random.randint(2,size=n) #0->environment, 1->skin
# site_name = ['environment' if i == 0 else 'skin' for i in site]
# e_coli = np.random.uniform(0,1,n) * 0.2
# c_acnes = e_coli * 3 + np.random.uniform(0,1,n) * 0.1 + np.random.uniform(0,1,n) * site * 0.5
# k_pneumoniae = np.random.uniform(0,1,n) * 0.2
# h_sapiens = k_pneumoniae * 1.5 + np.random.uniform(0,1,n) * 0.1 + np.random.uniform(0,1,n) * site * 0.5
# absolute_abundance = np.array([e_coli,c_acnes,k_pneumoniae,h_sapiens])
# relative_abundance = (absolute_abundance / absolute_abundance.sum(axis=0)).T
# relative_abundance_data = pd.DataFrame(relative_abundance, columns=['e_coli','c_acnes','k_pneumoniae','h_sapiens'])
# relative_abundance_data['sites'] = site_name
# relative_abundance_data.to_excel('./rel_abd.xlsx',index=False)
