# %% [markdown]
# # Tutorial on other types of figures
# %%
# Set up
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram

# %% [markdown]
# For this tutorial, we are using the iris test dataset, an example that is commonly used to demonstrate how python functions work. The dataset contains the sepal length	and width, petal length and width of various species. It is a great dataset to experiment with pandas.
# %% Data Setup
# Import iris dataset
sns.set_style("white")
iris = sns.load_dataset("iris")
iris_covariance = iris.cov() # Calculate the covariance
iris_correlation = iris.corr() # Calculate the correlation

# %% [markdown]
# # Heatmap
# Heatmap is a great way to illustrate 2D distance matrices or spatial data. One common usage is to represent covariance and correlation visually. In the following example, we explore how to plot a heatmap. Try running the commands sequentially to understand how commande like `ax.set_xticks` alter the graph.

# %% Plotting a Heatmap
# Method 1 - using matplotlib
fig = plt.figure() # Create figure
ax = fig.add_subplot(111) # Create axes
ax.imshow(iris_correlation, interpolation='nearest') # Makes 2D plot. 
ax.set_xticks(np.arange(len(iris_covariance.columns))) # 
ax.set_xticklabels(iris_covariance.columns)
ax.set_yticks(np.arange(len(iris_covariance.columns)))
ax.set_yticklabels(iris_covariance.columns)
ax.tick_params(axis='x',labelrotation=90) # Rotate X ticks by 90 degrees so that they don't overlap.
fig.show()

# Method 2 - using seaborn
# Notice that seaborn's function is a lot cleaner
sns.heatmap(iris_correlation)

# %% [markdown]
# # Dendrogram
#
# Once we obtain a distance-like measure between different variables (such as correlation and covariance), we perform hierarchical clustering to group variables that are closely related to each other. We use a dendrogram to represent that cluster. You can learn more about dendrograms from this [link](https://joernhees.de/blog/2015/08/26/scipy-hierarchical-clustering-and-dendrogram-tutorial/).

# %% Plotting a Dendogram
# Method 1 - matplotlib
distance_vec = pdist(iris_correlation) # Generate distance vector
dist = linkage(distance_vec, 'single') # Perform hierarchical clustering. 
# Outputs clusters defined by the distance (dist) between subclusters (idx1, idx2) 
# and the number of samples in the cluster.
# [idx1, idx2, dist, samples]
fig = plt.figure()
ax = fig.add_subplot(111)
dn = dendrogram(dist, labels=iris_correlation.columns, ax=ax) # Plot dendrogram
ax.tick_params(axis='x',labelrotation=90) 

# %% Plotting a Heatmap with a Dendogram
# Method 2 - seaborn
#sns.set_theme(color_codes=True)
iris = sns.load_dataset("iris")
species = iris.pop("species")
g = sns.clustermap(iris)

# %% [markdown]
# # 3D plot
# While matplotlib has a way to plot 3D figures, it lacks user interaction and therefore not very useful for data exploration. Here, we use the plotly package to plot a 3D scatter plot instead. Downloading the plotly package and try the following code.

# %% Plotting a 3D plot
import plotly.express as px
# px is the plotly equivalent of plt from matplotlib.
iris = sns.load_dataset("iris") # Loads the same iris data
fig = px.scatter_3d(iris, x='sepal_length', y='sepal_width', z='petal_width',
              color='species') # Run plotly scatter. 
fig.show()

# %% [markdown]
# # Summary
# These 2 tutorials are meant to showcase the types of figures you can plot on python. As we have learnt, there are usually multiple ways to plot the same plot, some easier than others, some more customizatble than others, some that you can interact with etc. Depending on the libraries you use, you may need to preprocess the data such as transforming the data from wide to long form. There is so many parameters in a figure that you can customize (like annotation, legends position, font, font sizes, colors, figure size etc.) that it is easy to get lost while striving for perfection. Until you are familiar with matplotlib and its children libraries like seaborn, always look for and adapt existing code that are available online as a first attempt. Unless you are creating a new type of visualization, most types of chart have already been coded for. 
#
# Good luck!
