# %% [markdown]
# # Tutorial2 - Plotting correlations, statistics, errors etc.
# # Matplotlib
# %%
# Set up
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram

# %% Data Setup
sns.set_style("white")
iris = sns.load_dataset("iris")
iris_covariance = iris.cov()
iris_correlation = iris.corr()

# %% Plotting a Heatmap
# Method 1
fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow(iris_covariance, interpolation='nearest', extent=[-1,1,-1,1])
ax.set_xticks([-0.75,-0.25,0.25,0.75])
ax.set_xticklabels(iris_covariance.columns)
ax.set_yticks([-0.75,-0.25,0.25,0.75])
ax.set_yticklabels(iris_covariance.columns)
ax.tick_params(axis='x',labelrotation=90)
fig.show()

# Method 2
sns.heatmap(iris_correlation)



# %% Plotting a Dendogram
# [](https://joernhees.de/blog/2015/08/26/scipy-hierarchical-clustering-and-dendrogram-tutorial/)
distance_vec = pdist(iris_correlation)
dist = linkage(distance_vec, 'single')
# [idx1, idx2, dist, samples]
fig = plt.figure()
ax = fig.add_subplot(111)
dn = dendrogram(dist, labels=iris_correlation.columns, ax=ax)
ax.tick_params(axis='x',labelrotation=90)

# %% Plotting a Heatmap with a Dendogram
#sns.set_theme(color_codes=True)
iris = sns.load_dataset("iris")
species = iris.pop("species")
g = sns.clustermap(iris)

# %% Plotting a 3D plot
import plotly.express as px
#df = px.data.iris()
iris = sns.load_dataset("iris")
fig = px.scatter_3d(iris, x='sepal_length', y='sepal_width', z='petal_width',
              color='species')
fig.show()
