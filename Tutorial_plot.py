# %% md
# # Tutorial - Plotting with matplotlib and seaborn
# # Matplotlib

# %%
# Set up
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

# %% Plot a line
x = np.arange(-10,10,1)
y = x**2
fig = plt.plot(x,y)

# %% Plot multiple lines
x = np.arange(-10,10,1)
y1 = x
y2 = x**2
y3 = x**3
plt.plot(x,y1)
plt.plot(x,y2)
plt.plot(x,y3)

# %% Figure legend
plt.close('all')
plt.plot(x,y1,label=r'$x$')
plt.plot(x,y2,label=r'$x^2$')
plt.plot(x,y3,label=r'$x^3$')
plt.legend()

# %% Axis labels and title
plt.close('all')
plt.plot(x,y1,label=r'$x$')
plt.plot(x,y2,label=r'$x^2$')
plt.plot(x,y3,label=r'$x^3$')
plt.legend()
plt.xlabel('x_axis')
plt.ylabel('y_axis')
plt.title('title')

# %% Axis range
plt.close('all')
plt.plot(x,y1,label=r'$x$')
plt.plot(x,y2,label=r'$x^2$')
plt.plot(x,y3,label=r'$x^3$')
plt.legend()
plt.xlim([-5,5])
plt.ylim([-100,100])

# %% Color
plt.close('all')
plt.plot(x,y1,label=r'$x$',color='red')
plt.plot(x,y2,label=r'$x^2$', color='blue')
plt.plot(x,y3,label=r'$x^3$',color='#000000')
plt.legend()

# %% Markers and linestyle
# linestyles - https://matplotlib.org/3.1.0/gallery/lines_bars_and_markers/linestyles.html
# markers - https://matplotlib.org/3.3.1/api/markers_api.html
plt.close('all')
plt.plot(x,y1,label=r'$x$',marker='.')
plt.plot(x,y2,label=r'$x^2$',marker='*')
plt.plot(x,y3,label=r'$x^3$',linestyle='-.')
plt.legend()

# %% Save plots
plt.close('all')
plt.plot(x,y1,label=r'$x$',marker='.')
plt.plot(x,y2,label=r'$x^2$',marker='*')
plt.plot(x,y3,label=r'$x^3$',linestyle='-')
plt.legend()
plt.savefig('./tutorial_fig.png')

# %% Subplots
plt.close('all')
ax1 = plt.subplot(1,3,1) #nrows, ncols, index
ax2 = plt.subplot(1,3,2)
ax3 = plt.subplot(1,3,3)
ax1.plot(x,y1,label=r'$x$',marker='.')
ax2.plot(x,y2,label=r'$x^2$',marker='*')
ax3.plot(x,y3,label=r'$x^3$',linestyle='-')
plt.legend()


# %% Object oriented approach
# Create canvas
## Axes are contained within figures
fig = plt.figure() # Figure object
ax1 = fig.add_subplot(2,2,1) # Axes object
ax2 = fig.add_subplot(2,2,2)
ax3 = fig.add_subplot(2,2,3)
ax4 = fig.add_subplot(2,2,4)
fig.show()

# Add a line object
lin1, = ax1.plot(x,y1) # Line object
fig.show()

# Change object properties
lin1.set_color('red') # Change line color to red
ax1.set_xlabel('x_label') # Change xaxis label in axes 1
ax1.set_title('Subplot 1') # Set title for axes 1
fig.suptitle('Figure title') # Set title for figure

# Miscellaneous commands
#plt.show() # show plots
#plt.close('all') # close all plots

# %% md
# # Seaborn

# %% Seaborn
sns.set()
plt.close('all')
plt.plot(x,y1,label=r'$x$',marker='.')
plt.plot(x,y2,label=r'$x^2$',marker='*')
plt.plot(x,y3,label=r'$x^3$',linestyle='-.')
plt.legend()

# %% md
# # Where to look for more information
# https://seaborn.pydata.org/tutorial.html
# https://matplotlib.org/3.3.1/tutorials/index.html#
# # Galleries
# https://seaborn.pydata.org/examples/index.html
# https://matplotlib.org/3.1.1/gallery/index.html
# # Other python plotting libraries
# 1. plotly
# 2. ggplot2
# # The grammar of graphics
# https://towardsdatascience.com/a-comprehensive-guide-to-the-grammar-of-graphics-for-effective-visualization-of-multi-dimensional-1f92b4ed4149
