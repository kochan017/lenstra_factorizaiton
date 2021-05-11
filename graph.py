'''
Purpose: Visualization of elliptic curves
@author: Ko Yat Chan
@last updated: 5/11/2021
'''

import numpy as np
import matplotlib.pyplot as plot
# set range of x and y
y, x = np.ogrid[-10:10:1000j, -10:10:100j]
# equation 1
eq1 = pow(y, 2) - pow(x, 3) - 1 # y^2 = x^3 + 1
# equation 2
eq2 = pow(y, 2) - pow(x, 3) + 1 # y^2 = x^3 - 1
# equation 3
eq3 = pow(y, 2) - pow(x, 3) + x # y^2 = x^3 - x
# equation 4
eq4 = pow(y, 2) - pow(x, 3) + 4 * x  # y^2 = x^3 - 4x
# set the line of the equation
plot.contour(x.ravel(), y.ravel(), eq4, [0])
# set the grid of the graph
plot.grid(alpha=.5, linestyle='--')
# set title based on which equation we use in plot.contour
# plot.title("y^2 = x^3 + 1")
# plot.title("y^2 = x^3 - 1")
# plot.title("y^2 = x^3 - x")
plot.title("y^2 = x^3 - 4x")
# plot x,y axis label
plot.xlabel('x')
plot.ylabel('y')
# plot legend
plot.legend()
# show the plot
plot.show()