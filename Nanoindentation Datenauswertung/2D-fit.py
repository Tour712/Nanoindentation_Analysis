# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 09:53:09 2022
fit a plane through a set of x,y,z data points

@author: eim01
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# These constants are to create random data for the sake of this example
N_POINTS = 100
TARGET_X_SLOPE = 0.3
TARGET_y_SLOPE = 0.6
TARGET_OFFSET  = 5
EXTENTS = 5
NOISE = 0.5



# Create random data.
# In your solution, you would provide your own xs, ys, and zs data.
xs = [np.random.uniform(2*EXTENTS)-EXTENTS for i in range(N_POINTS)]
ys = [np.random.uniform(2*EXTENTS)-EXTENTS for i in range(N_POINTS)]
zs = []
for i in range(N_POINTS):
    zs.append(xs[i]*TARGET_X_SLOPE + \
              ys[i]*TARGET_y_SLOPE + \
              TARGET_OFFSET + np.random.normal(scale=NOISE))

# #data from modulus mapping
# xr = [0, 5, 10, 15, 20, 0, 5, 10, 15, 20] #[um]
# yr = [0, 0, 0, 0, 0, 2, 2, 2, 2, 2]   #[um]
# zr = [5.193168, 5.191017, 5.189020, 5.187237, 5.185195,	5.183112,	5.183038,	5.183049,	5.181165, 5.179137	]   #[um]


#data from modulus mapping (23.06.22)
xs = [0, 0, 0, 10, 10, 10, 20, 20, 20] #[um]
ys = [0, 10, 20, 20, 10, 0, 0, 10, 20]   #[um]
zs = [3.815748, 3.805186, 3.794780, 3.730393, 3.740052, 3.751712, 3.687332, 3.675042, 3.668683]   #[um]

#data from modulus mapping on structured SI probe (30.06.22)
xr_2 = [0, 0, 0, 0, 0, 15, 15, 15, 15, 15, 30, 30, 30, 30, 30, 45, 45, 45, 45, 45, 60, 60, 60, 60, 60] #[um]
yr_2 = [0, 15, 30, 45, 60, 60, 45, 30, 15, 0, 0, 15, 30, 45, 60, 60, 45, 30, 15, 0, 0, 15, 30, 45, 60]   #[um]
zr_2 = [3.301858, 2.855460, 2.845203, 3.272744, 2.826360, 2.712046, 3.159698, 2.735282, 2.746988, 3.194620, 3.084131, 
        3.071779, 3.063480, 3.055256, 2.609020, 2.766842, 2.944491, 2.956198, 2.532009, 2.977644, 2.865404,
        2.855142, 2.408864, 2.834606, 2.388485]    #[um]

# plot raw data
plt.figure()
ax = plt.subplot(111, projection='3d')
ax.scatter(xr_2, yr_2, zr_2, color='b')
ax.set_xlabel('x [um]\n')
ax.set_ylabel('y [um]')
ax.set_zlabel('z [um]')
#%%
# do fit
tmp_A = []
tmp_b = []
for i in range(len(xs)):
    tmp_A.append([xs[i], ys[i], 1])
    tmp_b.append(zs[i])
b = np.matrix(tmp_b).T
A = np.matrix(tmp_A)

# Manual solution
fit = (A.T * A).I * A.T * b
errors = b - A * fit
residual = np.linalg.norm(errors)

# Or use Scipy
# from scipy.linalg import lstsq
# fit, residual, rnk, s = lstsq(A, b)

print("solution: %f x + %f y + %f = z" % (fit[0], fit[1], fit[2]))
#print("errors: \n", errors)
#print("residual:", residual)


#compute slope of plane (gradient)
Z_grad = fit[0:2]
max_grad = np.linalg.norm(Z_grad)
max_deg = np.rad2deg(np.arctan(max_grad))
print('Ebenensteigung []:', max_deg)

# plot plane
xlim = ax.get_xlim()
ylim = ax.get_ylim()
X,Y = np.meshgrid(np.arange(xlim[0], xlim[1]),
                  np.arange(ylim[0], ylim[1]))
Z = np.zeros(X.shape)
for r in range(X.shape[0]):
    for c in range(X.shape[1]):
        Z[r,c] = fit[0] * X[r,c] + fit[1] * Y[r,c] + fit[2]
ax.plot_wireframe(X,Y,Z, color='k')
#tex = "Ebenengleichung:" +str(fit[0])+'*x'+str(fit[1])+'*y'+str(fit[2])+'=z'
ax.set_xlabel('x [um]\n')
ax.set_ylabel('y [um]')
ax.set_zlabel('z [um]')
ax.set_title('Untersuchung der Ausrichtung\n'+'Ebenensteigung: %f°' % (max_deg))
tex = "Ebenengleichung:" +str(fit[0])+'*x'+str(fit[1])+'*y'+str(fit[2])+'=z'
#ax.text(18, 14, 3.85, tex, fontsize=11, va='bottom')
plt.show()

#print(Z)








