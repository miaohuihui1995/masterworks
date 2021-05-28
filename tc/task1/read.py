# -*- coding:utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

rho_group = np.load('rho_group.npy')
T = 1401
tau = 0.01

for i in range(93):
	if(rho_group[i][T - 1] > 0.01):
		print("{}: {}".format(i, rho_group[i][T - 1]))

x_mesh = np.linspace(0, (T - 1) * tau, T)
for i in range(93):
	plt.plot(x_mesh, rho_group[i])

plt.ylim(0., 1)
plt.xlabel("$\\omega t$")
plt.ylabel("$\\rho(t)$")
plt.grid()
plt.show()
