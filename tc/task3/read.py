# -*- coding:utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

rho_group = np.load('rho_group.npy')
T = 801
tau = 0.02

for i in range(2235):
	if(rho_group[i][T - 1] > 0.005):
#	if(rho_group[i][T - 1] < 0.005 and rho_group[i][T - 1] > 0.0005):
		print("[{}]: {}".format(i, rho_group[i][T - 1]))

x_mesh = np.linspace(0, (T - 1) * tau, T)

#plt.plot(x_mesh, rho_group[0])
for i in range(2235):
	plt.plot(x_mesh, rho_group[i])
#	if(rho_group[i][T - 1] > 0.005):
#	if(rho_group[i][T - 1] < 0.005 and rho_group[i][T - 1] > 0.0005):
#		plt.plot(x_mesh, rho_group[i])
plt.ylim(0., 0.05)
plt.xlabel("$\\omega t$")
plt.ylabel("$\\rho(t)$")
#plt.title("$G=4, g=3, g_{tun}=2, g_{spin}=1, \\gamma_{\\omega}=1, \\gamma_{\\Omega}=1, \\gamma_{\\omega_{spin}}=1$")
plt.grid()
plt.show()

