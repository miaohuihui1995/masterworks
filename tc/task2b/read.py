# -*- coding:utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

rho_group = np.load('rho_group.npy')
T = 701
tau = 0.02

sum1 = 0
for i in range(5542):
	sum1 += rho_group[i][T - 1]
	if(rho_group[i][T - 1] > 0.005):
		print("{}: {}".format(i, rho_group[i][T - 1]))
print(sum1)

x_mesh = np.linspace(0, (T - 1) * tau, T)
for i in range(5542):
	if(rho_group[i][T - 1] > 0.005):
		plt.plot(x_mesh, rho_group[i])
"""	if(i == 0):
		plt.plot(x_mesh, rho_group[i], label = "$|00000>|-1,-1,1>|-1,-1,1>$")
	elif(i == 15):
		plt.plot(x_mesh, rho_group[i], label = "$|00000>|2,-1,0>|2,-1,0>$")
	elif(i == 19):
		plt.plot(x_mesh, rho_group[i], label = "$|00000>|2,-1,0>|-1,2,0>$")
	elif(i == 23):
		plt.plot(x_mesh, rho_group[i], label = "$|00000>|2,-1,0>|-1,-1,-1>$")
	elif(i == 27):
		plt.plot(x_mesh, rho_group[i], label = "$|00000>|-1,2,0>|2,-1,0>$")
	elif(i == 30):
		plt.plot(x_mesh, rho_group[i], label = "$|00000>|-1,2,0>|-1,2,0>$")
	elif(i == 34):
		plt.plot(x_mesh, rho_group[i], label = "$|00000>|-1,2,0>|-1,-1,-1>$")
	elif(i == 38):
		plt.plot(x_mesh, rho_group[i], label = "$|00000>|-1,-1,-1>|2,-1,0>$")
	elif(i == 42):
		plt.plot(x_mesh, rho_group[i], label = "$|00000>|-1,-1,-1>|-1,2,0>$")
	elif(i == 49):
		plt.plot(x_mesh, rho_group[i], label = "$|00000>|-1,-1,-1>|-1,-1,-1>$")
	else:
		plt.plot(x_mesh, rho_group[i])
"""
plt.ylim(0, 0.4)
plt.xlabel("$\\omega t$")
plt.ylabel("$\\rho(t)$")
#plt.title("$InitialState:|00000>|-1,-1,1>|-1,-1,1>, g_{\\Omega\\uparrow}=g_{\\Omega\\downarrow}=4, g_{\\omega\\uparrow}=g_{\\omega\\downarrow}=3, g_{tun\\uparrow}=g_{tun\\downarrow}=1, g_{\\omega_{spin}}=2, \\gamma_{\\Omega\\uparrow}=\\gamma_{\\Omega\\downarrow}=\\gamma_{\\omega\\uparrow}=\\gamma_{\\omega\\downarrow}=\\gamma_{\\omega_{spin}}=1$")
plt.grid()
plt.show()

