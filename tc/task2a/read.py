# -*- coding:utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

rho_group = np.load('rho_group.npy')
T = 1401
tau = 0.01

sum1 = 0
for i in range(318):
	sum1 += rho_group[i][T - 1]
	if(rho_group[i][T - 1] > 0.005):
		print("{}: {}".format(i, rho_group[i][T - 1]))
print(sum1)

x_mesh = np.linspace(0, (T - 1) * tau, T)
for i in range(318):
	plt.plot(x_mesh, rho_group[i])
"""	if(i == 0):
		plt.plot(x_mesh, rho_group[i], label = "$|00000>|-1,-1,1>|-1,-1,1>$")
	if(i == 15):
		print("*")
#		plt.plot(x_mesh, rho_group[i], label = "$|00000>|-1,-1,-1>|-1,-1,1>$")
	elif(i == 19):
		print("*")
#		plt.plot(x_mesh, rho_group[i], label = "$|00000>|-1,-1,1>|-1,-1,-1>$")
	elif(i == 25):
		print("*")
#		plt.plot(x_mesh, rho_group[i], label = "$|00000>|-1,-1,-1>|-1,-1,-1>$")
	elif(i == 29):
		print("*")
#plt.plot(x_mesh, rho_group[i], label = "$|00000>|1,-1,-1>|-1,-1,-1>$")
	elif(i == 32):
		print("*")
#		plt.plot(x_mesh, rho_group[i], label = "$|00000>|-1,1,-1>|-1,-1,-1>$")
	elif(i == 38):
		plt.plot(x_mesh, rho_group[i], label = "$|00000>|-1,-1,-1>|1,-1,-1>$")
	elif(i == 44):
		print("*")
#		plt.plot(x_mesh, rho_group[i], label = "$|00000>|1,-1,-1>|-1,-1,-1>$")
	elif(i == 50):
		plt.plot(x_mesh, rho_group[i], label = "$|00000>|-1,1,-1>|-1,-1,-1>$")
	elif(i == 61):
		print("*")
#		plt.plot(x_mesh, rho_group[i], label = "$|00000>|-1,-1,-1>|1,-1,-1>$")
	elif(i == 69):
		print("*")
#		plt.plot(x_mesh, rho_group[i], label = "$|00000>|-1,-1,-1>|-1,-1,1>$")
	elif(i == 73):
		print("*")
#		plt.plot(x_mesh, rho_group[i], label = "$|00000>|-1,-1,1>|-1,-1,-1>$")
	elif(i == 83):
		print("*")
#		plt.plot(x_mesh, rho_group[i], label = "$|00000>|-1,-1,-1>|-1,-1,-1>$")
	elif(i == 93):
		print("*")
#		plt.plot(x_mesh, rho_group[i], label = "$|00000>|1,-1,-1>|-1,-1,-1>$")
	elif(i == 97):
		print("*")
#		plt.plot(x_mesh, rho_group[i], label = "$|00000>|-1,1,-1>|-1,-1,-1>$")
	elif(i == 107):
		print("*")
#		plt.plot(x_mesh, rho_group[i], label = "$|00000>|-1,-1,-1>|1,-1,-1>$")
	else:
		plt.plot(x_mesh, rho_group[i])
"""
#plt.legend()
plt.ylim(0, 0.4)
plt.xlabel("$\\omega t$")
plt.ylabel("$\\rho(t)$")
#plt.title("$InitialState:|00000>|-1,-1,1>|-1,-1,1>, g_{\\Omega\\uparrow}=g_{\\Omega\\downarrow}=4, g_{\\omega\\uparrow}=g_{\\omega\\downarrow}=3, g_{tun\\uparrow}=g_{tun\\downarrow}=1, g_{\\omega_{spin}}=2, \\gamma_{\\Omega\\uparrow}=\\gamma_{\\Omega\\downarrow}=\\gamma_{\\omega\\uparrow}=\\gamma_{\\omega\\downarrow}=\\gamma_{\\omega_{spin}}=1$")
plt.grid()
plt.show()

