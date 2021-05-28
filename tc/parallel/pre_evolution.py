# -*- coding:utf-8 -*-
import numpy as np
import scipy.linalg as scpl

# Программа для системы атомов с \Lambda-спектром
# Владелец: Мяо Хуэйхуэй
# Почта: miaohuihui3@gmail.com
# СКИ610, ВМК, МГУ
# Москва, Май 2021

# Входить начальное состояние и параметры, генерировать автоматически список базисных состояний, гамильтониан и декогерентые операторы, потом выполнять автоматически унитарную эволюцию
# Алгоритм возможно применяем для математического моделирования других сложных квантовых систем

# ---------- Начальное состояние и параметры----------

tau = 0.01 # Шаг времени
T = 1401 # Раз итерации

nt_ph = 5	# number of types of photon
g0 = 4.	# g_{\Omega\uparrow}	
g1 = 4.	# g_{\Omega\downarrow}
g2 = 2.	# g_{\omega\uparrow}
g3 = 2.	# g_{\omega\downarrow}
g4 = 1.	# g_{\omege_{spin}}
g5 = 1.	# g_{tun\uparrow}
g6 = 1.	# g_{tun\downarrow}
gamma0 = 1.	# \gamma_{\Omega\uparrow}
gamma1 = 1. # \gamma_{\Omega\downarrow}
gamma2 = 1. # \gamma_{\omega\uparrow}
gamma3 = 1. # \gamma_{\omega\downarrow}
gamma4 = 1. # \gamma_{\omega_{spin}}
initial_at = np.array([0, 1, 2])	# [orbit0, orbit1, orbit2]：состояние атома, у нас три орбиты; 0 - нет электрона на орбите, 1 - один электрон со спином вверх на орбите, -1 - один электрон со спином вниз на орбите и 2 - два электрона на орбите
initial_ph = np.array([0, 0, 0, 0, 0])	# [\Omega\uparrow, \Omega\downarrow, \omega\uparrow, \omega\downarrow, \omega_{spin}]: состояние фотонов, у нас итого пять типов фотонов
initial_bs = np.array([0, 0, 0, 0, 0, 0, 1, 2])	# Тензорное произведение состояний фотонов и атома, здесь это наше начальное состояние
list_bs = [] # Список базисных состояний
list_bs.append(initial_bs)

# ---------- Получение списка базисных состояний без утечки фотонов ----------

# ElectronMove:
# Функция описывает переход электрона между двумя разными орбитами(или местами)
# state: одно базисное состояние
# id1, id2: номер орбиты
# n1, n2: величина на орбите id1 или id2
# list_state: список базисных состояний
# Найти все возможные состояния, которые состояние state может переходить через переход электрона между орбитами id1 и id2, потом добавить их на список
def ElectronMove(state, seek, id1, id2, n1, n2, list_state):
	tmp = np.copy(state)
	list_state_tmp = []
	for i in range(len(list_state)):
		list_state_tmp.append(list_state[i])
	cnt = 0
	
	if(n1 != n2):
		if(n1 == 0 and n2 == 2):
			# Если нет электрона на орбите id1 и два электрона(один со спином вверх и другой со спином вниз) на орбите id2, то может быть, электрон со спином вверх переходит с орбиты id2 на орбиту id1, получить 1, -1
			tmp[id1] = 1
			tmp[id2] = -1
			if(seek != -1):
				tmp[seek] += 1
			for m in range(len(list_state_tmp)):
				if((tmp == list_state_tmp[m]).all()): # Сравнивать состояние получено со состояниями из списка, если это новое состояние, то добавить его на список
					cnt += 1
			if(cnt == 0):
				tmp1 = np.copy(tmp)
				list_state_tmp.append(tmp1)
			cnt = 0
			if(seek != -1):
				tmp[seek] -= 1

			# Или электрон со спином вниз переходит с орбиты id2 на орбиту id1, получить -1, 1
			tmp[id1] = -1
			tmp[id2] = 1
			if(seek != -1):
				tmp[seek + 1] += 1
			for m in range(len(list_state_tmp)):
				if((tmp == list_state_tmp[m]).all()):
					cnt += 1
			if(cnt == 0):
				tmp2 = np.copy(tmp)
				list_state_tmp.append(tmp2)
			cnt = 0

		elif(n1 == 2 and n2 == 0):
			# Если два электрона на орбите id1 и нет электрона на орбите id2, то может быть, электрон со спином вниз переходит с орбиты id1 на орбиту id2, получить 1, -1
			tmp[id1] = 1
			tmp[id2] = -1
			if(seek != -1):
				tmp[seek + 1] -= 1
			if(tmp[seek + 1] >= 0):
				for m in range(len(list_state_tmp)):
					if((tmp == list_state_tmp[m]).all()):
						cnt += 1
				if(cnt == 0):
					tmp1 = np.copy(tmp)
					list_state_tmp.append(tmp1)
				cnt = 0
			if(seek != -1):
				tmp[seek + 1] += 1

			# Или электрон со спином вверх переходит с орбиты id1 на орбиту id2, получить -1, 1
			tmp[id1] = -1
			tmp[id2] = 1
			if(seek != -1):
				tmp[seek] -= 1
			if(tmp[seek] >= 0):
				for m in range(len(list_state_tmp)):
					if((tmp == list_state_tmp[m]).all()):
						cnt += 1
				if(cnt == 0):
					tmp2 = np.copy(tmp)
					list_state_tmp.append(tmp2)
				cnt = 0

		elif(n1 == 1 and n2 == -1):
			# Если один электрон со спином вверх на орбите id1 и один электрон со спином вниз на орбите id2, то может быть, электрон со спином вверх переходит с орбиты id1 на орбиту id2, получить 0, 2
			tmp[id1] = 0
			tmp[id2] = 2
			if(seek != -1):
				tmp[seek] -= 1
			if(tmp[seek] >= 0):
				for m in range(len(list_state_tmp)):
					if((tmp == list_state_tmp[m]).all()):
						cnt += 1
				if(cnt == 0):
					tmp1 = np.copy(tmp)
					list_state_tmp.append(tmp1)
				cnt = 0
			if(seek != -1):
				tmp[seek] += 1

			# Или электрон со спином вниз переходит с орбиты id2 на орбиту id1, получить 2, 0
			tmp[id1] = 2
			tmp[id2] = 0
			if(seek != -1):
				tmp[seek + 1] += 1
			for m in range(len(list_state_tmp)):
				if((tmp == list_state_tmp[m]).all()):
					cnt += 1
			if(cnt == 0):
				tmp2 = np.copy(tmp)
				list_state_tmp.append(tmp2)
			cnt = 0

		elif(n1 == -1 and n2 == 1):
			# Если один электрон со спином вниз на орбите id1 и один электрон со спином вверх на орбите id2, то может быть, электрон со спином вниз переходит с орбиты id1 на орбиту id2, получить 0, 2
			tmp[id1] = 0
			tmp[id2] = 2
			if(seek != -1):
				tmp[seek + 1] -= 1
			if(tmp[seek + 1] >= 0):
				for m in range(len(list_state_tmp)):
					if((tmp == list_state_tmp[m]).all()):
						cnt += 1
				if(cnt == 0):
					tmp1 = np.copy(tmp)
					list_state_tmp.append(tmp1)
				cnt = 0
			if(seek != -1):
				tmp[seek + 1] += 1

			# Или электрон со спином вверх переходит с орбиты id2 на орбиту id1, получить 2, 0
			tmp[id1] = 2
			tmp[id2] = 0
			if(seek != -1):
				tmp[seek] += 1
			for m in range(len(list_state_tmp)):
				if((tmp == list_state_tmp[m]).all()):
					cnt += 1
			if(cnt == 0):
				tmp2 = np.copy(tmp)
				list_state_tmp.append(tmp2)
			cnt = 0

		elif((n1 == 0 and n2 == 1) or (n1 == -1 and n2 == 2)):
			# Здесь два случая: электрон со спином вверх переходит с орбиты id2 на орбиту id1, просто менять места цифр: 0, 1 --> 1, 0; -1, 2 --> 2, -1
			tmp[id1] = state[id2]
			tmp[id2] = state[id1]
			if(seek != -1):
				tmp[seek] += 1
			for m in range(len(list_state_tmp)):
				if((tmp == list_state_tmp[m]).all()):
					cnt += 1
			if(cnt == 0):
				tmp = np.copy(tmp)
				list_state_tmp.append(tmp)
			cnt = 0

		elif((n1 == 1 and n2 == 0) or (n1 == 2 and n2 == -1)):
			# Здесь два случая: электрон со спином вверх переходит с орбиты id1 на орбиту id2, просто менять места цифр: 1, 0 --> 0, 1; 2, -1 --> -1, 2
			tmp[id1] = state[id2]
			tmp[id2] = state[id1]
			if(seek != -1):
				tmp[seek] -= 1
			if(tmp[seek] >= 0):
				for m in range(len(list_state_tmp)):
					if((tmp == list_state_tmp[m]).all()):
						cnt += 1
				if(cnt == 0):
					tmp = np.copy(tmp)
					list_state_tmp.append(tmp)
				cnt = 0

		elif((n1 == 0 and n2 == -1) or (n1 == 1 and n2 == 2)):
			# Здесь два случая: электрон со спином вниз переходит с орбиты id2 на орбиту id1, просто менять места цифр: 0, -1 --> -1, 0; 1, 2 --> 2, 1
			tmp[id1] = state[id2]
			tmp[id2] = state[id1]
			if(seek != -1):
				tmp[seek + 1] += 1
			for m in range(len(list_state_tmp)):
				if((tmp == list_state_tmp[m]).all()):
					cnt += 1
			if(cnt == 0):
				tmp = np.copy(tmp)
				list_state_tmp.append(tmp)
			cnt = 0
	
		elif((n1 == -1 and n2 == 0) or (n1 == 2 and n2 == 1)):
			# Здесь два случая: электрон со спином вниз переходит с орбиты id1 на орбиту id2, просто менять места цифр: -1, 0 --> 0, -1; 2, 1 --> 1, 2
			tmp[id1] = state[id2]
			tmp[id2] = state[id1]
			if(seek != -1):
				tmp[seek + 1] -= 1
			if(tmp[seek + 1] >= 0):
				for m in range(len(list_state_tmp)):
					if((tmp == list_state_tmp[m]).all()):
						cnt += 1
				if(cnt == 0):
					tmp = np.copy(tmp)
					list_state_tmp.append(tmp)
				cnt = 0
	
	return list_state_tmp

# SpinMove:
# Функция описывает переход спинов для одного состояния
# state: одно базисное состояние
# list_state: список базисных состояний
# Переходить со спином вверх во спин вниз, испуская фотон; обратно, поглощая фотон
# Найти все возможные состояния, которые состояние state может переходить через переход спинов, потом добавить их на список
def SpinChange(nt_ph, state, list_state):
	tmp = np.copy(state)
	list_state_tmp = []
	for i in range(len(list_state)):
		list_state_tmp.append(list_state[i])
	cnt = 0
	
	for i in range(nt_ph, len(state)):
		if(tmp[i] == 1):
			tmp[i] *= -1 # Переходить со спином вверх во спин вниз
			tmp[nt_ph - 1] += 1
			for m in range(len(list_state_tmp)):
				if((tmp == list_state_tmp[m]).all()):
					cnt += 1
			if(cnt == 0):
				tmp1 = np.copy(tmp)
				list_state_tmp.append(tmp1)
			cnt = 0
			tmp[i] *= -1
			tmp[nt_ph - 1] -= 1

		elif(tmp[i] == -1):
			if(tmp[nt_ph - 1] > 0):
				tmp[i] *= -1 # Переходить со спином вниз во спин вверх
				tmp[nt_ph - 1] -= 1
				for m in range(len(list_state_tmp)):
					if((tmp == list_state_tmp[m]).all()):
						cnt += 1
				if(cnt == 0):
					tmp2 = np.copy(tmp)
					list_state_tmp.append(tmp2)
				cnt	= 0
				tmp[i] *= -1
				tmp[nt_ph - 1] += 1
			
	return list_state_tmp

# StateChange: ElectronMove + SpinChange
def StateChange(nt_ph, state, list_state):
	list_state = np.copy(ElectronMove(state, 0, nt_ph + 0, nt_ph + 2, state[nt_ph + 0], state[nt_ph + 2], list_state))
	list_state = np.copy(ElectronMove(state, 2, nt_ph + 1, nt_ph + 2, state[nt_ph + 1], state[nt_ph + 2], list_state))
	list_state = np.copy(SpinChange(nt_ph, state, list_state))
	return list_state

flag = True
length_last = 0

while(flag):
	list_bs1 = np.copy(list_bs)
	list_bs_compare = np.copy(list_bs)

	for i in range(length_last, len(list_bs)):
		list_bs1 = np.copy(StateChange(nt_ph, list_bs[i], list_bs1))

	if(len(list_bs1) == len(list_bs_compare)):
		flag = False
	else:
		length_last = len(list_bs)
		list_bs = np.copy(list_bs1)

print("List of basic states without consideration leak of photons: ")
print(list_bs)
print("Length of list: {}".format(len(list_bs)))

# Найти максимум количества свободных фотонов
max1 = 0
for i in range(len(list_bs)):
	sum1 = 0
	for j in range(nt_ph):
		sum1 += list_bs[i][j]
		if(sum1 > max1):
			max1 = sum1
print("Max number of free photon is {}".format(max1))

# ---------- Получение списка базисных состояний с утечкой фотонов ----------

list_sink = []
for i in range(len(list_bs)):
	list_sink.append(list_bs[i])
list_sink1 = []
cnts = [] # Содержать номер границ между разными подпространствами
cnts.append(0)
cnts.append(len(list_bs))

# Здесь использовать метод комбинаторики
# Например: у состояния фотонов |0, 1, 0, 1, 0> есть два фотона из двух разных типов, просто написать так |1, 1>
#			Генерировать все возможные случаи для двух мест, где сумма меньше или равна 2
#			Получить |2, 0>, |0, 2>, |1, 1>, |1, 0>, |0, 1>, |0, 0>
#			Потом |1, 1> вычитает их. Сохранять случаи, у результата которых нет отрицательного числа: |1, 1>, |1, 0>, |0, 1>, |0, 0>. |2, 0> несохранен, потому что |1, 1> - |2, 0> = |-1, 1>, там есть отрицательное число; аналогично, |0, 2> также несохранен
#			Т.е. получить |0, 1, 0, 1, 0>, |0, 1, 0, 0, 0>, |0, 0, 0, 1, 0>, |0, 0, 0, 0, 0>. Это все возможные состояния при утечке фотонв
for t in range(1, max1 + 1):
	cnt2 = 0	
	for i in range(len(list_bs)):
		sum1 = 0
		index = []
		cnt = 0
		for j in range(nt_ph):
			sum1 += list_bs[i][j]
			if(list_bs[i][j] > 0):
				index.append(j)
				cnt += 1

		if(sum1 > 0):
			arr = []
			for m0 in range(0, sum1 + 1):
				if(cnt > 1):
					for m1 in range(0, sum1 + 1):
						if(cnt > 2):
							for m2 in range(0, sum1 + 1):
								if(cnt > 3):
									for m3 in range(0, sum1 + 1):
										if(cnt > 4):
											for m4 in range(0, sum1 + 1):
												arr.append([m0, m1, m2, m3, m4])		
										else:
											arr.append([m0, m1, m2, m3])		
								else:
									arr.append([m0, m1, m2])		
						else:
							arr.append([m0, m1])		
				else:
					arr.append([m0])		

			arr1 = []
			for j in range(len(arr)):
				sum2 = 0
				for k in range(len(arr[j])):
					sum2 += arr[j][k]
				if(sum1 - sum2 == t):
					arr1.append(arr[j])

		if(sum1 > 0):
			for j in range(len(arr1)):
				bs_tmp = list_bs[i].copy()
				cnt1 = 0
				if(bs_tmp[index[0]] - arr1[j][0] >= 0):
					if(cnt > 1 and bs_tmp[index[1]] - arr1[j][1] >= 0):
						if(cnt > 2 and bs_tmp[index[2]] - arr1[j][2] >= 0):
							if(cnt > 3 and bs_tmp[index[3]] - arr1[j][3] >= 0):
								if(cnt > 4 and bs_tmp[index[4]] - arr1[j][4] >= 0):
									bs_tmp[index[4]] = arr1[j][4]
									bs_tmp[index[3]] = arr1[j][3]
									bs_tmp[index[2]] = arr1[j][2]
									bs_tmp[index[1]] = arr1[j][1]
									bs_tmp[index[0]] = arr1[j][0]
									cnt1 += 5
									if(cnt1 == len(index)):
										list_sink.append(bs_tmp)
										list_sink1 = np.copy(list_sink)
										list_sink = list(np.copy(list_sink1))
										cnt2 += 1
									continue
								else:
									bs_tmp[index[3]] = arr1[j][3]
									bs_tmp[index[2]] = arr1[j][2]
									bs_tmp[index[1]] = arr1[j][1]
									bs_tmp[index[0]] = arr1[j][0]
									cnt1 += 4
									if(cnt1 == len(index)):
										list_sink.append(bs_tmp)
										list_sink1 = np.copy(list_sink)
										list_sink = list(np.copy(list_sink1))
										cnt2 += 1
									continue
							else:
								bs_tmp[index[2]] = arr1[j][2]
								bs_tmp[index[1]] = arr1[j][1]
								bs_tmp[index[0]] = arr1[j][0]
								cnt1 += 3
								if(cnt1 == len(index)):
									list_sink.append(bs_tmp)
									list_sink1 = np.copy(list_sink)
									list_sink = list(np.copy(list_sink1))
									cnt2 += 1
								continue
						else:
							bs_tmp[index[1]] = arr1[j][1]
							bs_tmp[index[0]] = arr1[j][0]
							cnt1 += 2
							if(cnt1 == len(index)):
								list_sink.append(bs_tmp)
								list_sink1 = np.copy(list_sink)
								list_sink = list(np.copy(list_sink1))
								cnt2 += 1
							continue
					else:
						bs_tmp[index[0]] = arr1[j][0]
						cnt1 += 1
						if(cnt1 == len(index)):
							list_sink.append(bs_tmp)
							list_sink1 = np.copy(list_sink)
							list_sink = list(np.copy(list_sink1))
							cnt2 += 1
						continue
	cnts.append(cnt2)

for i in range(1, len(cnts)):
	cnts[i] += cnts[i - 1]

# Отбросить повторные состояния из списка базисных состояний с утечкой фотонов
cnt1 = 0
cnt2 = 0
cnt3 = 0
list_sink_unique = []
list_sink_unique.append(list_sink[0])
for i in range(len(list_sink)):
	cnt = 0
	for j in range(len(list_sink_unique)):
		if((list_sink[i] == list_sink_unique[j]).all()):
			if(i >= cnts[1] and i < cnts[2]):
				cnt1 += 1
			elif(i >= cnts[2] and i < cnts[3]):
				cnt2 += 1
			elif(i >= cnts[3] and i < cnts[4]):
				cnt3 += 1
			cnt += 1
	if(cnt == 0):
		list_sink_unique.append(list_sink[i])
cnts[2] -= cnt1
cnts[3] -= cnt1 + cnt2
cnts[4] -= cnt1 + cnt2 + cnt3

print("cnts: ")
print(cnts)
print("List of basic states with consideration leak of photons: ")
print(list_sink_unique)
print("Length of list: {}".format(len(list_sink_unique)))

"""
for i in range(len(list_sink)):
	cnt = 0
	for j in range(len(list_sink_unique)):
		if((list_sink[i] == list_sink_unique[j]).all()):
			cnt += 1
	if(cnt != 1):
		print("error")
"""

# ---------- Создание списка триад для построения гамильтониана ---------

# Здесь функции ParametersInElectronMove, ParametersInSpinChange, ParameterInStateChange похожи на функции ElectronMove, SpinChange, StateChange
# Идеа: выбрать одно базисное состояние, найти все возможные состояния, на которые оно может переходит через переход электрона или переход спинов. Потом сохранять его триаду
# Например: триада [i, j, g]: i - номер базисного состояния, j - номер состояния из списка может переходить, g - сила. Т.е. Величина g находится на месте (i, j) гамильтониана
# Сохранять все триады, чтобы построить гамильтониан

def ParametersInElectronMove(index, state, seek, id1, id2, n1, n2, g, gg, list_state):
	tmp = np.copy(state)
	parameters_H = []
	
	if(n1 != n2):
		if(n1 == 0 and n2 == 2):
			tmp[id1] = 1
			tmp[id2] = -1
			if(seek != -1):
				tmp[seek] += 1
			for m in range(len(list_state)):
				if((tmp == list_state[m]).all()):
					parameters_H.append(np.array([index, m, g]))
			if(seek != -1):
				tmp[seek] -= 1

			tmp[id1] = -1
			tmp[id2] = 1
			if(seek != -1):
				tmp[seek + 1] += 1
			for m in range(len(list_state)):
				if((tmp == list_state[m]).all()):
					parameters_H.append(np.array([index, m, gg]))

		elif(n1 == 2 and n2 == 0):
			tmp[id1] = 1
			tmp[id2] = -1
			if(seek != -1):
				tmp[seek + 1] -= 1
			if(tmp[seek + 1] >= 0):
				for m in range(len(list_state)):
					if((tmp == list_state[m]).all()):
						parameters_H.append(np.array([index, m, gg]))
			if(seek != -1):
				tmp[seek + 1] += 1

			tmp[id1] = -1
			tmp[id2] = 1
			if(seek != -1):
				tmp[seek] -= 1
			if(tmp[seek] >= 0):
				for m in range(len(list_state)):
					if((tmp == list_state[m]).all()):
						parameters_H.append(np.array([index, m, g]))

		elif(n1 == 1 and n2 == -1):
			tmp[id1] = 0
			tmp[id2] = 2
			if(seek != -1):
				tmp[seek] -= 1
			if(tmp[seek] >= 0):
				for m in range(len(list_state)):
					if((tmp == list_state[m]).all()):
						parameters_H.append(np.array([index, m, g]))
			if(seek != -1):
				tmp[seek] += 1

			tmp[id1] = 2
			tmp[id2] = 0
			if(seek != -1):
				tmp[seek + 1] += 1
			for m in range(len(list_state)):
				if((tmp == list_state[m]).all()):
					parameters_H.append(np.array([index, m, gg]))

		elif(n1 == -1 and n2 == 1):
			tmp[id1] = 0
			tmp[id2] = 2
			if(seek != -1):
				tmp[seek + 1] -= 1
			if(tmp[seek + 1] >= 0):
				for m in range(len(list_state)):
					if((tmp == list_state[m]).all()):
						parameters_H.append(np.array([index, m, gg]))
			if(seek != -1):
				tmp[seek + 1] += 1

			tmp[id1] = 2
			tmp[id2] = 0
			if(seek != -1):
				tmp[seek] += 1
			for m in range(len(list_state)):
				if((tmp == list_state[m]).all()):
					parameters_H.append(np.array([index, m, g]))

		elif((n1 == 0 and n2 == 1) or (n1 == -1 and n2 == 2)):
			tmp[id1] = state[id2]
			tmp[id2] = state[id1]
			if(seek != -1):
				tmp[seek] += 1
			for m in range(len(list_state)):
				if((tmp == list_state[m]).all()):
					parameters_H.append(np.array([index, m, g]))

		elif((n1 == 1 and n2 == 0) or (n1 == 2 and n2 == -1)):
			tmp[id1] = state[id2]
			tmp[id2] = state[id1]
			if(seek != -1):
				tmp[seek] -= 1
			if(tmp[seek] >= 0):
				for m in range(len(list_state)):
					if((tmp == list_state[m]).all()):
						parameters_H.append(np.array([index, m, g]))

		elif((n1 == 0 and n2 == -1) or (n1 == 1 and n2 == 2)):
			tmp[id1] = state[id2]
			tmp[id2] = state[id1]
			if(seek != -1):
				tmp[seek + 1] += 1
			for m in range(len(list_state)):
				if((tmp == list_state[m]).all()):
					parameters_H.append(np.array([index, m, gg]))
	
		elif((n1 == -1 and n2 == 0) or (n1 == 2 and n2 == 1)):
			tmp[id1] = state[id2]
			tmp[id2] = state[id1]
			if(seek != -1):
				tmp[seek + 1] -= 1
			if(tmp[seek + 1] >= 0):
				for m in range(len(list_state)):
					if((tmp == list_state[m]).all()):
						parameters_H.append(np.array([index, m, gg]))
	
	return parameters_H

def ParametersInSpinChange(index, nt_ph, state, g, list_state):
	tmp = np.copy(state)
	parameters_H = []
	
	for i in range(nt_ph, len(state)):
		if(tmp[i] == 1):
			tmp[i] *= -1
			tmp[nt_ph - 1] += 1
			for m in range(len(list_state)):
				if((tmp == list_state[m]).all()):
					parameters_H.append(np.array([index, m, g]))
			tmp[i] *= -1
			tmp[nt_ph - 1] -= 1

		elif(tmp[i] == -1):
			if(tmp[nt_ph - 1] > 0):
				tmp[i] *= -1
				tmp[nt_ph - 1] -= 1
				for m in range(len(list_state)):
					if((tmp == list_state[m]).all()):
						parameters_H.append(np.array([index, m, g]))
				tmp[i] *= -1
				tmp[nt_ph - 1] += 1
			
	return parameters_H

def ParametersInStateChange(index, nt_ph, state, list_state):
	parameters_H = []
	for i in range(len(ParametersInElectronMove(index, state, 0, nt_ph + 0, nt_ph + 2, state[nt_ph + 0], state[nt_ph + 2], g0, g1, list_state))):
		parameters_H.append(ParametersInElectronMove(index, state, 0, nt_ph + 0, nt_ph + 2, state[nt_ph + 0], state[nt_ph + 2], g0, g1, list_state)[i])
	for i in range(len(ParametersInElectronMove(index, state, 2, nt_ph + 1, nt_ph + 2, state[nt_ph + 1], state[nt_ph + 2], g2, g3, list_state))):
		parameters_H.append(ParametersInElectronMove(index, state, 2, nt_ph + 1, nt_ph + 2, state[nt_ph + 1], state[nt_ph + 2], g2, g3, list_state)[i])
	for i in range(len(ParametersInSpinChange(index, nt_ph, state, g4, list_state))):
		parameters_H.append(ParametersInSpinChange(index, nt_ph, state, g4, list_state)[i])
	return parameters_H

length_last = 0
parameters_H = []
for k in range(max1 + 1):
	list_tmp = []
	for i in range(cnts[k], cnts[k + 1]):
		list_tmp.append(list_sink_unique[i])
	for i in range(len(list_tmp)):
		for j in range(len(ParametersInStateChange(i, nt_ph, list_tmp[i], list_tmp))):
			parameters_H.append(ParametersInStateChange(i, nt_ph, list_tmp[i], list_tmp)[j])
	for i in range(length_last, len(parameters_H)):
		parameters_H[i][0] += cnts[k]
		parameters_H[i][1] += cnts[k]
		length_last = len(parameters_H)
		
# Отбросить повторные триады из списка
parameters_H_unique = []
paremeters_H_unique = np.copy(parameters_H[0])
for i in range(len(parameters_H)):
	cnt = 0
	for j in range(len(parameters_H_unique)):
		if((parameters_H[i] == parameters_H_unique[j]).all()):
			cnt += 1
	if(cnt == 0):
		if(parameters_H[i][0] < parameters_H[i][1]):
			parameters_H_unique.append(parameters_H[i])


print("List of triads for building H: ")
print(parameters_H_unique)
print("Length of list: {}".format(len(parameters_H_unique)))

# ---------- Создание списка триад для построения декогерентных операторов ---------

# Идеа: выбрать одно базисное состояние, найти все возможные состояния, на которые оно может переходит через утечку одного фотона. Потом сохранять его триаду
# Например: триада [i, j, gamma]: i - номер базисного состояния, j - номер состояния из списка может переходить. Т.е. получить декогерентный оператор A_{ij}
# Сохранять все триады, чтобы построить декогерентные операторы

parameters_L = []

for t in range(max1):
	list_tmp1 = []
	list_tmp2 = []
	for i in range(cnts[t], cnts[t + 1]):
		list_tmp1.append(list_sink_unique[i])
	for i in range(cnts[t + 1], cnts[t + 2]):
		list_tmp2.append(list_sink_unique[i])
	for i in range(len(list_tmp1)):
		sum1 = 0
		bs1 = list_tmp1[i]
		for m in range(nt_ph):
			sum1 += bs1[m]
		for j in range(len(list_tmp2)):
			cnt = 0
			index = -1
			sum2 = 0
			bs2 = list_tmp2[j]
			for n in range(nt_ph):
				sum2 += bs2[n]
			for k in range(nt_ph):
				if(bs1[k] != bs2[k]):
					cnt += 1
					index = k
			cnt1 = 0
			for k in range(nt_ph, len(initial_bs)):
				if(bs1[k] == bs2[k]):
					cnt1 += 1
			if(np.abs(sum1 - sum2) == 1 and cnt == 1 and cnt1 == len(initial_bs) - nt_ph):
				if(index == 0):
					parameters_L.append([i + cnts[t], j + cnts[t + 1], gamma0])
				elif(index == 1):
					parameters_L.append([i + cnts[t], j + cnts[t + 1], gamma1])
				elif(index == 2):
					parameters_L.append([i + cnts[t], j + cnts[t + 1], gamma2])
				elif(index == 3):
					parameters_L.append([i + cnts[t], j + cnts[t + 1], gamma3])
				elif(index == 4):
					parameters_L.append([i + cnts[t], j + cnts[t + 1], gamma4])

# Отбросить повторные триады из списка
parameters_L_unique = []
for i in range(len(parameters_L)):
	if(parameters_L[i][0] < parameters_L[i][1]):
		parameters_L_unique.append(parameters_L[i])

print("List of triads for building A: ")
print(parameters_L_unique)
print("Length of list: {}".format(len(parameters_L_unique)))


# Создать гамильтониан
H = np.eye(len(list_sink_unique))
for i in range(len(parameters_H_unique)):
	H[int(parameters_H_unique[i][0])][int(parameters_H_unique[i][1])] = parameters_H_unique[i][2]
H += H.T
for i in range(len(list_sink_unique)):
	H[i][i] = 1.

val, vec = np.linalg.eig(H)

#exp_H_p = scpl.expm(1j * H * tau)
#exp_H_m = scpl.expm(-1j * H * tau)

"""cnt1 = 0
cnt2 = 0
index_i_exp_H_p = []
index_j_exp_H_p = []
index_i_exp_H_m = []
index_j_exp_H_m = []
real_exp_H_p = []
imag_exp_H_p = []
real_exp_H_m = []
imag_exp_H_m = []
for j in range(len(list_sink_unique)):
	for i in range(len(list_sink_unique)):
		if(exp_H_p[i][j] != 0):
			cnt1 += 1
			index_i_exp_H_p.append(i)
			index_j_exp_H_p.append(j)
			real_exp_H_p.append(exp_H_p[i][j].real)
			imag_exp_H_p.append(exp_H_p[i][j].imag)
for i in range(len(list_sink_unique)):
	for j in range(len(list_sink_unique)):
		if(exp_H_m[i][j] != 0):
			cnt2 += 1
			index_i_exp_H_m.append(i)
			index_j_exp_H_m.append(j)
			real_exp_H_m.append(exp_H_m[i][j].real)
			imag_exp_H_m.append(exp_H_m[i][j].imag)
print("Length of triads for exp_H_p: {}".format(cnt1))
print("Length of triads for exp_H_m: {}".format(cnt2))"""

index_i_L = []
index_j_L = []
value_L = []
for i in range(len(parameters_L_unique)):
	index_i_L.append(parameters_L_unique[i][0])
	index_j_L.append(parameters_L_unique[i][1])
	value_L.append(parameters_L_unique[i][2])
	
index_i_H = []
index_j_H = []
value_H = []
for i in range(len(parameters_H_unique)):
	index_i_H.append(parameters_H_unique[i][0])
	index_j_H.append(parameters_H_unique[i][1])
	value_H.append(parameters_H_unique[i][2])

# Сохранять параметры из exp_H_p, exp_H_m, parameters_L_unique, list_sink_unique на файлы
#np.savetxt('index_i_exp_H_p.txt', index_i_exp_H_p, fmt='%d')
#np.savetxt('index_j_exp_H_p.txt', index_j_exp_H_p, fmt='%d')
#np.savetxt('real_exp_H_p.txt', real_exp_H_p, fmt='%0.18f')
#np.savetxt('imag_exp_H_p.txt', imag_exp_H_p, fmt='%0.18f')
#np.savetxt('index_i_exp_H_m.txt', index_i_exp_H_m, fmt='%d')
#np.savetxt('index_j_exp_H_m.txt', index_j_exp_H_m, fmt='%d')
#np.savetxt('real_exp_H_m.txt', real_exp_H_m, fmt='%0.18f')
#np.savetxt('imag_exp_H_m.txt', imag_exp_H_m, fmt='%0.18f')
np.savetxt('index_i_H.txt', index_i_H, fmt='%d')
np.savetxt('index_j_H.txt', index_j_H, fmt='%d')
np.savetxt('value_H.txt', value_H, fmt='%d')
np.savetxt('val.txt', val, fmt='%.18f')
np.savetxt('vec.txt', vec.T, fmt='%.18f')
np.savetxt('index_i_L.txt', index_i_L, fmt='%d')
np.savetxt('index_j_L.txt', index_j_L, fmt='%d')
np.savetxt('value_L.txt', value_L, fmt='%d')
np.savetxt('list_sink.txt', list_sink_unique, fmt='%d')
