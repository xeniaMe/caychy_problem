# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 15:42:29 2024

@author: Mi
"""

import numpy as np
import matplotlib.pyplot as plt

# Параметры системы
m = 1.0  # масса (кг)
L = 5.0  # длина стержня (м)
g0 = 9.81  # ускорение свободного падения (м/с^2)

# Начальные условия
x0 = 3  # начальное положение по x (м)
y0 = -4   # начальное положение по y (м)
vx0 = 1.0  # начальная скорость по x (м/с)
vy0 = 1.0  # начальная скорость по y (м/с)
T_0 = m * g0 * y0 - 0.5 * m * (vx0**2 + vy0**2)  #начальное значение T полученное из уравнений 

# Время моделирования
t_min = 0
t_max = 10
#dt = 0.01  # шаг по времени

#Для решения используется схема Розенброка
alpha = complex(1, 1) / 2 # параметр схемы Розенброка
#сетка
N = 500 #число интервалов
tau = (t_max - t_min) / N
t = np.linspace(t_min, t_max, N + 1) #время

# Правая часть системы уравнений
def F(u, g0, m, L, t):
    f = np.zeros(5)
    g = g0 + 0.05 * np.sin(2 * np.pi * t)
    f[0] = u[2]
    f[1] = u[3]
    f[2] = (u[4] * u[0]) / (m * L)
    f[3] = (u[4] * u[1]) / (m * L) - g
    f[4] = u[0]**2 + u[1]**2 - L**2
    return f

def D():
    d = np.zeros((5,5))
    for i in range(5):
        d[i,i] = 1
    return d

# Матрица Якоби
def F_u(u, mass):
    f_u = np.zeros((5,5))
    f_u[0,2] = 1
    f_u[1,3] = 1
    f_u[2,0] = u[4] / (mass * L)
    f_u[2,4] = u[0] / (mass * L)
    f_u[3,1] = u[4] / (mass * L)
    f_u[3,4] = u[1] / (mass * L)
    f_u[4,0] = 2 * u[0]
    f_u[4,1] = 2 * u[1]
    return f_u



u = np.zeros((N+1, 5))#массив для функции u
u[0,:] = [x0, y0, vx0, vy0, T_0] #первая строка массива с нач. усл.

#схема розенброка
for n in range(N):
    w1 = np.linalg.solve(D() - alpha * tau * F_u(u[n], m), F(u[n], g0, m, L, t[n]))
    u[n+1] = u[n] + tau * w1.real


# Построение графика

    
# Получаем значения x(t), y(t) и sqrt(x^2 + y^2)
x_values = u[:, 0]
y_values = u[:, 1]
r_values = np.sqrt(x_values**2 + y_values**2)

plt.figure(figsize=(20, 17))
plt.suptitle("интервал [0,100]", fontsize=20)
plt.subplot(3, 1, 1)
plt.plot(t, x_values, 'b-')
plt.title("График x(t)")
plt.xlabel("Время (с)")
plt.ylabel("x (м)")
plt.grid()

plt.subplot(3, 1, 2)
plt.plot(t, y_values, 'g-')
plt.title("График y(t)")
plt.xlabel("Время (с)")
plt.ylabel("y (м)")
plt.grid()

plt.subplot(3, 1, 3)
plt.plot(t, r_values, 'r-')
plt.title("График sqrt($x^2$ + $y^2$)(t)")
plt.xlabel("Время (с)")
plt.ylabel("sqrt($x^2$ + $y^2$) (м)")
plt.grid()
plt.show()
'''
plt.figure()
for n in range(N):
    #plt.plot((0, u[n,0]), (0, u[n,1]), 'b-')  # Линия стержня
    plt.plot(u[n,0], u[n,1], 'ro')  # Конец стержня
    
plt.tight_layout()
plt.show()
plt.title("Движение системы")
plt.xlabel("x (м)")
plt.ylabel("y (м)")
plt.axis('equal')
plt.grid()
'''

#Period = 2*(2*np.pi*np.sqrt(L/g(g0, t)))
#t = linespace(t_0, Period, N + 1)