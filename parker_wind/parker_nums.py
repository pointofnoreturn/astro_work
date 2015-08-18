import numpy as np
import math
import matplotlib.pyplot as plt

# v = u / c
# x = r / (GM/(2c^2))

def deriv(dY, x, y):
	derivative = eval(dY)
	return derivative

def RK4(dY, initX, initY, endX, h):
	XValues = np.arange(initX, endX + h, h)
	steps = int(math.ceil(XValues.size))
	YValues = np.zeros(steps)
	YValues[0] = initY
	for i in range(0, steps - 1):
		k1 = h * deriv(dY, XValues[i], YValues[i])
		k2 = h * deriv(dY, XValues[i] + 0.5 * h, YValues[i] + 0.5 * k1)
		k3 = h * deriv(dY, XValues[i] + 0.5 * h, YValues[i] + 0.5 * k2)
		k4 = h * deriv(dY, XValues[i] + h, YValues[i] + k3)
		YValues[i + 1] = YValues[i] + k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6
	return np.array([XValues, YValues])


delta = 0.0001	# small deviation from critical point

#y = u/c
#x = r/r_c

dudr = "(2 * y * (x - 1)) / (x * x * (y * y - 1))" # expression for du/dr
output1 = RK4(dudr, 1.0 + delta, 1.0 - delta, 2.5, 0.1) # integration forward
output2 = RK4(dudr, 1.0 - delta, 1.0 + delta, 0.5, -0.1) # integration backward
totalOutput = np.concatenate((output1, output2), axis = 1)
totalOutput = totalOutput[:, totalOutput[0,:].argsort()]

plt.plot(totalOutput[0,:], totalOutput[1,:])
plt.xlabel('r/r_c')
plt.ylabel('u/c')
plt.title('Stable manifold of Parker wind model')
plt.show()
