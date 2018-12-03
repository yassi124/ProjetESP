import math as m
PI = 3.14


def arange(x):
	T = []
	for i in range(x):
		T.append(i*1.0)
	return T

x	 = []
x = arange(256)
Vr = []

for i in x:
	Vr.append(m.sin(2*PI*50*(i/256)))

print(x)
print(Vr)
