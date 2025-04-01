import numpy as np

a = 778.479e9
mS = 1.98892e30
mJ = 1.89813e27

GM = 6.674e-11
mtot = mS + mJ
Omega = np.sqrt(GM*(mtot)/(a**3))
alpha = mJ/mtot
beta = mS/mtot

xS_prime = alpha*a
xJ_prime = beta*a

x01 = 2*a + xS_prime
x02 = 2*a + xS_prime
y0 = 0

v0x = -11000
v0y1 = 2000 - Omega*x01
v0y2 = 2000 - Omega*x02

print("Initial conditions:")
print("x0 = ", x01)
print("x0 = ", x02)
print("y0 = ", y0)
print("v0x = ", v0x)
print("v0y = ", v0y1)
print("v0y = ", v0y2)
