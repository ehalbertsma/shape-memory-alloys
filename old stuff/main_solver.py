from gekko import GEKKO
import numpy as np
import matplotlib.pyplot as plt

m = GEKKO()
m.time = np.linspace(0,20,100)
k = 10
f = m.Var(value=0.0)
y = m.Param(value=m.time)
z = m.Param(value=m.time)
t = m.Param(value=m.time)
m.Equation(f.dy()==f.dz().dt())
m.options.IMODE = 4
m.solve(disp=True)

plt.plot(m.time,f.value)
plt.xlabel('time')
plt.ylabel('f')
plt.show()