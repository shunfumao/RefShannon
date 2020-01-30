from gen_sg2 import *
import pdb
import matplotlib.pyplot as plt

x = [float(i)/100 for i in range(100)]
y_adjusty = [adjust_y_approx4(i) for i in x]
y = [i[0] for i in y_adjusty]
adjust_y = [i[1] for i in y_adjusty]

fig, ax = plt.subplots()
#plt.plot(x, y, marker='o')
plt.plot(x, adjust_y, marker='o')
ax.set_yscale('log')
plt.show()
pdb.set_trace()