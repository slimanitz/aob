from matplotlib import pyplot as plt
import numpy as np



x = [0,1,2,3,4,5,6]
y = [1,5.3,4.8,5.6,5.8,16.2,17.5]


plt.plot(x,y)
plt.xlabel("Optimisations")
plt.yticks(np.arange(min(y), max(y)+1, 1.0))
plt.ylabel("GFLOPS")
plt.title("Courbe GFLOPS par optimisation")
plt.grid()

plt.show()


