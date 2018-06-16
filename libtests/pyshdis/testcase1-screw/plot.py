import numpy as np
import matplotlib.pyplot as plt

maxstep = 1
mu = 1
nu = 0.333333
a = 1
b = 1 # 2.7412e-10
t = 1.5
Ngrid = 20
xmin= -20000
xmax=20000

sigma_eval_on_some_pts = np.load('./SingleLine_results/sigma_eval.npy')
print(sigma_eval_on_some_pts[5::9].shape)
plt.figure(figsize=(12, 9))
tau_ShE = sigma_eval_on_some_pts[5::9]
x3 = np.linspace(0, 3, 30)
plt.plot(x3, -tau_ShE, '^')
plt.legend(['ShElastic solution ($t/a='+str(t/a)+',l_{max}=20$)'])
plt.xlabel('z/t')
plt.ylabel(r'$\tau_{yz}a/\mu$')
plt.xlim(0, 3)
#    plt.ylim(0, 0.01)
plt.show()
