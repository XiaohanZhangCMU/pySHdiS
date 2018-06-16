from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from itertools import product, combinations


def paradis_draw(Ngrid, a, cell_min, cell_max, X, Y, Z):
  fig = plt.figure()
  ax = fig.gca(projection='3d')
  ax.set_aspect("equal")
  
  # draw cube
  r = [cell_min, cell_max]
  for s, e in combinations(np.array(list(product(r, r, r))), 2):
    if np.sum(np.abs(s-e)) == r[1]-r[0]:
      ax.plot3D(*zip(s, e), color="b")
      
  # draw sphere
  u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
  x = a*np.cos(u)*np.sin(v)
  y = a*np.sin(u)*np.sin(v)
  z = a*np.cos(v)
  ax.plot_wireframe(x, y, z, color="r")
  
  # draw paradis dislocation lines
  ax.scatter(X, Y, Z, color="g", s=3)
  
  plt.show()    



