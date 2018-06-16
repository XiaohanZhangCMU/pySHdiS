from ctypes import *
import numpy as np
import time, string, sys, getopt
#import subprocess, os, re, shutil, sys
#import commands
from Home_imp import *
#from glob import *
from numpy.ctypeslib import ndpointer
from void_dislocation_interaction_noplot import stress_soln_void_dislocation
from void_dislocation_interaction_noplot import void_mesh
#from display import paradis_draw

global PY3
import sys
if sys.version_info[0] < 3:
  PY3 = False
  print("Use Python2")
else:
  PY3 = True
  print("Use Python3")

if PY3:
  sys.path.append('../../../SHTOOLS-4.1')
else:
  sys.path.append('../../../SHTOOLS-4.0')

global useGPU, useX11, maxstep, mu, nu, a, b, t, Ngrid

#*********************************************************************************
#       Change Parameters Here
#*********************************************************************************
useGPU = False
useX11 = ".X11" #".noX11" 

maxstep = 2000000
mu = 54.6e9
nu = 0.324000
a = 4000
b = 1
t = 1.5
Ngrid = 20

if useGPU:
  from numba import vectorize, cuda
  
def init_paradis():
  global workdir, paradis_ctrlfile, paradis_datafile
  global paradis_result_dir, maxstep,  work_dir_length
  global shearModulus, nu
  global ParadisLibraryPath
  
  workdir = os.getcwd() 
  work_dir_length = len(workdir)
  taskname = sys.argv[1]
  paradis_ctrlfile =  workdir + '/'+taskname + '.ctrl'
  paradis_datafile =  workdir + '/'+taskname + '.data'
  ParadisLibraryPath = workdir+"/../../../lib/libparadisimp.so"+useX11

  if not os.path.exists(paradis_ctrlfile) or not os.path.exists(paradis_datafile):    
    print("Paradis input files missing.")
    exit(0)

# Export paradis functions
def paradis_bookkeeping():
  global paradislib, home, paradisinit, paradisstep, initializeParadisGPU, paradisAllSegmentStress
  global paradisGetNumPhysicNodes, paradisGetNodeList, paradis_SH_calc_stress
 
  argc = 2
  myargv = ctypes.c_char_p * argc
  if PY3:
    argv = myargv("paradis".encode('utf-8'),paradis_ctrlfile.encode('utf-8'))
  else:
    argv = myargv("paradis",paradis_ctrlfile)
  paradislib = cdll.LoadLibrary(ParadisLibraryPath)
  home = POINTER(Home_t)()

  #ParadisInit()
  paradisinit = paradislib.ParadisInit
  paradisinit.argtypes = [ ctypes.c_int,  POINTER(ctypes.c_char_p), POINTER(POINTER(Home_t))] 
  paradisinit(argc,byref(argv),byref(home))

  #AllSegmentStress()
  paradisAllSegmentStress = paradislib.AllSegmentStress
  paradisAllSegmentStress.argtypes = [ POINTER(Home_t), ctypes.c_double, ctypes.c_double, ctypes.c_double,ndpointer(ctypes.c_double, flags="C_CONTIGUOUS") ]
  paradisAllSegmentStress.restype = None

  #GetNumPysicNodes()
  paradisGetNumPhysicNodes = paradislib.GetNumPhysicNodes
  paradisGetNumPhysicNodes.argtypes = [ POINTER(Home_t) ] 
  paradisGetNumPhysicNodes.restype = ctypes.c_int

  #GetNodeList()
  paradisGetNodeList = paradislib.GetNodeList
  paradisGetNodeList.argtypes = [ POINTER(Home_t), POINTER(ctypes.c_double)]
  paradisGetNodeList.restype = None

  #SH_calc_stress()
  paradis_SH_calc_stress = paradislib.SH_calc_stress
  paradis_SH_calc_stress.argtypes = [ POINTER(Home_t), ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")] 
  paradis_SH_calc_stress.restype = None

  #InitializeParadisGPU()
  if useGPU:
    initializeParadisGPU = paradislib.InitializeParadisGPU
    initializeParadisGPU.argtypes = [POINTER(Home_t)]
    initializeParadisGPU(home)

  #ParadisStep()
  paradisstep = paradislib.ParadisStep
  paradisstep.argtypes = [ POINTER(Home_t)] 

  p_memSize = POINTER(c_int)
  p_param = home.contents.param
  home.contents.cycle = p_param.contents.cycleStart
  cycleEnd = p_param.contents.cycleStart + p_param.contents.maxstep
  initialDLBCycles = p_param.contents.numDLBCycles
  print(" Pyshdis init finished.")

'''
!!! Main Program Starts Here !!!

'''
def main(argv):

  init_paradis()
  paradis_bookkeeping()
  XX,YY,ZZ = void_mesh(Ngrid, mu, nu, a, b, t)

  assert(home.contents.param.contents.fmEnabled == 1)

  for tstep in range(0,maxstep):
    print('tstep = {0}/{1}'.format(tstep, maxstep))
    #t0 = time.time()  
    home.contents.cycle = tstep

    '''
    compute stress at each cavity point contributed from all segments

    '''
    tot_stress = np.zeros(XX.shape+(3,3))
  
    #@vectorize(['float32 (float32, float32)'], target = 'cuda')
    for i in np.ndindex(XX.shape):
      stress = np.zeros((3,3))
      paradisAllSegmentStress(home, XX[i], YY[i], ZZ[i], stress)
      tot_stress[i] = stress

    '''
    generate list of dislocation node positions to be evaluated

    '''
    nnodes = paradisGetNumPhysicNodes(home)
    XYZ = np.zeros((nnodes*3))
    paradisGetNodeList(home, XYZ.ctypes.data_as(POINTER(c_double)))

    '''
    tot_stress is converted to imaginary traction bc on cavity mesh XX, YY, ZZ
    call shelastic to evaluate stress at disl nodal positions X,Y,Z
    
    '''
    ptwz_sigma_eval = stress_soln_void_dislocation(Ngrid, mu,nu,a,b,t,XX, YY, ZZ, XYZ[::3],XYZ[1::3],XYZ[2::3],tot_stress)
    '''
    put pointwise stress back to node->SHstress, such that when  
    Compute.c:ComputeTF1SegSigbRem() is called, SHstress is added to each node
    WARNING: need to discuss nodal or segmental stress for this function?????
    
    '''
    sig_array = np.array(np.reshape(ptwz_sigma_eval,(ptwz_sigma_eval.size,1)), order='C')
    paradis_SH_calc_stress(home, sig_array)
    #np.save('sigma_eval.npy', sig_array)
    #t1 = time.time()

    paradisstep(home)
    #paradis_draw(Ngrid, a, xmin=-2000, xmax=2000, X, Y, Z)

if __name__ == "__main__":
   main(sys.argv[1:])
