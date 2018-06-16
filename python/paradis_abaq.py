import numpy as np
import time, string, sys, getopt
import subprocess, os, re, shutil
import commands
from Home import *
from ABQ import *
from glob import *

def xchdir(newdir):
  print '00000000000000000000--------', newdir
  print 'changing dir from .{0} to {1}'.format((os.getcwd())[work_dir_length:], newdir[work_dir_length:])
  os.chdir(newdir)

def clear_subdirectories():
  xchdir(ABQdir)
  commands.getstatusoutput('rm abaqus*')
  commands.getstatusoutput('rm .panfs*')
  commands.getstatusoutput('rm abq2paraInput.*')
  commands.getstatusoutput('rm abqJob.*')
  commands.getstatusoutput('rm .pan*')
  xchdir(ParaDiSdir)
  commands.getstatusoutput('rm -rf ABQ')
  commands.getstatusoutput('rm -rf Results')
  commands.getstatusoutput('rm abq2paraInput.*')
  xchdir(workdir)

def init_paradis_abq():
  global workdir, ParaDiSdir, ABQdir, paradis_ctrlfile, paradis_datafile
  global paradis_result_dir, maxstep,  work_dir_length
  global shearModulus, pois
  global abaqus_mesher, cae_pipe_file, submitjob_file, instructions_file
  global ParadisLibraryPath
  global results_folder
  
  workdir = os.getcwd() 
  work_dir_length = len(workdir)

#  paradis_ctrlfile = 'loop_test_hs.ctrl'
#  paradis_datafile = 'loop_test_hs.data'

#  paradis_ctrlfile = 'ABQfrs.ctrl'
#  paradis_datafile = 'ABQfrs.data'
  taskname = sys.argv[1]
#  taskname = 'frank_read_src_abaqus'
  paradis_ctrlfile =  workdir + '/'+taskname + '.ctrl'
  paradis_datafile =  workdir + '/'+taskname + '.data'
  abaqus_mesher =  workdir + '/'+taskname + '.py'
  cae_pipe_file = workdir + '/'+'cae_pipe.py'
  submitjob_file = workdir + '/'+'submitjob.py'

  results_folder =  workdir + '/'+taskname+'_results'

  commands.getstatusoutput('mkdir '+results_folder)
  
  ParaDiSdir =results_folder+'/ParaDiSdir/'
  ABQdir = results_folder+'/ABQdir/'
  instructions_file =  ABQdir + '/'+'instructions.py'

  commands.getstatusoutput('mkdir '+ParaDiSdir)
  commands.getstatusoutput('mkdir '+ABQdir)

  print 'ABQdir = ' , ABQdir 

  clear_subdirectories()

#  ParadisLibraryPath = "/home/xzhang11/Planet/Libs/ParaDiSDll/lib/libparadis.so"
  ParadisLibraryPath = workdir+"/../bin/libparadis.so"
  commands.getstatusoutput('mkdir '+ParaDiSdir+'/ABQ')

#  Check existence of paradis inputs and grep ParaDiS simulation constants for Abaqus
 
  xchdir(ParaDiSdir)
  print 'my current folder is ', os.getcwd(), 'and ', paradis_ctrlfile
  if not os.path.exists(paradis_ctrlfile) or not os.path.exists(paradis_datafile):    
    print "Paradis input files missing."
    exit(0)

  f = open( paradis_ctrlfile,'r')
  for line in f:
    if re.search('^dirname',line):
      paradis_result_dir = (line.split('=')[1])
      paradis_result_dir = re.sub('[!@#$ \n"]', '', paradis_result_dir)
      paradis_result_dir = paradis_result_dir + '/'
    if re.search('^maxstep',line):
      numbs = [int(s) for s in line.split() if s.isdigit()]
      maxstep = numbs[0]
    if re.search('^shearModulus',line):
      shearModulus = float(line.split('=')[1])
    if re.search('^pois',line):
      pois = float(line.split('=')[1])

  f.close()   
#  sys.exit()
  YoungsModulus = 2 * shearModulus *(1+pois)
  print "shearModulus = ", shearModulus
  print "pois = ", pois
  print "YoungsModulus = ", YoungsModulus
  xchdir(workdir)

#  Check necessary files for abaqus
  xchdir(ABQdir)  
  if not os.path.exists(cae_pipe_file)  or not os.path.exists(abaqus_mesher) or  not os.path.exists(submitjob_file):
    print "Abaqus scripts missing."
    exit(0)
  xchdir(workdir)

def paradis_bookkeeping():
  global paradislib, home, abqinit, paradisinit, paradisstep, abqstep

  xchdir(ParaDiSdir)
  argc = 2
  myargv = ctypes.c_char_p * argc
  argv = myargv("paradis",paradis_ctrlfile)
  paradislib = cdll.LoadLibrary(ParadisLibraryPath)
  p_home = POINTER(Home_t)
  p_abq    = POINTER(Abaqus_t)
  paradisinit = paradislib.ParadisInit
  paradisstep = paradislib.ParadisStep
  abqstep = paradislib.ABQ_Step

  abqinit = paradislib.ABQ_Init
  paradisinit.restype = p_home
  abqinit.restype = p_abq

  home = paradisinit(argc,byref(argv))
  
  p_memSize = POINTER(c_int)
  p_param = home.contents.param
  home.contents.cycle = p_param.contents.cycleStart
  cycleEnd = p_param.contents.cycleStart + p_param.contents.maxstep
  initialDLBCycles = p_param.contents.numDLBCycles

  xchdir(workdir)

'''
!!! Main Program Starts Here !!!

'''

def main(argv):

# sometimes I only want the program to clear up the folders - yes_clear = 1
  if len(argv) == 2:
    yes_clear = argv[0]
  else:
    yes_clear = 0

  init_paradis_abq()

#  Add necessary function names of paradis in the following function
  paradis_bookkeeping()
  
  if yes_clear:
    sys.exit()

#  create a pipe [insfile] and launch abaqus as engine wait
#  Load overhead.py to generate geometry and mesh in abaqus. 
    
  xchdir(ABQdir)
  abqlauncher = '/share/apps/abaqus/6.12/Commands/abaqus'
#  paradislauncher = '/home/xzhang11/Planet/Libs/ParaDiS/bin/paradisabq'
  if os.path.exists(instructions_file):
    print 'instruction exists, re-makefifo it.'
    commands.getstatusoutput('rm ' + instructions_file)
    commands.getstatusoutput('mkfifo '+ instructions_file)
#    commands.getstatusoutput('touch ' + instructions_file)
  
  xchdir(workdir)

  for tstep in range(0,maxstep):
    
    print 'tstep = {0}/{1}'.format(tstep, maxstep)
    t0 = time.time()  
    home.contents.cycle = tstep

    if tstep == 0:
      xchdir(ABQdir)
      p = subprocess.Popen( [abqlauncher,  'cae', 'noGUI=' + cae_pipe_file], cwd = ABQdir,stdout = subprocess.PIPE )

      ''' 
      xiaohan: here you can modify abaqus_mesher.py to create a new mesh
      
      '''
      commands.getstatusoutput('cat ' + abaqus_mesher + ' >> '+instructions_file)
      commands.getstatusoutput('cat ' + submitjob_file+ ' >> '+instructions_file)
      
      while not os.path.exists(ABQdir+'abq2paraInput.inp'):
        print 'wait for abaqus to finish abaqus_mesher'
        time.sleep(1)      
      p.kill()

    # move job1.inp to paradis dir
      xchdir(workdir)
      shutil.move(ABQdir + 'abq2paraInput.inp', ParaDiSdir + '/abq2paraInput.inp')
 

##############################
    # main loop starts
##############################

    xchdir(ParaDiSdir)
    abaqus = abqinit(home)  

    paradislib.SortNativeNodes(home)
    paradislib.CellCharge(home)

    paradislib.ABQ_stress_boundary(home,abaqus)
    paradislib.ABQ_calc_nodeF(home, abaqus)                                   
    paradislib.ABQ_writeInput(home, abaqus)


    while not os.path.exists(ParaDiSdir+'/ABQ/abqJob.inp'):
      print 'wait for paradis to prepare abqJob.inp'
      time.sleep(1)    

    # move job2.inp to abaqus dir
    shutil.move(ParaDiSdir + '/ABQ/abqJob.inp', ABQdir + 'abqJob.inp')
    
    # FEM starts to run and generates odb file
    p2 = subprocess.Popen( '{0} job=abqJob input=abqJob.inp interactive ask_delete=OFF'.format(abqlauncher), cwd = ABQdir, shell=True, stdout = subprocess.PIPE)
    trash = p2.communicate()
    # move odb file to paradis dir
    while not os.path.exists(ABQdir+'abqJob.sta'):
      print 'wait for abaqus to generate abqJob.odb'
      time.sleep(1)

    shutil.move(ABQdir+'/abqJob.odb', ParaDiSdir + 'ABQ/abqJob.odb') 
    shutil.move(ABQdir+'/abqJob.sta', ParaDiSdir + 'ABQ/abqJob.sta') 
    shutil.move(ABQdir+'/abqJob.dat', ParaDiSdir + 'ABQ/abqJob.dat') 

#    p2.kill()
#    break

    xchdir(ParaDiSdir)

    paradislib.readABQresult(home,abaqus);
    paradislib.ABQ_calc_node_stress(home, abaqus);

    paradislib.PrintStress(home,abaqus)
    paradislib.PrintForce(home,abaqus)

    xchdir(paradis_result_dir)
    paradisstep(home, abaqus)
    xchdir(workdir)

    s2=POINTER(POINTER(c_double))
    r =POINTER(c_double)
#    r[0] = 0
#    r[1] = 0
#    r[2] = 0
#    paradislib.ABQ_calc_stress(home, abaqus, r, s2)
    t1 = time.time()


if __name__ == "__main__":
   main(sys.argv[1:])
