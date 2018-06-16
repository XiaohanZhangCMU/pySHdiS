############################################################################
#
#    makefile: builds a version of paradis with new integrators and a 
#    new collision code
#
############################################################################
include ../makefile.sys
include ../makefile.setup

DEFS += -D_IMPLICIT
DEFS += -D_RETROCOLLISIONS

#
#	Needed for compiling shared library
#
DEFS += -fPIC
DEFS +=  -D_SHDIS
#DEFS += -D_GPU_SUBCYCLE
DEFS += -D_NOYOFFESTRESS

#enable thermally-activated cross-slip
#DEFS += -D_CROSSSLIP 
#print when cross slip events occur
#DEFS += -D_PRINT_CROSSSLIP  
#print when glide planes are violated
#DEFS += -D_DEBUG_GLIDE_VIOLATIONS

#prevent buffering of stdout output, so it continually flushes
#useful with mc2 when debugging
DEFS += -D_NOBUF_STDOUT

# pre-compute glide constraints
#DEFS += -D_PRECOMPUTE_GLIDE
# pre-compute FMM contribution during subcycling
#DEFS += -D_PRECOMPUTE_SUB_FMM

## The following flags for data output

#DEFS += -D_PRINTSTRESS
#DEFS += -D_WRITENODEFORCE

SRCDIR = ../src
INCDIR = ../include
BINDIR = ../bin
LIBDIR = ../lib

PARADISIMP_INCDIR = Include
#
#	The utilities use various source modules from the parallel
#       code located in the parent directory.  Maintain a list of
#	these files.
#
#	These modules are compiled in the parent directory with a
#	different set of preprocessor definitions than are needed
#	here, so we need to create links in this directory back to
#	the source modlues and create separate object modules for
#	these sources.
#

EXTERN_C_SRCS = CellCharge.c             \
      CommSendGhosts.c         \
      CommSendGhostPlanes.c    \
      CommSendMirrorNodes.c    \
      CommSendRemesh.c         \
      CommSendSecondaryGhosts.c \
      CommSendSegments.c       \
      CommSendVelocity.c        \
      CorrectionTable.c        \
      CrossSlip.c              \
      CrossSlipBCC.c           \
      DebugFunctions.c         \
      Decomp.c                 \
      DeltaPlasticStrain.c     \
      DeltaPlasticStrain_BCC.c \
      DeltaPlasticStrain_FCC.c \
      DisableUnneededParams.c  \
      DLBfreeOld.c             \
      deWitInteraction.c       \
      FixRemesh.c              \
      FMComm.c                 \
      FMSigma2.c               \
      FMSupport.c              \
      ForwardEulerIntegrator.c \
      FreeInitArrays.c         \
      GetDensityDelta.c        \
      GetNewGhostNode.c        \
      Gnuplot.c                \
      Heap.c                   \
      InitCellDomains.c        \
      InitCellNatives.c        \
      InitCellNeighbors.c      \
      InitHome.c               \
      InitRemoteDomains.c      \
      InitSendDomains.c        \
      InputSanity.c            \
      LoadCurve.c              \
      Matrix.c                 \
      Meminfo.c                \
      MemCheck.c               \
      Migrate.c                \
      MobilityLaw_BCC_0.c      \
      MobilityLaw_BCC_0b.c     \
      MobilityLaw_BCC_glide.c  \
      MobilityLaw_FCC_0b.c     \
      MobilityLaw_FCC_climb.c  \
      MobilityLaw_Relax.c      \
      NodeVelocity.c           \
      OsmoticForce.c           \
      ParadisInit.c            \
      ParadisFinish.c          \
      ParadisThread.c          \
      Parse.c                  \
      PickScrewGlidePlane.c    \
      Plot.c                   \
      ProximityCollision.c     \
      QueueOps.c               \
      ReadRestart.c            \
      ReadBinaryRestart.c      \
      ResetGlidePlanes.c       \
      RBDecomp.c               \
      RSDecomp.c               \
      RemapInitialTags.c       \
      Remesh.c                 \
      RemeshRule_3.c           \
      RemoteSegForces.c        \
      RemoveNode.c             \
      SemiInfiniteSegSegForce.c \
      SortNativeNodes.c        \
      SortNodesForCollision.c  \
      SplitSurfaceNodes.c      \
      Tecplot.c                \
      Util.c                   \
      WriteArms.c              \
      WriteAtomEye.c	       \
      WriteBinaryRestart.c     \
      WriteDensFlux.c          \
      WriteDensityField.c      \
      WriteForce.c             \
      WriteFragments.c         \
      WritePoleFig.c           \
      WritePovray.c            \
      WriteRestart.c           \
      WriteVelocity.c          \
      WriteVisit.c


EXTERN_CPP_SRCS = DisplayC.C       \
                  display.C 

EXTERN_HEADERS = Cell.h           \
                 Comm.h           \
                 DebugFunctions.h \
                 Decomp.h         \
                 DisplayC.h       \
                 FM.h             \
                 Force.h          \
                 InData.h         \
                 Init.h           \
                 Matrix.h         \
                 MirrorDomain.h   \
                 OpList.h         \
                 ParadisGen.h     \
                 ParadisThread.h  \
                 Parse.h          \
                 QueueOps.h       \
                 RBDecomp.h       \
                 RSDecomp.h       \
                 RemoteDomain.h   \
                 Restart.h        \
                 Tag.h            \
                 Topology.h       \
                 Util.h           \
                 WriteProp.h      \
                 display.h       

EXTERN_SRCS = $(EXTERN_C_SRCS) $(EXTERN_CPP_SRCS)
EXTERN_OBJS = $(EXTERN_C_SRCS:.c=.o) $(EXTERN_CPP_SRCS:.C=.o)

#
#       Define the sources in the partial/ directory
#

# Files that need to be merged with partial/ directory
#                   Initialize.c          (1 block)
#                   LocalSegForces.c      (3 blocks) 
#                   NodeForce.c           (6 blocks) 
#                   Param.c               (1 block) 
#                   Topology.c            (RemoveDoubleLinks, SplitNode, MergeNode) 

#
#	Define the exectutable, source and object modules for
#	the problem generator.
#


PARADISIMP     = paradisimp
PARADISIMP_BIN = $(BINDIR)/$(PARADISIMP)

LIBPARADISIMP = lib$(PARADISIMP).so
PARADISIMP_LIB = $(LIBDIR)/$(LIBPARADISIMP)

PARADISIMP_C_SRCS = IMP_Main.c                \
                    Collision.c               \
                    CommSendCoord.c           \
                    CommSendVelocitySub.c     \
                    CrossSlipFCC.c            \
                    FindPreciseGlidePlane.c   \
                    FixGlideViolations.c      \
                    GetNewNativeNode.c        \
                    GenerateOutput.c          \
                    ImplicitIntegrator.c      \
                    ImplicitIntegratorSub.c   \
                    ImplicitIntegratorSub0.c  \
                    Initialize.c			  \
                    LocalSegForces.c          \
                    MobilityLaw_FCC_0.c       \
                    MobilityLaw_Glide.c 	  \
                    NodeForce.c               \
                    NodeForceList.c           \
                    NodeForceList2.c          \
                    ParadisStep.c             \
                    Param.c                   \
                    PredictiveCollision.c     \
                    RemeshRule_2.c            \
                    RetroactiveCollision.c    \
                    RKFIntegrator.c           \
                    SegSegList.c              \
                    SortNodes.c               \
                    SparseMatrix.c            \
                    SubcycleIntegrator.c      \
                    Timer.c                   \
                    Topology.c                \
                    ThermalActivation.c       \
                    TrapezoidIntegrator.c     \
                    TrapezoidIntegratorSub.c  \
                    TrapezoidIntegratorSub1.c \
                    WriteParaview.c           \
                    WriteProp.c              \
										AllSegmentStress.c       \
										SH.c
                    
# CUDA objects
ifeq (D_GPU_SUBCYCLE,$(findstring D_GPU_SUBCYCLE,$(DEFS)))
PARADISIMP_CUDA_SRCS = SubcycleGPU.cu
endif


PARADISIMP_HEADERS = $(addprefix $(PARADISIMP_INCDIR)/, $(EXTERN_HEADERS))

PARADISIMP_INCS.gcc = -I$(HOME)/usr/include
PARADISIMP_LIBS.gcc = -lm

PARADISIMP_INCS.linux = -I/usr/include
PARADISIMP_LIBS.linux = -L/usr/lib/ -lm

PARADISIMP_INCS.mc2 = -I$(HOME)/usr/include
PARADISIMP_LIBS.mc2 = -L$(HOME)/usr/lib -lm

PARADISIMP_INCS.mac = -I$(HOME)/usr/include
PARADISIMP_LIBS.mac = -L$(HOME)/usr/lib -lm

PARADISIMP_INCS.linux.pc = -I/usr/include
PARADISIMP_LIBS.linux.pc = -L/usr/lib/ -lm


PARADISIMP_INCS.wcr = -I/opt/fftw-3.1.2/intel/include
PARADISIMP_LIBS.wcr = -L/opt/fftw-3.1.2/intel/lib/  -lm

PARADISIMP_INCS.vip = 
PARADISIMP_LIBS.vip =  -L/u/system/Power/libs/fftw-3.1.2/lib  -lm

PARADISIMP_INCS.su-ahpcrc = -I/lustre/home/mpotts/fftw/include
PARADISIMP_LIBS.su-ahpcrc = -L/lustre/home/mpotts/fftw/lib -llis -lm  

# For mjm
PARADISIMP_INCS.linux.opteron = -I/usr/cta/pet/MATH/include
PARADISIMP_LIBS.linux.opteron = -L/usr/cta/pet/MATH/lib/ -llis -lm

PARADISIMP_INCS.cygwin = 
PARADISIMP_LIBS.cygwin = -llis

PARADISIMP_INCS.mc-cc = 
PARADISIMP_LIBS.mc-cc = -llis

PARADISIMP_INCS.cray =-I/mnt/cfs/pkgs/PTOOLS/pkgs/fftw-2.1.5/include
PARADISIMP_LIBS.cray =-L/mnt/cfs/pkgs/PTOOLS/pkgs/fftw-2.1.5/lib/ -llis


PARADISIMP_INCS.Harold = -I/mnt/lustre/usrcta/pet/pkgs/fftw/3.2.2/include
PARADISIMP_LIBS.Harold = -L/mnt/lustre/usrcta/pet/pkgs/fftw/3.2.2/lib -llis


PARADISIMP_INCS = $(INCS_$(MODE).$(SYS)) $(XLIB_INCS) $(MPI_INCS) \
		  $(HDF_INCS)  $(PARADISIMP_INCS.$(SYS)) -I $(PARADISIMP_INCDIR)

PARADISIMP_LIBS = $(PARADISIMP_LIBS.$(SYS))

PARADISIMP_SRCS = $(PARADISIMP_C_SRCS) $(PARADISIMP_CPP_SRCS) $(PARADISIMP_CUDA_SRCS)
PARADISIMP_OBJS = $(PARADISIMP_C_SRCS:.c=.o) $(PARADISIMP_CPP_SRCS:.C=.o) $(PARADISIMP_CUDA_SRCS:.cu=.o)


###########################################################################
#
#	CUDA Compiler options
#	Do not forget to specify the target GPU architecture (-arch=...)
#	in CUDA_FLAGS for optimal performance.
#
###########################################################################
NV = nvcc
#NV_FLAGS = $(DEFS) -O3 $(CUDA_FLAG.$(SYS))
NV_FLAGS = -Xcompiler -fno-strict-aliasing -Xcompiler -D_IMPLICIT  -Xcompiler -D_RETROCOLLISIONS  -Xcompiler  -fPIC -Xcompiler  -D_GPU_SUBCYCLE  -Xcompiler -D_NOBUF_STDOUT -Xcompiler -shared -Xcompiler  -O3  $(CUDA_FLAG>$(SYS))
#--use_fast_math -arch=sm_20 -arch=sm_52 -g -G -Xcompiler -rdynamic -lineinfo -O3
NV_INCS = -I Include/ 
ifeq (D_GPU_SUBCYCLE,$(findstring D_GPU_SUBCYCLE,$(DEFS)))
NV_LIBS = $(CUDA_LIB.$(SYS))
endif


###########################################################################
#
#	Define a rule for converting .c files to object modules.
#	All modules are compile serially in this directory
#
###########################################################################

.c.o:		makefile ../makefile.sys ../makefile.setup
		- @ make $(PARADISIMP_HEADERS) > /dev/null 2>&1
		$(CC) $(OPT) $(CCFLAG) $(PARADISIMP_INCS) -c $<

.C.o:		makefile ../makefile.sys ../makefile.setup
		- @ make $(PARADISIMP_HEADERS) > /dev/null 2>&1
		$(CPP) $(OPT) $(CPPFLAG) $(PARADISIMP_INCS) -c $<

%.o:%.cu		makefile ../makefile.sys ../makefile.setup
		- @ make $(PARADISIMP_HEADERS) > /dev/null 2>&1
		$(NV) $(NV_FLAGS) $(NV_INCS) -c $<

###########################################################################
#
#	Define all targets and dependencies below
#
###########################################################################

all:		$(EXTERN_OBJS) $(PARADISIMP) 

clean:
		rm -f *.o $(EXTERN_SRCS) $(PARADISIMP_BIN) test_YoffeInfStress
		rm -f $(PARADISIMP_HEADERS)

depend:		 *.c $(SRCDIR)/*.c $(INCDIR)/*.h makefile
		makedepend -Y$(INCDIR) *.c  -fmakefile.dep

headers:
		echo $(PARADISIMP_HEADERS)

#
#	Create any necessary links in the current directory to source
#	modules located in the SRCDIR directory
#

$(EXTERN_SRCS): $(SRCDIR)/$@
		- @ ln -s  -f $(SRCDIR)/$@ ./$@ > /dev/null 2>&1

$(PARADISIMP_HEADERS): $(INCDIR)/$(notdir $@)
		- @ ln -s  -f ../$(INCDIR)/$(notdir $@) $@ > /dev/null 2>&1

# For vip
#$(EXTERN_SRCS): $(SRCDIR)/$@
#                ln -s  -f $(SRCDIR)/$@ ./$@ > /dev/null 2>&1


$(PARADISIMP):	$(PARADISIMP_BIN)
$(PARADISIMP_BIN): $(PARADISIMP_SRCS) $(PARADISIMP_OBJS) $(EXTERN_OBJS) $(HEADERS)
	echo $(PARADISIMP_OBJS)
		$(CPP) $(OPT) $(OPENMP_FLAG) $(PARADISIMP_OBJS) $(EXTERN_OBJS) -o $@  $(LIB) $(PARADISIMP_LIBS) $(NV_LIBS)

$(LIBPARADISIMP):	$(LIBDIR) $(PARADISIMP_LIB)
$(PARADISIMP_LIB): $(PARADISIMP_OBJS) $(EXTERN_OBJS) $(HEADERS)
		$(CPP) -shared $(OPT) $(OPENMP_FLAG) $(PARADISIMP_OBJS) $(EXTERN_OBJS) \
		-o $@ $(LIB) $(PARADISIMP_LIBS) $(NV_LIBS)


