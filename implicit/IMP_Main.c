/***************************************************************************
 *
 *  Function    : Main
 *  Description : main routine for ParaDiS simulation
 *
 **************************************************************************/
#include <stdio.h>
#include <sys/time.h>
#include <time.h>
#include "Home.h"
#include "Init.h"

#ifdef _LIS
#include "lis.h"
#endif

#ifdef PARALLEL
#include "mpi.h"
#endif

#ifdef _GPU_SUBCYCLE
#include "SubcycleGPU.h"
#endif

#ifdef FPES_ON
#include <fpcontrol.h>
#endif


int main (int argc, char *argv[])
{
        int     cycleEnd, memSize, initialDLBCycles;
        time_t  tp;
        Home_t  *home;
        Param_t *param;

/*
 *      On some systems, the getrusage() call made by Meminfo() to get
 *      the memory resident set size does not work properly.  In those
 *      cases, the function will try to return the current heap size 
 *      instead.  This initial call allows meminfo() to get a copy of
 *      the original heap pointer so subsequent calls can calculate the
 *      heap size by taking the diference of the original and current
 *      heap pointers.
 */
        Meminfo(&memSize);

/*
 *      on linux systems (e.g. MCR) if built to have floating point exceptions
 *      turned on, invoke macro to do so
 */
   
#ifdef FPES_ON
        unmask_std_fpes();
#endif


        ParadisInit(argc, argv, &home);
        
#ifdef _GPU_SUBCYCLE
        InitializeParadisGPU(home);
#endif
   
        home->cycle        = home->param->cycleStart;
        home->NumForceCall = 0;
        param              = home->param;
        cycleEnd           = param->cycleStart + param->maxstep;
        initialDLBCycles   = param->numDLBCycles;

/*
 *      Perform the needed number (if any) of load-balance-only
 *      steps before doing the main processing loop.  These steps
 *      perform only the minimal amount of stuff needed to
 *      estimate per-process load, move boundaries and migrate
 *      nodes among processsors to get a good initial balance.
 */
        gettimeofday(&home->clock_time_beg, NULL );
        TimerStart(home, INITIALIZE);

        if ((home->myDomain == 0) && (initialDLBCycles != 0)) {
            time(&tp);
            printf("  +++ Beginning %d load-balancing steps at %s",
                   initialDLBCycles, asctime(localtime(&tp)));
        }

        while (param->numDLBCycles > 0) {
            ParadisStep(home);
            home->cycle++;
            param->numDLBCycles--;
        }

        if ((home->myDomain == 0) && (initialDLBCycles != 0)) {
            time(&tp);
            printf("  +++ Completed load-balancing steps at %s",
                   asctime(localtime(&tp)));
        }

        TimerStop(home, INITIALIZE);

/*
 *      Any time spent doing the initial DLB-only steps should
 *      just be attributed to initialization time, so be sure to
 *      reset the other timers before going into the main
 *      computational loop
 */
        TimerInitDLBReset(home);

#ifdef _IMPLICIT
/*
 * 		If we are using the implicit integrator, initialize the linear
 * 		solver (lis).
 */
		if (strcmp(param->timestepIntegrator, "implicit") == 0) {
#ifdef _LIS
            lis_initialize(&argc, &argv);
#else
			Fatal("Linear Algebra Library (-D_LIS) is required for implicit integrator!");
#endif
        }
#endif

#ifdef _CROSSSLIP
/*
 *		Thermal activation requires that the random number generater is 
 *		initialized.
 */
		srand((long)time(NULL));
#endif

#ifdef _NOBUF_STDOUT
/*
 *	Force the stdout to flush.
 */
	setbuf(stdout, NULL);
#endif

/*
 *      The cycle number may have been incremented during the initial
 *      load-balance steps, so reset it to the proper starting
 *      value before entering the main processing loop.
 */
        home->cycle = home->param->cycleStart;

        while (home->cycle < cycleEnd) {
            ParadisStep(home);
            TimerClearAll(home);
        }

        ParadisFinish(home);

        exit(0);
}
