/*--------------------------------------------------------------------------
 *
 *	Implicit.h	Define the struct that holds all relevant data for the
 *		implicit solver.
 *
 *------------------------------------------------------------------------*/

#ifndef _Implicit_h
#define _Implicit_h

#include "Typedefs.h"

#define MAX_J_SIZE 1000 //temporary maximum Jacobian size to avoid malloc calls

struct _implicit {

		//real8		*Jval;
		//int			*Jrow;
		//int			*Jcol;
		//int			Ndof;	/*total number of degrees of freedom in the system*/
		//real8		Jval[MAX_J_SIZE];
		//int			Jrow[MAX_J_SIZE];
		//int			Jcol[MAX_J_SIZE];
		//int			noJuses;
		//int			Jsize;

};

typedef struct _implicit Implicit_t;

#endif
