/***************************************************************************
 *   
 *      Module:     DeltaPlasticStrain.
 *
 *      Description:  This contains a simple generic dispatch function that
 *                    will invoke the version of DeltaPlasticStrain*() 
 *                    appropriate to the type of material being simulated
 *
 ***************************************************************************/
#include "Home.h"
#include "Mobility.h"

void DeltaPlasticStrain(Home_t *home)
{
		int  i;

		switch(home->param->materialType) {

			case MAT_TYPE_BCC:
				DeltaPlasticStrain_BCC(home);
				break;

			case MAT_TYPE_FCC:
				DeltaPlasticStrain_FCC(home);
				break;

			default:
				DeltaPlasticStrain_BCC(home);
				break;
		}
		
/*********************************************
 *
 *		NOTE: PARADIS USED TO CALCULATE THE CHANGE OF ROTATION 
 *		      AXIS IN THE OPPOSITE DIRECTION. THE FOLLOWING LINES
 *		      FIX THIS PROBLEM.
 *
 *		      LAST EDIT: AMIN AGHAEI (07/22/2015)
 *********************************************/
		for (i = 0; i < 6; i++) {
            home->param->delpSpin[i] = -home->param->delpSpin[i];
        }

        return;
}
