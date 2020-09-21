/**************************************
 MODEL FILE
 **************************************/


/*
 This header file contains contains all the model dependet functions necessary for the program. 
 Note that model dependent parameters should go in g2parameters.h. 
 Note that all functions must be in program units which are given by the following rescallings (pr denotes quantity used in program). 
 Any other model functions should be added here (and may need to tweek g2function.cpp and g2output.cpp).
 
 
 B=mphi //Note that B may change depending on your model
 dt_(pr)=dt*B
 x_(pr)=x*B
 f_(pr)=f
 V_(pr)=1/B^2*V
 dV_(pr)/df_(pr)=1/B^2*dV/df
 
 
 Copyright (2013): 
 Kenyon College
 John T. Giblin, Jr
 Tate Deskins and Hillary Child
 Last Updated: 06.27.2013
*/

#include "g2header.h" //contains declerations for program functions.

// Model specific details about the run to be output to an information file
//change as needed
void modelinfo(FILE *info)
{
    // Name and description of model
    fprintf(info,"Chaotic (Quadratic) Inflation\n");
    fprintf(info,"V =1/2 m^2 phi^2 \n\n");
    
    // Model specific parameter values
    fprintf(info,"mass= %Le m_pl\n",mphi);
	fprintf(info,"g^2= %Le m_pl\n",gsq);
	
#if rand_init==1
    const char initType[6]="U";
#else
    const char initType[6]="Not u";
#endif
    fprintf(info,"%ssing random initial conditions",initType);
}


/***************************
 user defined model functions
 ***************************/
#define PHI field[s][0]
#define CHI field[s][1]
#define PHIDOT dfield[s][0]
#define CHIDOT dfield[s][1]

#define COUP gsq/mphi/mphi //coupling term

long double potential(int s, int i, int j, int k)//user defined potential
{
    return (0.5*PHI[i][j][k]*PHI[i][j][k] + 0.5*COUP*PHI[i][j][k]*PHI[i][j][k]*CHI[i][j][k]*CHI[i][j][k]);
}


long double dVdf(int s, int fld, int i, int j, int k)//user defined derivative of the potential

{
	switch (fld) {
		case 0://derivative with respect to phi
			return (PHI[i][j][k] + COUP*CHI[i][j][k]*CHI[i][j][k]*PHI[i][j][k]);
			break;
		case 1://derivative with respect to chi
			return (COUP*PHI[i][j][k]*PHI[i][j][k]*CHI[i][j][k]);
			break;
		default:
			break;
	}

}

inline long double effMass(int s, int fld)//the effective mass used for random inital conditions
{
    long double avemass=0.;  
    int i,j,k;
    switch (fld) {
        case 0:
			LOOP{
				avemass += (1. + COUP*CHI[i][j][k]*CHI[i][j][k]);
			}
			return avemass/gridsize;
	case 1:
			LOOP{
				avemass += (COUP*PHI[i][j][k]*PHI[i][j][k]);
			}
			return avemass/gridsize;
        default://sets mass as 1(rescaled)if there is no case structure
            return 1.;
            break;
    }
}
/*******************
 field initialization
 *******************/

void initfields()//here the user may decide how the fields will be initialized
{
    static int first=0,s=0;
    DECLARE_INDEX
    
    if(first==0)
    {
        fldLOOP//loops over fld i,j,k
        {
			field[s][fld][i][j][k]=f0[fld];//initialize each field as its initial value
            dfield[s][fld][i][j][k]=df0[fld];// initialize each field derivative as its initial value
        }
   	
    }
    
    if(first==1)
    {
#if rand_init==1
		for(fld=0; fld<nflds; fld++){
			randInit(field[s][fld],dfield[s][fld],effMass(s,fld));//adds random intial conditions ontop of mean value above
		}
        initDestroy();
		printf("Fields fluctuated\n");
#endif
		
//Any other model specific initialization can go here -- i.e. Bubbles, etc	
    }
	
    calcEnergy(0); //This is important -- needed for first step of evolution

    first++;
}
#undef PHI
#undef CHI
#undef PHIDOT
#undef CHIDOT







