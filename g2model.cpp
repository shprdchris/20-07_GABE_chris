/**************************************
 MODEL FILE
 **************************************/


/*
 This header file contains contains all the model dependet functions necessary for the program. 
 Note that model dependent parameters should go in g2parameters.h. 
 Note that all functions must be in program units which are given by the following rescallings (pr denotes quantity used in program). 
 Any other model functions should be added here (and may need to tweek g2function.cpp and g2output.cpp).
 
 
 Copyright (2013): 
 Kenyon College
 John T. Giblin, Jr
 Tate Deskins and Hillary Child
 Last Updated: 06.27.2013
*/

#include "g2header.h" //contains declerations for program functions.

/*long double mod(long double f1, long double f2){
	return sqrt(pow(f1,2), + pow(f2,2));
	}
long double mod2(long double f1, long double f2){
	return (pow(f1,2), + pow(f2,2));
	}*/

// Model specific details about the run to be output to an information file
//change as needed
void modelinfo(FILE *info)
{
    // Name and description of model
    fprintf(info,"R2HI inflation, just scalar sector\n");
  //  fprintf(info,"V =1/2 m^2 phi^2 \n\n");
    
    // Model specific parameter values
    fprintf(info,"beta = = %.10Le \n",b);

	
#if rand_init==1
    const char initType[6]="U"; //TODO what does this mean
#else
    const char initType[6]="Not u";
#endif
    fprintf(info,"%ssing random initial conditions",initType);
}


/***************************
 user defined model functions
 ***************************/
#define PHI field[s][0]
#define H1 field[s][1]
#define H2 field[s][2]
#define H3 field[s][3]
#define H4 field[s][4]
#define PHIDOT dfield[s][0]
#define H1DOT dfield[s][1]
#define H2DOT dfield[s][2]
#define H3DOT dfield[s][3]
#define H4DOT dfield[s][4]

//#define COUP gsq/mphi/mphi //coupling term

long double potential(int s, int i, int j, int k)//user defined potential 
{
		long double pw2h = pow(H1[i][j][k], 2) +  pow(H2[i][j][k],2) + pow(H3[i][j][k],2) + pow(H4[i][j][k],2);
			long double er23f = exp(sqrt23*PHI[i][j][k]);
	 return (
			(3*(b*l + pow(xi1,2))*(pw2h*pw2h*l + pow(1 -er23f + pw2h*xi1,2)/b))/  (4*er23f*er23f*l)	
	 );
}


long double dVdf(int s, int fld, int i, int j, int k)//user defined derivative of the potential 
{
	long double pw2h = pow(H1[i][j][k], 2) +  pow(H2[i][j][k],2) + pow(H3[i][j][k],2) + pow(H4[i][j][k],2);
	long double er6f = exp(sqrt(6)*PHI[i][j][k]);
	long double er23f = exp(sqrt23*PHI[i][j][k]);

	switch (fld) {
		case 0:	
			{
						return (
							(sqrt(1.5)*(b*l + pow(xi1,2))*(-1 + er23f - b*pw2h*pw2h*l + 
							   (-2 + er23f)*pw2h*xi1 - pw2h*pw2h*pow(xi1,2)))/
						   (b*er23f*er23f*l)				  
						   );	
							
			}
			break;

		case 1:
			{
				
			return (
				(3*H1[i][j][k]*(b*l + pow(xi1,2))*(b*pw2h*l - (-1 + er23f)*xi1 +    pw2h*pow(xi1,2)))/(b*er23f*er23f*l)
					);
			}
			break;

		case 2:
			{
		return (	
				(3*H2[i][j][k]*(b*l + pow(xi1,2))*(b*pw2h*l - (-1 + er23f)*xi1 +    pw2h*pow(xi1,2)))/(b*er23f*er23f*l)
					);	
			}
		break;

	case 3:
			{

		return (
			(3*H3[i][j][k]*(b*l + pow(xi1,2))*(b*pw2h*l - (-1 + er23f)*xi1 +    pw2h*pow(xi1,2)))/(b*er23f*er23f*l)
					);
			}
		break;

	case 4:
			{
		return (
				(3*H4[i][j][k]*(b*l + pow(xi1,2))*(b*pw2h*l - (-1 + er23f)*xi1 +    pw2h*pow(xi1,2)))/(b*er23f*er23f*l)
				);
			}
		break;

		default:
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
	long double mdiagV[nflds]; //first value. need to assign pointer to something
	long double evecV[nflds][nflds];
	long double rav2[nflds]; //average of random fluctuations
	long double rdav2[nflds];


	for(fld=0;fld<nflds;fld++){
		rav2[fld]=0;
		rdav2[fld]=0;
	}

   
    if(first==0)
    {
		
        fldLOOP//loops over fld i,j,k //TODO
        {
			field[s][fld][i][j][k]=f0[fld];//initialize each field as its initial value
            dfield[s][fld][i][j][k]=df0[fld];// initialize each field derivative as its initial value
        }
    }

    if(first==1)
    {
#if rand_init==1 

		mass_matrix_calc( 0, &evecV, &mdiagV, true); 
		
	//	printf("Initial mass eigenvalues are:\n %Le \n %Le \n %Le \n", mdiagV[0], mdiagV[1], mdiagV[2] ); //expect columns 1,2 to have swap symmetry. and elements (0,1), (0,2) to be identical.
		
		for(fld=0; fld<nflds; fld++){	//currently only initializing phi_L and heavy Higgs directions
			randInit(Rfield[s][fld],Rdfield[s][fld], mdiagV[fld] ); //initializes fluctuations in basis of mass eigenstatess 
		}
		
		fromDiag( Rfield ,Rdfield , &evecV ); //convert flucations to (phi, h_1, h_2)

		fldLOOP{
			field[s][fld][i][j][k] +=specMod* (*Rfield)[fld][i][j][k] ;  //add to background values
			dfield[s][fld][i][j][k] +=specMod*(*Rdfield)[fld][i][j][k];
		}
		
        initDestroy();
		//printf("Fields fluctuated\n");
	
#endif
		
//Any other model specific initialization can go here -- i.e. Bubbles, etc	
    }
  
    calcEnergy(0); //This is important -- needed for first step of evolution
//	printf("energ (model) = %Le \n", edrho[0]);
    first++;
}
#undef PHI
#undef H1
#undef H2
#undef H3
#undef H4
#undef PHIDOT
#undef HDOT
#undef H2DOT
#undef H3DOT
#undef H4DOT







