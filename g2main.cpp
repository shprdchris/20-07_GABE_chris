/****************************
 LatticeSolve Main File
 ****************************/

/*
 This is the .cpp file that contains the main function for the GABE2.0 program. 
 In this function the looping over time occurs aswell as the calls to the output functions and the initialize function. 
 Also all the header calls and macro definitions are made in this file.
 
 Copyright (2013): 
 Kenyon College
 John T. Giblin, Jr
 Tate Deskins and Hillary Child
 Last Updated: 06.27.2013
*/

/*********************
 GABE Header
 *********************/
 //contains all the parameters, changable and unchangeable
#include "g2header.h" //contains declerations for program functions and parameters file


long double t=starttime;//this is the variable that stores the time
long double (*field)[nflds][N][N][N];//this stores the field values for each step along the grid
long double (*dfield)[nflds][N][N][N];//this stores the derivative of the field for each step along the grid
long double a[2];//this stores the scale facator for each step
long double adot[2];// this stores the time derivative of the scale factor
long double edpot[2]; //this stores the average potential energy
long double edkin[2]; //this stores the average kinetic energy
long double edgrad[2]; //this stores the average gradient energy
long double edrho[2]; // this stores the avg. energy density over the box


int main()
{
    omp_set_nested(1);//allows for nested parallelization
    printf("\n\nLattice evolution program started\n\n");
    output_parameters(); //Outputs general run information (info.txt)
    printf("Info file made\n");
    
    field=(long double(*)[nflds][N][N][N]) malloc(sizeof(long double)*nflds*2*N*N*N);//allocates memory for the fields
    printf("'field' memory allocated\n");
    
    dfield=(long double(*)[nflds][N][N][N]) malloc(sizeof(long double)*nflds*2*N*N*N);//allocates memory for the fields' time derivatives
    printf("'dfield' memory allocated\n");
    
    printf("Memory allocated for all arrays\n\n");
    printf("Starting run...\n\n");
    
    
    for(t=starttime;t<=endtime;t+=dt)//This is the main for-loop over time. Note that the current t value is the time for the currently evaluated fields
    {
		
        if(t==starttime)
        {
            initfields();//initializes fields as defined in model.h
            printf("Fields initialized (no fluctuations)\n");
            initexpansion();//start expansion;
            printf("Expansion started\n");
			
            initfields();//do fluctuations
            printf("Model Specific Initialization Completed\n");

            printf("Time evolution begining\n");
            
			screenout();//outputs current pr time to screen (see g2output.cpp)
			outputslice();//outputs slice values to file

        }

		step();//evolves the fields one step
        
		screenout();//outputs current pr time to screen
        outputslice();//outputs slice values to file
        
	}
#if spec_output==1
    specClear();//clears fftw stuff needed for spectra at the end of the run
#endif
    printf("\nRun complete\n\n");
    output_parameters();//finalizes any necessary info for the run to info.txt)
    printf("GABE2.2 finished\n\n\n");
}
