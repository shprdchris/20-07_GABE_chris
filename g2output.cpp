/****************************
 OUTPUT FILE
 ****************************/

/*
 This header file contains all the functions for the output of the data. Note that the outputslice() function needs to be generalized for n-fields and right now only gives meaning full output for 1 field. The only information outputed is the value of the field for the slice, the times atwhich the slices were output and the info.dat file.
 
 Copyright (2013): 
 Kenyon College
 John T. Giblin, Jr
 Tate Deskins and Hillary Child
 Last Updated: 06.27.2013
*/

#include "g2header.h" //contains declerations for program functions.

//readable time variables

time_t tStart,tCurrent,tLast,tFinish; // Keep track of elapsed clock time

/*****************
 output stuffs
 *****************/

void outputfield(int first)//outputs the field values
{
    static FILE *slicefield;
    //static int fld,i,j,k;
    static char name[5000];
    
    sprintf(name,"./slices/slices_fields_%d.dat", first);
    //sprintf(name,"slices_field%d_%d.dat", fld, first);
    slicefield=fopen(name,"w");
    
#if field_outdim==3//outputs slice for 3dimensions
    for(int i=0;i<N;i+=field_sliceskip){
        for(int j=0;j<N;j+=field_sliceskip){
            for(int k=0;k<N;k+=field_sliceskip){
                for(fld=0;fld<nflds,fld++){
                      fprintf(slicefield,"%Le ", field[0][fld][i][j][k]);//,rescale_B*dfield[0][fld][i][j][k]);
                      fprintf(slicefield,"\n");
                } //loop over fields
				
                fprintf(slicefield,"\n");
            } //loop over k
            fprintf(slicefield,"\n");
        }
    }
    
#elif field_outdim==2//outputs slice for 2dimensions
    
    for(int j=0;j<N;j+=field_sliceskip)
    {
        for(int k=0;k<N;k+=field_sliceskip)
        {
            for(int fld=0;fld<nflds;fld++)
            {
                fprintf(slicefield,"%Le ", field[0][fld][0][j][k]);
            }
			fprintf(slicefield, "%Le ", sqrt(  pow(field[0][1][0][j][k],2)+pow(field[0][2][0][j][k],2))  ); //radial bit
            fprintf(slicefield, "\n");
        }
        fprintf(slicefield, "\n");
    }
    
    
#elif field_outdim==1//outputs slice for 1dimension
    
    for(int k=0;k<N;k+=field_sliceskip)
    {
        for(int fld=0;fld<nflds;fld++)
        {
            fprintf(slicefield,"%Le ", field[0][fld][0][j][k]);
        }
        fprintf(slicefield, "\n");
        
        
    }
    
#endif
    fclose(slicefield);
}



void meansvars()//calculates the mean and variance of each field
{
	char name_[5000];
	static FILE *meansvarsout_;
	int i, j, k, fld;
	long double av,av_sq,var;
	long double evecV[nflds][nflds], mdiagV[nflds], mass_sq[nflds];
	
	sprintf(name_,"./slices/meansvariances.dat");
	meansvarsout_=fopen(name_,"a");
	
	fprintf(meansvarsout_,"%Lf",t);

	for(fld=0;fld<nflds;fld++)
	{
		av=0.;
		var=0.;
		// Calculate field mean
		LOOP
		av += field[0][fld][i][j][k];
		av = av/(long double)gridsize; // Convert sum to average
		av_sq = av*av; // Calculate mean squared for finding variance
		// Calculate variance
		LOOP
		var += field[0][fld][i][j][k]*field[0][fld][i][j][k] - av_sq;

		var = var/(long double)gridsize; // Convert sum to variance
		// Output means and variances to files
		fprintf(meansvarsout_," %Le %Le",av,var);
		// Check for instability. See if the field has grown exponentially and become non-numerical at any point.
		if((av!=0. && av/av!=1.)) // These two separate checks between them work on all the compilers I've tested
		{
			printf("Unstable solution developed. Field %d not numerical at t=%Le\n",fld,t);
			output_parameters();
			fflush(meansvarsout_);
			exit(1);
		}
	} // End of loop over fields
	av=0;
	var=0;
	LOOP //average of radial component
		av+=sqrt(pow(field[0][1][i][j][k],2)+pow(field[0][2][i][j][k],2)+pow(field[0][3][i][j][k],2)+pow(field[0][4][i][j][k],2)		);
		av = av/(long double)gridsize;
		av_sq=av*av;
	LOOP
		var += pow(field[0][1][i][j][k],2)+pow(field[0][2][i][j][k],2)+pow(field[0][3][i][j][k],2)	+pow(field[0][4][i][j][k],2)	 - av_sq;
		var = var/(long double)gridsize;
	fprintf(meansvarsout_," %Le %Le",av,var);

	fprintf(meansvarsout_,"\n");
	fclose(meansvarsout_);
}

int slicewaitf()
{
    return (int) (endtime/dt/slicenumber);//calculates number of timesteps to wait for slicenumber slices by time endtime
}

void outputslice()//externally called function for outputing the data from the run
{
    static FILE  *slicetime;
    static int first=0;
    static int lastslice,lastspec;
    static int slicewaitm;
	static FILE *masstime;
	static FILE *conservationtime;
    
    if(first==0)
    {
        if(slicewait==0)//determins number of steps to wait between each slice output based off of g2parameters.h
            slicewaitm=slicewaitf();
        else
            slicewaitm=slicewait;
        
        lastslice=slicewaitm;
        lastspec=specnumber;
    }
    
    lastslice++;
    
    if(lastslice>=slicewaitm)//this statement delays the output every slicewait*dt
    {
        
        
        
#if field_outdim!=0
        outputfield(first);//outputs field slice
#endif
        
#if spec_output==1//outputs spectra if last spectra output was not too soon (see g2paramters.h)
        lastspec++;
        if (lastspec>=specnumber) {
            specOut(first/specnumber);
            lastspec=0;
        }
#endif
        
#if var_output!=0
        meansvars();//outputs mean and variance of fields //TODO move mass plotting to here so it doesn't have to only plot with spectra?
#endif
        
        long double avPot=0;
		long double avKin=0;
		calcE0 ( &avKin, &avPot);

		
        /*times file output */
        //this routine outputs time, scalefactor, scalfactor derivative, hubble constant, and energy componets at each slice output
        slicetime=fopen("./slices/slices_time.dat","a");
        //slicetime=fopen("slices_time.dat","a");
		fprintf(slicetime,"%Le %Le %Le %Le %Le %Le %Le %Le \n", t,edkin[0],edpot[0],edgrad[0],edrho[0], avKin, avPot,a[0] ); //some zeroes to not mess up mathematica notebooks

		fclose(slicetime);

		//static FILE *masstime;
		masstime = fopen("./slices/mass_vs_time.dat", "a");
		fprintf(masstime, "%Le %Le  \n", t, Vrrcalc(0) );//mass_sq[1], mass_sq[2]);
		fclose(masstime);

		conservationtime = fopen("./slices/energy_conservation.dat", "a");
		fprintf(conservationtime, "%Le %Le %Le %Le  \n", t, edrho[0], 3*pow(a[0]/adot[0],-2), (edrho[0] -  3*pow(a[0]/adot[0],-2) )/(3*pow(a[0]/adot[0],-2)	)); //t, total energy, energy expression in friedmann equation
		fclose(conservationtime);
        
        first++;
        lastslice=0;
    }
}



void output_parameters()//this creats info.dat which contains information about run parameters and statistics
{
    static FILE *info;
    
    static int first=1;
    if(first) // At beginning of run output run parameters
    {
        char expanType[13];
#if expansion_type==1
        sprintf(expanType,"a dot");
#elif expansion_type==2
        sprintf(expanType,"a double dot");
#else
        sprintf(expanType,"no");
#endif
        
        info=fopen("info.txt","a");
        
        fprintf(info,"--------------------------\n");
        fprintf(info,"Model Specific Information\n");
        fprintf(info,"--------------------------\n");
        modelinfo(info);
        
        fprintf(info,"\n--------------------------\n");
        fprintf(info,"General Program Information\n");
        fprintf(info,"-----------------------------\n");
        fprintf(info,"Grid size=%d^3\n",N);
        fprintf(info,"L=%Le\n",L);
        fprintf(info,"dt=%Le, dx=%Le\n",dt,dx);
        fprintf(info,"end time=%Le\n",endtime);
        fprintf(info,"B=%Le\n",rescale_B);
        fprintf(info, "\nUsing %s expansion\n",expanType);
        fprintf(info, "\nUsing %d cores\n",tot_num_thrds);
        fprintf(info, "%d momentum bins for spectra\n",((int)(sqrt(3.)*N/2+1)));
        time(&tStart);
        fprintf(info,"\nRun began at %s",ctime(&tStart)); // Output date in readable form
        first=0;
    }
    else // If not at beginning, record elapsed time for run
    {
        time(&tFinish);
        fprintf(info,"Run ended at %s",ctime(&tFinish)); // Output ending date
        fprintf(info,"\nRun from t=%Le to t=%Le took ",starttime,t);
        readable_time((int)(tFinish-tStart),info);
        fprintf(info,"\n");
       // fprintf(info, "\nFinal scale factor is %Le\n",a[0]);
    }
    
    fflush(info);
    return;
}

// Convert a time in seconds to a more readable form and print the results
void readable_time(int tt, FILE *info)
{
    int tminutes=60,thours=60*tminutes,tdays=24*thours;
    
    if(tt==0)
    {
        fprintf(info,"less than 1 second\n");
        return;
    }
    
    // Days
    if(tt>tdays)
    {
        fprintf(info,"%d days",tt/tdays);
        tt = tt%tdays;
        if(tt>0)
            fprintf(info,", ");
    }
    // Hours
    if(tt>thours)
    {
        fprintf(info,"%d hours",tt/thours);
        tt = tt%thours;
        if(tt>0)
            fprintf(info,", ");
    }
    // Minutes
    if(tt>tminutes)
    {
        fprintf(info,"%d minutes",tt/tminutes);
        tt = tt%tminutes;
        if(tt>0)
            fprintf(info,", ");
    }
    // Seconds
    if(tt>0)
        fprintf(info,"%d seconds",tt);
    fprintf(info,"\n");
    return;
}



void screenout()//this calculates the time ellapsed from last screen output before outputting current program time.
{
    
    time(&tCurrent);
    
    if(screentime==0)
        return;
    
    if(tCurrent-tLast>screentime)
    {
        printf("%Lf\n",t);
        time(&tLast);
    }
}
