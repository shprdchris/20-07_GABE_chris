/**********************************
 SPECTRA FILE
 **********************************/

/*
 This file contains the functions needed to output spectra.
 
 Copyright (2013): 
 Kenyon College
 John T. Giblin, Jr
 Tate Deskins and Hillary Child
 Last Updated: 06.27.2013
*/

#include "g2header.h" //contains declerations for program functions.

#if spec_output==1

static long double *indft, *dindft; //decleares input field for dft
static fftwl_complex *outdft, *doutdft; //decleares output field for dft
static fftwl_plan spec_plan, dspec_plan; //decleares dft plan for fftw

#define PHI field[s][0]
#define H1 field[s][1]
#define H2 field[s][2]
#define PHIDOT dfield[s][0]
#define H1DOT dfield[s][1]
#define H2DOT dfield[s][2]

void specOut(int first)
{   
    DECLARE_INDEX
    const long double spec_norm=pow(N,-6); //normalization for SQUARE of FT of fields (DFT is not normalized)
    int px,py,pz;//trackes real place in momentum grid
    long double pmagnitude;//stores the magnitude of the momentum vector for current point
    const int numbins=((int)(sqrt(3.)*N/2+1));//number of spectra bins based off of size of the grid
    long double f2_power[nflds][numbins]	, df2_power[nflds][numbins],ndensity[nflds][numbins],rhodensity[nflds][numbins],spec_numpts[nflds][numbins];// holds the power and number of points in each spectra bin
    float dp=2.*M_PI/L;
	long double omegasq[nflds][numbins],omega, omsq;
	long double mass_sq[nflds];
	long double extra_np= pow(rescale_B, -2)*pow(L,3); //see initialization
	long double extra_rhop=pow(rescale_B,2)*extra_np/(2*M_PI*M_PI);
	const long double norm_tot =spec_norm; 

	long double ff, dff;
	long double mdiagV[nflds];
	long double evecV[nflds][nflds];
	long double rhoTot =0;
	long double rhoRef = 0;

    if(first==0)
    {
        
        fftwl_plan_with_nthreads(tot_num_thrds);//lets fftw know the number of available threads (tot_num_thrds)
		fftwl_plan_with_nthreads(tot_num_thrds);//tells fftw that tot_num_thrds are available for use durring dft

        indft = (long double *) fftwl_malloc(sizeof(long double) * N * N * N);//allocates memory for the input array
		dindft = (long double *) fftwl_malloc(sizeof(long double) * N * N * N);//allocates memory for the input array

        outdft =  (fftwl_complex*) fftwl_malloc(sizeof(fftwl_complex) * N * N * (N/2+1));//allocates memory for the output array
		doutdft =  (fftwl_complex*) fftwl_malloc(sizeof(fftwl_complex) * N * N * (N/2+1));//allocates memory for the output array

        spec_plan = fftwl_plan_dft_r2c_3d(N, N, N, indft, outdft, FFTW_MEASURE);// defines the plan for fftw (moving from configuration space to momentum space
		dspec_plan = fftwl_plan_dft_r2c_3d(N, N, N, dindft, doutdft, FFTW_MEASURE);// defines the plan for fftw (moving from configuration space to momentum space
        
    }
    
	fldLOOP{ //extract field values 
		Sfield[0][fld][i][j][k] = field[0][fld][i][j][k];
		Sdfield[0][fld][i][j][k] = dfield[0][fld][i][j][k];
	}

	if (first==0)
	{
		mass_matrix_calc(0, &evecV, &mdiagV, true); //a^2 in a^2 m^2 is implmented inside this step
	}
	else
	{
		//printf("spectra, no rotation");
		mass_matrix_calc(0, &evecV, &mdiagV, true);
	}

	for (fld=0;fld<nflds;fld++){
		mass_sq[fld] = mdiagV[fld];
	}

	if(first==0)
	{
		toDiag (Sfield, Sdfield, &evecV, true); //convert (phi, h_1, h_2)->(phi,H_1,H_2)->mass eigenstates
	}
	else
	{
		toDiag  (Sfield, Sdfield, &evecV, true); //convert (phi, h_1, h_2)->(phi,H_1,H_2), no rotation
	}

    for(fld=0;fld<nflds;fld++){ 
        LOOP//this loop asigns the field values to the input array
        {
                indft[k + N*(j + N*i)]=Sfield[0][fld][i][j][k];//TODO: add fourier transform of radial
				dindft[k + N*(j + N*i)]=Sdfield[0][fld][i][j][k]; //CHRIS
        }
        
        fftwl_execute(spec_plan);//performs the dft 
		fftwl_execute(dspec_plan);
        
        for(i=0;i<numbins;i++){//initialize the power and number of points in each bin to 0
            f2_power[fld][i]=0.0;
			df2_power[fld][i]=0.0;
			ndensity[fld][i]=0.0;
			rhodensity[fld][i]=0.0;
            spec_numpts[fld][i]=0;
        }
        
        for(i=0;i<N;i++) {//this is the loop over momentum space
            px=(i<=N/2 ? i: i-N);//calculates the real momentum space position for the x direction
            for(j=0; j<N;j++) {
                py = (j<=N/2 ? j : j-N); //calculates the real momentum space position for the y direection
                for(k=1;k<N/2;k++) { //since the k=0 and k=n/2+1 are un matched we treat the seperately
                    pz=k;//calculates the real momentum space position for the z direction
                    pmagnitude = sqrt(pw2(px)+pw2(py)+pw2(pz));//calculates the magnitude of the momentum of the point
					ff= 2.*(pw2(outdft[k + (N/2+1)*(j + N*i)][0])+pw2(outdft[k + (N/2+1)*(j + N*i)][1]));
				
					dff = 2.*(pw2(doutdft[k + (N/2+1)*(j + N*i)][0])+pw2(doutdft[k + (N/2+1)*(j + N*i)][1]));
                    f2_power[fld][(int) pmagnitude]+=ff; //adds the power in this mode to the proper bin times two 
					
					

					df2_power[fld][(int) pmagnitude]+=dff; //adds the power in this mode to the proper bin times two 
					omsq =  fabs(pw2(dp*pmagnitude) + mass_sq[fld]); //omega^2 //already done rescaling a^2 m^2
					omegasq[fld][(int) pmagnitude] =omsq;
					if(omsq > 0.) // Guard against floating point errors
					  {
						omega = sqrt(omsq); // Omega = Sqrt(p^2 + m^2) 
						ndensity[fld][(int)pmagnitude] += 2.*extra_np*(ff*omega + dff/omega); // Add this mode's contribution to number density
						
						
					  }
					//  rhodensity[fld][(int)pmagnitude] += 2.*extra_rhop*(ff*omsq + dff);
					  //TODO add to another rhodensity that doesn't get divided by nbins. normalize to compare with rho_tot
					  if (i != 0){
						 rhoTot += 0.5*norm_tot*(ff*omsq + dff); //TODO: doesn't include dp properly
					  }
				
                    //since the symmetry in the dft we only need half of the z direction box (see fftw documentation)
                    spec_numpts[fld][(int) pmagnitude]+=2;//adds two points for the bin
                }
                //this is the same procedure as seen above with out the multiplication by 2 since the zero mode and the k=N/2 modes are unmatched
                k=0;
                pz=k;
                pmagnitude = sqrt(pw2(px)+pw2(py)+pw2(pz));
				ff=(pw2(outdft[k + (N/2+1)*(j + N*i)][0])+pw2(outdft[k + (N/2+1)*(j + N*i)][1]));;
                f2_power[fld][(int) pmagnitude]+= ff;//[0]:real, [1]:imaginary
				dff = (pw2(doutdft[k + (N/2+1)*(j + N*i)][0])+pw2(doutdft[k + (N/2+1)*(j + N*i)][1])); //[0]:real, [1]:imaginary
				df2_power[fld][(int) pmagnitude]+=dff;

				omsq = fabs(pw2(dp*pmagnitude) + mass_sq[fld]);
				omegasq[fld][(int) pmagnitude]=omsq;
					if(omsq > 0.) // Guard against floating point errors
					  {
						omega = sqrt(omsq); // Omega = Sqrt(p^2 + m^2)
						ndensity[fld][(int)pmagnitude] +=2.*extra_np*(ff*omega + dff/omega); //TODO: is the factor of 2 in front right? Think it is - already divided |f|^2 elsewhere
						
					  }
					 // rhodensity[fld][(int)pmagnitude] += 2.*extra_rhop*(ff*omega*omega + dff);
					  if (i!=0){ //TODO check boolean syntax
						 rhoTot += 0.5*norm_tot*(ff*omsq + dff); //TODO: add to total rhodensity w/o factor of 2
					  }

                spec_numpts[fld][(int) pmagnitude]+=1;
                k=N/2;
                pz=k;
                pmagnitude = sqrt(pw2(px)+pw2(py)+pw2(pz));
              	ff=(pw2(outdft[k + (N/2+1)*(j + N*i)][0])+pw2(outdft[k + (N/2+1)*(j + N*i)][1]));;
                f2_power[fld][(int) pmagnitude]+= ff;//[0]:real, [1]:imaginary
				dff = (pw2(doutdft[k + (N/2+1)*(j + N*i)][0])+pw2(doutdft[k + (N/2+1)*(j + N*i)][1])); //[0]:real, [1]:imaginary
				df2_power[fld][(int) pmagnitude]+=dff;
				//if ( (int) pmagnitude ==1){printf("df2_power (growing) = %Le \n", df2_power[fld])};
				omsq = fabs(pw2(dp*pmagnitude) + mass_sq[fld]);
				omegasq[fld][(int) pmagnitude]=omsq;
					if(omsq > 0.) // Guard against floating point errors
					  {
						omega = sqrt(omsq); // Omega = Sqrt(p^2 + m^2)
						ndensity[fld][(int)pmagnitude] += 2.*extra_np*(ff*omega + dff/omega); //TODO: is the factor of 2 in front right? Think it is - already divided |f|^2 elsewhere
					  }
					//  rhodensity[fld][(int)pmagnitude] +=2.*extra_rhop*(ff*omsq + dff);
					  if (i!=0){ //TODO check boolean syntax
						 rhoTot += 0.5*norm_tot*(ff*omsq + dff); //TODO: add to total rhodensity w/o factor of 2
					  }

                spec_numpts[fld][(int) pmagnitude]+=1;
                
            }
        }
        
        for(i=0;i<numbins;i++){ //convert the sums in each bin to averages
            if(spec_numpts[fld][i]>0){
			
					f2_power[fld][i]=f2_power[fld][i]*spec_norm/((long double) spec_numpts[fld][i]);
					
					
					df2_power[fld][i]=df2_power[fld][i]*spec_norm/((long double) spec_numpts[fld][i]);
					ndensity[fld][i]= pow(a[0],3)*0.25 *ndensity[fld][i]*spec_norm/((long double) spec_numpts[fld][i]); //0.5 for general case, another 0.5 because includes forward and backward
					rhodensity[fld][i]=pow(rescale_B,2)/(2*M_PI*M_PI)*omega*ndensity[fld][i];//rhodensity[fld][i]*spec_norm/((long double)/spec_numpts[fld][i]);
			}	
            else
                f2_power[fld][i]=0.;
        }	
	} //ends loop over fields

	long double V11p=0; //to print //TODO make neater
	long double V22p=0;
	long double V12p=0;
	long double Vrrp = 0;
	//Viicalc(0, &V11p, &V22p, &V12p, &Vrrp);
	
    //Print the spectra to a file
	static FILE *inhomE; //inhomogeneous energy density via fourier transform
	inhomE = fopen("./slices/inhomogeneous_energy.dat", "a");
	fprintf(inhomE,"%Le %Le \n", t, rhoTot);
	fclose(inhomE);

	/*static FILE *masstime;
	masstime = fopen("./slices/mass_vs_time.dat", "a");
	fprintf(masstime, "%Le %Le %Le %Le %Le %Le %Le \n", t,mass_sq[0], V11p, V22p, V12p, Vrrp, gamW);//mass_sq[1], mass_sq[2]);
	fclose(masstime);*/
    
    static FILE *slicespectra;
    static char name[500];
    
    sprintf(name,"./slices/slices_spectra_%d.dat", first);
    slicespectra=fopen(name,"w");

	
	

    for(i=0;i<numbins;i++)
    {
        fprintf(slicespectra,"%Le", 2*M_PI*i/L);//this prints the conformal mode k
        for(fld=0;fld<nflds;fld++){
		//	printf("massSq[2] = %Le \n", mass_sq[2]);
				fprintf(slicespectra," %Le %Le %Le %Le %Le %Le %Le %Le", omegasq[fld][i]*f2_power[fld][i], df2_power[fld][i],mass_sq[fld],omegasq[fld][i],ndensity[fld][i],rhodensity[fld][i], pow(2*M_PI*i/L, 2)*ndensity[fld][i], pow(2*M_PI*i/L, 2)*rhodensity[fld][i] ); //index i is proportional to rounded |p|
		}
        fprintf(slicespectra,"\n");
    }
	for (int fld=0;fld<nflds;fld++){
		for (int rr=0;rr<nflds;rr++){
				fprintf(slicespectra,"%Le ", evecV[rr][fld]);
		} //loop over compoments
	} //loop over eigenvalue index
	fprintf(slicespectra, "\n");

    fclose(slicespectra);
	}


void specClear(){//clears the memory used in the dft's
    fftwl_destroy_plan(spec_plan);
    fftwl_free(indft);
    fftwl_free(outdft);

	fftwl_destroy_plan(dspec_plan);
    fftwl_free(dindft);
    fftwl_free(doutdft);
}

#endif
