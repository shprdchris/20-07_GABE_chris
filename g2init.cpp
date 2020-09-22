/**********************************
 Initialization FILE
 **********************************/

/*
 This header file contains all the functions which are indeendent of the model needed to initialize the fields and random initial conditions.
 
 Copyright (2013): 
 Kenyon College
 John T. Giblin, Jr
 Tate Deskins and Hillary Child
 Last Updated: 06.27.2013
*/

#include "g2header.h" //contains declerations for program functions.

/** expansion initialization **/
void initexpansion()
{
#if expansion_type==0
    //for power law
    a[0]=1;
    adot[0]=0;
    a[1]=1;
    adot[1]=0;
    calcEnergy(0);
#elif expansion_type==1
	a[0]=1;
    calcEnergy(0);
    adot[0]=adf(0); 
#elif expansion_type==2
//for some other user defined expansion type
#endif
}

/** global parameters needed for fftwl**/
#if rand_init==1

long double *inic;//holds the random initial conditions generated in configuration space
fftwl_complex *fkf, *fkd; //holds the generated initial conditions in momentum space for the field (fkf) and its derivative (fkd)
fftwl_plan picf, picd;// stores the plans for the dft to move from momentum space to configuration space


void dftMemAlloc()
{
    fftwl_plan_with_nthreads(tot_num_thrds);//tells fftw that tot_num_thrds are available for use durring dft
	fftwl_plan_with_nthreads(tot_num_thrds);//tells fftw that tot_num_thrds are available for use durring dft
    
    inic = (long double *) fftwl_malloc(sizeof(long double) * N * N * N);//allocates memory for inic
    
    fkf   = (fftwl_complex*) fftwl_malloc(sizeof(fftwl_complex) * N * N * (N/2+1));//allocates memory for fkf
    fkd   = (fftwl_complex*) fftwl_malloc(sizeof(fftwl_complex) * N * N * (N/2+1));//allocates memory for fkd
    
    picf = fftwl_plan_dft_c2r_3d(N, N, N, fkf, inic, FFTW_MEASURE);// defines the fftw plan for the field value dft
    picd = fftwl_plan_dft_c2r_3d(N, N, N, fkd, inic, FFTW_MEASURE);// defines the fftw plan for the derivative of the field value dft
}

/**function to initialize random initial conditions**/

void randInit(long double Rf[][N][N],long double Rdf[][N][N],long double d2vdf2)// function that implements the random (or pseudo-random initial conditions for the field values 
	//assumes we're working in (L,H) eigenstates throughout
	
{
    static int first=0;//counter for if statements
    int i,j,k,tester;//iterators for the for loop
    long double i2,j2,k2, wk, randVar, randAngle;//scaled momentum space position squared(i2,j2,k2) omega (wk) and random variable and random angle
    long double dk2 = 4.*M_PI*M_PI/L/L;//the square of the distance of gridpoints in momentum space (delta k)^2
    long double fkRI1[2], fkRI2[2];//these hold part of the initial conditions durring calculation
    
    
    long double H0=sqrt(8.*M_PI*grav/3.*edrho[0]);//the hubble parameter
    if(first==0)
    {
        dftMemAlloc();//allocates memory for the dft
        first++;
        srand(randseed);//seeds the random number generator
    }
    
    /* first we calculate phi(0, x) and phidot(0, x)	*/
    
//this directive opens tot_num_thrds of threads each with the private variables j,k,i2,j2,k2
//since this is a for directive it splits up amoungs various i values (so i is default private)
#pragma omp parallel for private (j,k,i2,j2,k2,tester,wk,randVar,randAngle, fkRI1,fkRI2) num_threads (tot_num_thrds)
    for(i=0; i<N; i++){
        i2 = (i<(N/2+1) ? (i*i) : ((N-i)*(N-i)));//define square i position
        
        for(j=0; j<N; j++){
            j2 = (j<(N/2+1) ? (j*j) : ((N-j)*(N-j)));//define square j position
            
            for(k=0; k<N/2+1; k++){
                k2 = (k*k);//define square k position
                tester=i+j+k;//note only 0 when i,j and k are zero
                switch (tester) {//this initializes the zero mode to 0
                    case 0:
                        
                        fkf[k + (N/2+1)*(j + N*i)][0] = 0.;
                        fkf[k + (N/2+1)*(j + N*i)][1] = 0.;
                        
                        
                        fkd[k + (N/2+1)*(j + N*i)][0] = 0.;
                        fkd[k + (N/2+1)*(j + N*i)][1]= 0.;
                        break;
                        
                    default://if tester is not zero the following occurs                     
#if field_full_rand==1
                        
                        wk = sqrt(dk2*((long double)i2 + (long double)j2 + (long double)k2) + d2vdf2);//dk*i\doti+d^2v/df^2 
                                               
# ifdef spec_cut_off
                        /* the tanh part at the end is the mode cut off, this effectively damps out high frequency modes to 0*/
                        
                        randVar = (long double)rand()/(long double)RAND_MAX;//this is a random variable generated for each mode
                        randAngle = 2.*M_PI*(long double)rand()/(long double)RAND_MAX;//random angle to add arbitrary phase for each mode
                        
                        fkRI1[0] = rescale_B*cos(randAngle)*sqrt(-log(randVar)/(2.*wk*L*L*L))*(1.-tanh(spec_smooth*(sqrt(i2+j2+k2)-N*spec_cut_off)))/2.;//long double
                        //"forward" moving real wave
                        
                        fkRI1[1] = rescale_B*sin(randAngle)*sqrt(-log(randVar)/(2.*wk*L*L*L))*(1.-tanh(spec_smooth*(sqrt(i2+j2+k2)-N*spec_cut_off)))/2.;//imaginary
                        //"forward" moving imaginary wave    
                        
                                            
                        randVar = (long double)rand()/(long double)RAND_MAX;//unlike lattice easy (becuse that was an error
                        randAngle = 2.*M_PI*(long double)rand()/(long double)RAND_MAX;


                        //"backward" moving waves
                        fkRI2[0] = rescale_B*cos(randAngle)*sqrt(-log(randVar)/(2.*wk*L*L*L))*(1.-tanh(spec_smooth*(sqrt(i2+j2+k2)-N*spec_cut_off)))/2.;//long double //last bit is window function
                        fkRI2[1] = rescale_B*sin(randAngle)*sqrt(-log(randVar)/(2.*wk*L*L*L))*(1.-tanh(spec_smooth*(sqrt(i2+j2+k2)-N*spec_cut_off)))/2.;//imaginary
                        
# else
                        /*This is the same as the above with out any mode damping*/
                        randVar = (long double)rand()/(long double)RAND_MAX;
                        randAngle = 2.*M_PI*(long double)rand()/(long double)RAND_MAX;
                        
                        fkRI1[0] = rescale_B*cos(randAngle)*sqrt(-log(randVar)/(2.*wk*L*L*L));//long double
                        fkRI1[1] = rescale_B*sin(randAngle)*sqrt(-log(randVar)/(2.*wk*L*L*L));//imaginary
                         
                        randVar = (long double)rand()/(long double)RAND_MAX;//unlike lattice easy (becuse that was an error
                        randAngle = 2.*M_PI*(long double)rand()/(long double)RAND_MAX;

                        fkRI2[0] = rescale_B*cos(randAngle)*sqrt(-log(randVar)/(2.*wk*L*L*L));//long double
                        fkRI2[1] = rescale_B*sin(randAngle)*sqrt(-log(randVar)/(2.*wk*L*L*L));//imaginary
# endif
                        
                        /*This sums the generated forward and backwards waves to make the real and imaginary parts of the momentum space field and field derivatives*/
                        fkf[k + (N/2+1)*(j + N*i)][0] = (fkRI1[0] + fkRI2[0])/sqrt(2.);//long double
							
                        fkf[k + (N/2+1)*(j + N*i)][1] = (fkRI1[1] + fkRI2[1])/sqrt(2.);//imaginary

                        fkd[k + (N/2+1)*(j + N*i)][0] = (fkRI2[1] - fkRI1[1])*wk/sqrt(2.) - H0*fkf[k + (N/2+1)*(j + N*i)][0];//real  
                        fkd[k + (N/2+1)*(j + N*i)][1]= (fkRI1[0] - fkRI2[0])*wk/sqrt(2.) -   H0*fkf[k + (N/2+1)*(j + N*i)][1];//imaginary


#endif
#if field_full_rand==0
                      //this is the same procedure as above but instead of randomly genreating the modes it uses the mean of the distribution making a pseudo random distribution in configuration space  
                        wk = sqrt(dk2*((long double)i2 + (long double)j2 + (long double)k2) + d2vdf2);//dk*i\doti+d^2v/df^2
                        
                        
                        fkf[k + (N/2+1)*(j + N*i)][0] = rescale_B*sqrt(1/(2.*wk*L*L*L))*(1.-tanh(spec_smooth*(sqrt(i2+j2+k2)-N*spec_cut_off)))/2.;//long double
                        fkf[k + (N/2+1)*(j + N*i)][1] = 0.;//(fkRI1[1])*sqrt(2.);//imaginary
                        
                        
                        fkd[k + (N/2+1)*(j + N*i)][0] = - H0*fkf[k + (N/2+1)*(j + N*i)][0];//long double
                        fkd[k + (N/2+1)*(j + N*i)][1] = 0.;//(fkRI1[0])*wk*sqrt(2.) -   H0*fkf[k + (N/2+1)*(j + N*i)][1];//imaginary //commented out in original file
                        
#endif
                        
                        
                        break;
                }
            }
        }
    }
	 
    fftwl_execute(picf);//this envokes the dft for the field going from momentum to configuration
	
#pragma omp parallel for private (j,k) num_threads (tot_num_thrds) 
    LOOP
    { 
        Rf[i][j][k] += inic[k + N*(j + N*i)];// this is the addition of the random initial conditions to that dictated in the model file
        
    }
    
    fftwl_execute(picd);//this envokes the dft for the field going from momentum to configuration

#pragma omp parallel for private (j,k) num_threads (tot_num_thrds)
		
    LOOP
    {
        Rdf[i][j][k] += inic[k + N*(j + N*i)];// this is the addition of the random initial conditions to that dictated in the model file 
    }
	
}


/**function to free up memory used for the dft**/

void initDestroy()
{
	fftwl_destroy_plan(picf);
    fftwl_destroy_plan(picd);
    fftwl_free(inic);
    fftwl_free(fkf);
    fftwl_free(fkd);
}
#endif
