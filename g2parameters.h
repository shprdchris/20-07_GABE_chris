/****************************
 PARAMETERS FILE
 ****************************/

/*
 This headr file contains the adjustable and non adjustable global parameters for the GABE2.0 program.
 
 Copyright (2013): 
 Kenyon College
 John T. Giblin, Jr
 Tate Deskins and Hillary Child
 Last Updated: 06.27.2013
*/

/***************
 model constants
 ***************/
 #define num_flds 2// number of fields
const long double mphi=1.e-6;//mass of phi field
const long double phi0=0.193;//initial avg phi field value
const long double gsq=2.5e-7;//g^2 value for phi chi coupling
const long double f0[2]={phi0,0.};//array storing intial phi and chi field values
const long double df0[2]={-0.142231,0.};//array storing intial phi and chi field derivative values

const long double grav=1.; //the gravitational constant
const long double c=1.; //speed of light, shouldn't change, if you do, fix all equations


/***************************
 model independent parameters
 ***************************/
#define parallelize 1// for parallelization set to 1 and set other variables set to 0 for no parallelization
#define tot_num_thrds 4//total (max) number of threads to run durring program
const int randseed=44463132;//seed for rand number generator
const int N=128;//number of points along one side of grid
const long double L=20.;// length of one side of box in prgm units
const long double starttime=0.;//start time of simulation
const long double endtime=100.;//end time of simulations
const long double dt=0.01;//time step size
#define expansion_type 1//(0 for no expansion 1 for evolving from adot 2 for user defined expansion 
//(will need to adjust functions file (adot and such) and type two evolution in the step() function fnd g2init.cpp initexpansion() for user defined expansion )

/*************************
 model dependent parameters
 *************************/

const long double rescale_B=mphi;//rescallings

#define rand_init 1//1 to have random initialization 0 to not (see model file)
#define field_full_rand 1// 1 to have full random 0 to have symmetric kspace initilaization
// #define spec_cut_off sqrt(3.)/8.// do not define for no cut off or smoothing
// #define spec_smooth 0.5// determines the sharpness of the tanh window function

/****************
 output parameters
 *****************/
const long double screentime=60;// in seconds how frequently output prgm time to screen
const int slicewait=10;//how many dt's to wait between outputs (1 for no waiting) if 0 then slicenumber will be used.
const int slicenumber=20;//approx number of slices to output (only used if slicewati=0)
const int field_sliceskip=2;//howmany points to print in field profile (1 is every, 2 every two, 3 every three...)
const int specnumber=1; //how many spectra to out put (1= every output slice 2 every two....)
#define field_outdim 2// number of dimensions of output in field profile (0 for no output)
#define spec_output 1// 1 to ouptu spectra zero for no spectra output
#define var_output 1// 1 to output mean and varraince zero for no varriance output


/*********************************
 These are important DO NOT CHANGE
 *********************************/
const int nflds=num_flds; //stores number of fields for looping
const long double dx=L/((long double) N);//stores the change in x from point to point
const long double gridsize=N*N*N;//stores size of grid for averaging


#if parallelize!=1
#undef tot_num_thrds
#define tot_num_thrds 1
#endif

