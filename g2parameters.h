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
 #define num_flds 5// number of fields
//const long double mphi=1.e-6;//mass of phi field

const long double AA=1; //shift of CMB normalization
const long double l = 0.01; //lambda

//ALL THE PARAMETERS I WANT TO PASTE IN
const long double b= 1.06*pow(10,9);
const long double xi1 =sqrt( (2*pow(10,9)-b)*l	); 
const long double f0[num_flds]={0.176169,0.00487158, 0};
const long double df0[num_flds]={-0.362985 ,-0.00538003 ,0}; 
const long double kmax = 1000;
const long double kcut =125;;
const long double dt=   2*M_PI/kmax/40;
const int N=64;
const long double L= 2*M_PI*N/kmax;
const long double specMod= 1;

const bool doGamma=false;
//const int gamSkip=10;
//TODO: gamma implementation looks wrong: decay rate too slow cf 1 angular mode
//TODO: remove gamma: change all ifgammas to false, always

//OTHER PARAMETERS TO MODIFY
const long double starttime=0 ;//start time of simulationo
const long double endtime=4000;//end time of simulations //tHubble=20
#define expansion_type 1//(0 for no expansion 1 for evolving from adot 2 for user defined expansion

//const long double xi1 =sqrt(	(AA*AA*2*pow(10,9)-b)*l	);  //xi
const long double mm2 = l/( 3*(xi1*xi1 + l*b) ); //m2 of inflaton for phi_0 >0: set this timescale to unity
const long double mm2inv = 1/mm2; //for convenience
const long double g = 0.560499; // electroweak g
const long double aW =0.025; //alpha_weak


const double sqrt23 = 0.81649658092; //sqrt(2/3)
const double sqrt6 =2.44948974278; //sqrt(6)

const long double grav=1.; //the gravitational constant
const long double c=1.; //speed of light, shouldn't change, if you do, fix all equations


/***************************
 model independent parameters
 ***************************/
#define parallelize 1// for parallelization set to 1 and set other variables set to 0 for no parallelization
#define tot_num_thrds 8//total (max) number of threads to run durring program
const int randseed=44463132;//seed for rand number generator

 
//(will need to adjust functions file (adot and such) and type two evolution in the step() function fnd g2init.cpp initexpansion() for user defined expansion )

/*************************
 model dependent parameters
 *************************/

const long double rescale_B= sqrt(mm2)/**pow(f0[0],-1.+beta/2.)*/;//rescallings

//const bool iftheta =true; //1: do theta properly. 0: set theta=0 artificially

#define rand_init 1//1 to have random initialization 0 to not (see model file)
#define field_full_rand 1// 1 to have full random 0 to have symmetric kspace initilaization

#define spec_cut_off kcut/kmax// do not define for no cut off or smoothing
//#define spec_frac_width 1
//#defin    spec_cut_off*(1-spec_frac_width);
 #define spec_smooth 100// determines the sharpness of the tanh window function

/****************
 output parameters
 *****************/
const long double sliceOutTime = .1;
const long double specOutTime = 1;


const long double screentime=10;// in seconds how frequently output prgm time to screen
const int slicewait= (int) (sliceOutTime / dt);   ;//how many dt's to wait between outputs (1 for no waiting) if 0 then slicenumber will be used.
const int slicenumber=1;//approx number of slices to output (only used if slicewati=0)
const int field_sliceskip=1;//howmany points to print in field profile (1 is every, 2 every two, 3 every three...) =(
const int specnumber= (int)(specOutTime / (slicewait*dt)  )	; //how many spectra to out put (1= every output slice 2 every two....)

#define field_outdim 0// number of dimensions of output in field profile (0 for no output)
#define spec_output 1// 1 to output spectra, zero for no spectra output
#define var_output 1// 1 to output mean and varraince, zero for no varriance output

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

