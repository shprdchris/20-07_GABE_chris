/*********************************
 FUNCTIONS FILE
 **********************************/

/*
 Tis header file contains all the functions which are independent of the model needed to evolve the fields by the Runge-Kutta Second order method (for first order finite derivatives). The incr and decr commands set periodic boundary conditions on the lattice durring evolution.
 
 Copyright (2013): 
 Kenyon College
 John T. Giblin, Jr
 Tate Deskins and Hillary Child
 Last Updated: 06.27.2013
*/

#include "g2header.h" //contains declerations for program functions.

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

long double pw2(long double x)//squares long doubles
{
    return x*x;
}

/*long double mod(long double f1, long double f2)
{
	return pow( (pow(f1,2), + pow(f2,2)),0.5 );
}
long double mod2(long double f1, long double f2)
{
	return (pow(f1,2), + pow(f2,2));
}*/

inline int incr(int i)//for periodic boundaries
{
    return i==N-1? 0: i+1;
}

inline int decr(int i)//for periodic boundaries
{
    return i==0? N-1: i-1;
}

/** Laplacian Functions **/


long double laplacian(long double f[][N][N], int i, int j, int k)//this calculates the seven point laplacian 
{
    
    return (f[incr(i)][j][k]+f[decr(i)][j][k]
            +f[i][incr(j)][k]+f[i][decr(j)][k]
            +f[i][j][incr(k)]+f[i][j][decr(k)]
            -6.*f[i][j][k])/(dx*dx);
    
}

/** Spatial Derivative Functions **/
/*
The following fucntions are used to calculate the gradient energy of the field only; but they can also be implemented if derivative couplings are desired*/

//If you know which partial derivative you need

long double dfdi(long double f[][N][N], int i, int j, int k)// spatial derivative of the field f in i (x) direction direction
{
    return (f[incr(i)][j][k]-f[decr(i)][j][k])/2./dx;
}

long double dfdj(long double f[][N][N], int i, int j, int k)// spatial derivative of the field f in j (y) direction
{
    return (f[i][incr(j)][k]-f[i][decr(j)][k])/2./dx;
}

long double dfdk(long double f[][N][N], int i, int j, int k)// spatial derivative of the field f in k (z) direction
{
    return (f[i][j][incr(k)]-f[i][j][decr(k)])/2./dx;
}

//If you want to loop over spatial derivatives this form is somewhat more convienent
long double dfdx(long double f[][N][N], int x, int i, int j, int k)//spatial derivative of the field f in the "x" direction.
{
    switch (x)
    {
        case 0:
            return dfdi(f,i,j,k);
        case 1:
            return dfdj(f,i,j,k);
        case 2:
            return dfdk(f,i,j,k);
        default:
            return 0;
    }
}

/**Functions needed for self-consistent expansion**/

long double gradF2(long double f[][N][N],int i,int j,int k){
    
    return  dfdi(f,i,j,k)*dfdi(f,i,j,k)+dfdj(f,i,j,k)*dfdj(f,i,j,k)+dfdk(f,i,j,k)*dfdk(f,i,j,k);//this is the unscaled gradient fo the field at the point i,j,k
    
}

long double avgGrad(int s) //Find the average gradient energy
{
    long double grad=0;
    DECLARE_INDEX
    for(fld=0; fld<nflds; fld++){
#pragma omp parallel for private (j,k) reduction(+:grad) num_threads (tot_num_thrds)
        LOOP//loops over i,j,k
        {
			if (fld==0){
					grad-=field[s][fld][i][j][k]*laplacian(field[s][fld],i,j,k);//sums the gradient energy at each point
			}
			else{
				grad-=field[s][fld][i][j][k]*laplacian(field[s][fld],i,j,k);//*exp(-sqrt23*field[s][0][i][j][k])	;//sums the gradient energy at each point
			}
        }
    }
    return grad/gridsize/2./a[s]/a[s];//divides by the gridsize (to normalize) and 1/(2a^2) to get the gradient energy density
}

long double avgPot(int s) //Find the average potential energy
{
    long double pot=0;
    DECLARE_INDEX
    #pragma omp parallel for private (j,k) reduction(+:pot) num_threads (tot_num_thrds)
        LOOP//loops over i,j,k
        {
            pot+=potential(s,i,j,k);//sums the potential at every point
        }
    
    return pot/gridsize;//averages over the grid
}


long double avgKin(int s) //Find the average kinetic energy
{
    long double kin=0;
    DECLARE_INDEX
    for(fld=0; fld<nflds; fld++){
#pragma omp parallel for private (j,k) reduction(+:kin) num_threads (tot_num_thrds)
        LOOP//loops over i,j,k
		{
        {if (fld==0){
				 kin+=dfield[s][fld][i][j][k]*dfield[s][fld][i][j][k];//sums the square field derivative at every point
		}
		else{
				kin+=dfield[s][fld][i][j][k]*dfield[s][fld][i][j][k];//*exp(-sqrt23*field[s][0][i][j][k]);
			}
        }
	 }
	 } //field loop
    return kin/gridsize/2.;//divide by the grid size to get the average and 2
}

/*long double avgKinRho(int s) //Find the average kinetic energy
{
    long double kin=0;
    DECLARE_INDEX
	#pragma omp parallel for private (j,k) reduction(+:kin) num_threads (tot_num_thrds)
        LOOP//loops over i,j,k
		{
	    kin+=0.5*exp(-sqrt23*field[0][0][i][j][k])*pow(	(field[0][1][i][j][k]*dfield[0][1][i][j][k]+field[0][2][i][j][k]*dfield[0][2][i][j][k]+field[0][3][i][j][k]*dfield[0][3][i][j][k] +field[0][4][i][j][k]*dfield[0][4][i][j][k]  )/sqrt( pow(field[0][1][i][j][k],2)+pow(field[0][2][i][j][k],2)+pow(field[0][3][i][j][k],2) +pow(field[0][4][i][j][k],2)         ) , 2);
	 }
    return kin/gridsize;//divide by the grid size to get the average and 2
}*/

/*long double avgKinH(int s) //Find the average kinetic energy
{
    long double kin=0;
    DECLARE_INDEX
    for(fld=1; fld<nflds; fld++){
#pragma omp parallel for private (j,k) reduction(+:kin) num_threads (tot_num_thrds)
        LOOP//loops over i,j,k
		{
				kin+=dfield[s][fld][i][j][k]*dfield[s][fld][i][j][k]*exp(-sqrt23*field[s][0][i][j][k]);	
          }
	 } //field loop
    return kin/gridsize/2.;//divide by the grid size to get the average and 2
}*/

void calcEnergy(int s) //Calculate the total energy
{
	edkin[s]=avgKin(s);
    edpot[s]=avgPot(s);
    edgrad[s]=avgGrad(s);
    edrho[s]=edkin[s]+edpot[s]+edgrad[s];
	//edkinH[s]=avgKinH(s);
	//edkinRho[s]=avgKinRho(s);
}

long double calcdda(int s) //calculates contribution to rk2
{
	return -( adot[s]*adot[s]/a[s] + a[s]*(edkin[s]-edgrad[s]/3-edpot[s]) )/2;
}

void calcE0 ( long double *kin, long double* pot)
{
	DECLARE_INDEX
	long double av;
	long double dav;
	long double avArray[nflds];
	*kin = 0;
	*pot = 0;
	for (int fld=0;fld<nflds;fld++){
		av=0;
		dav=0;
		LOOP{
			av+=field[0][fld][i][j][k];
			dav+=dfield[0][fld][i][j][k];
		} //points
		av/=(long double)gridsize;
		dav/=(long double) gridsize;
		//printf("dav^2(%i) = %Le\n", fld, dav*dav);
		avArray[fld] = av;
		if (fld==0){
			*kin +=0.5*dav*dav;
		}
		else if (fld > 0){
			*kin+= 0.5*dav*dav*exp(-sqrt23*avArray[0]);
		}
	} //loop over fields
	*pot=3*(pow(xi1,2) + l*b)/pow(xi1,2)* (3*pow(xi1,2)*(pow(  ( pow(avArray[1],2) + pow(avArray[2],2) + pow(avArray[3],2) +  + pow(avArray[4],2) ) ,2  )*l + pow(-1 + exp(sqrt23*avArray[0]) - ( pow(avArray[1],2) + pow(avArray[2],2) + pow(avArray[3],2) +  + pow(avArray[4],2) )*xi1,2)/b))/(4.*exp(2*sqrt23*avArray[0])*l);
}

void calckin0H ( long double *kin)
{
	DECLARE_INDEX
	long double av;
	long double dav1;
	long double dav2;
	long double dav3;
	long double dav4;

	av =0;
	dav1=0;
	dav2=0;
	dav3=0;
	dav4=0;

	*kin = 0;
		LOOP{
			av+=field[0][0][i][j][k];
			dav1+=dfield[0][1][i][j][k];
			dav2+=dfield[0][2][i][j][k];
			dav3+=dfield[0][3][i][j][k];
			dav4+=dfield[0][4][i][j][k];
		} //points
		av/=(long double)gridsize;
		dav1/=(long double) gridsize;
		dav2/=(long double) gridsize;
		dav3/=(long double) gridsize;
		dav4/=(long double) gridsize;

		*kin = 0.5*exp(-sqrt23*av)*( dav1*dav1 + dav2*dav2 +dav3*dav3 + dav4*dav4);
}

long double adf(int s)//the friedman equation
{
    return sqrt(/*8.*M_PI* grav*/1/3.*edrho[s])*a[s]; //using the reduced planck mass throughout
}


long double ddfield( int s, int fld, int i, int j, int k)//evaluates the double time derivative of the field fld (s) at i,j,k. //TODO wrong in curved space. can't just uncomment: would pick up extra bits from non-canonical.
{
		switch (fld)
		{
			case 0:
				return (laplacian(field[s][fld],i,j,k)/a[s]/a[s] - dVdf(s,fld,i,j,k) - 3.*adot[s]/a[s]*dfield[s][fld][i][j][k]
							/*-1/sqrt(6)*exp( -sqrt23*field[s][0][i][j][k] )*(	pow(dfield[s][1][i][j][k],2)		+ pow(dfield[s][2][i][j][k],2) +  pow(dfield[s][3][i][j][k],2) + + pow(dfield[s][4][i][j][k],2) 
																				- pow(dfdi(field[s][1],i,j,k),2)/pow(a[s],2)	- pow(dfdi(field[s][2],i,j,k),2)/pow(a[s],2) -pow(dfdi(field[s][3],i,j,k),2)/pow(a[s],2) -pow(dfdi(field[s][4],i,j,k),2)/pow(a[s],2) 
																				- pow(dfdj(field[s][1],i,j,k),2)/pow(a[s],2)	- pow(dfdj(field[s][2],i,j,k),2)/pow(a[s],2) -pow(dfdj(field[s][3],i,j,k),2)/pow(a[s],2) -pow(dfdj(field[s][4],i,j,k),2)/pow(a[s],2)
																				- pow(dfdk(field[s][1],i,j,k),2)/pow(a[s],2)	- pow(dfdk(field[s][2],i,j,k),2)/pow(a[s],2) -pow(dfdk(field[s][3],i,j,k),2)/pow(a[s],2) -pow(dfdk(field[s][4],i,j,k),2)/pow(a[s],2)
																				)*/
																				
						); 
			case 1:
				return ( laplacian(field[s][fld],i,j,k)/a[s]/a[s] - /*exp(sqrt23*field[s][0][i][j][k])* */dVdf(s,fld,i,j,k) - 3.*adot[s]/a[s]*dfield[s][fld][i][j][k]
						/*+ sqrt23*(dfield[s][0][i][j][k]*dfield[s][fld][i][j][k]	
									-dfdi(field[s][0],i,j,k)*dfdi(field[s][fld],i,j,k)/pow(a[s],2)
									-dfdj(field[s][0],i,j,k)*dfdj(field[s][fld],i,j,k)/pow(a[s],2)
									-dfdk(field[s][0],i,j,k)*dfdk(field[s][fld],i,j,k)/pow(a[s],2) 
									)*/
									
						);			
			case 2:
			return ( laplacian(field[s][fld],i,j,k)/a[s]/a[s] - /*exp(sqrt23*field[s][0][i][j][k])* */dVdf(s,fld,i,j,k) - 3.*adot[s]/a[s]*dfield[s][fld][i][j][k]
						/*+ sqrt23*(dfield[s][0][i][j][k]*dfield[s][fld][i][j][k]	
									-dfdi(field[s][0],i,j,k)*dfdi(field[s][fld],i,j,k)/pow(a[s],2)
									-dfdj(field[s][0],i,j,k)*dfdj(field[s][fld],i,j,k)/pow(a[s],2)
									-dfdk(field[s][0],i,j,k)*dfdk(field[s][fld],i,j,k)/pow(a[s],2) 
									)*/
									
						);		
			case 3:
				return ( laplacian(field[s][fld],i,j,k)/a[s]/a[s] - /*exp(sqrt23*field[s][0][i][j][k])* */dVdf(s,fld,i,j,k) - 3.*adot[s]/a[s]*dfield[s][fld][i][j][k]
						/*+ sqrt23*(dfield[s][0][i][j][k]*dfield[s][fld][i][j][k]	
									-dfdi(field[s][0],i,j,k)*dfdi(field[s][fld],i,j,k)/pow(a[s],2)
									-dfdj(field[s][0],i,j,k)*dfdj(field[s][fld],i,j,k)/pow(a[s],2)
									-dfdk(field[s][0],i,j,k)*dfdk(field[s][fld],i,j,k)/pow(a[s],2) 
									)*/
									
						);		
			case 4:
				return ( laplacian(field[s][fld],i,j,k)/a[s]/a[s] - /*exp(sqrt23*field[s][0][i][j][k])* */dVdf(s,fld,i,j,k) - 3.*adot[s]/a[s]*dfield[s][fld][i][j][k]
						/*+ sqrt23*(dfield[s][0][i][j][k]*dfield[s][fld][i][j][k]	
									-dfdi(field[s][0],i,j,k)*dfdi(field[s][fld],i,j,k)/pow(a[s],2)
									-dfdj(field[s][0],i,j,k)*dfdj(field[s][fld],i,j,k)/pow(a[s],2)
									-dfdk(field[s][0],i,j,k)*dfdk(field[s][fld],i,j,k)/pow(a[s],2) 
									)*/
									
						);		
		}
}


/** RK2 function **/

void step()//this steps (integrates) the field and it time derivative via the rk2 meathod.
{
    
    DECLARE_INDEX

	
#if expansion_type==0
    //no expansion note that this only caclulates the energies at the end of the full step

    for(fld=0;fld<nflds;fld++)//the first part of the RK2 step 
    {
    //paralleleizes over the index i
#pragma omp parallel for private (j,k) num_threads (tot_num_thrds) 
        LOOP//loops over fld i,j,k
        {
            field[1][fld][i][j][k]=field[0][fld][i][j][k]+.5*dt*dfield[0][fld][i][j][k];
            dfield[1][fld][i][j][k]=dfield[0][fld][i][j][k]+.5*dt*ddfield(0,fld,i,j,k);
        }
    }
	
	for(fld=0;fld<nflds;fld++)//the second part of the RK2 step 
    {
    //paralleleizes over the index i
#pragma omp parallel for private (j,k) num_threads (tot_num_thrds)
        LOOP//This returns the actuall value of the field and derivative at t
        {
            field[0][fld][i][j][k]=field[0][fld][i][j][k]+dt*dfield[1][fld][i][j][k];
            dfield[0][fld][i][j][k]=dfield[0][fld][i][j][k]+dt*ddfield(1,fld,i,j,k);
        }
    }	
    calcEnergy(0);//calcualtes the energy at the end of the step.
	
	
	
#elif expansion_type==1 //already know edkin[0] etc from energy calc at end
	
    for(fld=0;fld<nflds;fld++)//first step of the Rk2 integration
    {
#pragma omp parallel for private (j,k) num_threads (tot_num_thrds)
        LOOP//loops over fld i,j,k
        {
            field[1][fld][i][j][k]=field[0][fld][i][j][k]+.5*dt*dfield[0][fld][i][j][k];
            dfield[1][fld][i][j][k]=dfield[0][fld][i][j][k]+.5*dt*ddfield(0,fld,i,j,k);
        }
    }
	
    a[1]=a[0]+.5*dt*adot[0];//this does the first step of the RK2 for the scale factor
	adot[1]=adot[0]+.5*dt*calcdda(0);		//evolving scale factor via raych. equation
    calcEnergy(1);//this calculates the energy based on this half step
	//adot[1]=adf(1);		//friedmann directly
    
    	for(fld=0;fld<nflds;fld++)//second step of the Rk2 integration
    {
#pragma omp parallel for private (j,k) num_threads (tot_num_thrds)
        LOOP//This returns the actuall value of the field and derivative at t
        {
            field[0][fld][i][j][k]=field[0][fld][i][j][k]+dt*dfield[1][fld][i][j][k];
            dfield[0][fld][i][j][k]=dfield[0][fld][i][j][k]+dt*ddfield(1,fld,i,j,k);
        }
    }
			
    
	a[0]=a[0]+dt*adot[1];//this calclates the full step scale factor
	adot[0]=adot[0]+dt*calcdda(1);
    calcEnergy(0);//calculates the energy at the full step
    //adot[0]=adf(0);//then calculates adot based off of the full step
    
    
    
    
 /*   #elif expansion_type==2
	
    for(fld=0;fld<nflds;fld++)//first step of the Rk2 integration
    {
#pragma omp parallel for private (j,k) num_threads (tot_num_thrds)
        LOOP//loops over fld i,j,k
        {
            field[1][fld][i][j][k]=field[0][fld][i][j][k]+.5*dt*dfield[0][fld][i][j][k];
            dfield[1][fld][i][j][k]=dfield[0][fld][i][j][k]+.5*dt*ddfield(0,fld,i,j,k);
        }
    }
	a[1]=a[0]=+.5*dt*adot[0];
	adot[1]=adot[0]+.5*dt*calcdda(0);
	
    /* this may need to chage based off of user defined expansion*/
   // a[1]=a[0]+.5*dt*adot[0];//this does the first step of the RK2 for the scale factor
 /*   calcEnergy(1);//this calculates the energy based on this half step
    //adot[1]=adf(1);//this updates adot based off of the energy at this step
    
    	for(fld=0;fld<nflds;fld++)//second step of the Rk2 integration
    {
#pragma omp parallel for private (j,k) num_threads (tot_num_thrds)
        LOOP//This returns the actuall value of the field and derivative at t
        {
            field[0][fld][i][j][k]=field[0][fld][i][j][k]+dt*dfield[1][fld][i][j][k]; //just RK when converting N 2nd-order into 2N 1st-order
            dfield[0][fld][i][j][k]=dfield[0][fld][i][j][k]+dt*ddfield(1,fld,i,j,k);
        }
    }	
    a[0]=a[0]=+dt*adot[1];
	adot[0]=adot[0]+dt*calcdda(1);
	//a[0]=a[0]+dt*adot[1];//this calclates the full step scale factor
    calcEnergy(0);//calculates the energy at the full step
   // adot[0]=adf(0);//then calculates adot based off of the full step*/
#endif
    
}

 
 //void mass_matrix_calc(int s, int fld, long double& mdiag, long double ev0, long double ev1, long double ev2	) //currently accepts small (h_1, h_2): good! Want to output to Sfield
  void mass_matrix_calc(int s, long double(* evec)[nflds][nflds], long double(* mdiag)[nflds], bool ifdiag) //currently accepts small (h_1, h_2): good! Want to output to Sfield
 {
	int i,j,k;
	long double V00=0;
	long double V11=0;
	long double V22=0;
	long double V33=0;
	long double V44=0;

	long double  V10=0;
	long double  V20=0;
	long double  V30=0;
	long double  V40=0;

	long double  V12=0;
	long double V13=0;
	long double V14=0;
	long double V23=0;
	long double V24=0;
	long double V34=0;

	long double pw2h; // |h|^2
	long double er6f; //exp(sqrt(6)*phi)
	long double er23f; // \exp(sqrt(2/3) * phi)

	Matrix5L m2_matrix; //object of matrix class
//	Matrix3L m2_matrix_3d; //

		for (int ii=0; ii<nflds; ii++){
			for (int jj=0; jj<nflds; jj++){
					m2_matrix(ii,jj)=0;
			}
		}

	LOOP
    {		
		pw2h = pow(H1[i][j][k],2) + pow(H2[i][j][k],2) +  pow(H3[i][j][k],2) +  pow(H4[i][j][k],2);
		er6f = exp(sqrt(6)*PHI[i][j][k]);
		er23f = exp(sqrt23*PHI[i][j][k]);

			V00 += -((-2 + er23f - pw2h*xi1)*(b*l + pow(xi1,2)))/(b*pow(er23f,2)*l);
				

			 V11 += (3*(b*l + pow(xi1,2))*(b*l*(2*pow(H1[i][j][k],2) + pw2h) + 
						 xi1*(1 - er23f + 2*pow(H1[i][j][k],2)*xi1 + pw2h*xi1)))/(b*er23f*l);


			V22 += (3*(b*l + pow(xi1,2))*(b*l*(2*pow(H2[i][j][k],2) + pw2h) + 
						 xi1*(1 - er23f + 2*pow(H2[i][j][k],2)*xi1 + pw2h*xi1)))/(b*er23f*l);

			V33 += (3*(b*l + pow(xi1,2))*(b*l*(2*pow(H3[i][j][k],2) + pw2h) + 
						 xi1*(1 - er23f + 2*pow(H3[i][j][k],2)*xi1 + pw2h*xi1)))/(b*er23f*l);

			V44 +=  (3*(b*l + pow(xi1,2))*(b*l*(2*pow(H4[i][j][k],2) + pw2h) + 
						 xi1*(1 - er23f + 2*pow(H4[i][j][k],2)*xi1 + pw2h*xi1)))/(b*er23f*l);


			  V10 +=-(sqrt(6)*H1[i][j][k]*xi1*(b*l + pow(xi1,2))/er23f  )/(b*l);			  
			  V20 +=-(sqrt(6)*H2[i][j][k]*xi1*(b*l + pow(xi1,2))/er23f)/(b*l);
			  V30 +=-(sqrt(6)*H3[i][j][k]*xi1*(b*l + pow(xi1,2))/er23f)/(b*l);
			  V40 +=-(sqrt(6)*H4[i][j][k]*xi1*(b*l + pow(xi1,2))/er23f)/(b*l);

			  V12 += (6*H1[i][j][k]*H2[i][j][k]*pow(b*l + pow(xi1,2),2))/(b*er23f*l);
			  V13 += (6*H1[i][j][k]*H3[i][j][k]*pow(b*l + pow(xi1,2),2))/(b*er23f*l);
			  V14 += (6*H1[i][j][k]*H4[i][j][k]*pow(b*l + pow(xi1,2),2))/(b*er23f*l);
			  V23 += (6*H2[i][j][k]*H3[i][j][k]*pow(b*l + pow(xi1,2),2))/(b*er23f*l);
			  V24 += (6*H2[i][j][k]*H4[i][j][k]*pow(b*l + pow(xi1,2),2))/(b*er23f*l);
			  V34 += (6*H3[i][j][k]*H4[i][j][k]*pow(b*l + pow(xi1,2),2))/(b*er23f*l);
			  
    }
    V00 /=  gridsize;
	V11 /=  gridsize;
	V22 /=  gridsize;
	V33 /=  gridsize;
	V44 /=  gridsize;

	V10 /=  gridsize;
	V20 /=  gridsize;
	V30 /=  gridsize;
	V40 /=  gridsize;

	V12 /=  gridsize;
	V13 /=  gridsize;
	V14 /=  gridsize;
	V23 /=  gridsize;
	V24 /=  gridsize;
	V34 /=  gridsize;

	m2_matrix(0,0) = V00;
	m2_matrix(1,1) = V11;
	m2_matrix(2,2) = V22;
	m2_matrix(3,3) = V33;
	m2_matrix(4,4) = V44;

	m2_matrix(0,1) = V10;
	m2_matrix(1,0) = V10;
	m2_matrix(0,2) = V20;
	m2_matrix(2,0) = V20;
	m2_matrix(0,3) = V30;
	m2_matrix(3,0) = V30;
	m2_matrix(0,4) = V40;
	m2_matrix(4,0) = V40;

	m2_matrix(1,2) = V12;
	m2_matrix(2,1) = V12;
	m2_matrix(1,3) = V13;
	m2_matrix(3,1) = V13;
	m2_matrix(1,4) = V14;
	m2_matrix(4,1) = V14;
	m2_matrix(2,3) = V23;
	m2_matrix(3,2) = V23;
	m2_matrix(2,4) = V24;
	m2_matrix(4,2) = V24;
	m2_matrix(3,4) = V34;
	m2_matrix(4,3) = V34;

	
	//std::cout<<"m2_matrix is: \n"<<m2_matrix<<"\n";

	if (ifdiag==true) //finds eigenvalues and eigenvectors. not proper basis after t=0
	{
		Eigen::EigenSolver<Matrix5L> es(m2_matrix);

		for (int fld=0;fld<nflds;fld++)
			{	
				(*mdiag)[fld] = es.eigenvalues()[fld].real(); 
				}	
	

	
		for (int fld=0;fld<nflds;fld++)
				{for (int rr=0;rr<nflds;rr++){
					
						(*evec)[rr][fld]=(long double) es.eigenvectors().real().col(fld).array()[rr];
					//	printf("evec[%i][%i] is: %Le \n", rr, fld, (*evec)[rr][fld]);
					}
				
			}
		
		
	}

	else //just use the diagonal components of the mass matrix 
	{
		for (int fld=0;fld<nflds;fld++)
			{	
				(*mdiag)[fld] = m2_matrix(fld,fld); 
				}		
	
		for (int fld=1;fld<nflds;fld++)
				{for (int rr=0;rr<nflds;rr++){
					if (rr == fld){
						(*evec)[rr][fld]= 1;	
					}
					else{
						(*evec)[rr][fld]=0;
					}
			}
		}
	}
 }
	
void fromDiag ( long double (*Sfield)[nflds][N][N][N],long double (*Sdfield)[nflds][N][N][N] , long double(*evec)[nflds][nflds]  ){ //TODO are the dimensions of these right
		DECLARE_INDEX
			Matrix5L uT; // each column: one eigenvector
			Vector5L fVec, dfVec;

			for (int cc=0;cc<nflds;cc++){ //loop over columns
				for(int rr=0;rr<nflds;rr++){
					uT(rr,cc) =	(*evec)[rr][cc]	;	//fill the matrix
				}
			}
			
			LOOP{	//for all points
					for (fld=0;fld<nflds;fld++){
						fVec(fld) = (*Sfield)[fld][i][j][k]; //vector of field values
						dfVec(fld) = (*Sdfield)[fld][i][j][k];
					} //field loop
				//	std::cout<<"fVec = "<<fVec<<"\n";
					if (i==1 & j==1 & k==1){
				//		std::cout<<"fVec before fromDiag is: \n"<<fVec<<"\n";
					}
				fVec = uT*fVec; //matrix multiplication 
				dfVec = uT*dfVec;
				
						
				(*Sfield)[0][i][j][k]= fVec(0); //overwrite Sfield[i][j][k] with new values
				(*Sdfield)[0][i][j][k]=dfVec(0);

				for(int ii=1;ii< (nflds) ;ii++)
				{
					(*Sfield)[ii][i][j][k]=fVec(ii) * exp(fVec(0) / sqrt(6)	);
					(*Sdfield)[ii][i][j][k]=dfVec(ii) * exp(fVec(0) / sqrt(6))	+ fVec(ii)*dfVec(0)/sqrt(6) * exp(fVec(0) / sqrt(6)	);
				}
			}//points loop
}

void toDiag ( long double (*Sfield)[nflds][N][N][N],long double(* Sdfield)[nflds][N][N][N] , long double(*evec)[nflds][nflds], bool ifdiag){ 
			DECLARE_INDEX
			Matrix5L uT, uTinv; // each column: one eigenvector
			Vector5L fVec, dfVec;

		if (ifdiag==true) //rotate into something like eigenstates. doesn't catch the angular mode though
		{
			for (int cc=0;cc<nflds;cc++){ //loop over columns
				for(int rr=0;rr<nflds;rr++){
					uT(rr,cc) =	(*evec)[rr][cc]	;	//fill the matrix
				}
			}
		}
			LOOP{	//for all points
					fVec(0) = field[0][0][i][j][k]; //vector of field values
					dfVec(0) = dfield[0][0][i][j][k]; //vector of field values
					
					for ( int ii=1; ii<nflds ;ii++)
					{
						fVec(ii) = field[0][ii][i][j][k] / exp(field[0][0][i][j][k] / sqrt(6)	);
						dfVec(ii) = dfield[0][ii][i][j][k] / exp(field[0][0][i][j][k] / sqrt(6))	- (1/sqrt(6) )*dfield[0][0][i][j][k]*field[0][ii][i][j][k] / exp(field[0][0][i][j][k] / sqrt(6))	;
					}
			
		if (ifdiag==true)
		{
				uTinv = uT.inverse().eval();	//invert the matrix
				//std::cout<<"inverse matrix is: \n"<<uTinv<<"\n";

				fVec = uTinv*fVec; //matrix multiplication
				dfVec = uTinv*dfVec; //disabled diagonalization
		}

				for (fld=0;fld<nflds;fld++){ //why do I have to close the field loop and do again?
				(*Sfield)[fld][i][j][k] = fVec(fld); //save diagonalized canonical fields
				(*Sdfield)[fld][i][j][k] = dfVec(fld);
				}
			}//positions loop		
}


void Viicalc ( int s, long double *V11, long double *V22, long double *V12, long double *Vrr){ //function to spit out V11, V12, V22
	DECLARE_INDEX
	long double V11temp=0;
	long double V22temp=0;
	long double V12temp=0;
	long double Vrrtemp=0;

	long double pw2h;
	long double er6f;
	long double er23f;

	LOOP{
		pw2h = pow(H1[i][j][k],2) + pow(H2[i][j][k],2) +  pow(H3[i][j][k],2) +  pow(H4[i][j][k],2);
		er6f = exp(sqrt(6)*PHI[i][j][k]);
		er23f = exp(sqrt23*PHI[i][j][k]);


		  V11temp +=(3*(b*l + pow(xi1,2))*(b*l*(2*pow(H1[i][j][k],2) + pw2h) + 
					   xi1*(1 - er23f + 2*pow(H1[i][j][k],2)*pw2h*xi1*xi1)))/(b*pow(er23f,2)*l);
				  

		   V22temp +=(3*(b*l + pow(xi1,2))*(b*l*(2*pow(H2[i][j][k],2) + pw2h) + 
					xi1*(1 - er23f + 2*pow(H2[i][j][k],2)*pw2h*xi1*xi1)))/(b*pow(er23f,2)*l);

		   V12temp +=(6*H1[i][j][k]*H2[i][j][k]*pow(b*l + pow(xi1,2),2))/(b*er23f*l);

		   Vrrtemp +=(3*(b*l + pow(xi1,2))*(3*b*l*pw2h + xi1*(1 - er23f + 3*pw2h*xi1)))/(b*pow(er23f,2)*l);
	}
	V11temp/=gridsize;
	V22temp/=gridsize;
	V12temp/=gridsize;
	Vrrtemp/=gridsize;

	*V11 = V11temp;
	*V22 = V22temp;
	*V12 = V12temp;
	*Vrr = Vrrtemp;
}

long double Vrrcalc ( int s){ //function to spit out V11, V12, V22
	DECLARE_INDEX
	
	long double Vrrtemp=0;

	long double pw2h;
	long double er6f;
	long double er23f;

	LOOP{
		pw2h = pow(H1[i][j][k],2) + pow(H2[i][j][k],2) +  pow(H3[i][j][k],2) +  pow(H4[i][j][k],2);
		er6f = exp(sqrt(6)*PHI[i][j][k]);
		er23f = exp(sqrt23*PHI[i][j][k]);

		   Vrrtemp +=(3*(b*l + pow(xi1,2))*(3*b*l*pw2h + xi1*(1 - er23f + 3*pw2h*xi1)))/(b*pow(er23f,2)*l);
	}
	Vrrtemp/=gridsize;

	return Vrrtemp;
}