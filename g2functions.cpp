/**********************************
 FUNCTIONS FILE
 **********************************/

/*
 This header file contains all the functions which are independent of the model needed to evolve the fields by the Runge-Kutta Second order method (for first order finite derivatives). The incr and decr commands set periodic boundary conditions on the lattice durring evolution.
 
 Copyright (2013): 
 Kenyon College
 John T. Giblin, Jr
 Tate Deskins and Hillary Child
 Last Updated: 06.27.2013
*/

#include "g2header.h" //contains declerations for program functions.

long double pw2(long double x)//squares long doubles
{
    return x*x;
}

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
            grad-=field[s][fld][i][j][k]*laplacian(field[s][fld],i,j,k);//sums the gradient energy at each point
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
            kin+=dfield[s][fld][i][j][k]*dfield[s][fld][i][j][k];//sums the square field derivative at every point
        }
    }
    return kin/gridsize/2.;//divide by the grid size to get the average and 2
}


void calcEnergy(int s) //Calculate the total energy
{
    edkin[s]=avgKin(s);
    edpot[s]=avgPot(s);
    edgrad[s]=avgGrad(s);
    edrho[s]=edkin[s]+edpot[s]+edgrad[s];
}

long double adf(int s)//the friedman equation
{
    return sqrt(8.*M_PI*grav/3.*edrho[s])*a[s];
}

long double ddfield( int s, int fld, int i, int j, int k)//evaluates the double time derivative of the field fld (s) at i,j,k.
{
   return (laplacian(field[s][fld],i,j,k)/a[s]/a[s] - dVdf(s,fld,i,j,k) - 3.*adot[s]/a[s]*dfield[s][fld][i][j][k]);
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
	
	
	
#elif expansion_type==1
	
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
    calcEnergy(1);//this calculates the energy based on this half step
    adot[1]=adf(1);//this updates adot based off of the energy at this step
    
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
    calcEnergy(0);//calculates the energy at the full step
    adot[0]=adf(0);//then calculates adot based off of the full step
    
    
    
    
    #elif expansion_type==2
	
    for(fld=0;fld<nflds;fld++)//first step of the Rk2 integration
    {
#pragma omp parallel for private (j,k) num_threads (tot_num_thrds)
        LOOP//loops over fld i,j,k
        {
            field[1][fld][i][j][k]=field[0][fld][i][j][k]+.5*dt*dfield[0][fld][i][j][k];
            dfield[1][fld][i][j][k]=dfield[0][fld][i][j][k]+.5*dt*ddfield(0,fld,i,j,k);
        }
    }
	
    /* this may need to chage based off of user defined expansion*/
    a[1]=a[0]+.5*dt*adot[0];//this does the first step of the RK2 for the scale factor
    calcEnergy(1);//this calculates the energy based on this half step
    adot[1]=adf(1);//this updates adot based off of the energy at this step
    
    	for(fld=0;fld<nflds;fld++)//second step of the Rk2 integration
    {
#pragma omp parallel for private (j,k) num_threads (tot_num_thrds)
        LOOP//This returns the actuall value of the field and derivative at t
        {
            field[0][fld][i][j][k]=field[0][fld][i][j][k]+dt*dfield[1][fld][i][j][k];
            dfield[0][fld][i][j][k]=dfield[0][fld][i][j][k]+dt*ddfield(1,fld,i,j,k);
        }
    }	
    /* this may need to change based off of user defined expansion*/
	a[0]=a[0]+dt*adot[1];//this calclates the full step scale factor
    calcEnergy(0);//calculates the energy at the full step
    adot[0]=adf(0);//then calculates adot based off of the full step
#endif
    
}
