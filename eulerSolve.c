//David Ruffner
//Computational Physics
//11-10-09
//1D Euler solver
//---------------------------------------------------------------------------
//Used some template functions from Jonathan Zrake's advection code

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h> 
//#include "mpi.h"
#include "euler.h"
#include "config.h"
#include "myconfig.h"

#define STRING_LENGTH 128
#define gamma 1.5
#define g .1


// Emulate numpy.linspace behavior
// -----------------------------------------------------------------------
void linspace2nd(double *vals, double a, double b, int N, double *dx)
{
  *dx = (b-a) / (N);//Important, I found that if I used the factor N instead
  //of N-1, then I wouldn't get an over count as my analytical solution went
  // around the boundary
  int i;
  for (i=0; i<(N+4); ++i)//For the two ghost cell on each side
    {
      vals[i] = a + (i-2) * (*dx);
    }
}


// Apply a function to a whole array
// -----------------------------------------------------------------------
void apply_over_array(double (*f)(double), double *y, double *x, int N)
{
  int i;
  for (i=0; i<N; ++i)
    {
      y[i] = f(x[i]);
    }
}

// Apply a function to a whole array of Conserved Structs
// -----------------------------------------------------------------------
void apply_over_arrayCons2(struct Conserved (*f)(double x),
			   struct Conserved *y, double *x, int N)
{
  int i;
  for (i=0; i<N; ++i)
    {
      y[i].rho = f(x[i]).rho;
      y[i].E = f(x[i]).E;
      y[i].px = f(x[i]).px;
    }
}


// Initial condition functions
// -----------------------------------------------------------------------
//Shock Tube with interface at x=0
struct Conserved shock(double x, struct Conserved Ur, struct Conserved Ul)
{
  if ( x < 0)
    {
      return Ul;
    }
  else
    {
      return Ur;
    }
}

//sets values of U array for a shock tube
//---------------------------------------------------------------------------
void shockInitial(struct Conserved *y, struct Conserved Ur, 
		  struct Conserved Ul, double *x, int N)
{
  int i;
  for (i=0; i<N; ++i)
    {
      y[i] = shock(x[i],Ur, Ul);
      //printf("%f \n", y[i].rho);
    }
}


// Gaussian
double initial1(double x)
{
  return exp(-x*x / 0.1);
}
// Square
double initial2(double x)
{
  if(x<= 0.1 && x>= -0.1)
    {
      return 1;
    }
  else
    {
      return 0;
    }
}









//Sign function
//--------------------
double sgn(double x)
{
  if(x==0)
    {
      return 1.0;
    }
  else
    {
      return x/fabs(x);
    }
}

//Min Mod Function
//-----------------------------------------------------------------
double minmod(double x, double y, double z)
{
  return .25*fabs( sgn(x)+sgn(y) )*(sgn(x)+sgn(y))*Min2(fabs(x),Min2(fabs(y),fabs(z) ) );
}

//PLM interpolation
//------------------------------------------------------------------

//Reconstruct (Assumes two ghost points at the end of the array)
//-------------------------------------------------------------
void reconstruct(struct Conserved *U, struct Conserved *dUdx, int N)
{
  double theta = 1.0;
  int i;
  for(i=0; i< N+2 ; ++i)
    {
      dUdx[i].rho = minmod(theta*(U[i+1].rho-U[i].rho),
			   .5*(U[i+2].rho-U[i].rho),
			   theta*(U[i+2].rho-U[i+1].rho) );
      dUdx[i].px = minmod(theta*(U[i+1].px-U[i].px),
			   .5*(U[i+2].px-U[i].px),
			  theta*(U[i+2].px-U[i+1].px) );
      dUdx[i].E = minmod(theta*(U[i+1].E-U[i].E),
			   .5*(U[i+2].E-U[i].E),
			 theta*(U[i+2].E-U[i+1].E) );
    }
}

//Flux at Fiph in secondOrder
//----------------------------------------------------------------------
void Flux_at_iph2nd(struct Conserved (*flux)(struct Conserved, 
	  struct Conserved), 
	  struct Conserved *U, struct Conserved *Fiph, double N)
{
  struct Conserved *dUdx = (struct Conserved*) calloc(N+2, 
						     sizeof(struct Conserved));
  reconstruct(U,dUdx, N);
  
  int i;
  struct Conserved Ul;
  struct Conserved Ur;
  for(i=0; i<N+1; ++i)
    {
      Ul.rho = U[i+1].rho + .5*dUdx[i].rho;
      Ul.px = U[i+1].px + .5*dUdx[i].px;
      Ul.E = U[i+1].E + .5*dUdx[i].E;

      Ur.rho = U[i+2].rho - .5*dUdx[i+1].rho;
      Ur.px = U[i+2].px - .5*dUdx[i+1].px;
      Ur.E = U[i+2].E - .5*dUdx[i+1].E;

      Fiph[i] = flux(Ur, Ul);
    }
}
      

//Numerical Approx of dx*dU/dt  2nd Order
//----------------------------------------------------------------
void L2nd(struct Conserved *U, struct Conserved *LU, int N, double dx)
{
  struct Conserved *Fiph = (struct Conserved*) calloc(N+1, 
					      sizeof(struct Conserved));
  //N+1 Fiph cells because flux is needed on both sides of each U
  
  Flux_at_iph2nd(F_HLL, U, Fiph, N);

  int j;
  for(j=1; j < N+1; ++j)
    {
      LU[j-1].rho = - Fiph[j].rho + Fiph[j-1].rho;
      LU[j-1].px = - Fiph[j].px + Fiph[j-1].px - dx*U[j-1].rho*g;
      LU[j-1].E = - Fiph[j].E + Fiph[j-1].E - dx*U[j-1].px*g;
      //printf("LU.px = %f \n", LU[j-1].px);
    } 
  //printf("Fiph[0].px = %f, Fiph[1].px = %f \n", Fiph[0].px, Fiph[1].px);   
  free(Fiph);
  
}



//Boundary Conditions
//----------------------------------------------------------------------
//First Order:assuming U is N+2 long

//Outflow
void BC_1stO_out(struct Conserved *U, int N)
{
  U[0].rho = U[1].rho;
  U[0].px  = U[1].px;
  U[0].E   = U[1].E;

  U[N+1].rho = U[N].rho;
  U[N+1].px  = U[N].px;
  U[N+1].E   = U[N].E;
}
//Reflecting
void BC_1stO_rflt(struct Conserved *U, int N)
{
  U[0].rho = U[1].rho;
  U[0].px  = -U[1].px;
  U[0].E   = U[1].E;

  U[N+1].rho = U[N].rho;
  U[N+1].px  = -U[N].px;
  U[N+1].E   = U[N].E;
}

//Reflecting left out right
void BC_1stO_Rout_Lrflt(struct Conserved *U, int N)
{
  U[0].rho = U[1].rho;
  U[0].px  = -U[1].px;
  U[0].E   = U[1].E;

  U[N+1].rho = U[N].rho;
  U[N+1].px  = U[N].px;
  U[N+1].E   = U[N].E;
}

//Second Order Outflow
//-----------------------------------------
void BC_2ndO_out(struct Conserved *U, int N)
{
  U[0] = U[2];
  U[1] = U[2];

  U[N+2] = U[N+1];
  U[N+3] = U[N+1];
  
}
//Second Order Reflecting
//-----------------------------------------
void BC_2ndO_rfl(struct Conserved *U, int N)
{
  U[0] = U[2];
  U[1] = U[2];
  U[0].px = -U[0].px;
  U[1].px = -U[1].px;  

  U[N+2] = U[N+1];
  U[N+3] = U[N+1];
  U[N+2].px = -U[N+2].px;
  U[N+3].px = -U[N+3].px;  
  
}


//Applys the appropriate BCs based on BCnum
//------------------------------------------
void BC(int BCnum, struct Conserved *U, int N)
{
  if(BCnum == 1)
    {
      BC_2ndO_out(U,N);
    }
  else if(BCnum == 2)
    {
      BC_2ndO_rfl(U,N);
    }
   
}
  



//Time step of Euler Equations 2nd order, using F_HLL
//----------------------------------------------------------------------
double evolveF_HLL_2ndO(struct Conserved *U, int N, double dx, double CFL,
			double BCnum)
{
  struct Conserved *LU = (struct Conserved*) calloc(N, 
			    sizeof(struct Conserved));
  struct Conserved *Utemp = (struct Conserved*) calloc(N+4, 
			    sizeof(struct Conserved));
  double lamdaMax;
  
  double dt;
  BC(BCnum, U, N); 
  lamdaMax = MaxEigenVal(U, N);
  

  

  dt = CFL*dx/lamdaMax;//Satisfies the Courant Condition
  
  L2nd(U, LU, N, dx);//This creates the dx*dU/dx for all Ui
 
  //third order Runge-Kutta Time Step
  //------------------------------------------------------------------------
  
  int j;
  for(j=2; j < N+2; ++j)
    {
      //with gravity
      Utemp[j].rho = U[j].rho + (dt/dx)*LU[j-2].rho;
      Utemp[j].px = U[j].px + (dt/dx)*LU[j-2].px ;
      Utemp[j].E = U[j].E + (dt/dx)*LU[j-2].E ;
    }




  BC(BCnum, Utemp, N);

  L2nd(Utemp, LU, N, dx);
  
  
  for(j=2; j < N+2; ++j)
    {
      Utemp[j].rho = .75*U[j].rho + 0.25*Utemp[j].rho+ 0.25*(dt/dx)*LU[j-2].rho;
      Utemp[j].px  = .75*U[j].px  + 0.25*Utemp[j].px + 0.25*(dt/dx)*LU[j-2].px;
      Utemp[j].E   = .75*U[j].E   + 0.25*Utemp[j].E  + 0.25*(dt/dx)*LU[j-2].E;
    }

  BC(BCnum, Utemp, N);
  L2nd(Utemp, LU, N, dx);
  
  
  for(j=2; j < N+2; ++j)
    {
      U[j].rho = (1/3.0)*U[j].rho + (2/3.0)*Utemp[j].rho+ (2/3.0)*(dt/dx)*LU[j-2].rho;
      U[j].px  = (1/3.0)*U[j].px  + (2/3.0)*Utemp[j].px + (2/3.0)*(dt/dx)*LU[j-2].px;
      U[j].E   = (1/3.0)*U[j].E   + (2/3.0)*Utemp[j].E  + (2/3.0)*(dt/dx)*LU[j-2].E;
    }
  
  //--------------------------------------------------------------------------
  
  free(Utemp);
  free(LU);
  
  return dt;
}



int main(int argc, char **argv)
{
  
  
  
  //Get Parameters
  //----------------------------------------------------
  int    Nx    =  get_int_param("param.cfg", "Nx", 10);
  int    BCnum2 = get_int_param("param.cfg", "BCnum2", 1);
  int    InitialType = get_int_param("param.cfg", "InitialType", 1);
  double runTime   = get_double_param("param.cfg", "runTime", 0);
  double time   = get_double_param("time.cfg", "t", 0);
  double x0    = -3.0;
  double x1    = 5.0;

  //put in initial values of Ur Ul
  struct Conserved Ur;
  struct Conserved Ul;
  double Mrpre = get_double_param("param.cfg", "Mr.pre", 0);
  double Mrrho = get_double_param("param.cfg", "Mr.rho", 0);
  double Mrv   = get_double_param("param.cfg", "Mr.v", 0);
  
  double Mlpre = get_double_param("param.cfg", "Ml.pre", 0);
  double Mlrho = get_double_param("param.cfg", "Ml.rho", 0);
  double Mlv   = get_double_param("param.cfg", "Ml.v", 0);
  Ur = mixed_to_cons(Mrpre, Mrrho, Mrv );
  Ul = mixed_to_cons(Mlpre, Mlrho, Mlv );
  //------------------------------------------

  
  
  
  //Gets the dx
   
  double dx2nd;
  double *x2nd = (double*) calloc(Nx+4, sizeof(double));
  linspace2nd(x2nd, x0, x1, Nx, &dx2nd); 
  


  //Setting the Initial Conditions
  //-----------------------------------------------------------------------
  struct Conserved *Uinitial2nd = (struct Conserved*) 
    calloc(Nx+4, sizeof(struct Conserved));//Size N+4 for the four ghost cells
  struct Primitive *Pinitial2nd = (struct Primitive*) 
    calloc(Nx+4, sizeof(struct Primitive));//Size N+4 for the four ghost cells
  
  
  
  if(InitialType == 1)
    {
      
      shockInitial(Uinitial2nd,Ur,Ul,x2nd,Nx+4);
           
     
      //Writing initial Conditions to a file.
      
      write_ascii_tableCons("run/initialIn2nd.cfg",x2nd, 
			    Uinitial2nd ,Nx, "w");

    }
  else if(InitialType == 2)
    {
      
      read_ascii_tableCons("run/testFinalCons2nd.cfg",x2nd, Uinitial2nd,Nx);
    }
  else if(InitialType == 3)
    {
      
      read_ascii_tablePrim("run/initialIn2nd.cfg",x2nd, Pinitial2nd,Nx);
      int i;
      
      for(i=0;i<Nx+4;++i)
	{
	  Uinitial2nd[i] = prim_to_cons(Pinitial2nd[i]);
	}
      
    }
  else if(InitialType == 4)
    {
      
      read_ascii_tableCons("run/initialIn2nd.cfg",x2nd, Uinitial2nd,Nx);
      
    }
  //-----------------------------------------------------------------------


  //Time Evolution
  //------------------------------------------------------------------------------
  double CFL = 0.5;
  
  
  double timeLocal = 0;  
  int j=0;
     
  
  
  while(timeLocal<= runTime)
    {
      j++;
      timeLocal = timeLocal + evolveF_HLL_2ndO(Uinitial2nd, Nx, dx2nd, CFL, 
				     BCnum2);
    }
  time = time + timeLocal;
  printf("For Second Order the time after the %d th  step is: %f \n" ,j, time);
  
  //output time
  FILE *outfile = fopen("time.cfg", "w");
  
  fprintf(outfile, "t = %f\n", time);
  
  fclose(outfile);


  
  int i;  
  for(i=0;i<Nx+4;++i)
     {
       Pinitial2nd[i] = cons_to_prim(Uinitial2nd[i]);
     }
  //----------------------------------------------------------------------------




  //output results
  //-------------------------------------------------------------
  
  
  write_ascii_tableCons("run/testFinalCons2nd.cfg",x2nd, Uinitial2nd,Nx, "w");
  write_ascii_tablePrim("run/testFinalPrim2nd.cfg",x2nd, Pinitial2nd,Nx, "w");
  

  //--------------------------------------------------------   


  free(x2nd);
  
  free(Pinitial2nd);
  free(Uinitial2nd);
  
  
  return 0;
}
