//David Ruffner
//Computational Physics
//11-10-09
//2D Euler solver
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

void print2DarrayPrimPre(struct Primitive **array, int nrows, int ncolumns);

double MaxEigenVal(struct Conserved **U, int Nx, int Ny)//should I go through
//ghost cells as well
{
 
  int i,j;
  struct Primitive **Uprim;
  Uprim    = malloc((Ny+4)* sizeof(struct Primitive *));
     
  if(Uprim == NULL)
    {
      fprintf(stderr, "out of memory\n");
      return 0;
    }
  for(i = 0; i < (Ny+4); i++)
    {
      Uprim[i]    = malloc((Nx+4) * sizeof(struct Primitive));
      
      if(Uprim[i] == NULL)
	{
	  fprintf(stderr, "out of memory\n");
	  return 0;
	}
    }

  //printf("in maxEigenVal \n");

  double lamda1X  = 0;
  double lamda2X  = 0;
  double lamda1Y  = 0;
  double lamda2Y  = 0;
  double lamdaMax = 0;
  double tempMax = 0;

  

  for(i=2; i<Ny+2; ++i)
    {
      for(j=2; j<Nx+2; j++)
	{
	  Uprim[i][j] = cons_to_prim(U[i][j]);
	  EigenvaluesX(Uprim[i][j], &lamda1X, &lamda2X);
	  EigenvaluesY(Uprim[i][j], &lamda1Y, &lamda2Y);
	  tempMax = Max2(Max2(fabs(lamda1X), fabs(lamda2X)),
			 Max2(fabs(lamda1Y), fabs(lamda2Y)));
          //maybe abs value
	  lamdaMax = Max2(tempMax, lamdaMax);
	}
    }
  //printf("mid maxEigenVal \n");
  //printf("Ny %d\n", Ny);
  //print2DarrayPrimPre(Uprim, Ny+4,Nx+4);
  for(i=0; i<Ny+4; i++)
    {
      
      free(Uprim[i]);
    }
  free(Uprim);
  
  //printf("end maxEigenVal \n");
  
  //double lamdaMax = 1.0;
  //printf("lamdaMax %f\n",lamdaMax);
  return lamdaMax;
}

//Output the time in simulation
//---------------------------------------------------
void write_ascii_time(const char *outname,
			   double time, char* perm)
{
  FILE *outfile = fopen(outname, perm);
  
  fprintf(outfile, "t = %f\n", time);
  fclose(outfile);
}

// Output for plotting
// -----------------------------------------------------------------------
void write_ascii_table(const char *outname,
			   double *vals, int N, char* perm)
{
  FILE *outfile = fopen(outname, perm);
  int i;
  for (i=2; i<N+2; ++i)//for 2nd euler
    {
      fprintf(outfile, "%f\n", vals[i]);
	
    }
  fclose(outfile);
}










// makes a 2D array of Rho values
// -----------------------------------------------------------------------
void write_ascii_tableRhoXY(const char *outname,
			   struct Conserved **vals, int Nx, int Ny, char* perm)
{
  FILE *outfile = fopen(outname, perm);
  int i,j;
  for (i=2; i<Ny+2; ++i)//for 2nd euler
    {
      for(j=2; j<Nx+2; ++j)
	{
	  fprintf(outfile, "%f ", vals[i][j].rho);
	}
      fprintf(outfile, "\n");
    }
  fclose(outfile);
}

// makes a 2D array of Rho values
// -----------------------------------------------------------------------
void write_ascii_tablepxXY(const char *outname,
			   struct Conserved **vals, int Nx, int Ny, char* perm)
{
  FILE *outfile = fopen(outname, perm);
  int i,j;
  for (i=2; i<Ny+2; ++i)//for 2nd euler
    {
      for(j=2; j<Nx+2; ++j)
	{
	  fprintf(outfile, "%f ", vals[i][j].px);
	}
      fprintf(outfile, "\n");
    }
  fclose(outfile);
}

// makes a 2D array of py values
// -----------------------------------------------------------------------
void write_ascii_tablepyXY(const char *outname,
			   struct Conserved **vals, int Nx, int Ny, char* perm)
{
  FILE *outfile = fopen(outname, perm);
  int i,j;
  for (i=2; i<Ny+2; ++i)//for 2nd euler
    {
      for(j=2; j<Nx+2; ++j)
	{
	  fprintf(outfile, "%f ", vals[i][j].py);
	}
      fprintf(outfile, "\n");
    }
  fclose(outfile);
}

// makes a 2D array of E values
// -----------------------------------------------------------------------
void write_ascii_tableEXY(const char *outname,
			   struct Conserved **vals, int Nx, int Ny, char* perm)
{
  FILE *outfile = fopen(outname, perm);
  int i,j;
  for (i=2; i<Ny+2; ++i)//for 2nd euler
    {
      for(j=2; j<Nx+2; ++j)
	{
	  fprintf(outfile, "%f ", vals[i][j].E);
	}
      fprintf(outfile, "\n");
    }
  fclose(outfile);
}

//Make 2Darrays from Mixed values
//--------------------------------------------------------
void write_ascii_tableMixRhoXY(const char *outname, 
			       struct Mixed **vals, int Nx,
			       int Ny, char* perm)
{
  FILE *outfile = fopen(outname, perm);
  int i,j;
  for (i=2; i<Ny+2; ++i)//for 2nd euler
    {
      for(j=2; j<Nx+2; ++j)
	{
	  fprintf(outfile, "%f ", vals[i][j].rho);
	}
      fprintf(outfile, "\n");
    }
  fclose(outfile);
}

void write_ascii_tableMixvxXY(const char *outname, 
			       struct Mixed **vals, int Nx,
			       int Ny, char* perm)
{
  FILE *outfile = fopen(outname, perm);
  int i,j;
  for (i=2; i<Ny+2; ++i)//for 2nd euler
    {
      for(j=2; j<Nx+2; ++j)
	{
	  fprintf(outfile, "%f ", vals[i][j].vx);
	}
      fprintf(outfile, "\n");
    }
  fclose(outfile);
}

void write_ascii_tableMixvyXY(const char *outname, 
			       struct Mixed **vals, int Nx,
			       int Ny, char* perm)
{
  FILE *outfile = fopen(outname, perm);
  int i,j;
  for (i=2; i<Ny+2; ++i)//for 2nd euler
    {
      for(j=2; j<Nx+2; ++j)
	{
	  fprintf(outfile, "%f ", vals[i][j].vy);
	}
      fprintf(outfile, "\n");
    }
  fclose(outfile);
}

void write_ascii_tableMixpreXY(const char *outname, 
			       struct Mixed **vals, int Nx,
			       int Ny, char* perm)
{
  FILE *outfile = fopen(outname, perm);
  int i,j;
  for (i=2; i<Ny+2; ++i)//for 2nd euler
    {
      for(j=2; j<Nx+2; ++j)
	{
	  fprintf(outfile, "%f ", vals[i][j].pre);
	}
      fprintf(outfile, "\n");
    }
  fclose(outfile);
}

void read_ascii_X(const char *inname, double *x, int N)
//This N is the actual N it reads
{
  
  FILE *infile = fopen(inname, "r");
  int i;

  
  for (i=0; i<N; ++i)
    {
      fscanf(infile, "%lf\n", &x[i]);
    }
  fclose(infile);
}

void read_ascii_tableConsRho(const char *inname, 
			  struct Conserved **vals, int Nx, int Ny)

{
  
  FILE *infile = fopen(inname, "r");
  printf("scanning\n");
  int i, j;

  
  for (i=0; i<Ny; ++i)
    {
      for(j= 0; j<Nx; ++j)
	{
	  fscanf(infile, "%lf", &vals[i+2][j+2].rho);//accounts for ghost cells
	}
    }
  fclose(infile);
}

void read_ascii_tableConspx(const char *inname, 
			  struct Conserved **vals, int Nx, int Ny)

{
  
  FILE *infile = fopen(inname, "r");
  printf("scanning\n");
  int i, j;

  
  for (i=0; i<Ny; ++i)
    {
      for(j= 0; j<Nx; ++j)
	{
	  fscanf(infile, "%lf", &vals[i+2][j+2].px);//accounts for ghost cells
	}
    }
  fclose(infile);
}

void read_ascii_tableConspy(const char *inname, 
			  struct Conserved **vals, int Nx, int Ny)

{
  
  FILE *infile = fopen(inname, "r");
  printf("scanning\n");
  int i, j;

  
  for (i=0; i<Ny; ++i)
    {
      for(j= 0; j<Nx; ++j)
	{
	  fscanf(infile, "%lf", &vals[i+2][j+2].py);//accounts for ghost cells
	}
    }
  fclose(infile);
}


void read_ascii_tableConsE(const char *inname, 
			  struct Conserved **vals, int Nx, int Ny)

{
  
  FILE *infile = fopen(inname, "r");
  printf("scanning\n");
  int i, j;

  
  for (i=0; i<Ny; ++i)
    {
      for(j= 0; j<Nx; ++j)
	{
	  fscanf(infile, "%lf", &vals[i+2][j+2].E);//accounts for ghost cells
	}
    }
  fclose(infile);
}

void read_ascii_tablePrimpre(const char *inname, 
			  struct Primitive **vals, int Nx, int Ny)

{
  
  FILE *infile = fopen(inname, "r");
  printf("scanning\n");
  int i, j;

  
  for (i=0; i<Ny; ++i)
    {
      for(j= 0; j<Nx; ++j)
	{
	  fscanf(infile, "%lf", &vals[i+2][j+2].pre);//accounts for ghost cells
	}
    }
  fclose(infile);
}

void read_ascii_tablePrimvx(const char *inname, 
			  struct Primitive **vals, int Nx, int Ny)

{
  
  FILE *infile = fopen(inname, "r");
  printf("scanning\n");
  int i, j;

  
  for (i=0; i<Ny; ++i)
    {
      for(j= 0; j<Nx; ++j)
	{
	  fscanf(infile, "%lf", &vals[i+2][j+2].vx);//accounts for ghost cells
	}
    }
  fclose(infile);
}

void read_ascii_tablePrimvy(const char *inname, 
			  struct Primitive **vals, int Nx, int Ny)

{
  
  FILE *infile = fopen(inname, "r");
  printf("scanning\n");
  int i, j;

  
  for (i=0; i<Ny; ++i)
    {
      for(j= 0; j<Nx; ++j)
	{
	  fscanf(infile, "%lf", &vals[i+2][j+2].vy);//accounts for ghost cells
	}
    }
  fclose(infile);
}


void read_ascii_tablePrime(const char *inname, 
			  struct Primitive **vals, int Nx, int Ny)

{
  
  FILE *infile = fopen(inname, "r");
  printf("scanning\n");
  int i, j;

  
  for (i=0; i<Ny; ++i)
    {
      for(j= 0; j<Nx; ++j)
	{
	  fscanf(infile, "%lf", &vals[i+2][j+2].e);//accounts for ghost cells
	}
    }
  fclose(infile);
}


void read_ascii_tableMixedrho(const char *inname, 
			  struct Mixed **vals, int Nx, int Ny)

{
  
  FILE *infile = fopen(inname, "r");
  printf("scanning\n");
  int i, j;

  
  for (i=0; i<Ny; ++i)
    {
      for(j= 0; j<Nx; ++j)
	{
	  fscanf(infile, "%lf", &vals[i+2][j+2].rho);
	  //accounts for ghost cells
	}
    }
  fclose(infile);
}

void read_ascii_tableMixedvx(const char *inname, 
			  struct Mixed **vals, int Nx, int Ny)

{
  
  FILE *infile = fopen(inname, "r");
  printf("scanning\n");
  int i, j;

  
  for (i=0; i<Ny; ++i)
    {
      for(j= 0; j<Nx; ++j)
	{
	  fscanf(infile, "%lf", &vals[i+2][j+2].vx);
	  //accounts for ghost cells
	}
    }
  fclose(infile);
}

void read_ascii_tableMixedvy(const char *inname, 
			  struct Mixed **vals, int Nx, int Ny)

{
  
  FILE *infile = fopen(inname, "r");
  printf("scanning\n");
  int i, j;

  
  for (i=0; i<Ny; ++i)
    {
      for(j= 0; j<Nx; ++j)
	{
	  fscanf(infile, "%lf", &vals[i+2][j+2].vy);
	  //accounts for ghost cells
	}
    }
  fclose(infile);
}

void read_ascii_tableMixedpre(const char *inname, 
			  struct Mixed **vals, int Nx, int Ny)

{
  
  FILE *infile = fopen(inname, "r");
  printf("scanning\n");
  int i, j;

  
  for (i=0; i<Ny; ++i)
    {
      for(j= 0; j<Nx; ++j)
	{
	  fscanf(infile, "%lf", &vals[i+2][j+2].pre);
	  //accounts for ghost cells
	}
    }
  fclose(infile);
}

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

void print2DarrayConsPre(struct Primitive **array, int nrows, int ncolumns)
	{
	int i, j;
	printf("now print 2D array \n");
	printf("The Pre values");
	for(i = 0; i < nrows; i++)
		{
		  printf("\n");
		  for(j = 0; j < ncolumns; j++)
		    {
			
		      printf("%2.2f \t", array[i][j].pre);
		    }
		}
	printf("\n");
	}

void print2DarrayConsRho(struct Conserved **array, int nrows, int ncolumns)
	{
	int i, j;
	printf("now print 2D array \n");
	printf("The Rho values");
	for(i = 0; i < nrows; i++)
		{
		  printf("\n");
		  for(j = 0; j < ncolumns; j++)
		    {
			
		      printf("%2.2f \t", array[i][j].rho);
		    }
		}
	printf("\n");
	}

void print2DarrayConsPx(struct Conserved **array, int nrows, int ncolumns)
	{
	int i, j;
	printf("now print 2D array \n");
	printf("The Px values");
	for(i = 0; i < nrows; i++)
		{
		  printf("\n");
		  for(j = 0; j < ncolumns; j++)
		    {
			
		      printf("%2.2f \t", array[i][j].px);
		    }
		}
	printf("\n");
	}

void print2DarrayConsPy(struct Conserved **array, int nrows, int ncolumns)
	{
	int i, j;
	printf("now print 2D array \n");
	printf("The Py values");
	for(i = 0; i < nrows; i++)
		{
		  printf("\n");
		  for(j = 0; j < ncolumns; j++)
		    {
			
		      printf("%2.2f \t", array[i][j].py);
		    }
		}
	printf("\n");
	}

void print2DarrayConsE(struct Conserved **array, int nrows, int ncolumns)
	{
	int i, j;
	printf("now print 2D array \n");
	printf("The E values");
	for(i = 0; i < nrows; i++)
		{
		  printf("\n");
		  for(j = 0; j < ncolumns; j++)
		    {
			
		      printf("%2.2f \t", array[i][j].E);
		    }
		}
	printf("\n");
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







void reconstructX(struct Conserved **U, struct Conserved **dUdx, int Nx,
		  int Ny)
{
  
  double theta = 1.5;
  int i, j;
  //printf("In  reconstructX \n");
  //printf("the dUdx rho values: \n");
  //print2DarrayConsRho(dUdx,Ny,Nx+2);
  //printf("the U rho values: \n");
  //print2DarrayConsRho(U,Ny+4,Nx+4);
  for(i = 2; i< Ny+2; i++)
    {
      for(j = 0; j< Nx+2; ++j)
	{
	  dUdx[i-2][j].rho = minmod(theta*(U[i][j+1].rho-U[i][j].rho),
			   .5*(U[i][j+2].rho-U[i][j].rho),
			   theta*(U[i][j+2].rho-U[i][j+1].rho) );
	  //printf("dUdx[%d][%d].rho: %f \n", i-2,j, dUdx[i-2][j].rho);
	  dUdx[i-2][j].px = minmod(theta*(U[i][j+1].px-U[i][j].px),
			   .5*(U[i][j+2].px-U[i][j].px),
			  theta*(U[i][j+2].px-U[i][j+1].px) );
	  dUdx[i-2][j].py = minmod(theta*(U[i][j+1].py-U[i][j].py),
			   .5*(U[i][j+2].py-U[i][j].py),
			  theta*(U[i][j+2].py-U[i][j+1].py) );
	  dUdx[i-2][j].E = minmod(theta*(U[i][j+1].E-U[i][j].E),
			   .5*(U[i][j+2].E-U[i][j].E),
			 theta*(U[i][j+2].E-U[i][j+1].E) );
	}
    }
  //printf("After the dUdx rho values: \n");
  //print2DarrayConsRho(dUdx,Ny,Nx+2);
  //printf("about to leave  reconstructX \n");
}

void reconstructY(struct Conserved **U, struct Conserved **dUdy, int Nx,
		  int Ny)
{
  double theta = 1.5;
  int i, j;
  
  for(i = 0; i< Ny+2; i++)
    {
      for(j = 2; j< Nx+2; ++j)
	{
	  
	  dUdy[i][j-2].rho = minmod(theta*(U[i+1][j].rho-U[i][j].rho),
			   .5*(U[i+2][j].rho-U[i][j].rho),
			   theta*(U[i+2][j].rho-U[i+1][j].rho) );
	  dUdy[i][j-2].px = minmod(theta*(U[i+1][j].px-U[i][j].px),
			   .5*(U[i+2][j].px-U[i][j].px),
			   theta*(U[i+2][j].px-U[i+1][j].px) );
	  dUdy[i][j-2].py = minmod(theta*(U[i+1][j].py-U[i][j].py),
			   .5*(U[i+2][j].py-U[i][j].py),
			   theta*(U[i+2][j].py-U[i+1][j].py) );
	  dUdy[i][j-2].E = minmod(theta*(U[i+1][j].E-U[i][j].E),
			   .5*(U[i+2][j].E-U[i][j].E),
			   theta*(U[i+2][j].E-U[i+1][j].E) );
	  
	  


	  
	}
    }
  
}

int Flux_at_iph2ndX(struct Conserved **U,struct Conserved **Fiph,int  Nx,
		    int Ny)
{
  //printf("In  Fiphx \n");
  
  int i,j;
  struct Conserved **dUdx;
  dUdx = malloc((Ny)* sizeof(struct Conserved *));
  if( dUdx == NULL)
    {
      fprintf(stderr, "out of memory\n");
      return 0;
    }
  for(i = 0; i < (Ny); i++)
    {
      dUdx[i] = malloc((Nx+2) * sizeof(struct Conserved));
      if(dUdx[i] == NULL )
	{
	  fprintf(stderr, "out of memory\n");
	  return 0;
	}
    }

  reconstructX(U, dUdx, Nx, Ny);
  
  struct Conserved Ul;
  struct Conserved Ur;
  //printf("In  Fiphx after dUdx is assigned \n");
  for(i=2; i< Ny+2; i++)
    {
      for(j = 0; j< Nx+1; ++j)
	{
	  Ul.rho = U[i][j+1].rho + .5*dUdx[i-2][j].rho;
	  Ul.px = U[i][j+1].px + .5*dUdx[i-2][j].px;
	  Ul.py = U[i][j+1].py + .5*dUdx[i-2][j].py;
	  Ul.E = U[i][j+1].E + .5*dUdx[i-2][j].E;

	  Ur.rho = U[i][j+2].rho - .5*dUdx[i-2][j+1].rho;
	  Ur.px = U[i][j+2].px - .5*dUdx[i-2][j+1].px;
	  Ur.py = U[i][j+2].py - .5*dUdx[i-2][j+1].py;
	  Ur.E = U[i][j+2].E - .5*dUdx[i-2][j+1].E;
	  //printf("i %d, j %d\n", i,j);
	  //printf("Ur.rho %f\n", Ur.rho);
	  //printf("Ul.rho %f\n", Ul.rho);
	  Fiph[i-2][j] = F_HLLX(Ur, Ul);
	}
    }
  //printf("After the Fiph rho values: \n");
  //print2DarrayConsRho(Fiph,Ny,Nx+1);
  //printf("After the dUdx rho values: \n");
  //print2DarrayConsRho(dUdx,Ny,Nx+2);
  //printf("about to leave  fiphX \n");
  for(i = 0; i < Ny; i++)
    {
      free(dUdx[i]);
	
    }
  free(dUdx);
  //printf("leaving fiphX\n");
  return 0;
  //need to free stuff
}
    
int Flux_at_iph2ndY(struct Conserved **U,struct Conserved **Giph,int  Nx,
		    int Ny)
{
  
  //printf("In  FiphY \n");
  int i,j;
  struct Conserved **dUdy;
  dUdy = malloc((Ny+2)* sizeof(struct Conserved *));
  if( dUdy == NULL)
    {
      fprintf(stderr, "out of memory\n");
      return 0;
    }
  for(i = 0; i < (Ny+2); i++)
    {
      dUdy[i] = malloc((Nx) * sizeof(struct Conserved));
      if(dUdy[i] == NULL )
	{
	  fprintf(stderr, "out of memory\n");
	  return 0;
	}
    }
  
  
  reconstructY(U, dUdy, Nx, Ny);


  //printf("In  FiphY after reconstruct \n");
  struct Conserved Ul;
  struct Conserved Ur;

  //printf("in FiphY the U rho values: \n");
  //print2DarrayConsRho(U,Ny+4,Nx+4);
  for(i=0; i< Ny+1; i++)
    {
      for(j = 0; j< Nx; ++j)
	{
	  Ul.rho = U[i+1][j+2].rho + .5*dUdy[i][j].rho;
	  Ul.px = U[i+1][j+2].px + .5*dUdy[i][j].px;
	  Ul.py = U[i+1][j+2].py + .5*dUdy[i][j].py;
	  Ul.E = U[i+1][j+2].E + .5*dUdy[i][j].E;

	  Ur.rho = U[i+2][j+2].rho - .5*dUdy[i+1][j].rho;
	  Ur.px = U[i+2][j+2].px - .5*dUdy[i+1][j].px;
	  Ur.py = U[i+2][j+2].py - .5*dUdy[i+1][j].py;
	  Ur.E = U[i+2][j+2].E - .5*dUdy[i+1][j].E;

	  

	  Giph[i][j] = F_HLLY(Ur, Ul);
	}
    }
  
  
  for(i = 0; i < Ny+2; i++)
    {
      
      free(dUdy[i]);	
    }
  free(dUdy);
  
  //printf("Leaving  FiphY \n");
  
  return 0;
  //need to free stuff
}

int L(struct Conserved **U, struct Conserved **LU, int Nx, 
	  int Ny, double dx, double dy)
{
  //printf("In L \n");
  int i,j;
  struct Conserved **Fiph;
  struct Conserved **Giph;
  Fiph = malloc((Ny)* sizeof(struct Conserved *));
  Giph = malloc((Ny+1)* sizeof(struct Conserved *));
  if(Fiph == NULL)
    {
      fprintf(stderr, "out of memory\n");
      return 0;
    }
  for(i = 0; i < (Ny); i++)
    {
      Fiph[i] = malloc((Nx+1) * sizeof(struct Conserved));
      if(Fiph[i] == NULL)
	{
	  fprintf(stderr, "out of memory\n");
	  return 0;
	}
    }
  if(Giph == NULL )
    {
      fprintf(stderr, "out of memory\n");
      return 0;
    }
  for(i = 0; i < (Ny+1); i++)
    {
      Giph[i] = malloc((Nx) * sizeof(struct Conserved));
      if(Giph[i] == NULL )
	{
	  fprintf(stderr, "out of memory\n");
	  return 0;
	}
    }
  
  Flux_at_iph2ndX(U, Fiph, Nx, Ny);
  Flux_at_iph2ndY(U, Giph, Nx, Ny);
  
  
  for(j=1; j<Nx+1; ++j)
    {
      for(i=1; i<Ny+1; ++i)
	{
	  LU[i-1][j-1].rho = (- Fiph[i-1][j].rho + Fiph[i-1][j-1].rho)/dx + 
	                (- Giph[i][j-1].rho + Giph[i-1][j-1].rho)/dy;
	  LU[i-1][j-1].px = (- Fiph[i-1][j].px + Fiph[i-1][j-1].px)/dx + 
	                (- Giph[i][j-1].px + Giph[i-1][j-1].px)/dy;
	  LU[i-1][j-1].py = (- Fiph[i-1][j].py + Fiph[i-1][j-1].py)/dx + 
	                    (- Giph[i][j-1].py + Giph[i-1][j-1].py)/dy
	  -U[i+1][j+1].rho*g;
	  LU[i-1][j-1].E = (- Fiph[i-1][j].E + Fiph[i-1][j-1].E)/dx + 
	                   (- Giph[i][j-1].E + Giph[i-1][j-1].E)/dy
	    -U[i+1][j+1].py*g;
	}
    }
  //printf("After the LU rho values: \n");
  //print2DarrayConsRho(LU,Ny,Nx);
  //printf("about to leave  L \n");
  //Need to free Fiph and Giph

  for(i = 0; i < Ny; i++)
    {
	free(Fiph[i]);
	
    }
  free(Fiph);
  for(i = 0; i < Ny+1; i++)
    {
	free(Giph[i]);
	
    }
  free(Giph);


  return 0;
}

//BCs
//------------------------------------------------------
void BC_out(struct Conserved **U, int Nx, int Ny);

void BC_rfl(struct Conserved **U, int Nx, int Ny);

void BC_per(struct Conserved **U, int Nx, int Ny);
void BC_RfP(struct Conserved **U, int Nx, int Ny);
void BC_RfO(struct Conserved **U, int Nx, int Ny);//reflecting outflow

void BC(int BCnum, struct Conserved **U,int  Nx, int Ny)
{
  if(BCnum == 1)
    {
      BC_out(U, Nx, Ny);
    }
  if(BCnum == 2)
    {
      BC_rfl(U, Nx, Ny);
    }
  if(BCnum == 3)
    {
      BC_per(U, Nx, Ny);
    }
  if(BCnum == 4)
    {
      BC_RfP(U, Nx, Ny);
    }
  if(BCnum == 5)
    {
      BC_RfO(U, Nx, Ny);
    }
}

//Outflow BCs
//---------------------------------------------------
void BC_out(struct Conserved **U, int Nx, int Ny)
{
  int i,j;
  //left
  //right
  for(i=2; i<Ny+2; i++)
    {
      U[i][0] = U[i][2];
      U[i][1] = U[i][2];
      U[i][Nx+2] = U[i][Nx+1];
      U[i][Nx+3] = U[i][Nx+1];
      
    }  
  //top and bottom
  for(j=2; j<Nx+2; j++)
    {
      U[0][j] = U[2][j];
      U[1][j] = U[2][j];
      U[Ny+2][j] = U[Ny+1][j];
      U[Ny+3][j] = U[Ny+1][j];
    }  
  
}

//Reflecting BCs
//---------------------------------------------------
void BC_rfl(struct Conserved **U, int Nx, int Ny)
{
  int i,j;
  //left
  //right the xcomp is opposit
  for(i=2; i<Ny+2; i++)
    {
      U[i][0] = U[i][2];
      U[i][1] = U[i][2];
      U[i][0].px = -U[i][2].px;
      U[i][1].px = -U[i][2].px;
      U[i][Nx+2] = U[i][Nx+1];
      U[i][Nx+3] = U[i][Nx+1];
      U[i][Nx+2].px = -U[i][Nx+1].px;
      U[i][Nx+3].px = -U[i][Nx+1].px;
      
    }  
  //top and bottom the ycomp is opposite
  for(j=2; j<Nx+2; j++)
    {
      U[0][j] = U[3][j];
      U[1][j] = U[2][j];
      U[0][j].py = -U[3][j].py;//should these be both U[2][j]?
      U[1][j].py = -U[2][j].py;
      U[Ny+2][j] = U[Ny][j];
      U[Ny+3][j] = U[Ny+1][j];
      U[Ny+2][j].py = -U[Ny][j].py;
      U[Ny+3][j].py = -U[Ny+1][j].py;
    }  
  
}

//Periodic BCs
//---------------------------------------------------
void BC_per(struct Conserved **U, int Nx, int Ny)
{
  int i,j;
  //left
  //right the xcomp is opposit
  for(i=2; i<Ny+2; i++)
    {
      U[i][0] = U[i][Nx];
      U[i][1] = U[i][Nx+1];
      U[i][Nx+2] = U[i][2];
      U[i][Nx+3] = U[i][3];
            
    }  
  //top and bottom the ycomp is opposite
  for(j=2; j<Nx+2; j++)
    {
      U[0][j] = U[Ny][j];
      U[1][j] = U[Ny+1][j];
      U[Ny+2][j] = U[2][j];
      U[Ny+3][j] = U[3][j];
    }  
  
}

//Reflect bottom top, Periodic sides BCs
//---------------------------------------------------
void BC_RfP(struct Conserved **U, int Nx, int Ny)
{
  int i,j;
  //left
  //right the xcomp is opposit
  for(i=2; i<Ny+2; i++)
    {
      U[i][0] = U[i][Nx];
      U[i][1] = U[i][Nx+1];
      U[i][Nx+2] = U[i][2];
      U[i][Nx+3] = U[i][3];
            
    }  

  //top and bottom the ycomp is opposite
  for(j=2; j<Nx+2; j++)
    {
      U[0][j] = U[2][j];
      U[1][j] = U[2][j];
      U[0][j].py = -U[2][j].py;
      U[1][j].py = -U[2][j].py;
      U[Ny+2][j] = U[Ny+1][j];
      U[Ny+3][j] = U[Ny+1][j];
      U[Ny+2][j].py = -U[Ny+1][j].py;
      U[Ny+3][j].py = -U[Ny+1][j].py;
    }  
  
  
}
void BC_RfO(struct Conserved **U, int Nx, int Ny)
{
  int i,j;
  
  //left
  //right
  for(i=2; i<Ny+2; i++)
    {
      U[i][0] = U[i][2];
      U[i][1] = U[i][2];
      U[i][Nx+2] = U[i][Nx+1];
      U[i][Nx+3] = U[i][Nx+1];
      
    }  
  //top and bottom the ycomp is opposite
  for(j=2; j<Nx+2; j++)
    {
      U[0][j] = U[2][j];
      U[1][j] = U[2][j];
      U[0][j].py = -U[2][j].py;
      U[1][j].py = -U[2][j].py;
      U[Ny+2][j] = U[Ny+1][j];
      U[Ny+3][j] = U[Ny+1][j];
      U[Ny+2][j].py = -U[Ny+1][j].py;
      U[Ny+3][j].py = -U[Ny+1][j].py;
    }  
  
}

double evolveF_HLL_2ndO(struct Conserved **U, int Nx, int Ny, double dx,
			double dy, double CFL, double BCnum)
{
  //printf("In time evolve \n");
  int i,j;
  struct Conserved **LU;
  struct Conserved **Utemp;
  LU    = malloc((Ny)* sizeof(struct Conserved *));//there are no ghost cells
                                                   //for LU
  Utemp = malloc((Ny+4)* sizeof(struct Conserved *));
  if(LU == NULL || Utemp == NULL)
    {
      fprintf(stderr, "out of memory\n");
      return 0;
    }
  for(i = 0; i < (Ny+4); i++)
    {
      Utemp[i] = malloc((Nx+4) * sizeof(struct Conserved));
      if(Utemp[i] == NULL )
	{
	  fprintf(stderr, "out of memory\n");
	  return 0;
	}
    }
  for(i = 0; i < (Ny); i++)
    {
      LU[i]    = malloc((Nx) * sizeof(struct Conserved));
      if(LU[i] == NULL)
	{
	  fprintf(stderr, "out of memory\n");
	  return 0;
	}
    }

  double lamdaMax;
  double dt;
  BC(BCnum, U, Nx, Ny);
  lamdaMax = MaxEigenVal(U, Nx, Ny);
  //printf("lamdaMax %f \n", lamdaMax);

  dt = CFL*Min2(dx,dy)/(sqrt(2)*lamdaMax);//Courant Condition
  //printf("dt %f\n", dt);
  //printf("In time evolve after BCs \n");
  L(U, LU, Nx, Ny, dx, dy);//creats dx*dU/dx for all Uij






 
  
  //printf("the LU rho vals in 1st part\n");
  //print2DarrayConsRho(LU, Nx, Ny);






  //printf("dt %f\n", dt);
  
  //third order Runge-Kutta time step
  //-------------------------------------------------------------
  for(i=2; i<Ny+2; i++)
    {
      for(j = 2; j<Nx+2; j++)
	{
	  //Need to switch U to Utemp to get back 3rd order
	  Utemp[i][j].rho = U[i][j].rho + (dt)*LU[i-2][j-2].rho;
	  Utemp[i][j].px = U[i][j].px + (dt)*LU[i-2][j-2].px ;
	  //testint to see how important these factors are
	  // should be
	  //Utemp[i][j].py = U[i][j].py + (dt)*LU[i-2][j-2].py ;
	  Utemp[i][j].py = U[i][j].py + (dt)*LU[i-2][j-2].py ;
	  Utemp[i][j].E = U[i][j].E + (dt)*LU[i-2][j-2].E ;
	}
    }
  
  //printf("In time evolve after part1  first time step \n");
  
  //print2DarrayConsRho(Utemp,Ny+4,Nx+4);
  //print2DarrayConsPy(Utemp,Ny+4,Nx+4);
  //print2DarrayConsE(Utemp,Ny+4,Nx+4);
  
  
  BC(BCnum, Utemp, Nx, Ny);
  L(Utemp, LU, Nx, Ny, dx, dy);
  
  for(i=2; i<Ny+2; i++)
    {
      for(j = 2; j<Nx+2; j++)
	{
	  Utemp[i][j].rho = .75*U[i][j].rho + .25*Utemp[i][j].rho +
	    .25*(dt)*(LU[i-2][j-2].rho);
	  Utemp[i][j].px = .75*U[i][j].px + .25*Utemp[i][j].px +
	                   .25*(dt)*LU[i-2][j-2].px;
	  Utemp[i][j].py = .75*U[i][j].py + .25*Utemp[i][j].py +
	                   .25*(dt)*LU[i-2][j-2].py;
	  Utemp[i][j].E = .75*U[i][j].E + .25*Utemp[i][j].E +
	                   .25*(dt)*LU[i-2][j-2].E;
	}
    }
  
  
  BC(BCnum, Utemp, Nx, Ny);
  L(Utemp, LU, Nx, Ny, dx, dy);

  for(i=2; i<Ny+2; i++)
    {
      for(j = 2; j<Nx+2; j++)
	{
	  U[i][j].rho = (1/3.0)*U[i][j].rho + (2/3.0)*Utemp[i][j].rho +
	                   (2/3.0)*(dt)*LU[i-2][j-2].rho;
	  U[i][j].px = (1/3.0)*U[i][j].px + (2/3.0)*Utemp[i][j].px +
	                   (2/3.0)*(dt)*LU[i-2][j-2].px;
	  U[i][j].py = (1/3.0)*U[i][j].py + (2/3.0)*Utemp[i][j].py +
	                   (2/3.0)*(dt)*LU[i-2][j-2].py;
	  U[i][j].E = (1/3.0)*U[i][j].E + (2/3.0)*Utemp[i][j].E +
	                   (2/3.0)*(dt)*LU[i-2][j-2].E;
	}
    }

  
  //-------------------------------------------
  //printf("In time evolve after part3 first time step \n");
  //print2DarrayConsPy(U,Ny+4,Nx+4);
  //free stuff
  for(i = 0; i < Ny; i++)
    {
	free(LU[i]);
	
    }
  free(LU);
  for(i = 0; i < Ny+4; i++)
    {
	free(Utemp[i]);
	
    }
  free(Utemp);
  
  return dt;
}



int main(int argc, char **argv)
{
  //Get Parameters
  //----------------------------------------------------
  int    Nx    =  get_int_param("param.cfg", "Nx", 10);
  int    Ny    =  get_int_param("param.cfg", "Ny", 10);
  int    BCnum = get_int_param("param.cfg", "BCnum", 1);
  int    InitialType = get_int_param("param.cfg", "InitialType", 1);
  double runTime   = get_double_param("param.cfg", "runTime", 0);
  double time   = get_double_param("time.cfg", "t", 0);
  double x0    = -1.0;
  double x1    = 1.0;
  double y0    = -1.00;
  double y1    = 1.00;

  //put in initial values of Ur Ul
  struct Conserved Ur;
  struct Conserved Ul;
  struct Mixed Mr,Ml;
  Mr.pre = get_double_param("param.cfg", "Mr.pre", 0);
  Mr.rho = get_double_param("param.cfg", "Mr.rho", 0);
  Mr.vx   = get_double_param("param.cfg", "Mr.vx", 0);
  Mr.vy   = get_double_param("param.cfg", "Mr.vy", 0);
  
  Ml.pre = get_double_param("param.cfg", "Ml.pre", 0);
  Ml.rho = get_double_param("param.cfg", "Ml.rho", 0);
  Ml.vx   = get_double_param("param.cfg", "Ml.vx", 0);
  Ml.vy   = get_double_param("param.cfg", "Ml.vy", 0);
  
  Ur = mixed_to_cons(Mr);
  Ul = mixed_to_cons(Ml);
  //------------------------------------------

  //Gets the dx
  double dx;
  double *x = (double*) calloc(Nx+4, sizeof(double));
  linspace2nd(x, x0, x1, Nx, &dx); 
  //Gets the dx
  double dy;
  double *y = (double*) calloc(Ny+4, sizeof(double));
  linspace2nd(y, y0, y1, Ny, &dy);

  //x and y values 


  //Setting the Initial Conditions
  //allocate 2D arrays of structs //Now x y are in the struct
  //-----------------------------------------------------------------------
  int i,j;
  struct Conserved **U;
  struct Primitive **P;
  struct Mixed     **M;
  U = malloc((Ny+4)* sizeof(struct Conserved *));
  P = malloc((Ny+4)* sizeof(struct Conserved *));
  M = malloc((Ny+4)* sizeof(struct Conserved *));
  if(U == NULL || P == NULL || M == NULL)
    {
      fprintf(stderr, "out of memory\n");
      return 0;
    }
  for(i = 0; i < (Ny+4); i++)
    {
      U[i] = malloc((Nx+4) * sizeof(struct Conserved));
      P[i] = malloc((Nx+4) * sizeof(struct Conserved));
      M[i] = malloc((Nx+4) * sizeof(struct Conserved));
      if(U[i] == NULL || P[i] == NULL  || M[i] == NULL )
	{
	  fprintf(stderr, "out of memory\n");
	  return 0;
	}
    }


  

  //Initial Conditions
  //------------------------------------------------------------------
  if(InitialType == 0)
    {
      read_ascii_tableConsRho("run/FinalRho.cfg", U, Nx, Ny);
      read_ascii_tableConspx("run/Finalpx.cfg", U, Nx, Ny);
      read_ascii_tableConspy("run/Finalpy.cfg", U, Nx, Ny);
      read_ascii_tableConsE("run/FinalE.cfg", U, Nx, Ny);
      
      
    }
  else if(InitialType == 1)
    {
      read_ascii_tableMixedrho("run/Initialrho.cfg", M, Nx, Ny);
      read_ascii_tableMixedvx("run/Initialvx.cfg", M, Nx, Ny);
      read_ascii_tableMixedvy("run/Initialvy.cfg", M, Nx, Ny);
      read_ascii_tableMixedpre("run/Initialpre.cfg", M, Nx, Ny);
      for(i = 2 ; i < (Ny+2); i++)
	{
	  for(j = 2; j<(Nx+2); j++)
	    {
	      U[i][j] = mixed_to_cons(M[i][j]);
	    }
	}	    
      
    }
  else if(InitialType == 2)
    {
      //Shock: high pressure top, low bottom
      for(i = 2 ; i < (Ny+2); i++)
	{
	  for(j = 2; j<(Nx+2); j++)
	    {
	      if(y[i]>= 0)
		{
		  U[i][j]=Ul;
		}
	      else
		{
		  U[i][j] = Ur;
		}
	  	  
	    }
	}
    }
  else if(InitialType == 3)
    {
      //Shock Diagonal: high pressure bottom right, top left
      for(i = 2 ; i < (Ny+2); i++)
	{
	  for(j = 2; j<(Nx+2); j++)
	    {
	      if(y[i]>= -x[j] )
		{
		  U[i][j]=Ur;
		}
	      else
		{
		  U[i][j] = Ul;
		}
	  	  
	    }
	}
    }
  else if(InitialType == 4)
    {
      //Kelvin-Helmholtz
      for(i = 2 ; i < (Ny+2); i++)
	{
	  for(j = 2; j<(Nx+2); j++)
	    {
	      if(y[i]>= -.25 && y[i] <= .25 )// inner part
		{
		  
		  U[i][j]=Ur;
		}
	      else
		{
		  U[i][j] = Ul;
		}
	      if(fabs(x[j])<=.2)
		{
		  U[i][j].py = 0.01;
		}
	  	  
	    }
	}
    }
  else if(InitialType == 5)
    {
      //Density advection
      for(i = 2 ; i < (Ny+2); i++)
	{
	  for(j = 2; j<(Nx+2); j++)
	    {
	      /*double r;
	      r = sqrt(pow(x[j],2)+pow(y[i],2));
	      if(r <= 1.0 )// inner part
		{		  
		  U[i][j]=Ur;
		}
	      else
		{
		  U[i][j] = Ul;
		}
	      */

	      if(fabs(x[j])<=1 && fabs(y[i])<=1 )// inner part
		{		  
		  U[i][j]=Ur;
		}
	      else
		{
		  U[i][j] = Ul;
		}
	    }
	}
    }
  else if(InitialType == 6)
    {
      //Circular Pressure Wave //need to set v to zero, density similar, and 
      //just do slightly different pressures
      for(i = 2 ; i < (Ny+2); i++)
	{
	  for(j = 2; j<(Nx+2); j++)
	    {
	      double r;
	      r = sqrt(pow(x[j],2)+pow(y[i],2));
	      if(r <= .25 )// inner part	        
		{		  
		  U[i][j]=Ur;
		}
	      else
		{
		  U[i][j] = Ul;
		}
	    }
	}
    }
  else if(InitialType == 7 || InitialType == 8)
    {
      //For Gravity 
      //Two Gases on top of each other
      //With RayleighTaylor Perturbation
      for(i = 2 ; i < (Ny+2); i++)
	{
	  for(j = 2; j<(Nx+2); j++)
	    {	
	      double rho = Mr.rho;
	      double Ptop = Mr.pre;
	      double rho2 = Ml.rho;
	      if(y[i]>= 0)
		{
		  
		  M[i][j] = Mr;
		  
		  M[i][j].pre = Ptop +rho*g*(y1- y[i]);
		  //Add perturbation
		  if(y[i]<.05 && fabs(x[j])<.025 && InitialType == 8)
		    {
		      M[i][j].vy = -.1;
		    }
		  U[i][j] = mixed_to_cons(M[i][j]);
		}
	      else
		{ 
		  
		  M[i][j] = Ml;
		  
		  M[i][j].pre = Ptop +rho*g*y1+ rho2*g*(-y[i]);
		  
		  U[i][j] = mixed_to_cons(M[i][j]);
		}
	      
	    }
	}
    }
    else if(InitialType == 9)
    {
      
      for(i = 2 ; i < (Ny+2); i++)
	{
	  for(j = 2; j<(Nx+2); j++)
	    {
	      M[i][j] = Mr;
	      double r;
	      r = sqrt(pow(x[j],2)+pow(y[i],2));
	      M[i][j].rho = Mr.rho*exp(-16*pow(r,2))+Ml.rho; 
	      U[i][j] = mixed_to_cons(M[i][j]);
	      /*if(r <= .25 )// inner part	        
		{		  
		  U[i][j]=Ur;
		}
	      else
		{
		  U[i][j] = Ul;
		}
	      */
	    }
	}
    }
    else if(InitialType == 10)
      {
      //For Gravity 
      //Two Gases on top of each other
      //With RayleighTaylor Perturbation
      for(i = 2 ; i < (Ny+2); i++)
	{
	  for(j = 2; j<(Nx+2); j++)
	    {	
	      double rho = Mr.rho;
	      double Ptop = Mr.pre;
	      double rho2 = Ml.rho;
	      double vyPerp = .2;
	      double h=.5;
	      double w = .5;
	      
	      if(y[i]<h && x[j]<(w+.25) && x[j]>.25 )
		{
		    {
		      M[i][j].rho = rho2;
		      M[i][j].pre = Ptop +rho*g*(y1-h)+ rho2*g*(h-y[i]);
		      M[i][j].vx = Ml.vx;
		      M[i][j].vy = vyPerp;
		    }
		  U[i][j] = mixed_to_cons(M[i][j]);
		}
	      else if(y[i]>= 0)
		{
		  
		  M[i][j] = Mr;
		  
		  M[i][j].pre = Ptop +rho*g*(y1- y[i]);
		  //Add perturbation in density
		   U[i][j] = mixed_to_cons(M[i][j]);
		}
	      else
		{ 
		  
		  M[i][j] = Ml;
		  
		  M[i][j].pre = Ptop +rho*g*y1+ rho2*g*(-y[i]);
		  if(x[j]<(w+.25) && x[j]>.25 )
		    {
		      M[i][j].vy = vyPerp;
		    }
		  
		  U[i][j] = mixed_to_cons(M[i][j]);
		}
	      
	    }
	}
      }
      
  //printf("here are the initial conditions: \n");
  //print2DarrayConsRho(U,Ny+4,Nx+4);
  //print2DarrayConsPx(U,Ny+4,Nx+4);
  //print2DarrayConsE(U,Ny+4,Nx+4);
  //----------------------------------------------------------

  //Make sure there is an M array of mixed values
  for(i = 2 ; i < (Ny+2); i++)
	{
	  for(j = 2; j<(Nx+2); j++)
	    {
	      M[i][j] = cons_to_mixed(U[i][j]);
	    }
	}

  //Output initial conditions
  //--------------------------------------------------------------
  write_ascii_table("run/Initialx.cfg", x, Nx,"w");
  write_ascii_table("run/Initialy.cfg", y, Ny,"w");
  
  
  
  write_ascii_tableRhoXY("run/InitialRho.cfg",U, Nx, Ny, "w");
  write_ascii_tablepxXY("run/Initialpx.cfg",U, Nx, Ny, "w");
  write_ascii_tablepyXY("run/Initialpy.cfg",U, Nx, Ny, "w");
  write_ascii_tableEXY("run/InitialE.cfg",U, Nx, Ny, "w");

  write_ascii_tableMixpreXY("run/Initialpre.cfg",M, Nx,Ny,"w");
  write_ascii_tableMixvxXY( "run/Initialvx.cfg",M, Nx,Ny,"w");
  write_ascii_tableMixvyXY( "run/Initialvy.cfg",M, Nx,Ny,"w");



  printf("dx = %f, dy = %f\n", dx,dy);

  //Time Evolution
  //-----------------------------------------------------------------
  double CFL = 0.5;
  double timeLocal =0;
  int k = 0;
  printf("time %f\n", time);
  while(timeLocal<= runTime)
    {
      //printf("In while loop \n");
      k++;
      //printf("step %d \n", k);
      //timeLocal=timeLocal+.1;
      timeLocal = timeLocal + evolveF_HLL_2ndO(U, Nx, Ny, dx, dy, CFL, BCnum);
    }
  time = time + timeLocal;
  printf("Time after the %d th step is: %f \n", k, time);
  write_ascii_time("time.cfg", time,"w");
  //printf("\n");
  //printf("here are the final conditions: \n");
  //print2DarrayConsRho(U,Ny+4,Nx+4);
  //print2DarrayConsPx(U,Ny+4,Nx+4);
  //print2DarrayConsE(U,Ny+4,Nx+4);

  for(i = 2 ; i < (Ny+2); i++)
	{
	  for(j = 2; j<(Nx+2); j++)
	    {
	      M[i][j] = cons_to_mixed(U[i][j]);
	    }
	}
  //Output results
  //--------------------------------------------------------------
  write_ascii_table("run/Finalx.cfg", x, Nx,"w");
  write_ascii_table("run/Finaly.cfg", y, Ny,"w");
  
  
  
  write_ascii_tableRhoXY("run/FinalRho.cfg",U, Nx, Ny, "w");
  write_ascii_tablepxXY("run/Finalpx.cfg",U, Nx, Ny, "w");
  write_ascii_tablepyXY("run/Finalpy.cfg",U, Nx, Ny, "w");
  write_ascii_tableEXY("run/FinalE.cfg",U, Nx, Ny, "w"); 

  write_ascii_tableMixpreXY("run/Finalpre.cfg",M, Nx,Ny,"w");
  write_ascii_tableMixvxXY( "run/Finalvx.cfg",M, Nx,Ny,"w");
  write_ascii_tableMixvyXY( "run/Finalvy.cfg",M, Nx,Ny,"w");

  

  printf("hello \n");
  
  for(i = 0; i < Ny+4; i++)
    {
	free(U[i]);
	free(P[i]);
	free(M[i]);
    }
  free(U);
  free(P);
  free(M);
  free(x);
  free(y);



  return 0;
}
