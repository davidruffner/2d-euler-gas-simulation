#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h> 

#define STRING_LENGTH 128


#include "euler.h"




// Output for plotting
// -----------------------------------------------------------------------
void write_ascii_tableCons(const char *outname,
			   double *x, struct Conserved *vals, int N, char* perm)
{
  FILE *outfile = fopen(outname, perm);
  int i;
  for (i=2; i<N+2; ++i)//for 2nd euler
    {
      fprintf(outfile, "%f %f %f %f\n", x[i], vals[i].rho, vals[i].E, 
	      vals[i].px);
    }
  fclose(outfile);
}

// -----------------------------------------------------------------------
void write_ascii_tablePrim(const char *outname,
			   double *x, struct Primitive *vals, 
			   int N, char* perm)
{
  FILE *outfile = fopen(outname, perm);
  int i;
  for (i=2; i<N+2; ++i)//for 2nd euler
    {
      fprintf(outfile, "%f %f %f %f\n", x[i], vals[i].pre, vals[i].e, 
	      vals[i].vx);
    }
  fclose(outfile);
}

// -----------------------------------------------------------------------
void read_ascii_tableCons(const char *inname, double *x, 
			  struct Conserved *vals, int N)
//This N is the actual N it reads
{
  
  FILE *infile = fopen(inname, "r");
  int i;

  
  for (i=2; i<N+2; ++i)
    {
      fscanf(infile, "%lf %lf %lf %lf\n", &x[i], &vals[i].rho, &vals[i].E, 
	      &vals[i].px);
    }
  fclose(infile);
}

// -----------------------------------------------------------------------
void read_ascii_tablePrim(const char *inname, double *x, 
			  struct Primitive *vals, int N)
//This N is the actual N it reads
{
  
  FILE *infile = fopen(inname, "r");
  int i;

  
  for (i=2; i<N+2; ++i)
    {
      fscanf(infile, "%lf %lf %lf %lf\n", &x[i], &vals[i].pre, &vals[i].e, 
	      &vals[i].vx);
    }
  fclose(infile);
}
