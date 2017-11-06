
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h> 

#define STRING_LENGTH 128
#define gamma 1.5
#define g .1

#include "euler.h"

//Functions for dealing with conserved and primitive quantities
//---------------------------------------------------------------------
/*struct Conserved
{
  double rho, E, px, py;
};

struct Primitive
{
  double pre, e, vx, vy;
}; */

//Converts from input for shock tube to Conserved
//---------------------------------------------------------------
struct Conserved mixed_to_cons(struct Mixed mix)
{
  struct Conserved U;
  U.rho = mix.rho;
  U.px  = mix.vx*mix.rho;
  U.py  = mix.vy*mix.rho;
  U.E   = 0.5*mix.rho*(pow(mix.vx,2)+pow(mix.vy,2)) + mix.pre/(gamma-1);
  return U;
}

//---------------------------------------------------------------
struct Mixed cons_to_mixed(struct Conserved U)
{
  struct Mixed M;
  M.rho = U.rho;
  M.vx  = U.px/U.rho;
  M.vy  = U.py/U.rho;
  M.pre = (gamma-1)*(U.E - (pow(U.px,2)+pow(U.py,2))/(2*U.rho));
  return M;
}

//--------------------------------------------------------
struct Conserved prim_to_cons(struct Primitive inprim)
{
  double P;
  double e;
  double vx;
  double vy;
  //double gamma;

  
  struct Conserved outcons;

  P = inprim.pre;
  e = inprim.e;
  vx = inprim.vx;
  vy = inprim.vy;
   
  
  outcons.E = P*(0.5*(pow(vx,2) + pow(vy,2))+e)/((gamma-1)*e);
  outcons.px = P*vx/((gamma-1)*e);
  outcons.py = P*vy/((gamma-1)*e);
  outcons.rho = P/((gamma-1)*e);

  return outcons;
}

//---------------------------------------------------
struct Primitive cons_to_prim(struct Conserved incons)
{
  double rho;
  double px;
  double py;
  double E;
  //double gamma;

  

  struct Primitive outprim;

  rho = incons.rho;
  px = incons.px;
  py = incons.py;
  E = incons.E;
   
  
  outprim.pre = (gamma-1)*(E-0.5*(pow(px,2)+pow(py,2))/rho);
  outprim.vx = px/rho;
  outprim.vy = py/rho;
  outprim.e = E/rho - 0.5*(pow(px,2)+pow(py,2))/pow(rho,2);

  return outprim;
}

struct Conserved equalNotxy(struct Conserved Uin, struct Conserved Uout)
{
  Uout.rho = Uin.rho;
  Uout.px = Uin.px;
  Uout.py = Uin.py;
  Uout.E = Uin.E;

  return Uout;
}
//Flux at for a given value of U in x direction
//--------------------------------------------
struct Conserved fluxX(struct Primitive prim)
{
  double P;
  double e;
  double vx;
  double vy;

  struct Conserved cons;
  struct Conserved outFlux;
  double rho;
  double px;
  double py;
  double E;
  //double gamma;

  //gamma = 1.5;//for a 3/2 gas

  P = prim.pre;
  e = prim.e;
  vx = prim.vx;
  vy = prim.vy;
  
  cons = prim_to_cons(prim);

  rho = cons.rho;
  E   = cons.E;
  px  = cons.px;
  py  = cons.py;

  //printf("%f %f %f \n" , P, E, vx);

  outFlux.rho = rho*vx;
  outFlux.px   = rho*pow(vx,2) + P;
  outFlux.py   = rho*vx*vy;
  outFlux.E  = (E+P)*vx;

  return outFlux;
}

//Flux at for a given value of U in y direction
//--------------------------------------------
struct Conserved fluxY(struct Primitive prim)
{
  double P;
  double e;
  double vx;
  double vy;

  struct Conserved cons;
  struct Conserved outFlux;
  double rho;
  double px;
  double py;
  double E;
  //double gamma;

  //gamma = 1.5;//for a 3/2 gas

  P = prim.pre;
  e = prim.e;
  vx = prim.vx;
  vy = prim.vy;
  
  cons = prim_to_cons(prim);

  rho = cons.rho;
  E   = cons.E;
  px  = cons.px;
  py  = cons.py;

  //printf("%f %f %f \n" , P, E, vx);

  outFlux.rho = rho*vy;
  outFlux.px   = rho*vx*vy;
  outFlux.py   = rho*pow(vy,2) + P;
  outFlux.E  = (E+P)*vy;

  return outFlux;
}

//Eigenvalues of the flux jacobian, effective wave speeds
//---------------------------------------------------------------
void EigenvaluesX(struct Primitive inprim, double* lamdap, double* lamdam)
{
  
  double vx;
  double vy;
  double P;
  double rho;
  double cs;

  vx = inprim.vx;
  vy = inprim.vy;
  P  = inprim.pre;
  rho = prim_to_cons(inprim).rho;

  //double gamma;
  //gamma = 1.5;//for a 3/2 gas
 
  cs = pow(gamma*P/rho,0.5);
  //v2= pow(vx,2)+pow(vy,2);
  //v = sqrt(v2);
  
  
  //should the eigenvalues use the absolute value of velocity?
  //or just the component
  //I think it should be just the component because this slows down 
  // the waves that have to go against the flow
  *lamdap = vx + cs;
  *lamdam = vx - cs;
}

void EigenvaluesY(struct Primitive inprim, double* lamdap, double* lamdam)
{
  
  double vx;
  double vy;
  double P;
  double rho;
  double cs;

  vx = inprim.vx;
  vy = inprim.vy;
  P  = inprim.pre;
  rho = prim_to_cons(inprim).rho;

 
  cs = pow(gamma*P/rho,0.5);
  
  
  *lamdap = vy + cs;
  *lamdam = vy - cs;
}



// Minimum of two numbers
//---------------------------------
double Min2(double num1, double num2)
{
  if(num1>num2)
    {
      return num2;
    }
  else
    {
      return num1;
    }
}
// Max of two numbers
//---------------------------------
double Max2(double num1, double num2)
{
  if(num1<num2)
    {
      return num2;
    }
  else
    {
      return num1;
    }
}

//Finds Alphas needed for HLL function
//-------------------------------------------
double Alphavalue(double lamdal, double lamdar)
{
  double alpha;
  alpha = 0;
  alpha = Max2(0 , Max2(lamdal,lamdar) );
  return alpha; 
}



//The HLL flux X
//------------------------------------------------
struct Conserved F_HLLX(struct Conserved Ur, struct Conserved Ul)
{
  struct Conserved fluxHLL;

  struct Conserved fluxr;
  struct Conserved fluxl;
  struct Primitive Pr;
  struct Primitive Pl;

  double lamdarP = 0;
  double lamdarM = 0;
  double lamdalP = 0;
  double lamdalM = 0;

  double alphaP = 0;
  double alphaM = 0;
  
  Pr = cons_to_prim(Ur);
  Pl = cons_to_prim(Ul);

  fluxr = fluxX(Pr);
  fluxl = fluxX(Pl);

  EigenvaluesX(Pr, &lamdarP, &lamdarM);
  EigenvaluesX(Pl, &lamdalP, &lamdalM);

  //printf("The eigenvales \n");
  //printf("lamdarP: %f \n", lamdarP);
  //printf("lamdarM: %f \n", lamdarM);
  //printf("lamdalP: %f \n", lamdalP);
  //printf("lamdalM: %f \n", lamdalM);

  alphaP = Alphavalue(lamdalP, lamdarP);
  alphaM = Alphavalue(-lamdalM, -lamdarM);

  //printf("fluxr.px = %f \n" , fluxr.px);
  //printf("fluxl.px = %f \n" , fluxl.px);

  fluxHLL.rho = (alphaP*fluxl.rho + alphaM*fluxr.rho - 
		 alphaP*alphaM*(Ur.rho-Ul.rho) )/(alphaP+alphaM);
  //printf("fluxl.rho %f \n", fluxl.rho);
  //printf("fluxr.rho %f \n", fluxr.rho);
  //printf("alphaP %f \n", alphaP);
  //printf("alphaM %f \n", alphaM);
  //printf("%f \n", fluxHLL.rho);

  fluxHLL.E = (alphaP*fluxl.E + alphaM*fluxr.E - 
		 alphaP*alphaM*(Ur.E-Ul.E) )/(alphaP+alphaM);

  fluxHLL.px = (alphaP*fluxl.px + alphaM*fluxr.px - 
		 alphaP*alphaM*(Ur.px-Ul.px) )/(alphaP+alphaM);

  //printf("fluxl.px %f \n", fluxl.px);
  //printf("fluxr.px %f \n", fluxr.px);
  //printf("%f \n", fluxHLL.px);
  fluxHLL.py = (alphaP*fluxl.py + alphaM*fluxr.py - 
		 alphaP*alphaM*(Ur.py-Ul.py) )/(alphaP+alphaM);

  return fluxHLL;
}

//The HLL flux Y
//--------------------------------------------------------------------
struct Conserved F_HLLY(struct Conserved Ur, struct Conserved Ul)
{
  struct Conserved fluxHLL;

  struct Conserved fluxr;
  struct Conserved fluxl;
  struct Primitive Pr;
  struct Primitive Pl;

  double lamdarP = 0;
  double lamdarM = 0;
  double lamdalP = 0;
  double lamdalM = 0;

  double alphaP = 0;
  double alphaM = 0;
  
  Pr = cons_to_prim(Ur);
  Pl = cons_to_prim(Ul);

  fluxr = fluxY(Pr);
  fluxl = fluxY(Pl);

  EigenvaluesY(Pr, &lamdarP, &lamdarM);
  EigenvaluesY(Pl, &lamdalP, &lamdalM);

  //printf("The eigenvales \n");
  //printf("lamdarP: %f \n", lamdarP);
  //printf("lamdarM: %f \n", lamdarM);
  //printf("lamdalP: %f \n", lamdalP);
  //printf("lamdalM: %f \n", lamdalM);

  alphaP = Alphavalue(lamdalP, lamdarP);
  alphaM = Alphavalue(-lamdalM, -lamdarM);

  //printf("fluxr.px = %f \n" , fluxr.px);
  //printf("fluxl.px = %f \n" , fluxl.px);

  fluxHLL.rho = (alphaP*fluxl.rho + alphaM*fluxr.rho - 
		 alphaP*alphaM*(Ur.rho-Ul.rho) )/(alphaP+alphaM);
  //printf("fluxl.rho %f \n", fluxl.rho);
  //printf("fluxr.rho %f \n", fluxr.rho);
  //printf("alphaP %f \n", alphaP);
  //printf("alphaM %f \n", alphaM);
  //printf("%f \n", fluxHLL.rho);

  fluxHLL.E = (alphaP*fluxl.E + alphaM*fluxr.E - 
		 alphaP*alphaM*(Ur.E-Ul.E) )/(alphaP+alphaM);

  fluxHLL.px = (alphaP*fluxl.px + alphaM*fluxr.px - 
		 alphaP*alphaM*(Ur.px-Ul.px) )/(alphaP+alphaM);

  //printf("fluxl.px %f \n", fluxl.px);
  //printf("fluxr.px %f \n", fluxr.px);
  //printf("%f \n", fluxHLL.px);
  fluxHLL.py = (alphaP*fluxl.py + alphaM*fluxr.py - 
		 alphaP*alphaM*(Ur.py-Ul.py) )/(alphaP+alphaM);

  return fluxHLL;
}
