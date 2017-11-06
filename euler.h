#ifndef EULER_H
#define EULER_H



struct Conserved
{
  double rho, E, px, py;
};
struct Primitive
{
  double pre, e, vx, vy;
};
struct Mixed
{
  double rho, vx, vy, pre;
};
struct Conserved mixed_to_cons(struct Mixed);
struct Mixed cons_to_mixed(struct Conserved);

struct Conserved prim_to_cons(struct Primitive);
struct Primitive cons_to_prim(struct Conserved);

struct Conserved equalNotxy(struct Conserved Uin, struct Conserved Uout);

void EigenvaluesX(struct Primitive, double*, double*);//(inprim,lamdap,lamdam)
void EigenvaluesY(struct Primitive, double*, double*);



double Alphavalue(double, double );//(lamdal,lamdar)


double Min2(double , double);
double Max2(double , double);

struct Conserved fluxX(struct Primitive);
struct Conserved fluxY(struct Primitive);

struct Conserved F_HLLX(struct Conserved, struct Conserved);
struct Conserved F_HLLY(struct Conserved, struct Conserved);
                      //(Ur,Ul)

#endif // EULER_H
