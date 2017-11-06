#ifndef MYCONFIG_H
#define MYCONFIG_H
void write_ascii_tableCons(const char *outname,
			   double *x, struct Conserved *vals, 
			   int N, char* perm);
void write_ascii_tablePrim(const char *outname,
			   double *x, struct Primitive *vals, int N, char* perm);

void read_ascii_tableCons(const char *inname, double *x, 
			  struct Conserved *vals, int N);
void read_ascii_tablePrim(const char *inname, double *x, 
			  struct Primitive *vals, int N);
#endif
