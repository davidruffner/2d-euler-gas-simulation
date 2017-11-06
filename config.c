
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define STRING_LENGTH 128


// Prototypes for private functions
// ---------------------------------------------------------------------
static int read_par(const char *fname, const char *parnam, char *parval);


// Public function definitions
// ---------------------------------------------------------------------
double get_double_param(const char *fname, const char *parnam, double defval)
{
  char parval[STRING_LENGTH];

  if (read_par(fname, parnam, parval))
    {
      printf("Parameter %s not found in %s. Using default: %f.\n", parnam, fname, defval);
      return defval;
    }
  else
    {
      return atof(parval);
    }
}
int get_int_param(const char *fname, const char *parnam, int defval)
{
  char parval[STRING_LENGTH];

  if (read_par(fname, parnam, parval))
    {
      printf("Parameter %s not found in %s. Using default: %d.\n", parnam, fname, defval);
      return defval;
    }
  else
    {
      return atoi(parval);
    }
}


// Private function definitions
// ---------------------------------------------------------------------
static int read_par(const char *fname, const char *parnam, char *parval)
{

  FILE *pfile = fopen(fname, "r");
  char *charp;
  char  line[STRING_LENGTH];
  char  lhs [STRING_LENGTH];

  while ((charp = fgets(line, STRING_LENGTH, pfile)))
    {

      // skip lines beginning with #
      if ((charp = strchr(line, '#')))
	*charp = '\0';

      // break the line into variable name and value 
      if (strchr(line, '='))
	{
	  int i;
	  char *cp0, *cp1;
	  cp0 = strtok(line, "=");
	  cp1 = strtok(NULL, "=");

	  for (i=0; i<STRING_LENGTH; ++i)
	    {
	      if ((charp=strchr(cp0, '\t'))) *charp = ' ';
	      else break;
	    }

	  cp0 = strtok(cp0, " ");

	  for (i=0; i<STRING_LENGTH; ++i)
	    {
	      if ((charp = strchr(cp1, '\t'))) *charp = ' ';
	      else break;
	    }
	  for (i=0; i<STRING_LENGTH; ++i)
	    {
	      if ((charp = strchr(cp1, '"'))) *charp = ' ';
	      else break;
	    }
	  for (i=0; i<STRING_LENGTH; ++i)
	    {
	      if ((charp = strchr(cp1, '\n'))) *charp = ' ';
	      else break;
	    }

	  cp1 = strtok(cp1, " ");
	  strcpy(lhs, cp0);

	  if (strcmp(lhs, parnam) == 0)
	    {
	      strcpy(parval, cp1);
	      fclose(pfile);
	      return 0;
	    }
	}
    }

  fclose(pfile);
  return 1;
}
