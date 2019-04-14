//---------------------------------------------------------------------------

#ifndef UtilsH
#define UtilsH
//---------------------------------------------------------------------------

#include <string.h>

#define STRING_LEN 256
#define NUM_EXPS 47

// discrete exponents for scrollers
float utilsExp(int i);
int   utilsLog(float f);

// parse line of the form "var = val"
void utilsParseLine(char *s, char *var, char *val);

#endif
