#include "Utils.h"
#include <string.h>
#include "Vector3d.h"

//---------------------------------------------------------------------------

static float series[NUM_EXPS] =
{0,1,2,3,4,5,6,7,8,9,
   10,20,30,40,50,60,70,80,90,
   100,200,300,400,500,600,700,800,900,
   1000,2000,3000,4000,5000,6000,7000,8000,9000,
   10000,20000,30000,40000,50000,60000,70000,80000,90000,
   100000};

//---------------------------------------------------------------------------
float utilsExp(int i)
//---------------------------------------------------------------------------
{
  if (i < 0) i = 0;
  if (i >= NUM_EXPS) i = NUM_EXPS-1;
  return series[i];
}

//---------------------------------------------------------------------------
int utilsLog(float f)
//---------------------------------------------------------------------------
{
  int i = 0;
  while (i < NUM_EXPS && f > (series[i] + 0.5)) i++;
  if (i >= NUM_EXPS) i = NUM_EXPS - 1;
  return i;
}


//--------------------------------------------------------------------
bool isSeparator(char ch)
//--------------------------------------------------------------------
{
  return (ch <= ' ') || (ch == '=');
}


//--------------------------------------------------------------------
void utilsParseLine(char *s, char *var, char *val)
//--------------------------------------------------------------------
{
  int len = strlen(s);
  int i = 0;
  while ((i < len) && isSeparator(s[i])) i++;
  int j = 0;
  while ((i < len) && !isSeparator(s[i])) {
    var[j] = s[i]; i++; j++;
  }
  var[j] = 0;
  while ((i < len) && isSeparator(s[i])) i++;
  j = 0;
  if ((i < len) && (s[i] == '"')) {
    i++;
    while ((i < len) && (s[i] != '"')) {
      val[j] = s[i]; i++; j++;
    }
  }
  else if ((i < len) && (s[i] == '(')) {
    i++;
    while ((i < len) && (s[i] != ')')) {
      val[j] = s[i]; i++; j++;
    }
  }
  else {
    while ((i < len) && !isSeparator(s[i])) {
      val[j] = s[i]; i++; j++;
    }
  }
  val[j] = 0;
}

