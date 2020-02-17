#include <stdlib.h>
#include "logsum10.h"

// statically allocate the lookup table
float flogsum10_lookup[p7_LOGSUMTEN_TBL]; 

// initialize
int esl_flogsum10_init(void)
{
  static int firsttime = TRUE;
  if (!firsttime) {
      return eslOK;
  }
  firsttime = FALSE;

  int i;
  for (i = 0; i < p7_LOGSUMTEN_TBL; i++) {
    flogsum10_lookup[i] = log10(1. + pow(10., (double) -i / p7_LOGSUMTEN_SCALE));
  }
  return eslOK;
}

// exact calculation that this package approximates
float exact_logsum10(float a, float b)
{
    return log10(pow(10, a) + pow(10, b));
}

// calculate the error between the exact and approximate functions
float esl_flogsum_10_error(float a, float b)
{
  float approx = esl_flogsum10(a,b);
  float exact = exact_logsum10(a, b);
  return (pow(10, approx) - pow(10, exact));
}

//
// Code below is conditionally compiled when this symbol is defined
//
#ifdef LOGSUM10_EXAMPLE
/* gcc -o example -g -O2 -I. -L. -DLOGSUM10_EXAMPLE logsum.c -lm
 * ./example -0.5 -0.5
 */

int
main(int argc, char **argv)
{
  float a = atof(argv[1]);
  float b = atof(argv[2]);
  float result;

  esl_flogsum10_init();
  result = esl_flogsum10(a, b);
  printf("esl_flogsum10(%f,%f) = %f\n", a, b, result);

  result = exact_logsum10(a, b);
  printf("log10(10^%f + 10^%f) = %f\n", a, b, result);

  printf("Absolute error in probability: %f\n", esl_flogsum_10_error(a,b));
  return eslOK;
}
#endif /*p7LOGSUM10_EXAMPLE*/

