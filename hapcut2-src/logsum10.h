//
// logsum10 -- a port of Sean Eddy's fast table-driven log sum, converted to log base10
// This code was originally part of HMMER. This version is used with
// Sean Eddy's permission as public domain code.
//
#ifndef LOGSUMTEN_H
#define LOGSUMTEN_H
#include <assert.h>
#include <stdio.h>
#include <math.h>

/* p7_LOGSUM10_SCALE defines the precision of the calculation; the
 * default of 1000.0 means rounding differences to the nearest 0.001
 * nat. p7_LOGSUM10_TBL defines the size of the lookup table; the
 * default of 16000 means entries are calculated for differences of 0
 * to 16.000 nats (when p7_LOGSUM10_SCALE is 1000.0).  e^{-p7_LOGSUM10_TBL /
 * p7_LOGSUM10_SCALE} should be on the order of the machine FLT_EPSILON,
 * typically 1.2e-7.
 */
#define p7_LOGSUMTEN_TBL   16000
#define p7_LOGSUMTEN_SCALE 1000.f
#define ESL_MAX(a,b)    (((a)>(b))?(a):(b))
#define ESL_MIN(a,b)    (((a)<(b))?(a):(b))
#define eslINFINITY     INFINITY
#define TRUE            1
#define FALSE           0
#define eslOK           1

// initialize the lookup table
// this must be called before any calls to esl_flogsum10
int esl_flogsum10_init(void);

// approximation to log(10^a + 10^b)
static inline float esl_flogsum10(float a, float b)
{
  extern float flogsum10_lookup[p7_LOGSUMTEN_TBL]; 
  const float max = ESL_MAX(a, b);
  const float min = ESL_MIN(a, b);
  return (min == -eslINFINITY || (max-min) >= 15.7f) ? max : max + flogsum10_lookup[(int)((max-min)*p7_LOGSUMTEN_SCALE)];
}

#endif
