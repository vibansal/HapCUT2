#include <stdarg.h>
#include <stdio.h>
#include<stdlib.h>
#include <time.h>
#include <errno.h>
#include<math.h>
#include "common.h"
#include<string.h>

void check_input_0_or_1(char* x){
    if (!(strcmp(x, "0") == 0 || strcmp(x, "1") == 0)){
        fprintf(stderr, "\nERROR: Invalid input \"%s\" for <0/1> option flag.\n",x);
        exit(1);
    }
}

/** Adapted from http://stackoverflow.com/questions/3673226/how-to-print-time-in-format-2009-08-10-181754-811 */
int fprintf_time(FILE *stream, const char *format, ...)
{
	time_t timer;
	char buffer[26];
	struct tm* tm_info;
	va_list arg;
	int done;

	time(&timer);
	tm_info = localtime(&timer);
	strftime(buffer, 26, "%Y:%m:%d %H:%M:%S", tm_info);

	fprintf(stream, "[%s] ", buffer);

	va_start(arg, format);
	done = vfprintf(stream, format, arg);
	va_end(arg);

	return done;

}

float phred(float x){
	if (x > 0 && x <= 1.05){
		if (x > 1.0){
			x = 1.0;
		}
		return -10.0 * log10(x);
	}else{
		fprintf(stderr, "Invalid prob-to-phred conversion: %f\n", x);
		exit(1);
	}
}

float unphred(float x){
	if (x > 0){
		return pow(10, x / -10.0);
	}else{
		fprintf(stderr, "Invalid phred-to-prob conversion: %f\n", x);
		exit(1);
	}
}


char* concatStrings(char** var_list,int n,char sep)
{
	// concatenate string GT + ':' + DP + ':' from an array 
	int i=0,l=0,j=0,k=0;
	for (i=0;i<n;i++) l += strlen(var_list[i]) +1; 
	char* concat = calloc(l+1,sizeof(char)); 
	for (i=0;i<n;i++)
	{
		for (j=0;j<strlen(var_list[i]);j++) concat[k++] = var_list[i][j];
		concat[k++] = sep; 	
	}	
	concat[k] = '\0';
	return concat; 
}


// split a string using a single separator, '\n' and '\0' are also delimitors at end of string
int splitString(char* input,char sep,char** var_list)
{
	int n=0,i=0,s=0,j=0;
	while (1)
	{
		if (input[i] == sep || input[i] == '\n' || input[i] == '\0')
		{
			if (i-s > 0) {
				var_list[n] = malloc(i-s+1); 
				for (j = s; j < i; j++) var_list[n][j - s] = input[j]; 
				var_list[n][j-s] = '\0';
				n +=1;
			}
			s = i+1; // start of next string
		}
		if (input[i] =='\n' || input[i] == '\0') break;
		i++;
	}
	return n;	
}


// split a string using a single separator, all memory is allocated in the function, not working
int splitString_full(char* input,char sep,char** out)
{
	int n=0,i=0,s=0,j=0;
	while (input[i] != '\n' && input[i] != '\0')
	{
		if (input[i] == sep) n++;
		i++;
	}
	int strings = n+1;
	out = (char**)malloc(sizeof(char*)*strings);
	i=0; n=0;
	while (1)
	{
		if (input[i] == sep || input[i] == '\n' || input[i] == '\0')
		{
			if (i-s > 0) {
				out[n] = malloc(i-s+1); 
				for (j = s; j < i; j++) out[n][j - s] = input[j]; 
				out[n][j-s] = '\0';
				n +=1;
			}
			s = i+1; // start of next string
		}
		if (input[i] =='\n' || input[i] == '\0') break;
		i++;
	}
	return strings;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
