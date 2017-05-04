#include <stdarg.h>
#include <stdio.h>
#include<stdlib.h>
#include <time.h>
#include <errno.h>
#include<math.h>
#include "common.h"


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
