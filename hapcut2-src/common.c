#include <stdarg.h>
#include <stdio.h>
#include <time.h>
#include <errno.h>
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
