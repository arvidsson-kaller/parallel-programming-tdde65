/*
  File: thresfilter.h
  Declaration of pixel structure and thresfilter function.
 */

#ifndef _THRESFILTERMPI_H_
#define _THRESFILTERMPI_H_

#include "pixel.h"
unsigned averageMPI(const unsigned size, pixel *src);

void thresfilterMPI(const unsigned size, const unsigned average, pixel *src);

#endif
