/*
  File: blurfilter.h
  Declaration of pixel structure and blurfilter function.
 */

#ifndef _BLURFILTERPTHREADS_H_
#define _BLURFILTERPTHREADS_H_

#include "pixel.h"

void blurfilterPthreads(const int xsize, const int ysize, pixel* src, const int radius, const double *w);

#endif
