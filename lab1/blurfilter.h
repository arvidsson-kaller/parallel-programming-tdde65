/*
  File: blurfilter.h
  Declaration of pixel structure and blurfilter function.
 */

#ifndef _BLURFILTER_H_
#define _BLURFILTER_H_

#include "pixel.h"

void blurfilter(const int xsize, const int ysize, pixel* src, const int radius, const double *w);

#endif
