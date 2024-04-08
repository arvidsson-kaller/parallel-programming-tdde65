# Lab1

## A. Ptheads

Every pixel is an indepedent task. Therefore if for example the image is `64x64` pixels, and we have 4 threads, then each thread needs to compute `16x64` tasks.
Every pixel/task requires the same computation (no load balancing).