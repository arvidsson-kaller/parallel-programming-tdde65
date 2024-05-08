#include <math.h>

int is_power_of_two(unsigned int x)
{
	// A power of two has only one bit set. So, if x is a power of two,
	// only one bit will be set, and x - 1 will have all the bits set
	// to the right of the original set bit. For example:
	// 8 (1000) -> 7 (0111)
	// 16 (10000) -> 15 (01111)
	// If x is 0, it will return false.
	return (x && !(x & (x - 1)));
}

void calculate_dimensions(int p, int *cols, int *rows)
{
	if (is_power_of_two(p) && p != 1)
	{
		unsigned diff = -1;
		for (int i = 1; i < p; i *= 2)
		{
			unsigned tmp = abs(i - p / i);
			if (tmp < diff)
			{
				diff = tmp;
				*cols = i;
				*rows = p / i;
			}
		}
	}
	else
	{
		double root = sqrt(p);
		unsigned flo = floor(root);
		unsigned cei = ceil(root);
		if (cei * flo == p)
		{
			*cols = cei;
			*rows = flo;
		}
		else
		{
			*cols = p;
			*rows = 1;
		}
	}
	if (*rows > *cols)
	{
		unsigned tmp = *cols;
		*cols = *rows;
		*rows = tmp;
	}
}