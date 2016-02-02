#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <inttypes.h>
#include <math.h>

#include "strmap.h"

struct hashAux
{
  FILE *outfile;
  int index;
  double value1;
};
