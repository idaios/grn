#include "binary_network.h"

#define MAXSTEPS 10000

#define BASE 3

#define CHAR_BIT 8

#ifdef _U128
typedef unsigned __int128 ui;

typedef unsigned __int128 uint128_t;
#else
typedef unsigned long long  ui;
#endif

typedef ui T;

typedef unsigned char uc;

unsigned int samplestep;

FILE *selFile;

FILE *genoutFile;

FILE *mutfile;

FILE *choicefile;

#ifdef _U128

/*      UINT64_MAX 18446744073709551615ULL */
#define P10_UINT64 10000000000000000000ULL   /* 19 zeroes */

#define E10_UINT64 19

#define STRINGIZER(x)   # x

#define TO_STRING(x)    STRINGIZER(x)

static int print_u128_u(uint128_t  u128)
{
  
  int rc;
  if (u128 > UINT64_MAX)
    {
      uint128_t leading  = u128 / P10_UINT64;
      uint64_t  trailing = u128 % P10_UINT64;
      rc = print_u128_u(leading);
      rc += printf("%." TO_STRING(E10_UINT64) PRIu64, trailing);
    }
  else
    {
      uint64_t u64 = u128;
      rc = printf("%" PRIu64, u64);
    }
  return rc;
}

static int fprint_u128_u(FILE *outfile, uint128_t  u128)
{
  int rc;
  if (u128 > UINT64_MAX)
    {
      uint128_t leading  = u128 / P10_UINT64;
      uint64_t  trailing = u128 % P10_UINT64;
      rc = fprint_u128_u(outfile, leading);
      rc += fprintf(outfile, "%." TO_STRING(E10_UINT64) PRIu64, trailing);
    }
  else
    {
      uint64_t u64 = u128;
      rc = fprintf(outfile, "%" PRIu64, u64);
    }
  return rc;
}



static int sfprint_u128_u(char str[], const unsigned int n, uint128_t  u128)
{
  int rc;
  if (u128 > UINT64_MAX)
    {
      uint128_t leading  = u128 / P10_UINT64;
      uint64_t  trailing = u128 % P10_UINT64;
      rc = sprint_u128_u(str, leading);

      if( rc >= 0 && rc < n )
      	rc += snprintf(str + rc, n-rc, "%." TO_STRING(E10_UINT64) PRIu64, trailing);
      else
	{
	  fprintf(stderr, "Cannot print such a long string....\n");
	  assert(0);
	}
    }
  else
    {
      uint64_t u64 = u128;
      rc = snprintf(str, n,  "%" PRIu64, u64);
    }
  return rc;
}

#else

int sprintf64(char str[], const unsigned int n, ui u)
{
  int rc;
  
  rc = snprintf( str, n, "%llu", u);

  assert( rc >= 0 );

  assert( rc < n);

  return rc;
}

#endif

int getKey( char str[], const unsigned int n, ui u )
{

  int rc;

#ifdef _U128
  rc = sfprint_u128_u(str, n, u);
#else
  rc = sprintf64(str, n, u);
#endif  
  
  return rc;

}


//uc *tempv;

static void iter(const char *key, const char *value, struct hashAux *obj)
{
  fprintf(obj->outfile, "%d\t%s\t%s\t%f\n", obj->index, key, value, obj->value1);
}


ui getRand( unsigned int  n )
{

  ui k = 1;

  unsigned int i = 0;

  double r = 0;

  for(i = 0; i < n; ++i)
    {
      r = (double)rand() / (RAND_MAX );

      if( r > 0.5 )
	k = k | ( (ui)1<<i);
	
    }

  return k;
}


unsigned int mypopcount(ui v)
{
  unsigned int c = 0;
  
  v = v - ((v >> 1) & (T)~(T)0/3);                           // temp
  v = (v & (T)~(T)0/15*3) + ( (v >> 2) & (T)~(T)0/15*3);      // temp
  v = (v + (v >> 4)) & (T)~(T)0/255*15;                      // temp
  c = (T)(v * ((T)~(T)0/255)) >> (sizeof(T) - 1) * CHAR_BIT; // count
  return c;

}


int binary_search(double A[], double key, int imin, int imax)
{
  int imid;

  
  if( key > A[imax] )
    assert( key <= A[imax]);

  if( key == A[imax]) 
    return imax;

  if(key <= A[imin])
    return imin;

  // continue searching while [imin,imax] is not empty
  while ( imax - imin > 1)
    {
      // calculate the midpoint for roughly equal partition

      imid = (imin+imax)/2;
      
      /* fprintf(stderr, "A mid: %f, min: %f, max: %f, imid: %d, imin: %d, imax: %d\n", A[imid], A[imin], A[imax], imid, imin, imax); */

      if(A[imid] == key)
        // key found at index imid
        return imid;
 
      // determine which subarray to search
      else if (A[imid] < key)
        // change min index to search upper subarray
        imin = imid ;
      else if( key < A[imid] )        
        // change max index to search lower subarray
        imax = imid;

      else
	{
	  fprintf(stderr, "Shouldn't be here... binary search... \n");
	  assert(0);
	}
    }

  return imax;

}


void vectorToInts ( ui m, ui *a, ui *b, ui n)
{
  
  ui i = 0, r = m;

  uc d;

  *a = *b = 0;

  if( m == 0 ) 
    {
      return;
    }    
  
  while( r > 0 )
    {
      
      d = r % 3;

      r = r / 3;
      
      if( d == 1 )
	
#ifdef _U128
	(*a) = (*a) | ((ui)1 << i); //base2[i];
#else
        (*a) = (*a) | (1LL << i); //base2[i];
#endif
        else if( d == 2)
#ifdef _U128
	  (*b) = (*b) | ((ui)1 << i); //base2[i];
#else
          (*b) = (*b) | (1LL << i); //base2[i];
#endif

	  ++i;
    }
  
  assert( i < n * n + 1);
  
}


ui multiplyVectorMatrix( ui v, ui m, ui n)
{
  ui vv, i, k = 0, r, r1 = 0, r2 = 0, r0= 0, f = 0, r1a = 0, mm = 0;

  for(i = 0; i <n; ++i)
    {
            
      vv = v << k;

      r = vv & m;
      
      r1 = mypopcount(r); //__builtin_popcountll(r); //popcount[r];

      r = vv & (~m);

      r2 = mypopcount(r); //popcount[r];

      r0 = (r1 > r2);

      r0 = r0 << i;

      f = f | r0;

      k += n;

    }

  return f;
}



ui multiplyVectorMatrix1a( ui v, ui m, ui n)
{

  ui vv, i, k = 0, r, r1 = 0, r2 = 0, r0= 0, f = 0;

  for(i = 0; i <n; ++i)
    {
      
      
      vv = v << k;

      r = vv & m;

      r1 = mypopcount(r); //popcount[r];

      r = vv & (~m);

      r2 = mypopcount(r); //popcount[r];

      r0 = (r1 > r2);

      r0 = r0 << i;

      f = f | r0;

      k += n;

    }

  return f;
}


/* ui getPopDif( ui k ) */
/* { */
/*   ui l = (ui)1;   */
/* } */



ui multiplyVectorMatrix2( ui v, ui a, ui b, ui n)
{
  ui vv, i, k = 0, r, r1 = 0, r2 = 0, r0= 0, f = 0;

  for(i = 0; i <n; ++i)
    {
      
      vv = v << k;

      r = vv & a;

      r1 = mypopcount(r); //popcount[r]

      r = vv & (b);

      //printf("i: %llu, r2: %llu, vv: %llu\n", i, r, vv);

      r2 = mypopcount(r); //popcount[r];

      r0 = (r1 > r2);

      r0 = r0 << i;

      f = f | r0;

      k += n;

      //printf("f: %llu, r: %llu, r1: %llu, popcount: %llu, half: %llu\n", f, r, r1, popcount[r], n/2);
      
    }

  return f;
}

ui getDif( ui a, ui b)
{

  
  
  return mypopcount(a^b); //popcount[a^b];

}


void evolveNeutrality( ui v, ui N, ui m, ui g, ui opt, float mu, ui n, ui *pop, ui *popM,  ui *lastPositions, ui *noPeriod, ui *maxPeriod, ui totalVectors, ui *successes, ui d)
{
  ui i,j, p, mdim = n*n, a = 0, b = 0, stepsToEquilibrium = MAXSTEPS, inprevious, steps, period = 0, k;

  float r;

  /* initialization for all members of the population*/
  for( i = 0; i < N; ++i)
    {
      /* set each individual's matrix and vector */
      pop[i] = v;
      
      popM[i] = m;
      
    }
    
  for( i = 0; i < g; ++i)
    {

      for( p = 0; p < N; ++p)
	pop[p] = v;
      
      for( j = 0; j < totalVectors; ++j)
	{
	  successes[ j ] = 0;
	  //lastPositions[ j ] = MAXSTEPS;
	}

      noPeriod[i] = 0;

      maxPeriod[i] = 0;
      
      /* mutation */
      for( p = 0; p < N; ++p)
	{

	  
	  
	  r = rand()/(1. + RAND_MAX);

	  /* mutate according to a mutation probability mu 
	     mutation rate should equal the total mutation rate for the 
	     network 
	  */
	  
	  //printf("matrix prin: %llu , r: %f", popM[p], r); 

	  if( r < mu)
	    {

	      
	      
	      /* choose a random location on the matrix to mutate */
	      k = rand() % mdim;

	      /* effectively mutation causes a toggle between 0 and 1 */
#ifdef _U128
popM[p] ^= ( (ui)1 << k ); //base2[k];
#else
popM[p] ^= ( 1LL << k ); //base2[k];
#endif

}


	  //printf("  matrix meta: %llu\n", popM[p]); 

	  a = b = 0;

	  vectorToInts( popM[p], &a, &b, n);
	  
	  /* initialize the period for each individual */
	  period = 0;
	  
	  /* reset the lastPositions vector. This will 
	     help to detect the period 
	  */
	  for( j = 0; j < totalVectors; ++j)
	    lastPositions[j] = MAXSTEPS;
	  
	  lastPositions[ pop[ p ] ] = 0;
	  
	  inprevious = pop[ p ];
	  
	  /* then do the development */
	  for(steps = 0; steps < MAXSTEPS; ++steps)
	    {
	      
	      //printf("input: %llu, matrix: %llu, output: ", pop[p], popM[p]);
	      if( d == 3)
		pop[p] = multiplyVectorMatrix2(pop[p], a, b, n);
	      
	      else if( d == 2)
		pop[p] = multiplyVectorMatrix( pop[p], popM[p], n);
	      
	      //printf("%llu\n", pop[p]);
	      
	      if( pop[p] == inprevious )
		{
		  stepsToEquilibrium = steps;
		  period = 1;
		  break;
		}
	      
	      if( lastPositions[ pop[p] ] < MAXSTEPS )
		
		{
		  period = steps + 1 - lastPositions[ pop[p] ];
		  stepsToEquilibrium = steps;
		  break;
		}
	      
	      
	      lastPositions[ pop[p] ] = steps + 1;
	      
	      inprevious = pop[p];
	      
	    }
	  
	  if( period == 0 )
	    noPeriod[i]++;

	  
	  if( stepsToEquilibrium < MAXSTEPS && period == 1)
	    successes[ pop[p] ]++;
	  else if( period > 1 )
	    {
	      p--;
	      continue;
	    }
	  	  
	}



      	  
      if( period > maxPeriod[i])
	maxPeriod[i] = period;	  

#ifdef _U128
      printf("------ Generation : ");
      print_u128_u(i);
      printf(" ------\n");
#else
      printf("------ Generation : %llu -----\n", i);

       for( j = 0; j < totalVectors; ++j) 
	 printf("Successes of vector  %llu: %llu\n", j, successes[j]); 
#endif
      
    }
  
}



void evolveDriftNeutrality( ui v, ui N, ui m, ui g, ui opt, float mu, ui n, ui *pop, ui *popM,  ui *lastPositions, ui *noPeriod, ui *maxPeriod, ui totalVectors, ui *successes, ui d, ui *tempM)
{
  ui i, ii, j, p, mdim = n*n, a = 0, b = 0, stepsToEquilibrium = MAXSTEPS, inprevious, steps, period = 0, k;

  float r;


  /* initialization for all members of the population*/
  for( i = 0; i < N; ++i)
    {
      /* set each individual's matrix and vector */
      pop[i] = v;
      
      popM[i] = m;
      
    }
    
  for( i = 0; i < g; ++i)
    {

      for( p = 0; p < N; ++p)
	pop[p] = v;
      
      for( j = 0; j < totalVectors; ++j)
	{
	  successes[ j ] = 0;
	  //lastPositions[ j ] = MAXSTEPS;
	}

      noPeriod[i] = 0;

      maxPeriod[i] = 0;
      
      /* mutation */
      for( p = 0; p < N; ++p)
	{
	  
	  r = rand()/(1. + RAND_MAX);

	  /* mutate according to a mutation probability mu 
	     mutation rate should equal the total mutation rate for the 
	     network 
	  */
	  
	  //printf("matrix prin: %llu , r: %f", popM[p], r); 

	  if( r < mu)
	    {
	      
	      /* choose a random location on the matrix to mutate */
	      k = rand() % mdim;

	      /* effectively mutation causes a toggle between 0 and 1 */
#ifdef _U128
	      popM[p] ^= ( (ui)1 << k ); //base2[k];
#else
	      popM[p] ^= ( 1LL << k ); //base2[k];
#endif
	    }


	  //printf("  matrix meta: %llu\n", popM[p]); 

	  a = b = 0;

	  vectorToInts( popM[p], &a, &b, n);
	  
	  /* initialize the period for each individual */
	  period = 0;
	  
	  /* reset the lastPositions vector. This will 
	     help to detect the period 
	  */
	  for( j = 0; j < totalVectors; ++j)
	    lastPositions[j] = MAXSTEPS;
	  
	  lastPositions[ pop[ p ] ] = 0;
	  
	  inprevious = pop[ p ];
	  
	  /* then do the development */
	  for(steps = 0; steps < MAXSTEPS; ++steps)
	    {
	      
	      //printf("input: %llu, matrix: %llu, output: ", pop[p], popM[p]);
	      if( d == 3)
		pop[p] = multiplyVectorMatrix2(pop[p], a, b, n );
	      
	      else if( d == 2)
		pop[p] = multiplyVectorMatrix( pop[p], popM[p], n);
	      
	      //printf("%llu\n", pop[p]);
	      
	      if( pop[p] == inprevious )
		{
		  stepsToEquilibrium = steps;
		  period = 1;
		  break;
		}
	      
	      if( lastPositions[ pop[p] ] < MAXSTEPS )
		
		{
		  period = steps + 1 - lastPositions[ pop[p] ];
		  stepsToEquilibrium = steps;
		  break;
		}
	      
	      
	      lastPositions[ pop[p] ] = steps + 1;
	      
	      inprevious = pop[p];
	      
	    }
	  
	  if( period == 0 )
	    noPeriod[i]++;

	  /* if period == 1 then the 
	     genotype pop[p] i.e. of individual p
	     has succeeded
	  */
	  if( stepsToEquilibrium < MAXSTEPS && period == 1)
	    successes[ pop[p] ]++;
	  else if( period > 1 )
	    {
	      /* we can choose an individual based on its
		 distance from the optimum
		 even if it has a periodic equilibrium
	      */

	      p--;
	      continue;
	    }
	  	  
	} // for each member of the population

      /* choose parents for the next generation */
      for( p = 0; p < N; ++p )
	{
	  
	  ii = rand() % N;

	  tempM[p] = popM[ii];
	  
	}

      popM = tempM;
      	  
      if( period > maxPeriod[i])
	maxPeriod[i] = period;	  

#ifdef _U128
#else
      printf("------ Drift Generation : %llu -----\n", i);
      
      for( j = 0; j < totalVectors; ++j)
	printf("Successes of vector  %llu: %llu\n", j, successes[j]);
#endif
      
    }
  
}


void evolveSelection( ui v, ui N, ui m, ui g, ui opt, float mu, ui n, ui *pop, ui *popM,  ui *lastPositions, ui *noPeriod, ui *maxPeriod, ui totalVectors, ui *successes, ui d, ui totalMatrices, ui *tempM)
{
  
  ui i, j, ii, jj, p, mdim = n*n, a = 0, b = 0, stepsToEquilibrium = MAXSTEPS, inprevious, steps, period = 0, k, similarity[N], totalSimilarity = 0, maxSimilarity = 0, matrixSum = 0;
  
  float r;
  
  maxSimilarity = n * N;
  
  ui choiceVector[ maxSimilarity ];
  
  /* initialization for all members of the population
     there are N individuals
  */
  for( i = 0; i < N; ++i)
    {
      /* set each individual's matrix and vector to the same initial value */

      /* pop[i] is the "expression" vector for each individual */
      pop[i] = v;
      
      /* popM is the regulatory matrix for each individual */
      popM[i] = m;
      
    }
  

  /* do the following steps for every one generation. There are 'g' generations */
  for( i = 0; i < g; ++i)
    {
      
      /* initialize matrixSum to 0 */
      matrixSum = 0;
      
      /* reset the vector for each individual */
      for( p = 0; p < N; ++p)
	pop[p] = v;

      
      for( j = 0; j < totalVectors; ++j)
	{
	  successes[ j ] = 0;
	  //lastPositions[ j ] = MAXSTEPS;
	}
      
      
      /* noPeriod[i] = 0; */
      
      /* maxPeriod[i] = 0; */
      
      /* initialize total similarity */
      totalSimilarity = ii = 0;      
      
      for( p = 0; p < N; ++p)
	{
	  
	  /* get a random number */
	  r = rand()/(1. + RAND_MAX);

	  /* mutate according to a mutation probability mu 
	   mutation rate should equal the total mutation rate for the 
	   network 
	  */

	  //printf("matrix prin: %llu , r: %f", popM[p], r); 
	  
	  /* mutation: if random number
	   is smaller than then mutation probability then mutate. Thus,
	  we mutate with probability mu */
	  if( r < mu)
	    {
	      
	      /* choose a random location on the matrix to mutate */
	      k = rand() % mdim;
	      
	      /* effectively mutation causes a toggle between 0 and 1 */
#ifdef _U128
	      popM[p] ^= (ui)( (ui)1 << k); //base2[k];
#else
	      popM[p] ^= (ui)(1LL << k); //base2[k];
#endif
	    }

	  //printf("  matrix meta: %llu\n", popM[p]); 

	  a = b = 0;
	  if( d == 3)
	    vectorToInts( popM[p], &a, &b, n);
	  
	  /* initialize the period for each individual */
	  period = 0;
	  
	  /* reset the lastPositions vector. This will 
	     help to detect the period 
	  */
	  for( j = 0; j < totalVectors; ++j)
	    lastPositions[j] = MAXSTEPS;
	  
	  lastPositions[ pop[ p ] ] = 0;

	  inprevious = pop[ p ];
	  
	  /* then do the development */
	  for(steps = 0; steps < MAXSTEPS; ++steps)
	    {
	      
	      //printf("input: %llu, matrix: %llu, output: ", pop[p], popM[p]);
	      
	      /* allowing 3 expression levels */
	      if( d == 3)
		pop[p] = multiplyVectorMatrix2(pop[p], a, b, n );
	      
	      /* allowing 2 expression levels */
	      else if( d == 2)
		pop[p] = multiplyVectorMatrix( pop[p], popM[p], n);
	      
	      //printf("%llu\n", pop[p]);
	      
	      /* equilibrium POINT is when the current state is 
		 the same as the previous
	      */
	      if( pop[p] == inprevious )
		{
		  stepsToEquilibrium = steps;
		  period = 1;
		  break;
		}
	      
	      /* equilibrium PERIOD */
	      if( lastPositions[ pop[p] ] < MAXSTEPS )
		
		{
		  period = steps + 1 - lastPositions[ pop[p] ];
		  stepsToEquilibrium = steps;
		  break;
		}
	      	      
	      lastPositions[ pop[p] ] = steps + 1;
	      
	      inprevious = pop[p];
	      
	    }

	  /* print_u128_u(steps); */
	  /* printf("\n"); */

	  
	  
	  /* if( period == 0 ) */
	  /*   noPeriod[i]++; */

	  if( stepsToEquilibrium < MAXSTEPS && period == 1)
	    {
	      matrixSum += stepsToEquilibrium;

	      //successes[ pop[p] ]++;
	    }
	  else if( period > 1 )
	    {

	      /* here we ignore individuals who have 
		 cyclic equilibria 
	      */
	      p--;
	      continue;
	    }

	  /* 
	     similarity from the optimum is the 
	     common set bits between the optimum and the 
	     equilibrium expression state
	  */
	     
	  similarity[p] = ( n - mypopcount( opt^pop[p] ) ); //popcount[ opt^pop[p] ];
	  

	  /* similarity of the total population; like the average fitness */
	  totalSimilarity += similarity[p];

	  
	  /* this is a trick to choose the parent more efficiently 
	   each individual has a weight equivalent to its similarity 
	  from the optimum */
	  for( jj = 0; jj < similarity[p]; ++jj)
	    {
	      choiceVector[jj+ii] = p;
	    }

	  ii += similarity[p];
	  	  
	}

#ifdef _U128

      if( i == 0)
	{
	  printf("RelativeSimilarity\tMeanSteps\n");
	}
      
      printf("%f\t%f\n", (double)totalSimilarity/(double)maxSimilarity, (double)matrixSum/N);
#else
      if( i % samplestep == 0 )
	{
	  printf("SelectionTotalSimilarity:\t%llu\tmax:\t%llu\tperc:\t%f\tsteps:\t%llu\n", totalSimilarity, maxSimilarity, (double)totalSimilarity/maxSimilarity, matrixSum);
	}
      
#endif

      /* selection */
      for( p = 0; p < N; ++p)
	{
	  
	  ii = rand() % totalSimilarity;
	  
	  ui pChoice = choiceVector[ii];
	  
	  //printf("pchoice: %llu, %llu, %llu\n", pChoice, ii, totalSimilarity);

	  tempM[p] = popM[ pChoice ];
	}

      popM = tempM;
      
            	  
      /* if( period > maxPeriod[i]) */
      /* 	maxPeriod[i] = period;	   */

#ifdef _U128
      
#else
      printf("------ Generation : %llu -----\n", i);

      for( j = 0; j < totalVectors; ++j)
	printf("Successes of vector  %llu: %llu\n", j, successes[j]);
#endif
      
    }
  
}

double getFitness(double dist, unsigned int fitnesstype)
{

  double fitness = 0;

  if( fitnesstype == 0 )
    fitness = 1.;
  
  else if( fitnesstype == 1)
    {
      fitness = 1. - dist;
	
    }

  else if( fitnesstype == 100 )
    fitness =  exp( -(dist*dist)/2. );

  else if( fitnesstype > 1 && fitnesstype < 100)
    fitness =  pow(1. - dist, fitnesstype);

  else if ( fitnesstype > 100 && fitnesstype < 200)
    fitness =  1. - pow(dist, fitnesstype - 100);

  assert(fitness >= 0 );

  return fitness;
}

void evolveSelection2( ui v, ui N, ui m, ui g, ui opt, float mu, ui n, ui *pop, ui *popM,  ui *lastPositions, ui *noPeriod, ui *maxPeriod, ui totalVectors, ui *successes, ui d, ui totalMatrices, ui *tempM, ui cyclic, unsigned int changeOptimum, unsigned int changeOptGen, unsigned int mutStop, unsigned int fitnessType )
{
  
  ui a = 0, b = 0;

  unsigned int cyclicDevelopment = 0, i, j, matrixSum = 0, stepsToEquilibrium = MAXSTEPS, p = 0, numberOfMutations = 0; 
  
  double  popCountA = 0., popCountB = 0., popCountVarA = 0., popCountVarB = 0.;
  
  ui mdim = n*n, inprevious, steps, k, similarity[N], totalSimilarity = 0, maxSimilarity = 0, devStates[MAXSTEPS], cycle;

  float inMu = mu;

  /** I'd like to store the number of each present haplotype */
  unsigned int K = 1024;
  char key[K];

  /* 
     This section is temporary
       
  */
    
  /* 
     this version of the function uses 
     a double weight vector to choose individuals 
     for the next  generation
  */
  

  unsigned int period = 0;
  
  /* for( i = 0; i < N; ++i) */
  /*   similarity[p] = 1; */
  
  float r;
  
  maxSimilarity = n * N;

  //ui choiceVector[ maxSimilarity ];
  
  double weightsVector[N], fitnessVector[N], tempDist = 0, sumWeights = 0;
    
  /* initialization for all members of the population
     there are N individuals
  */
  for( i = 0; i < N; ++i)
    {
      /* set each individual's matrix and vector to the same initial value */

      /* pop[i] is the "expression" vector for each individual */
      pop[i] = v;
      
      /* popM is the regulatory matrix for each individual */
      popM[i] = m;

      //fprintf(stderr, "i: %d, pop: %d, popM: %d\n", i, pop[i], popM[i]);

    }

  int rc = getKey(key, K, m);
  
  StrMap *sm;
  
  struct hashAux *aux = malloc(sizeof(struct hashAux));
  aux->outfile = genoutFile;
  aux->index = 0;

  /* do the following steps for every one generation. There are 'g' generations */
  for( i = 0; i < g; ++i)
    {

      numberOfMutations = 0;

      aux->index = i;
      sm = sm_new(N);

      if( sm == NULL)
	{
	  fprintf(stderr, "Coulnd;t create a hash table of size %u\n", N);
	  assert(sm != NULL);
	}

      if(changeOptimum && (i % changeOptGen) == 0)
	{
	  opt = getRand((unsigned int)n);
	  fprintf(stderr, "opt: %llu\n", opt);
	}

      if( mutStop > 0 && i == mutStop)
	mu = 0;


      popCountVarB = popCountVarA = popCountA = popCountB = 0.;

      cyclicDevelopment = 0;
      
      /* initialize matrixSum to 0 */
      matrixSum = 0;
      
      /* reset the vector for each individual */
      for( p = 0; p < N; ++p)
	pop[p] = v;

      
      for( j = 0; j < totalVectors; ++j)
	{
	  successes[ j ] = 0;
	  //lastPositions[ j ] = MAXSTEPS;
	}
            
      /* noPeriod[i] = 0; */
      
      /* maxPeriod[i] = 0; */
      
      /* initialize total similarity */
      //totalSimilarity = ii = 0;      

      sumWeights = 0;
      

      /* p: every person
	 N: all persons in the population
      */
      for( p = 0; p < N; ++p)
	{

	  unsigned int personClass = 0;

	  /* initialization to 1. which is the maximum distance */
	  tempDist = 1.;

	  /* initialization to 0, this is a scaled version of the fitness */
	  fitnessVector[p] = 0;
	  
	  /* cumulative vector */
	  weightsVector[p] = 0;
	  
	  
	  
	  /* get a random number */
	  r = rand()/(1. + RAND_MAX);

	  /* mutate according to a mutation probability mu 
	     mutation rate should equal the total mutation rate for the 
	     network 
	  */
	  
	  //printf("matrix prin: %llu , r: %f", popM[p], r); 
	  
	  /* mutation: if random number
	     is smaller than then mutation probability then mutate. Thus,
	     we mutate with probability mu 
	  */
	  if( mu>0 && r < mu)
	    {
	      
	      numberOfMutations++;

	      /* choose a random location on the matrix to mutate */
	      k = rand() % mdim;
	      
	      /* effectively mutation causes a toggle between 0 and 1 */
#ifdef _U128
	      popM[p] ^= (ui)( (ui)1 << k); //base2[k];
#else
	      popM[p] ^= (ui)(1LL << k); //base2[k];
#endif
	    }

	  //printf("  matrix meta: %llu\n", popM[p]); 

	  a = b = 0;

	  if( d == 3)
	    {
	      
	      vectorToInts( popM[p], &a, &b, n);
	      
	    }
	  else if( d == 2 )
	    {

	    }
	  
	  /* initialize the period for each individual */
	  period = 0;
	  
	  /* reset the lastPositions vector. This will 
	     help to detect the period 
	  */
	  for( j = 0; j < totalVectors; ++j)
	    lastPositions[j] = MAXSTEPS;
	  
	  fprintf( stderr, "p: %u, pop[p]: %u, totalVectors: %u\n", p, pop[p], totalVectors);
	  lastPositions[ pop[ p ] ] = 0;

	  inprevious = pop[ p ];
	  
	  /* then do the development */
	  for(steps = 0; steps < MAXSTEPS; ++steps)
	    {
	      
	      //printf("input: %llu, matrix: %llu, output: ", pop[p], popM[p]);
	      
	      /* allowing 3 expression levels */
	      if( d == 3)
		pop[p] = multiplyVectorMatrix2(pop[p], a, b, n );
	      
	      /* allowing 2 expression levels */
	      else if( d == 2)
		pop[p] = multiplyVectorMatrix( pop[p], popM[p], n);

	      devStates[steps] = pop[p];
	      
	      //printf("%llu\n", pop[p]);
	      
	      /* equilibrium POINT is when the current state is 
		 the same as the previous
	      */
	      if( pop[p] == inprevious )
		{
		  stepsToEquilibrium = steps;
		  period = 1;
		  break;
		}
	      
	      /* equilibrium PERIOD */
	      if( lastPositions[ pop[p] ] < MAXSTEPS )
		{
		  period = steps + 1 - lastPositions[ pop[p] ];
		  
		  stepsToEquilibrium = steps;
		  
		  break;
		}
	      
	      lastPositions[ pop[p] ] = steps + 1;
	      
	      inprevious = pop[p];
	      
	    } // steps
	  
	  rc = getKey(key, K, popM[p]);

	  if( sm_exists(sm, key) == 0 )
	    {
	      int sm_i = sm_put(sm, key, "1");

	      assert( sm_i == 1);
	    }
	  else
	    {
	      char buf[1024];
	      buf[0] = 0;
	      int sm_i = sm_get(sm, key, buf, sizeof(buf));
	      //fprintf(genoutFile, "%u\t%s\t%s\n", i, key, buf);
	      int bufInt = atoi(buf);
	      bufInt++;
	      sprintf( buf, "%d", bufInt);
	      sm_i = sm_put(sm, key, buf);
	      assert( sm_i == 1);
	      
	    }

	  
	  /* print_u128_u(steps); */
	  /* printf("\n"); */
	  	  
	  /* if( period == 0 ) */
	  /*   noPeriod[i]++; */

	  //fprintf(stderr, "period\t%u\n", period);

	  if( stepsToEquilibrium < MAXSTEPS && period == 1)
	    {
	      matrixSum += stepsToEquilibrium;

	      /* average number of differences between the optimum and the currect vector of the popualtions */
	      tempDist = (double)mypopcount( opt^pop[p] ) / (double)n; 

	      personClass = 1;
	      
	      //successes[ pop[p] ]++;
	    }
	  else if( period > 1 )
	    {

	      if( ! cyclic )
		{
		  /* 
		     here we ignore individuals who have 
		     cyclic equilibria 
		  */
		  fprintf(stderr, "1. no constant equilibrium achieved for %u in generation %u. vector %u(final), %u(initial), matrix %u ....\n", p, i, pop[p], v,  popM[p]);

		  tempDist = 1.;
		  
		  personClass = 2;
		}
	      else
		{
		  /* here we accept cyclic equilibria 
		     distance is measured based on every 
		     step of the equilibrium cycle
		  */
		  tempDist = 0.;

		  personClass = 3;
	      
		  for( cycle = 0; cycle < period; ++cycle)
		    {
		      tempDist += mypopcount(opt^devStates[ stepsToEquilibrium - cycle ] );
		    }
		  /* average distance when cycles > 1 */
		  tempDist /= (double)( n * period );
		  
		  cyclicDevelopment++;
		  
		}
	    }
	  else
	    {
	      /* here no equilibrium is achieved at all */

	      fprintf(stderr, "2. no equilibrium achieved....\n");

	      personClass = 4;
	      
	      tempDist = 1.;
	      	      
	    }
	
	  
	  if( d == 3)
	    {
	      
	      popCountA += (double)mypopcount(a);

	      popCountB += (double)mypopcount(b);

	      popCountVarA += ( (double)mypopcount(a) * mypopcount(a) );

	      popCountVarB += ( (double)mypopcount(b) * mypopcount(b) );
	      
	    }
	  else if( d == 2 )
	    {
	      unsigned int tmpPopCount = mypopcount( popM[p]);

	      popCountA += tmpPopCount;
	      
	      popCountVarA += (tmpPopCount * tmpPopCount);

	      //fprintf(stderr, "popcount: %f\t%f\t%u\t%u\t%f\t%f\t%f\n", popCountA, popCountVarA, tmpPopCount, p+1, popCountVarA/(p+1), (popCountA*popCountA)/((p+1)*(p+1)) , popCountVarA/(p+1) - (popCountA*popCountA)/((p+1)*(p+1))  );
	    }

	  
	  /* calculate the fitness for each person p, maximum is 1, minimum is 0 */
	  
	  fitnessVector[p] = getFitness(tempDist, fitnessType);

	  //fprintf(stderr, "fitness for %u in generation %u is %f (distance: %f)\n", p, i, fitnessVector[p], tempDist);

	  sumWeights += fitnessVector[p];

	  if( i > 593 && i < 600)
	    fprintf(stderr, "sumweights : %f, generation: %u, p: %u, fitness: %f, tempDist: %f, personClass: %u\n", sumWeights, i, p, fitnessVector[p], tempDist, personClass);

	  weightsVector[p] = sumWeights;

	  if( p > 0 )
	    assert( weightsVector[p] >= weightsVector[p-1]);
	  	  	  
	  /* 
	     similarity from the optimum is the 
	     common set bits between the optimum and the 
	     equilibrium expression state
	  */
	  
	  /* similarity of the total population; like the average fitness */
	  //totalSimilarity += similarity[p];
	  
	  /* this is a trick to choose the parent more efficiently 
	   each individual has a weight equivalent to its similarity 
	  from the optimum */
	  /* for( jj = 0; jj < similarity[p]; ++jj) */
	  /*   { */
	  /*     choiceVector[jj+ii] = p; */
	  /*   } */

	  /* ii += similarity[p]; */
	  	  
	}

      aux->index = i;

      aux->value1 = sumWeights;

      sm_enum(sm, iter, aux);

      fprintf(stderr, "%f\t%f\t%f\t%f\n", popCountVarA, popCountVarB, popCountA*popCountA, popCountB*popCountB);

      popCountA /= N;

      popCountB /= N;

      popCountVarA /= N;

      popCountVarB /= N;

      fprintf(stderr, "%f\t%f\t%f\t%f\n", popCountVarA, popCountVarB, popCountA*popCountA, popCountB*popCountB);

      popCountVarA -= (popCountA*popCountA);

      popCountVarB -= (popCountB*popCountB);
      

      if(i == 0 )
	{
	  fprintf(selFile, "generation\tmeanA\tmeanB\tvarA\tvarB\tfitness\tsteps\n");
	}

#ifdef _U128
      fprintf(selFile, "%u\t%f\t%f\t%f\t%f\t%f\t%u\n", i, popCountA, popCountB, popCountVarA, popCountVarB, sumWeights, matrixSum);
#else
      fprintf(selFile, "%u\t%f\t%f\t%f\t%f\t%f\t%u\n", i, popCountA, popCountB, popCountVarA, popCountVarB, sumWeights, matrixSum);
#endif
      

#ifdef _U128

      if( i == 0)
	{
	  printf("RelativeSimilarity\tMeanSteps\n");
	}
      
      printf("%f\t", sumWeights);
      printf("%u\t", matrixSum);
      printf("\n");
#else
      if( i % samplestep == 0 )
	{
	  
	  if( i == 0 )
	    printf("RelativeSimilarity\tMeanSteps\n");

	  printf("%f\t%u\n", sumWeights, matrixSum);
	  
	}
      
#endif


      fprintf(choicefile, "%u", i);
      
      /* selection */
      for( p = 0; p < N; ++p)
	{
	  
	  /* weightsVector[N-1] : sum of all fitnesses 
	     randValue is a random number between 0 and sum of fitnesses 
	  */
	  double randValue = ( (double)rand()/(RAND_MAX+1.) ) * ( weightsVector[N-1] ); 
	  
	  ui pChoice = binary_search( weightsVector, randValue, 0, N-1);


	  unsigned int ii;
	  if( i > 594 && i < 600 && p == 0 )
	    for( ii = 0; ii < N; ++ii)
	      fprintf(stderr, "randValue\t%u\t%u\t%f\n", i, ii, weightsVector[ii]);

	  if(i >= 593 && i <=597)
	    fprintf(stderr, "randValue: %u\t\t%f\tmin: %f\tmax: %f\tchoice: %u\n", i,randValue, weightsVector[0], weightsVector[N-1], pChoice);  
	  
	  /* fprintf(stderr, "pchoice: %d\n", pChoice);  */
	  assert(pChoice < N);
	  
	  tempM[p] = popM[ pChoice ];

	  fprintf(choicefile, "\t%u", pChoice);
	  
	  /* if( pChoice > totalVectors - 1) */
	  /*   { */
	  /*     fprintf(stderr, "pChoice should be smaller than totalVectors: %u < %u\n", pChoice, totalVectors ); */
	      
	  /*     assert( pChoice < totalVectors ); */

	  /*   } */

	  successes[ pop[pChoice] ]++;
	  
	}

      fprintf(choicefile, "\n");
      

      popM = tempM;
                  	  
      /* if( period > maxPeriod[i]) */
      /* 	zmaxPeriod[i] = period;	   */

#ifdef _U128
      
      /* fprintf(selFile, "------ Generation : %u -----\n", i); */
      
      /* for( j = 0; j < totalVectors; ++j) */
      /* 	{ */
      /* 	  fprintf(selFile, "Successes of vector  "); */
      /* 	  fprint_u128_u(selFile, j); */
      /* 	  fprintf(selFile, ": "); */
      /* 	  fprint_u128_u(selFile, successes[j]); */
      /* 	  fprintf(selFile, "\n"); */
      /* 	} */

#else
      
      /* fprintf(selFile, "------ Generation : %u -----\n", i); */
      
      /* for( j = 0; j < totalVectors; ++j) */
      /* 	fprintf(selFile, "Successes of vector  %u: %llu\n", j, successes[j]); */
#endif

      sm_delete(sm);

      fprintf(mutfile, "%u\t%u\n", i, numberOfMutations); 
      
    }


  free(aux);

    
}


ui development(ui input, ui matrix, ui states, ui n, ui *lastPositions,  ui *period, ui *stepsToEquilibrium)
{
  ui a = 0, b = 0, inprevious = input, in = input, steps;

  *stepsToEquilibrium = MAXSTEPS;

  if( states == 3 )
    vectorToInts(matrix, &a, &b, n);

  lastPositions[ in ] = 0;


  for(steps = 0; steps  < MAXSTEPS; ++steps)
    {
      
      if( states == 3 )
	in = multiplyVectorMatrix2( in, a, b, n );
      
      else if( states == 2)
	in = multiplyVectorMatrix( in, matrix, n );
      
      else
	{
	  fprintf(stderr, "Only states 2 or 3 are supported for the matrix.... will exit now\n");
	  assert(0);
	}

      if( in ==  inprevious )
	{
	  *stepsToEquilibrium = steps;
	  *period = 1;
	  break;
	}

      if( lastPositions[ in ] < MAXSTEPS )
	{
	  *period = steps + 1 - lastPositions[ in ];
	  *stepsToEquilibrium = steps;
	  break;
	}

      lastPositions[ in ] = steps + 1;

      inprevious = in;

    }
  
  return in;

}


void runExhaustiveDevelopment(ui n, ui d, unsigned long long  totalMatrices, unsigned long long totalVectors, ui *lastPositions, ui *matrixSum, ui *matrixMax, ui *maxPeriod, ui *noPeriod, ui *successes)
{

  ui a = 0, b = 0,  m, i, j, period, in, stepsToEquilibrium = MAXSTEPS;
        
  FILE *resfile = fopen("binaryNetworkOutput.csv", "w");
      
  FILE *matresfile = fopen("matrixOutput.csv", "w");
  
  fprintf(resfile, "input\tmatrix\teqsteps\teqvector\tperiod\n");
      
  for( m = 0; m < totalMatrices; ++m )
    {
	  
      a = b = 0;
	  
      vectorToInts( m, &a, &b, n);
	  
      *matrixSum = 0; 
	  
      *matrixMax = 0;
	  
      for( i = 0; i < totalVectors; ++i)
	{
	      
	  period = 0;
	      
	  /* reset the lastPositions vector. This will 
	     help to detect the period 
	  */
	  for( j = 0; j < totalVectors; ++j)
	    lastPositions[j] = MAXSTEPS;
	      
	  in = i;
	      
	  in = development( in, m, d, n, lastPositions, &period, &stepsToEquilibrium);
	      
	  if( period == 0 )
	    (*noPeriod)++;
	      
	  if( period > *maxPeriod)
	    *maxPeriod = period;	  
	      
	  (*matrixSum) += stepsToEquilibrium;
	      
	  if( stepsToEquilibrium > *matrixMax )
	    *matrixMax = stepsToEquilibrium;					      
	  
#ifdef _U128
#else
	  fprintf( resfile, "%llu\t%llu\t%llu\t%llu\t%llu\n", i, m, stepsToEquilibrium, in, period);
#endif
	      
	  if( stepsToEquilibrium < MAXSTEPS && period == 1)
	    successes[ in ]++;
	      
	}
	
#ifdef _U128
#else
	fprintf( matresfile, "%llu\t%f\t%llu\n", m, (*matrixSum)/(float)totalVectors, *matrixMax);
#endif
	
    }
    
  fclose(resfile);
      
  fclose(matresfile);

  fclose(selFile);

  fclose(genoutFile);  

  
      
}


int main(int argc, char ** argv)
 {


   selFile = stdout;

   mutfile = stdout;
   
   choicefile = stdout;
   
   unsigned long long i, j, n = 1;

   char seloutFileName[1000], genoutFileName[1000], mutoutFileName[1000], choiceFileName[1000]; 
   
   seloutFileName[0] = genoutFileName[0] = mutoutFileName[0] = choiceFileName[0] = 0;
  
  ui in, m, d = 2, a = 0, b = 0,
    stepsToEquilibrium = MAXSTEPS, *lastPositions, period = 0, maxPeriod = 0, noPeriod = 0,
    developAll = 0, matrixSum, generations = 1000, populationSize = 1000, cyclic = 0, startVector = 0, startMatrix = 0, optimum = 0;

  ui totalVectors = 1, totalMatrices = 1;

  unsigned long long seed = 0;

  ui  *successes;

  double mutationRate = 0;

  unsigned int changeOptGen = 0, changeOptimum = 0, mutStop = 0, startMatrixSet = 0, startVectorSet = 0, optimumSet = 0, fitnessType = 1;

  samplestep = 1;

  
  /* read in the command line options */
  for(i = 1; i < argc; ++i)
    {

      if( strcmp(argv[i], "-sv") == 0 )
	{
	  startVector = strtoull(argv[++i], NULL, 10);
	  startVectorSet = 1;
	  continue;
	}

      if( strcmp( argv[i], "-sm") == 0 )
	{
	  startMatrixSet = 1;
	  startMatrix = strtoull( argv[++i], NULL, 10);
	  continue;
	}

      if( strcmp( argv[i], "-opt") == 0 )
	{
	  optimumSet = 1;
	  optimum = strtoull( argv[++i], NULL, 10);
	  continue;
	}

      /* how many genes exist */
      if( strcmp( argv[i], "-n") == 0)
	{
	  n = atoi( argv[++i] );
	  continue;
	}

      /* how many regulatory modes exist	
	 d = 2 : either positive or negative
	 d = 3 : positive, negative or does not affect
      */
      if( strcmp( argv[i], "-d") == 0)
	{
	  d = atoi( argv[++i] );
	  continue;
	}

      /* seed for random number generation */
      if( strcmp( argv[i], "-seed") == 0 )
	{
	  seed = atol( argv[++i]);
	  continue;
	}

      if( strcmp( argv[i], "-changeOptimum") == 0 )
	{
	  changeOptimum = 1;
	  changeOptGen = atoi( argv[++i]);
	  continue;
	}
      /* mutation rate */
      if( strcmp( argv[i], "-mu") == 0 )
	{
	  mutationRate = atof( argv[++i]);
	  continue;
	}

      /* how many generations evolution should take place */
      if( strcmp( argv[i], "-generations") == 0 )
	{
	  generations = atoi( argv[++i]); 
	  continue;
	}

      /* the size of the population */
      if( strcmp( argv[i], "-N") == 0)
	{
	  populationSize = atoi( argv[++i]);
	  continue;
	}

      /* sample step ??? */
      if( strcmp( argv[i], "-ss") == 0 )
	{
	  samplestep = atoi( argv[++i]);
	  continue;
	}

      if( strcmp( argv[i], "-mutStop") == 0 )
	{
	  mutStop = atoi(argv[++i]);
	  assert( mutStop >= 0 );
	  continue;
	}

      /* develop all possible combinations?? */
      if( strcmp( argv[i], "-developAll") == 0 )
	{
	  developAll = 1;
	  continue;
	}

      /* allow or not cyclic equilibria */
      if( strcmp(argv[i], "--cyclic") == 0 )
	{
	  cyclic = 1;
	  continue;
	}

      if( strcmp( argv[i], "-selout") == 0)
	{
	  strcpy(seloutFileName, argv[++i]);
	  continue;
	}

      if( strcmp( argv[i], "-mutout") == 0 )
	{
	  strcpy( mutoutFileName, argv[++i]);
	  continue;
	}

      
      if( strcmp( argv[i], "-choiceout") == 0 )
	{
	  strcpy( choiceFileName, argv[++i]);
	  continue;
	}


      if( strcmp( argv[i], "-genout") == 0 )
	{
	  strcpy(genoutFileName, argv[++i]);
	  continue;
	}

      if( strcmp( argv[i], "-fitnessType") == 0 )
	{
	  fitnessType = atoi(argv[++i]);
	  continue;
	}


      fprintf( stderr, "Argument %s is invalid.... exit NOW!\n", argv[i]);
      
      exit(-1);

    }

  if(optimumSet == 0 )
    optimum = rand();
  
  if(startVectorSet == 0 )
    fprintf(stderr, "You have not provided an initial vector\n");
    //startVector = rand()/RAND_MAX;
  
   

  if( seloutFileName[0] != 0 )
    selFile = fopen(seloutFileName, "w");
  
  if( genoutFileName[0] != 0 )
    genoutFile = fopen(genoutFileName, "w");

  if( mutoutFileName[0] != 0 )
    mutfile = fopen( mutoutFileName, "w");
  
  if( choiceFileName[0] != 0 )
    choicefile = fopen( choiceFileName, "w");


  /* seed using the timer */
  if( seed == 0 )
    {
      seed = time(NULL);
      
      fprintf( stderr, "\nWarning: timer used as seed %llu\n", seed);
      
      srand( seed );
    }
  else
    /* seed using the seed */
    srand( seed );
    
  assert( n > 1);

  assert( d == 2 || d == 3);


  /* the total number of vectors should be 2^n */
  for( i = 0; i < n; ++i)
#ifdef _U128
    totalVectors = totalVectors << 1;
#else
  totalVectors = totalVectors << 1LL;
#endif

    
  /* the total number of matrices will be d^(n*n) 
     n*n is the size of the interaction matrix, for n genes (n*n possible interactions)
     If each regulation has d possible states (-1, 1, for d=2), (-1, 0, 1 for d = 3)  then
     there are d^(n*n) possible matrices
  */
  for( i = 0; i < n*n; ++i)
    totalMatrices = totalMatrices * d; 

  
  /* get a random starting matrix if not set */
  if(startMatrixSet == 0 )
    startMatrix = (ui)(rand()/RAND_MAX * totalMatrices);


  /* get a random starting vector if not set */
  if( startVectorSet == 0 )
    startVector = (ui)( rand()/RAND_MAX * totalVectors );
  
  
  lastPositions = calloc( totalVectors, sizeof(ui) );
  
  successes = calloc(totalVectors, sizeof(ui) );

  /* set the last time we saw a vector at last steps */
  for( i = 0; i < totalVectors; ++i)
    lastPositions[i] = MAXSTEPS;

  
#ifdef _U128
  fprintf(stderr, "totalMatrices: ");
  fprint_u128_u(stderr, totalMatrices);
  fprintf(stderr, "\n");
#else
  fprintf(stderr, "totalMatrices: %llu\n", totalMatrices);
#endif
  
  /* for( i = 0; i < totalMatrices; ++i) */
  /*   { */
  /*     popcount[i] = __builtin_popcountll(i); */
  /*   } */

  

#ifdef _U128
#else
  fprintf( stderr, "totalVectors: %llu, totalMatrices: %llu\n", totalVectors, totalMatrices );
#endif
  
  if(developAll)
    {
      
      FILE *resfile = fopen("binaryNetworkOutput.csv", "w");
      
      FILE *matresfile = fopen("matrixOutput.csv", "w");
      
      fprintf(resfile, "input\tmatrix\teqsteps\teqvector\tperiod\n");
      
      matrixSum = 0;
      
      ui matrixMax = 0;

      for( m = 0; m < totalMatrices; ++m )
	{
	  
	  a = b = 0;
	  
	  vectorToInts( m, &a, &b, n);
	  
	  matrixSum = 0; 
	  
	  matrixMax = 0;
	  
	  for( i = 0; i < totalVectors; ++i)
	    {
	      
	      period = 0;
	      
	      /* reset the lastPositions vector. This will 
		 help to detect the period 
	      */
	      for( j = 0; j < totalVectors; ++j)
		lastPositions[j] = MAXSTEPS;
	      
	      in = i;
	      
	      in = development( in, m, d, n, lastPositions,  &period, &stepsToEquilibrium);
	      
	      if( period == 0 )
		noPeriod++;
	      
	      if( period > maxPeriod)
		maxPeriod = period;	  
	      
	      matrixSum += stepsToEquilibrium;
	      
	      if( stepsToEquilibrium > matrixMax )
		matrixMax = stepsToEquilibrium;					      
	      
#ifdef _U128
#else
	      fprintf( resfile, "%llu\t%llu\t%llu\t%llu\t%llu\n", i, m, stepsToEquilibrium, in, period);
#endif
	      
	      if( stepsToEquilibrium < MAXSTEPS && period == 1)
		successes[ in ]++;
	      
	    }
	  
#ifdef _U128

#else
	  fprintf( matresfile, "%llu\t%f\t%llu\n", m, matrixSum/(float)totalVectors, matrixMax);
#endif
	    
	    }

      
#ifdef _U128

#else
	printf( "maxPeriod: %llu\n", maxPeriod);
	printf("noPeriod: %llu\n", noPeriod);
#endif
	

#ifdef _U128

#else
	for( i = 0; i < totalVectors; ++i)
	  printf("successes of %llu are %llu\n", i, successes[i]);
#endif
	fclose(resfile);
	
	fclose(matresfile);
	
        
	}

    ui *pop = calloc(populationSize, sizeof(ui) );

    ui *popM = calloc(populationSize, sizeof(ui));
    
    ui *tempM = calloc(populationSize, sizeof(ui));
    
    ui *noPeriodGen = calloc(generations, sizeof(ui) );
    
    ui *maxPeriodGen = calloc(generations, sizeof(ui) );
    
//    evolveNeutrality( 2, 100, 12, generations, 2, mutationRate, n, pop, popM,  lastPositions, noPeriodGen, maxPeriodGen, totalVectors, successes, d);

//  evolveDriftNeutrality( 2, 100, 12, generations, 2, mutationRate, n, pop, popM,  lastPositions, noPeriodGen, maxPeriodGen, totalVectors, successes, d, tempM);
    
    evolveSelection2( startVector, populationSize, startMatrix, generations, optimum, mutationRate, n, pop, popM,  lastPositions, noPeriodGen, maxPeriodGen, totalVectors, successes, d, totalMatrices, tempM, cyclic, changeOptimum, changeOptGen, mutStop, fitnessType);
    
    free(pop);
    
    free(popM);
    
    free(noPeriodGen);
    
    free(maxPeriodGen);
    
    free(lastPositions);
    
    free(tempM);
    
    free(successes);

    fclose( mutfile);

    fclose(choicefile);
    
    //  evolve( ui v, ui N, ui m, ui g, ui opt);
    
    return 1;
    
    }
