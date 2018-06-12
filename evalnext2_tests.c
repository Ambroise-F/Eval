
#include "int128g.c"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#include <math.h>

#define DEBUG 0
#define LONG long long int
#define ULONG unsigned long long int

#define REAL double
#define FMA(a,b,c) fma(a,b,c)
#define FLOOR(a) floor(a)

REAL pr;
REAL u;

LONG * array( int n ) { return (LONG *) malloc( (LONG) n*sizeof(LONG) ); }

ULONG seed, mult;

LONG rand64s(LONG p) {
    LONG x,y;
    seed = mult*seed;
    x = seed >> 32;
    seed = mult*seed;
    y = seed >> 32;
    x = (x<<31) | y;
    x = x % p;
    return(x);
}


/******************************************************************************************/
/*  Roman's assembler utilities                                                           */
/******************************************************************************************/


        /* z += a1:a0 */
        #define zadd(z,a0,a1) __asm__(\
        "       addq    %4, %0  \n\t" \
        "       adcq    %5, %1  \n\t" \
                : "=&r"(z[0]), "=&r"(z[1]) : "0"(z[0]), "1"(z[1]), "r"(a0), "r"(a1))

        /* z -= a1:a0 */
        #define zsub(z,a0,a1) __asm__(\
        "       subq    %4, %0  \n\t" \
        "       sbbq    %5, %1  \n\t" \
                : "=&r"(z[0]), "=&r"(z[1]) : "0"(z[0]), "1"(z[1]), "r"(a0), "r"(a1))

        /* z = a*b */
        #define zmul(z,a,b) __asm__(\
        "       mulq    %%rdx   \n\t" \
                : "=a"(z[0]), "=d"(z[1]) : "a"(a), "d"(b))

        /* z += a*b */
        #define zfma(z,a,b) do {        \
        unsigned long u,v;              \
        __asm__(                        \
        "       mulq    %%rdx           \n\t" \
        "       addq    %%rax, %0       \n\t" \
        "       adcq    %%rdx, %1       \n\t" \
                : "=&r"(z[0]), "=&r"(z[1]), "=a"(u), "=d"(v) : "0"(z[0]), "1"(z[1]), "a"(a), "d"(b));\
        } while (0)

        /* z -= a*b */
        #define zfms(z,a,b) do {        \
        unsigned long u,v;              \
        __asm__(                        \
        "       mulq    %%rdx           \n\t" \
        "       subq    %%rax, %0       \n\t" \
        "       sbbq    %%rdx, %1       \n\t" \
                : "=&r"(z[0]), "=&r"(z[1]), "=a"(u), "=d"(v) : "0"(z[0]), "1"(z[1]), "a"(a), "d"(b));\
        } while (0)

        /* z[0] = z % p */
        /* z[1] = z / p */
        /* quotient can overflow */
        #define zdiv(z,p) __asm__(\
        "       divq    %4      \n\t" \
                : "=a"(z[1]), "=d"(z[0]) : "a"(z[0]), "d"(z[1]), "r"(p))

        /* z = z % p safe */
        #define zmod(z,p) __asm__(\
        "       divq    %4      \n\t" \
        "       xorq    %0, %0  \n\t" \
                : "=a"(z[1]), "=d"(z[0]) : "a"(z[0]), "d"(z[1] < p ? z[1] : z[1] % p), "r"(p))

        /* z = z << s */
        #define zshl(z,s) __asm__(\
        "       shldq   %%cl, %0, %1    \n\t" \
        "       shlq    %%cl, %0        \n\t" \
                : "=&r"(z[0]), "=&r"(z[1]) : "0"(z[0]), "1"(z[1]), "c"(s))

/******************************************************************************************/
/*  Zp utilities                                                                          */
/******************************************************************************************/

static inline LONG add64s(LONG a, LONG b, LONG p) { LONG t; t = (a-p)+b; t += (t>>63) & p; return t; }
static inline LONG sub64s(LONG a, LONG b, LONG p) { LONG t; t = a-b; t += (t>>63) & p; return t; }
static inline LONG neg64s(LONG a, LONG p) { LONG t; t = -a; t += (t>>63) & p; return t; }
static inline LONG mul64s(LONG a, LONG b, LONG p) {
        LONG q, r;
        __asm__ __volatile__(           \
        "       mulq    %%rdx           \n\t" \
        "       divq    %4              \n\t" \
        : "=a"(q), "=d"(r) : "0"(a), "1"(b), "rm"(p));
        return r;
}

typedef struct {
	REAL p;	/* p */
        REAL u;	/* recriprocal of p */
} recreal;


REAL add50s(REAL a, REAL b){
  REAL s;
  s = a + b;
  return s>pr?s-pr:s;
}//pr global var

#pragma omp declare simd
REAL mul50s(REAL x, REAL y){
  REAL h,l,b,c,d,g,gs,gb;
  int ts,tb;
  h = x*y;
  l = FMA(x,y,-h);
  b = h*u;
  c = FLOOR(b);
  d = FMA(-c,pr,h);
  g = d + l;
  ts = g<0; //too small
  tb = g>=pr;//too big
  gs = g+pr;
  gb = g-pr;
  if (ts) { g=gs; }
  if (tb) { g=gb; }
  return g;    
}



/* a^n mod p assuming 0 <= a < p < 2^63 */
LONG powmod64s( LONG a, LONG n, LONG p )
{   LONG r,s;
    a += (a>>63) & p; // protect from bad input
    if( n==0 ) return 1;
    if( n==1 ) return a;
    for( r=1, s=a; n>0; n /= 2 ) { if( n & 1 ) r = mul64s(r,s,p); s = mul64s(s,s,p); }
    return r;
}


/******************************************************************************************/
/* Polynomial utilities                                                                   */
/******************************************************************************************/

/* print an array in form [a0,a1,...,an-1] */
void vecprint64s( LONG *A, int n )
{   int i;
    printf("[");
    for( i=0; i<n; i++ ) { printf("%lld",A[i]); if( i<n-1 ) printf(", "); }
    printf("]");
    return;
}

void polnsparseprint64s( LONG *C, LONG *X, int t, int n ) {
    int numbits,i,j;
    LONG m,d,mask,D[63];
    if( t==0 ) { printf("0;\n"); return; }
    numbits = 63/n;
    mask = (1LL<<numbits)-1;
    for( i=0; i<t; i++ ) {
        m = X[i];
        printf("%lld",C[i]);
        for( j=0; j<n; j++ ) {
            D[j] = m & mask;
            m = m >> numbits;
        }
        for( j=1; j<=n; j++ ) {
            d = D[n-j];
            if( d!=0 ) printf("*x%d^%lld",j,d);
        }
        if( i!=t-1 ) printf("+");
    }
    return;
}

void pol2insert( LONG *C, LONG *X, int n ) {
    // insert X[n] into X[0],X[1],...,X[n-1]
    LONG x,c;
    c = C[n];
    x = X[n];
    while( n>0 && X[n-1]>x ) { X[n] = X[n-1]; C[n] = C[n-1]; n--; }
    C[n] = c;
    X[n] = x;
    return;
}

void pol2sort64s( LONG *C, LONG *X, int n ) {
    // sort monomials in X into ascending order
    int i;
    for( i=1; i<n; i++ ) pol2insert(C, X, i); // insertion sort
    return;
}


/******************************************************************************************/
/* Polynomial routines                                                                    */
/******************************************************************************************/

void randpoly64s( int n, int *DEG, int t, LONG *A, LONG *X, LONG p )
{   // monic in x[n-1]
    int numbits,i,j;
    LONG x,c,y;
    numbits = 63/n;
    for( i=0; i<t; i++ ) {
        for( c=0; c==0; c = rand64s(p) );
        A[i] = c;
        x = rand64s(DEG[0]);
        for( j=1; j<n; j++ ) { y = rand64s(DEG[j]+1); x = (x << numbits) | y; }
        X[i] = x;
    }
    // make it monic in x[0] and not be divisible by x[0]
    A[t-1] = 1;
    X[t-1] = (LONG)DEG[0] << (numbits*(n-1));
    X[0  ] = X[0] & ( (1LL << numbits*(n-1)) - 1 );
}

LONG evalnext2( LONG *A, LONG *X, LONG *G, int na, 
                LONG *B, LONG *C, LONG *Y, LONG *Z, 
                int * nave,
                LONG p )
{   

    // a = sum( A[i] X[i], i=1..na ) where the monomials in X are sorted in ascending order
    //Compute A[i] = G[i] A[i] mod p and create B,Y = sum( A[i] X[i], i=1..na )
    //Compute A[i] = G[i] A[i] mod p and create C,Z = sum( A[i] X[i], i=1..na ) 
    // So that's two evaluations, one in B,Y and the second in C,Z

    int i,j,k,m,n;
    LONG c,d,t,x,ntotal,ncount;
    recint P;
    P = recip1(p);
    k = m = 0;
    ntotal = 0;
    ncount = 0;

    REAL tr,cr,dr,pr;
    REAL *Ar, *Gr;
    pr = (REAL)p;
    u = (REAL)1.0/u;
    //Ar = malloc(sizeof(REAL)*na); //aligned
    //Gr = malloc(sizeof(REAL)*na); //aligned
    
    posix_memalign((void **)&Ar, 32,na*sizeof(REAL));
    posix_memalign((void **)&Gr, 32,na*sizeof(REAL));
    

    // conversion
    for( i=0; i<na; i++ ) {  
      Ar[j] = (REAL)A[j]; Gr[j] = (REAL)G[j];
    }
    
    for( i=0; i<na; i=n ) {
        for( x=X[i],n=i+1; n<na && X[n]==x; n++ );
        ntotal += (n-i);
        ncount ++;

	// ( A * G * G ) % P

	// Mieux sur knighstlanding (deux boucles : 2% plus lent)
	
	/* c = d = 0; */
        /* for( j=i; j<n; j++ ) { */
        /*   t = mulrec64(A[j],G[j],P); c = add64s(t,c,p); */
        /*   t = mulrec64(t,G[j],P); A[j] = t; d = add64s(t,d,p); */
        /* } */

	// Mieux sur sfpn1 (une seule boucle : 13% plus lent)
	
//        for( c=0,j=i; j<n; j++ ) { t = mulrec64(A[j],G[j],P); A[j] = t; c = add64s(t,c,p); }
//        for( d=0,j=i; j<n; j++ ) { t = mulrec64(A[j],G[j],P); A[j] = t; d = add64s(t,d,p); }


	// En flottants
#pragma omp declare reduction (add50s:REAL:omp_out=add50s(omp_out,omp_in)) initializer (omp_priv=0.0)
	//sfpn1 1b -> ??? ms
	//sfpn1 2b -> 620 ms
	cr=0.0;
	dr=0.0;

	//1b
/* #pragma omp simd reduction (add50s:cr,dr) */
/* 	for( j=i; j<n; j++ ) { *\/ */
/* 	  t = mulrec64(A[j],G[j],P); c = add64s(t,c,p); */
/* 	  t = mulrec64(t,G[j],P); A[j] = t; d = add64s(t,d,p); */
/* 	} */	


	
	//2b
//#pragma omp simd reduction (add50s:cr)
 	for(j=i; j<n; j++ ) {
	  tr = mul50s(Ar[j],Gr[j]);
	  Ar[j] = tr;
	  cr = add50s(tr,cr);
	}
	
//#pragma omp simd reduction (add50s:dr)
	for(j=i; j<n; j++ ) {
	  tr = mul50s(Ar[j],Gr[j]);
	  Ar[j] = tr;
	  dr = add50s(tr,dr);
	}

	
    
        if( cr!=0.0 ) { B[k] = (LONG)c; Y[k] = x; k++; }
        if( dr!=0.0 ) { C[m] = (LONG)d; Z[m] = x; m++; }

        /* if( c!=0 ) { B[k] = c; Y[k] = x; k++; } */
        /* if( d!=0 ) { C[m] = d; Z[m] = x; m++; } */
    }

    if( ncount>0 ) nave[0] = ntotal/ncount;
    t = ((LONG) k << 32) | m;
    return t;
}


void evalinit( LONG *X, int t, int n, LONG *BETA, int m, LONG *G, LONG p )
{   // Evaluate the last m variables in the monoimials in X
    // Put the monomial evaluations in G and update the monomials in X
    int i,j,numbits,shift,d;
    LONG x,y,z,mask;
   
    numbits = 63/n;
    mask = (1LL<<numbits)-1;

    for( i=0; i<t; i++ ) {
       z = 1;
       x = X[i];
       for( j=m-1; j>=0; j-- ) {
           d = x & mask;
           z = mul64s(z,powmod64s(BETA[j],d,p),p);
           x = x >> numbits;
       }
       G[i] = z;
// for this test we will just put x back into X[i]
       y = 0;
       shift = 63/(n-m);
       for( j=0; j<n-m; j++ ) {
           y = y | ( (x&mask) << (j*shift) );
           x = x >> numbits;
       }
       X[i] = y;
    }
    return;
}

/******************************************************************************************/
/* Test routine                                                                           */
/******************************************************************************************/

void test( int numevals,
           int n, // #variables
           int m, // #variables evaluated = n-2
           int t, // #terms
           int d, // #max degree in each variable
           LONG p, // prime
           int print // print polynomials to check in Maple
         ) {
          
int i,tB,tC,D[32],dB,dC,nave;
LONG s, *A, *X, *B, *Y, *C, *Z, *G, *BETA, mask;
clock_t T1,T2;

   for( i=0; i<n; i++ ) D[i] = d;

   A = array(t);
   X = array(t);

      T1 = clock();
   randpoly64s( n, D, t, A, X, p );
   //printf("A := "); polnsparseprint64s( A, X, t, n ); printf(";\n");
      T2 = clock();
      if( !print ) printf("randpoly time = %lld ms\n", (T2-T1)/1000 );
   pol2sort64s( A, X, t );
      for( i=0; i<t-1; i++ ) if( X[i]>X[i+1] ) printf("NOT SORTED\n"); 
      T1 = clock();
      if( !print ) printf("sort time = %lld ms\n", (T1-T2)/1000 );


if( print ) {
   printf("p := %lld;\n",p);
   printf("A := "); polnsparseprint64s( A, X, t, n ); printf(";\n");
}

   BETA = array(m);
   for( i=0; i<m; i++ ) BETA[i] = rand64s(p);

   G = array(t);

      T1 = clock();
   evalinit( X, t, n, BETA, m, G, p );
      T2 = clock();
      if( !print ) printf("evalinit time = %lld ms\n", (T2-T1)/1000 );
   

   B = array(t); Y = array(t);
   C = array(t); Z = array(t);

      T1 = clock();
   for( i=0; i<numevals; i++ ) s = evalnext2( A, X, G, t, B, C, Y, Z, &nave, p );
      T2 = clock();
      if( !print ) printf("-----------------------------\nevalnext2 time = %lld ms\n-----------------------------\n", (T2-T1)/1000 );

if( !print )   printf("Average evaluation length = %d\n",nave);

   tB = s >> 32;
   mask = ((1LL) << 32) - 1;
   tC = s & mask;

if( print ) {
   printf("\nB := "); polnsparseprint64s( B, Y, tB, 2 ); printf(":\n");

   printf("BB := eval(A,{");
   for( i=3; i<=n; i++ ) if( i==3 ) printf("x%d=%lld",i,BETA[i-m]);
                         else       printf(",x%d=%lld",i,BETA[i-m]);
   printf("}) mod %lld:  # should equal B \n",p);

   printf("BB - B; # should be zero\n");

   printf("\nC := "); polnsparseprint64s( C, Z, tC, 2 ); printf(":\n");
   printf("CC := eval(A,{");
   for( i=3; i<=n; i++ ) {
       t = mul64s(BETA[i-m],BETA[i-m],p);
       if( i==3 ) printf("x%d=%lld",i,t); else printf(",x%d=%lld",i,t);
   }
   printf("}) mod %lld:  # should equal C \n",p);
   printf("CC - C; # should be zero\n");
}

   free(A); free(B); free(C); free(X); free(Y); free(Z); free(G); free(BETA);

}

int main() {

int d,n,m,t,N;
LONG p;

   seed = 1;
   mult = 6364136223846793003ll;

if( DEBUG ) {

   p = 11;  // prime
   t = 10;  // #terms
   n = 5;   // #variables
   m = n-2; // #variables evaluated -- must be n-2
   d = 3;   // degree
   N = 1;   // 2 evaluations

   test( 1, n, m, t, d, p, 1 );

}

   printf("\n  *** Timing test *** \n");
   p = 9223372036854775783LL; //2**63 - 25 (sur 63 bits)
   //p = 1125899906842597; //(sur 50 bits)
   t = 50000;   printf("  #terms = %d\n",t);
   n = 6;       printf("  #vars = %d\n",n);
   m = n-2;     printf("  #vars evaluated = %d\n",m);
   d = 10;      printf("  degree = %d\n",d);
   N = 5000;    printf("  #evaluations = %d\n\n", 2*N);

   test( N, n, m, t, d, p, 0 );
   return 1;

}
