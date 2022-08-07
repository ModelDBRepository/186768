: $Id: stats.mod,v 1.8 2006/06/30 19:03:36 billl Exp $
 
:* COMMENT
COMMENT
randwd   randomly chooses n bits to set to 1
hamming  v.hamming(v1) is hamming distance between 2 vecs
flipbits v.flipbits(scratch,num) flips num rand chosen bits
flipbalbits v.flipbalbits(scratch,num) balanced flipping
vpr      v.vpr prints out vector as 1 (x[i]>0) or 0 (x[i]<=0)
fac      not vec related - returns factorial
logfac   not vec related - returns log factorial
vseed    set some C level randomizer seeds
slope(num) does a linear regression to find the slope, assuming num=timestep of vector
vslope(num) does a linear regression to find the slope, assuming num=timestep of vector
stats(num) does a linear regression, assuming num=timestep of vector
vstats(v2) does a linear regression, using v2 as the x-coords
ENDCOMMENT

NEURON {
  SUFFIX nothing
  GLOBAL  STATS_INSTALLED
}

PARAMETER {
  : BVBASE = -1. : defined in vecst.mod
  STATS_INSTALLED=0
}

ASSIGNED { }

VERBATIM
#include <stdlib.h>
#include <math.h>
//#include <sys/time.h> 
extern double BVBASE;
extern double* hoc_pgetarg();
extern double hoc_call_func(Symbol*, int narg);
extern FILE* hoc_obj_file_arg(int narg);
extern Object** hoc_objgetarg();
extern void vector_resize();
extern int vector_instance_px();
extern void* vector_arg();
extern double* vector_vec();
extern double hoc_epsilon;
extern void set_seed();
extern int ivoc_list_count(Object*);
extern Object* ivoc_list_item(Object*, int);
static void hxe() { hoc_execerror("",0); }

typedef struct BVEC {
 int size;
 int bufsize;
 short *x;
 Object* o;
} bvec;
ENDVERBATIM
 
:* v1.slope(num) does a linear regression to find the slope, assuming num=timestep of vector

VERBATIM
static double slope(void* vv) {
	int i, n;
	double *x, *y;
        double timestep, sigxy, sigx, sigy, sigx2;
	/* how to get the instance data */
	n = vector_instance_px(vv, &y);

        if(ifarg(1)) { 
          timestep = *getarg(1); 
        } else { printf("You must supply a timestep\n"); return 0; }

        sigxy= sigx= sigy= sigx2=0; // initialize these

        x = (double *) malloc(sizeof(double)*n);
        for(i=0; i<n; i++) {
          x[i] = timestep*i;
          sigxy += x[i] * y[i];
          sigx  += x[i];
          sigy  += y[i];
          sigx2 += x[i]*x[i];
        }
        return (n*sigxy - sigx*sigy)/(n*sigx2 - sigx*sigx);
}
ENDVERBATIM
 
:* v1.vslope(v2) does a linear regression, using v2 as the x-coords

VERBATIM
static double vslope(void* vv) {
	int i, n;
	double *x, *y;
        double timestep, sigxy, sigx, sigy, sigx2;
	/* how to get the instance data */
	n = vector_instance_px(vv, &y);

        if(ifarg(1)) {
          if(vector_arg_px(1, &x) != n ) {
            hoc_execerror("Vector size doesn't match.", 0); 
          }
          sigxy= sigx= sigy= sigx2=0; // initialize these

          for(i=0; i<n; i++) {
            sigxy += x[i] * y[i];
            sigx  += x[i];
            sigy  += y[i];
            sigx2 += x[i]*x[i];
          }
        }         
        return (n*sigxy - sigx*sigy)/(n*sigx2 - sigx*sigx);
}
ENDVERBATIM
 
:* v1.stats(num) does a linear regression, assuming num=timestep of vector

VERBATIM
static double stats(void* vv) {
	int i, n;
	double *x, *y;
        double timestep, sigxy, sigx, sigy, sigx2, sigy2;
        double r, m, b;
	/* how to get the instance data */
	n = vector_instance_px(vv, &y);

        if(ifarg(1)) { 
          timestep = *getarg(1); 
        } else { printf("You must supply a timestep\n"); return 0; }

        sigxy= sigx= sigy= sigx2=sigy2= 0; // initialize these

        x = (double *) malloc(sizeof(double)*n);
        for(i=0; i<n; i++) {
          x[i] = timestep*i;
          sigxy += x[i] * y[i];
          sigx  += x[i];
          sigy  += y[i];
          sigx2 += x[i]*x[i];
          sigy2 += y[i]*y[i];
        }
        m = (n*sigxy - sigx*sigy)/(n*sigx2 - sigx*sigx);
        b = (sigy*sigx2 - sigx*sigxy)/(n*sigx2 - sigx*sigx);
        r = (n*sigxy - sigx*sigy)/(sqrt(n*sigx2-sigx*sigx) * sqrt(n*sigy2-sigy*sigy));

        printf("Examined %d data points\n", n);
        printf("slope     = %f\n", m);
        printf("intercept = %f\n", b);
        printf("R         = %f\n", r);
        printf("R-squared = %f\n", r*r);
        return 1;
}
ENDVERBATIM
 
:* v1.vstats(v2) does a linear regression, using v2 as the x-coords

VERBATIM
static double vstats(void* vv) {
	int i, n;
	double *x, *y;
        double timestep, sigxy, sigx, sigy, sigx2, sigy2;
        double r, m, b;
	/* how to get the instance data */
	n = vector_instance_px(vv, &y);

        if(ifarg(1)) {
          if(vector_arg_px(1, &x) != n ) {
            hoc_execerror("Vector size doesn't match.", 0); 
          }
          sigxy= sigx= sigy= sigx2=sigy2=0; // initialize these

          for(i=0; i<n; i++) {
            sigxy += x[i] * y[i];
            sigx  += x[i];
            sigy  += y[i];
            sigx2 += x[i]*x[i];
            sigy2 += y[i]*y[i];
          }
          m = (n*sigxy - sigx*sigy)/(n*sigx2 - sigx*sigx);
          b = (sigy*sigx2 - sigx*sigxy)/(n*sigx2 - sigx*sigx);
          r = (n*sigxy - sigx*sigy)/(sqrt(n*sigx2-sigx*sigx) * sqrt(n*sigy2-sigy*sigy));

          printf("Examined %d data points\n", n);
          printf("slope     = %f\n", m);
          printf("intercept = %f\n", b);
          printf("R         = %f\n", r);
          printf("R-squared = %f\n", r*r);
          return 1;
        } else {
          printf("You must supply an x vector!\n");
          return 0;
        }
}
ENDVERBATIM
 
:* v1.randwd(num[,v2]) will randomly flip num bits from BVBASE to 1
: does v1.fill(BVBASE); optionally fill v2 with the indices
VERBATIM
static double randwd(void* vv) {
	int i, ii, jj, nx, ny, flip, flag;
	double* x, *y;
	/* how to get the instance data */
	nx = vector_instance_px(vv, &x);
        flip = (int) *getarg(1);
        if (ifarg(2)) { /* write a diff vector to z */
          flag = 1; ny = vector_arg_px(2, &y);
          if (ny!=flip) { hoc_execerror("Opt vector must be size for # of flips", 0); }
        } else { flag = 0; }
        if (flip>=nx) { hoc_execerror("# of flips exceeds (or ==) vector size", 0); }
	for (i=0; i < nx; i++) { x[i] = BVBASE; }
	for (i=0,jj=0; i < flip; i++) { /* flip these bits */
	  ii = (int) ((nx+1)*((double)rand()/RAND_MAX)); // replaced drand48() with (double)rand()/RAND_MAX
	  if (x[ii]==BVBASE) {
	    x[ii] = 1.; 
            if (flag) { y[jj] = ii; jj++; }
	  } else {
	    i--;
	  }
	}
	return flip;
}
ENDVERBATIM
 
:* v1.hamming(v2[,v3]) compares v1 and v2 for matches, v3 gives diff vector
VERBATIM
static double hamming(void* vv) {
	int i, nx, ny, nz, prflag;
	double* x, *y, *z,sum;
	sum = 0.;
	nx = vector_instance_px(vv, &x);
	ny = vector_arg_px(1, &y);
        if (ifarg(2)) { /* write a diff vector to z */
          prflag = 1; nz = vector_arg_px(2, &z);
        } else { prflag = 0; }
	if (nx!=ny || (prflag && nx!=nz)) {
	  hoc_execerror("Vectors must be same size", 0);
	}
	for (i=0; i < nx; ++i) {
	  if (x[i] != y[i]) { sum++; 
            if (prflag) { z[i] = 1.; }
          } else if (prflag) { z[i] = 0.; }
        }
	return sum;
}
ENDVERBATIM

:* v1.flipbits(scratch,num) flips num bits
: uses scratch vector of same size as v1 to make sure doesn't flip same bit twice
VERBATIM
static double flipbits(void* vv) {	
	int i, nx, ny, flip, ii;
	double* x, *y;

	nx = vector_instance_px(vv, &x);
	ny = vector_arg_px(1, &y);
        flip = (int)*getarg(2);
	if (nx != ny) {
	  hoc_execerror("Scratch vector must be same size", 0);
	}
	for (i=0; i<nx; i++) { y[i]=x[i]; } /* copy */
	for (i=0; i < flip; i++) { /* flip these bits */
	  ii = (int) ((nx+1)*((double)rand()/RAND_MAX)); // replaced drand48() with (double)rand()/RAND_MAX
	  if (x[ii]==y[ii]) { /* hasn't been touched */
	    x[ii]=((x[ii]==1.)?BVBASE:1.);
	  } else {
	    i--; /* do it again */
	  }
	}
	return flip;
}
ENDVERBATIM
 
:* v1.flipbalbits(scratch,num) flips num bits making sure to balance every 1
: flip with a 0 flip to preserve initial power
: uses scratch vector of same size as v1 to make sure doesn't flip same bit twice
VERBATIM
static double flipbalbits(void* vv) {	
	int i, nx, ny, flip, ii, next;
	double* x, *y;

	nx = vector_instance_px(vv, &x);
	ny = vector_arg_px(1, &y);
        flip = (int)*getarg(2);
	if (nx != ny) {
	  hoc_execerror("Scratch vector must be same size", 0);
	}
	for (i=0; i<nx; i++) { y[i]=x[i]; } /* copy */
        next = 1; /* start with 1 */
	for (i=0; i < flip;) { /* flip these bits */
	  ii = (int) ((nx+1)*((double)rand()/RAND_MAX)); // replaced drand48() with (double)rand()/RAND_MAX
	  if (x[ii]==y[ii] && y[ii]==next) { /* hasn't been touched */
	    next=x[ii]=((x[ii]==1.)?BVBASE:1.);
            i++;
	  }
	}
	return flip;
}
ENDVERBATIM
 
:* v1.vpr() prints out neatly
VERBATIM
static double vpr(void* vv) {
  int i, nx;
  double* x;
  FILE* f;
  nx = vector_instance_px(vv, &x);
  if (ifarg(1)) { 
    f = hoc_obj_file_arg(1);
    for (i=0; i<nx; i++) {
      if (x[i]>BVBASE) { fprintf(f,"%d",1); 
      } else { fprintf(f,"%d",0); }
    }
    fprintf(f,"\n");
  } else {
    for (i=0; i<nx; i++) {
      if (x[i]>BVBASE) { printf("%d",1); 
      } else { printf("%d",0); }
    }
    printf("\n");
  }
  return 1.;
}
ENDVERBATIM

:* v1.bin(targ,invl) place counts for each interval
VERBATIM
static double bin(void* vv) {	
  int i, j, nx, ny, next, maxsz;
  double* x, *y, invl, loc;

  nx = vector_instance_px(vv, &x);
  ny = vector_arg_px(1, &y);
  invl = (int)*getarg(2);

  vv=vector_arg(1);
  maxsz=vector_buffer_size(vv);
  vector_resize(vv, maxsz);
  if (x[nx-1]/invl>(double)(maxsz-1)) {
    printf("Need size %d in target vector (%d)\n",(int)(x[nx-1]/invl+1),maxsz); 
    hoc_execerror("",0); }
  for (j=0; j<maxsz; j++) y[j]=0.;
  for (i=0,j=0,loc=invl; i<nx && j<maxsz; i++) {
    if (x[i]<loc) y[j]++; else {
      while (x[i]>loc) { loc+=invl; j++; }
      y[j]++;
    }
  }
  vector_resize(vv, j);
  return (double)j;
}
ENDVERBATIM

:* PROCEDURE install_stats()
PROCEDURE install_stats () {
  STATS_INSTALLED=1
VERBATIM
  install_vector_method("slope", slope);
  install_vector_method("vslope", vslope);
  install_vector_method("stats", stats);
  install_vector_method("vstats", vstats);
  install_vector_method("randwd", randwd);
  install_vector_method("hamming", hamming);
  install_vector_method("flipbits", flipbits);
  install_vector_method("flipbalbits", flipbalbits);
  install_vector_method("vpr", vpr);
  install_vector_method("bin", bin);
ENDVERBATIM
}

:* fac (n) 
: from numerical recipes p.214
FUNCTION fac (n) {
VERBATIM {
    static int ntop=4;
    static double a[101]={1.,1.,2.,6.,24.};
    static double cof[6]={76.18009173,-86.50532033,24.01409822,
      -1.231739516,0.120858003e-2,-0.536382e-5};
    int j,n;
    n = (int)_ln;
    if (n<0) { hoc_execerror("No negative numbers ", 0); }
    if (n>100) { /* gamma function */
      double x,tmp,ser;
      x = _ln;
      tmp=x+5.5;
      tmp -= (x+0.5)*log(tmp);
      ser=1.0;
      for (j=0;j<=5;j++) {
        x += 1.0;
        ser += cof[j]/x;
      }
      return exp(-tmp+log(2.50662827465*ser));
    } else {
      while (ntop<n) {
        j=ntop++;
        a[ntop]=a[j]*ntop;
      }
    return a[n];
    }
}
ENDVERBATIM
}
 
:* logfac (n)
: from numerical recipes p.214
FUNCTION logfac (n) {
VERBATIM {
    static int ntop=4;
    static double a[101]={1.,1.,2.,6.,24.};
    static double cof[6]={76.18009173,-86.50532033,24.01409822,
      -1.231739516,0.120858003e-2,-0.536382e-5};
    int j,n;
    n = (int)_ln;
    if (n<0) { hoc_execerror("No negative numbers ", 0); }
    if (n>100) { /* gamma function */
      double x,tmp,ser;
      x = _ln;
      tmp=x+5.5;
      tmp -= (x+0.5)*log(tmp);
      ser=1.0;
      for (j=0;j<=5;j++) {
        x += 1.0;
        ser += cof[j]/x;
      }
      return (-tmp+log(2.50662827465*ser));
    } else {
      while (ntop<n) {
        j=ntop++;
        a[ntop]=a[j]*ntop;
      }
    return log(a[n]);
    }
}
ENDVERBATIM
}

: unable to get the drand here to recognize the same fseed used in rand
FUNCTION vseed () {
  VERBATIM

#ifdef WIN32
  double seed;
  if (ifarg(1)) seed=*getarg(1); else {
    printf("TIME ACCESS NOT PRESENT IN WINDOWS\n");
    hxe();
  }
  srand((unsigned)seed); // changed from srand48 to srand - mjb
  set_seed(seed);
  return seed;
#else
  //struct  timeval tp;
  //struct  timezone tzp;
  double seed;
  if (ifarg(1)) seed=*getarg(1); else {
    //gettimeofday(&tp,&tzp); // tzp unknown storage size
    seed=2; //Wtp.tv_usec;
  }
  srand((unsigned)seed); // srand48
  set_seed(seed);
  srandom(seed);
  return seed;
#endif
  ENDVERBATIM
}

: from Numerical Recipes in C
FUNCTION gammln (xx) {
  VERBATIM {
    double x,tmp,ser;
    static double cof[6]={76.18009173,-86.50532033,24.01409822,-1.231739516,0.120858003e-2,-0.536382e-5};
    int j;
    x=_lxx-1.0;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.0;
    for (j=0;j<=5;j++) {
      x += 1.0;
      ser += cof[j]/x;
    }
    return -tmp+log(2.50662827465*ser);
  }
  ENDVERBATIM
}

FUNCTION betai(a,b,x) {
VERBATIM {
  double bt;
  if (_lx < 0.0 || _lx > 1.0) {printf("Bad x in routine BETAI\n"); hxe();}
  if (_lx == 0.0 || _lx == 1.0) bt=0.0;
  else
  bt=exp(gammln(_la+_lb)-gammln(_la)-gammln(_lb)+_la*log(_lx)+_lb*log(1.0-_lx));
  if (_lx < (_la+1.0)/(_la+_lb+2.0))
  return bt*betacf(_la,_lb,_lx)/_la;
  else
  return 1.0-bt*betacf(_lb,_la,1.0-_lx)/_lb;
 }
ENDVERBATIM
}

VERBATIM
#define ITMAX 100
#define EPS 3.0e-7
ENDVERBATIM

FUNCTION betacf(a,b,x) {
VERBATIM {
  double qap,qam,qab,em,tem,d;
  double bz,bm=1.0,bp,bpp;
  double az=1.0,am=1.0,ap,app,aold;
  int m;
  void nrerror();

  qab=_la+_lb;
  qap=_la+1.0;
  qam=_la-1.0;
  bz=1.0-qab*_lx/qap;
  for (m=1;m<=ITMAX;m++) {
    em=(double) m;
    tem=em+em;
    d=em*(_lb-em)*_lx/((qam+tem)*(_la+tem));
    ap=az+d*am;
    bp=bz+d*bm;
    d = -(_la+em)*(qab+em)*_lx/((qap+tem)*(_la+tem));
    app=ap+d*az;
    bpp=bp+d*bz;
    aold=az;
    am=ap/bpp;
    bm=bp/bpp;
    az=app/bpp;
    bz=1.0;
    if (fabs(az-aold) < (EPS*fabs(az))) return az;
  }
  printf("a or b too big, or ITMAX too small in BETACF"); return -1.;
}
ENDVERBATIM
}
