: $Id: vecst.mod,v 1.326 2006/06/24 23:38:36 billl Exp $
 
:* COMMENT
COMMENT
thresh   turns analog vec to BVBASE,1 vec separating at thresh (scalar or vec)
triplet  return location of a triplet in a vector
onoff    turns state vec on or off depending on values in other vecs
bpeval   used by backprop: vo=outp*(1-outp)*del
w        like where but sets chosen nums // modified from .wh in 1.23
dest.xing(src,tvec,thresh) determines where vector crosses in a positive direction 
src.updown(thresh,destlist) determines crossings in both directions+peaks
dest.snap(src,tvec,dt) interpolate src with tvec to prior dt step, saves only highest value
xzero    count 0 crossings by putting in a 1 or -1
negwrap  wrap negative values around to pos (with flag just set to 0)
indset(ind,val) sets spots indexed by ind to val
ismono([arg]) checks if a vector is monotonically increaseing (-1 decreasing)
count(num) count the number of nums in vec
scr.fewind(ind,veclist) // uses ind as index into other vecs
ind.findx(vecAlist,vecBlist) // uses ind as index into vecA's to load in vecB's
ind.sindx(vecAlist,vecBlist) // replace ind elements in vecA's with vecB's values
ind.sindv(vecAlist,valvec) // replace ind elements in vecA's with vals
scr.nind(ind,vec1,vec2[,vec3,vec4]) // uses ind to elim elems in other vecs
ind.keyind(key,vec1,vec2[,vec3,vec4]) // pick out bzw. values from vectors
ind.slct(key,args,veclist) // pick out bzw. values from vectors
ind.slor(key,args,veclist) // do OR rather than AND function of slct
vdest.intrp(flag) // interpolate numbers replacing numbers given as flag
v.insct(v1,v2)   // return v1 intersect v2
vdest.cull(vsrc,key)   // remove values found in key
vdest.redundout(vsrc[,INDFLAG])  // remove repeat values, with flag return indices
ind.mredundout(veclistA[,INDFLAG,veclistB]) // remove repeats from parallel vectors
v.cvlv(v1,v2)   // convolve v1 with v2
v2d      copy into a double block
d2v      copy out of a double block NB: same as vec.from_double()
smgs     rewrite of sumgauss form ivovect.cpp
smsy     sum syn potentials off a tvec
iwr      write integers
ird      read integers
ident     give pointer addresses and sizes
lcat     concatentate all vectors from a list
fread2   like vec.fread but uses flag==6 for unsigned int
vfill    fill vdest with multiple instances of vsrc until reach size
inv    return inverse of number multiplied by optional num
slone(src,val) select indices where ==VAL from sorted vector SRC
piva.join(pivb,veclista,veclistb) copies values from one set of vecs to other
Non-vector routines
isojt    compare 2 objects to see if they're of the same type
eqojt    compare 2 object pointers to see if they point to same thing
sumabs   return sum of absolute values
rnd    round off to nearest integer
ENDCOMMENT

NEURON {
  SUFFIX nothing
  : BVBASE is bit vector base number (typically 0 or -1)
  GLOBAL BVBASE, RES, VECST_INSTALLED, SHM_UPDOWN, DEBUG_VECST
}

PARAMETER {
  BVBASE = -1.
  VECST_INSTALLED=0
  DEBUG_VECST=0
  : misc
  ERR=-1.3479e121

  : 0 args
  ALL=-1.3479e120
  NEG=-1.3478e120
  POS=-1.3477e120
  CHK=-1.3476e120
  NOZ=-1.3475e120
  : 1 arg
  GTH=-1.3474e120
  GTE=-1.3473e120
  LTH=-1.3472e120
  LTE=-1.3471e120
  EQU=-1.3470e120
  EQV=-1.3469e120 : value equal to same row value in parallel vector
  EQW=-1.3468e120 : value found in other vector
  NEQ=-1.3467e120
  SEQ=-1.3466e120
  RXP=-1.3465e120
  : 2 args
  IBE=-1.3464e120
  EBI=-1.3463e120
  IBI=-1.3462e120
  EBE=-1.3461e120

  SHM_UPDOWN=4   : used in updown() for measuring sharpness
}

ASSIGNED { RES }

VERBATIM
#ifndef NRN_VERSION_GTEQ_8_2_0
#include <stdlib.h>
#include <math.h>
#include <sys/time.h> 
extern double* hoc_pgetarg();
extern double hoc_call_func(Symbol*, int narg);
extern FILE* hoc_obj_file_arg(int narg);
extern Object** hoc_objgetarg();
extern void vector_resize();
extern int vector_instance_px();
extern void* vector_arg();
extern double* vector_vec();
extern double hoc_epsilon;
extern double chkarg();
extern void set_seed();
extern int ivoc_list_count(Object*);
extern Object* ivoc_list_item(Object*, int);
extern int hoc_is_double_arg(int narg);
extern char* hoc_object_name(Object*);
char ** hoc_pgargstr();
#endif
static int ismono1 (double *x, int n, int flag);
int list_vector_px(Object *ob, int i, double** px);
int list_vector_px2 (Object *ob, int i, double** px, void** vv);
int list_vector_resize (Object *ob, int i, int sz);
static double sc[6];
static void hxe() { hoc_execerror("",0); }

typedef struct BVEC {
 int size;
 int bufsize;
 short *x;
 Object* o;
} bvec;

#define VRRY 50

#define BYTEHEADER int _II__;  char *_IN__; char _OUT__[16]; int BYTESWAP_FLAG=0;
#define BYTESWAP(_X__,_TYPE__) \
    if (BYTESWAP_FLAG == 1) { \
	_IN__ = (char *) &(_X__); \
	for (_II__=0;_II__<sizeof(_TYPE__);_II__++) { \
		_OUT__[_II__] = _IN__[sizeof(_TYPE__)-_II__-1]; } \
	(_X__) = *((_TYPE__ *) &_OUT__); \
    }

#define UNCODE(_X_,_J_,_Y_) {(_Y_)=floor((_X_)/sc[(_J_)])/sc[4]; \
                             (_Y_)=floor(sc[4]*((_Y_)-floor(_Y_))+0.5);}
ENDVERBATIM

:* v1.ident() gives addresses and sizes
VERBATIM
static double ident(void* vv) {
  int nx,bsz; double* x;
  nx = vector_instance_px(vv, &x);
  bsz=vector_buffer_size(vv);
  printf("Obj*%x Dbl*%x Size: %d Bufsize: %d\n",vv,x,nx,bsz);
  return 0.0;
}
ENDVERBATIM
 
:* v1.indset(ind,x[,y]) sets indexed values to x and other values to optional y
VERBATIM
static double indset(void* vv) {
	int i, nx, ny, nz, flag;
	double* x, *y, *z, val, val2 ;
	nx = vector_instance_px(vv, &x);
	ny = vector_arg_px(1, &y);
        if (hoc_is_object_arg(2)) { 
          flag=1;
          nz = vector_arg_px(2, &z); 
          if (ny!=nz) { hoc_execerror("v.indset: Vector sizes don't match.", 0); }
        } else { flag=0; val=*getarg(2); }
        if (ifarg(3)) { 
          val2 = *getarg(3); 
          for (i=0; i<nx; i++) { x[i]=val2; }
        }
	for (i=0; i<ny; i++) {
	  if (y[i] > nx) { hoc_execerror("v.indset: Index exceeds vector size", 0); }
          if (flag) x[(int)y[i]]=z[i]; else x[(int)y[i]]=val;
        }
	return (double)i;
}
ENDVERBATIM
 
VERBATIM
// Maintain parallel int vector to avoid slowness of repeated casts 
static unsigned int scrsz=0;
static unsigned int bufsz=0;
static unsigned int *scr;
static double dcr[100]; // scratch area for doubles
ENDVERBATIM

:* tmp.fewind(ind,veclist)
: picks out numbers from multiple vectors using index ind
VERBATIM
static double fewind (void* vv) {
  int i, j, nx, ni, nv[VRRY], num;
  Object* ob;
  double *x, *ind, *vvo[VRRY];
  nx = vector_instance_px(vv, &x);
  ni = vector_arg_px(1, &ind);
  ob = *hoc_objgetarg(2);
  num = ivoc_list_count(ob);
  if (num>VRRY) hoc_execerror("ERR: fewind can only handle VRRY vectors", 0);
  for (i=0;i<num;i++) { 
    nv[i] = list_vector_px(ob, i, &vvo[i]);
    if (nx!=nv[i]) { printf("fewind ERR %d %d %d\n",i,nx,nv[i]);
      hoc_execerror("Vectors must all be same size: ", 0); }
  }
  if (ni>scrsz) { 
    if (scrsz>0) { free(scr); scr=(unsigned int *)NULL; }
    scrsz=ni+10000;
    scr=(unsigned int *)ecalloc(scrsz, sizeof(int));
  }
  for (i=0;i<ni;i++) { scr[i]=(int)ind[i]; // copy into integer array 
    if (scr[i]>=nx || scr[i]<0) { printf("fewind ERR1 %d %d\n",scr[i],nx);
    hoc_execerror("Index vector out-of-bounds", 0); }
  }
  for (j=0;j<num;j++) {
    for (i=0;i<ni;i++) x[i]=vvo[j][scr[i]];    
    for (i=0;i<ni;i++) vvo[j][i]=x[i];   
    list_vector_resize(ob, j, ni);
  }
  return (double)ni;
}
ENDVERBATIM


:* ind.findx(vecAlist,vecBlist) 
: uses ind as index into vecA's to load in vecB's (for select); a nondestructive fewind()
VERBATIM
static double findx (void* vv) {
  int i, j, ni, nx, av[VRRY], bv[VRRY], num;
  Object *ob1, *ob2;
  double *ind, *avo[VRRY], *bvo[VRRY];
  ni = vector_instance_px(vv, &ind);
  ob1 = *hoc_objgetarg(1);
  ob2 = *hoc_objgetarg(2);
  num = ivoc_list_count(ob1);
  i = ivoc_list_count(ob2);
  if (i!=num) hoc_execerror("findx ****ERRA****: lists have different counts", 0);
  if (num>VRRY) hoc_execerror("findx ****ERRB****: can only handle VRRY vectors", 0);
  for (i=0;i<num;i++) { 
    av[i]=list_vector_px(ob1, i, &avo[i]); // source vectors
    if (av[0]!=av[i]) { printf("findx ****ERRC**** %d %d %d\n",i,av[0],av[i]);
      hoc_execerror("Src vectors must all be same size: ", 0); }
  }
  nx=av[0]; // size of source vecs
  for (i=0;i<num;i++) { 
    bv[i]=list_vector_px2(ob2, i, &bvo[i], &vv); // dest vectors 
    if (vector_buffer_size(vv)<ni) { 
      printf("findx ****ERRD**** arg#%d need:%d sz:%d\n",num+i+1,ni,vector_buffer_size(vv));
      hoc_execerror("Destination vector with insufficient size: ", 0); 
    } else {
      vector_resize(vv, ni);
    }
  }
  if (ni>scrsz) { 
    if (scrsz>0) { free(scr); scr=(unsigned int *)NULL; }
    scrsz=ni+10000;
    scr=(unsigned int *)ecalloc(scrsz, sizeof(int));
  }
  for (i=0;i<ni;i++) { scr[i]=(int)ind[i]; // copy into integer array
    if (scr[i]>=nx || scr[i]<0) { printf("findx ****ERRE**** **** IND:%d SZ:%d\n",scr[i],nx);
      hoc_execerror("Index vector out-of-bounds", 0); }
  }
  for (j=0,i=0;j<num;j++) for (i=0;i<ni;i++) bvo[j][i]=avo[j][scr[i]];    
  return (double)ni;
}
ENDVERBATIM

:* ind.sindx(vecAlist,vecBlist)
: uses ind as index into vecA's to replace with elements from vecB's
VERBATIM
static double sindx (void* vv) {
  int i, j, ni, nx, av[VRRY], bv[VRRY], num;
  Object *ob1, *ob2;
  double *ind, *avo[VRRY], *bvo[VRRY];
  ni = vector_instance_px(vv, &ind);
  ob1 = *hoc_objgetarg(1);
  ob2 = *hoc_objgetarg(2);
  num = ivoc_list_count(ob1);
  i = ivoc_list_count(ob2);
  if (num!=i) hoc_execerror("sindx ****ERRA****: two vec lists have different counts", 0);
  if (num>VRRY) hoc_execerror("sindx ****ERRB****: can only handle VRRY vectors", 0);
  for (i=0;i<num;i++) { 
    av[i]=list_vector_px(ob1, i, &avo[i]); // dest vectors
    if (av[0]!=av[i]) { printf("sindx ****ERRC**** %d %d %d\n",i,av[0],av[i]);
      hoc_execerror("Dest. vectors must all be same size: ", 0); }
  }
  nx=av[0]; // size of dest vecs
  for (i=0;i<num;i++) { 
    bv[i]=list_vector_px(ob2, i, &bvo[i]); // source vectors 
  if (bv[i]!=ni) { 
    printf("sindx ****ERRD**** arg#%d does note match ind length %d vs %d\n",num+i+1,ni,bv[i]);
    hoc_execerror("Source vector with insufficient size: ", 0); }
  }
  if (ni>scrsz) { 
    if (scrsz>0) { free(scr); scr=(unsigned int *)NULL; }
    scrsz=ni+10000;
    scr=(unsigned int *)ecalloc(scrsz, sizeof(int));
  }
  for (i=0;i<ni;i++) { scr[i]=(int)ind[i]; // copy into integer array
    if (scr[i]>=nx || scr[i]<0) { 
      printf("sindx ****ERRE**** IND:%d SZ:%d\n",scr[i],nx);
      hoc_execerror("Index vector out-of-bounds", 0); }
  }
  for (j=0,i=0;j<num;j++) for (i=0;i<ni;i++) avo[j][scr[i]]=bvo[j][i];
  return (double)ni;
}
ENDVERBATIM

:* ind.sindv(vecAlist,valvec)
: uses ind as index into vecA's to replace with values in valvec
VERBATIM
static double sindv(void* vv) {
  int i, j, ni, nx, av[VRRY], bv, num;
  Object* ob;
  double *ind, *avo[VRRY], *bvo;
  ni = vector_instance_px(vv, &ind);
  ob = *hoc_objgetarg(1);
  bv=vector_arg_px(2, &bvo); // source vector
  num = ivoc_list_count(ob);
  if (num>VRRY) hoc_execerror("sindv ****ERRA****: can only handle VRRY vectors", 0);
  for (i=0;i<num;i++) { 
    av[i]=list_vector_px(ob, i, &avo[i]); // dest vectors
    if (av[0]!=av[i]) { printf("sindv ****ERRC**** %d %d %d\n",i,av[0],av[i]);
      hoc_execerror("Dest. vectors must all be same size: ", 0); }
  }
  nx=av[0]; // size of source vecs
  if (bv!=num) { 
    printf("sindv ****ERRD**** Vector arg does note match list count %d vs %d\n",num,bv);
    hoc_execerror("Source vector is wrong size: ", 0); }
  if (ni>scrsz) { 
    if (scrsz>0) { free(scr); scr=(unsigned int *)NULL; }
    scrsz=ni+10000;
    scr=(unsigned int *)ecalloc(scrsz, sizeof(int));
  }
  for (i=0;i<ni;i++) { scr[i]=(int)ind[i]; // copy into integer array
    if (scr[i]>=nx || scr[i]<0) { 
      printf("sindv ****ERRE**** IND:%d SZ:%d\n",scr[i],nx);
      hoc_execerror("Index vector out-of-bounds", 0); }
  }
  for (j=0,i=0;j<num;j++) for (i=0;i<ni;i++) avo[j][scr[i]]=bvo[j];
  return (double)ni;
}
ENDVERBATIM

:* ind.slct(key,args,veclist)
: picks out indices of numbers in key from multiple vectors
VERBATIM
static double slct (void* vv) {
  int i, j, k, m, n, p, ni, nk, na, nv[VRRY], num, fl, lc, field[VRRY];
  Object* lob;
  double *ind, *key, *arg, *vvo[VRRY], val;

  ni = vector_instance_px(vv, &ind); // vv is ind
  nk = vector_arg_px(1, &key);
  na = vector_arg_px(2, &arg);
  lob = *hoc_objgetarg(3);
  num = ivoc_list_count(lob);
  if (num>VRRY) hoc_execerror("ERR: vecst::slct can only handle VRRY vectors", 0);

  for (i=0,j=0,k=0;i<num;i++,j++) { 
    nv[i] = list_vector_px(lob, i, &vvo[i]);
    if (ni!=nv[i] && (key[j]!=EQW || k!=1)) {
      printf("vecst::slct ERR %d %d %d %d %d\n",i,j,k,ni,nv[i]);
      hoc_execerror("index and searched vectors must all be same size: ", 0); 
    }
    if (key[j]==EQW || key[j]==EQV) if (k==0){j--;k++;} else k=0; // EQW,EQV take 2 vector args 
  }
  for (j=0;j<nk;j++) { // look for fields in a mkcode() coded double
    field[j]=-1;
    if (key[j]<=EBE && key[j]>=ALL) { field[j]=0;
    } else for (m=1;m<=5;m++) {
      if (key[j]<=EBE*(m+1) && key[j]>=ALL*(m+1)) { // m is is field key 1-5
        key[j]/=(m+1);
        field[j]=m;
      }
    }
    if (field[j]==-1) {printf("vecst::slct ERRF %d %g\n",j,key[j]); hxe(); }
  }
  if (2*nk!=na) { printf("vecst::slct ERR3 %d %d\n",nk,na); 
    hoc_execerror("Arg vector must be double key length",0); }
  for (i=0,n=0;i<nk;i++) if (key[i]==EQV || key[i]==EQW) n++; // special cases take 2 vec args
  if (nk+n!=num) { 
    printf("vecst::slct ERR2 %d(keys)+%d(EQV/W)!=%d(vecs)\n",nk,n,num); 
    hoc_execerror("Key length must be number of vecs + num of EQV/W",0); }
  for (j=0,k=0,m=0;j<ni;j++) { // j steps through elements of vectors
    for (i=0,m=0,n=0,fl=1;i<num;i++,n++,m+=2) { // i steps thru key, m thru args
      if (field[n]==0) val=vvo[i][j]; else UNCODE(vvo[i][j],field[n],val);
      if (key[n]==ALL) continue; // OK - do nothing
      if (key[n]==NOZ)        { if (val==0.) {fl=0; break;} else continue; 
      } else if (key[n]==POS) { if (val<=0.) {fl=0; break;} else continue; 
      } else if (key[n]==NEG) { if (val>=0.) {fl=0; break;} else continue; 
      } else if (key[n]==GTH) { if (val<=arg[m]) {fl=0; break;} else continue; 
      } else if (key[n]==GTE) { if (val< arg[m]) {fl=0; break;} else continue; 
      } else if (key[n]==LTH) { if (val>=arg[m]) {fl=0; break;} else continue; 
      } else if (key[n]==LTE) { if (val> arg[m]) {fl=0; break;} else continue; 
      } else if (key[n]==EQU) { if (val!=arg[m]) {fl=0; break;} else continue; 
      } else if (key[n]==EQV) { if (val!=vvo[i+1][j]) { 
          fl=0; break;} else { i++; continue; }
      } else if (key[n]==EQW) {  // check value against values in following vec
        fl=0; // assume it's not going to match
        for (p=0;p<nv[i+1];p++) if (val==vvo[i+1][p]) {fl=1; break;}
        if (fl==0) break; else { i++; continue; }
      } else if (key[n]==NEQ) { if (val==arg[m]) {fl=0; break;} else continue; 
      } else if (key[n]==IBE) { if ((val< arg[m])||(val>=arg[m+1])) {
          fl=0; break; } else continue;  // IBE="[)" include-bracket-exclude
      } else if (key[n]==EBI) { if ((val<=arg[m])||(val> arg[m+1])) { 
          fl=0; break; } else continue;   // "(]" : exclude-bracket-include
      } else if (key[n]==IBI) { if ((val< arg[m])||(val> arg[m+1])) {
          fl=0; break; } else continue;   // "[]" : include-bracket-include
      } else if (key[n]==EBE) { if ((val<=arg[m])||(val>=arg[m+1])) {
          fl=0; break; } else continue;   // "()" : exclude-bracket-exclude
      } else {printf("vecst::slct ERR4 %g\n",key[n]); hoc_execerror("Unknown key",0);}
    }
    if (fl) ind[k++]=j; // all equal
  }
  vector_resize(vv, k);
  return (double)k;
}
ENDVERBATIM

:* ind.slor(key,args,vec1,vec2[,vec3,vec4,...])
: picks out indices of numbers in key from multiple vectors
VERBATIM
static double slor(void* vv) {
  int i, j, k, m, n, p, ni, nk, na, nv[VRRY], num, fl, field[VRRY];
  Object* lob;
  double *ind, *key, *arg, *vvo[VRRY], val;

  ni = vector_instance_px(vv, &ind); // vv is ind
  nk = vector_arg_px(1, &key);
  na = vector_arg_px(2, &arg);
  lob = *hoc_objgetarg(3);
  num = ivoc_list_count(lob);
  if (num>VRRY) hoc_execerror("ERR: vecst::slor can only handle VRRY vectors", 0);

  for (i=0,j=0,k=0;i<num;i++,j++) { 
    nv[i] = list_vector_px(lob, i, &vvo[i]);
    if (ni!=nv[i] && (key[j]!=EQW || k!=1)) {
      printf("vecst::slct ERR %d %d %d %d %d\n",i,j,k,ni,nv[i]);
      hoc_execerror("index and searched vectors must all be same size: ", 0); 
    }
    if (key[j]==EQW || key[j]==EQV) if (k==0){j--;k++;} else k=0; // k counts 2 vecs
  }
  for (j=0;j<num;j++) { // look for fields
    field[j]=-1;
    if (key[j]<=EBE && key[j]>=ALL) { field[j]=0;
    } else for (m=1;m<=5;m++) {
      if (key[j]<=EBE*(m+1) && key[j]>=ALL*(m+1)) { // m is is field key 1-5
        key[j]/=(m+1);
        field[j]=m;
      }
    }
    if (field[j]==-1) {printf("vecst::slct ERRF %g\n",key[j]); hxe(); }
  }
  if (2*nk!=na) { printf("vecst::slor ERR3 %d %d\n",nk,na); 
    hoc_execerror("Arg vector must be double key length",0); }
  for (i=0,n=0;i<nk;i++) if (key[i]==EQV || key[i]==EQW) n++; // special case takes 2 vec args
  if (nk+n!=num) { 
    printf("vecst::slor ERR2 %d(keys)+%d(EQV)!=%d(vecs)\n",nk,n,num); 
    hoc_execerror("Key length must be number of vecs + num of EQV",0); }
  for (j=0,k=0,m=0;j<ni;j++) { // j steps through elements of vectors
    for (i=0,m=0,n=0,fl=0;i<num;i++,n++,m+=2) { // i steps thru key, m thru args
      if (field[n]==0) val=vvo[i][j]; else UNCODE(vvo[i][j],field[n],val);
      if (key[n]==ALL) {fl=1; break;} // OK - do nothing
      if (key[n]==NOZ)        { if (val==0.) continue; else {fl=1; break;} 
      } else if (key[n]==POS) { if (val<=0.) continue; else {fl=1; break;} 
      } else if (key[n]==NEG) { if (val>=0.) continue; else {fl=1; break;} 
      } else if (key[n]==GTH) { if (val<=arg[m]) continue; else {fl=1; break;} 
      } else if (key[n]==GTE) { if (val< arg[m]) continue; else {fl=1; break;} 
      } else if (key[n]==LTH) { if (val>=arg[m]) continue; else {fl=1; break;} 
      } else if (key[n]==LTE) { if (val> arg[m]) continue; else {fl=1; break;} 
      } else if (key[n]==EQU) { if (val!=arg[m]) continue; else {fl=1; break;} 
      } else if (key[n]==EQV) { if (val!=vvo[i+1][j]) continue; else { 
          i++; fl=1; break; }
      } else if (key[n]==EQW) {  // check value against values in following vec
        fl=0; // assume it's not going to match
        for (p=0;p<nv[i+1];p++) if (val==vvo[i+1][p]) {fl=1; break;}
        if (fl==1) break; else { i++; continue; }
      } else if (key[n]==NEQ) { if (val==arg[m]) continue; else {fl=1; break;} 
      } else if (key[n]==IBE) { if ((val< arg[m])||(val>=arg[m+1])) {
          continue; } else {fl=1; break;}  // IBE="[)" include-bracket-exclude
      } else if (key[n]==EBI) { if ((val<=arg[m])||(val> arg[m+1])) { 
          continue; } else {fl=1; break;}   // "(]" : exclude-bracket-include
      } else if (key[n]==IBI) { if ((val< arg[m])||(val> arg[m+1])) {
          continue; } else {fl=1; break;}   // "[]" : include-bracket-include
      } else if (key[n]==EBE) { if ((val<=arg[m])||(val>=arg[m+1])) {
          continue; } else {fl=1; break;}   // "()" : exclude-bracket-exclude
      } else {printf("vecst::slor ERR4 %g\n",key[n]); hoc_execerror("Unknown key",0);}
    }
    if (fl) ind[k++]=j; // all equal
  }
  vector_resize(vv, k);
  return (double)k;
}
ENDVERBATIM

:* v.iwr()
VERBATIM
static double iwr(void* vv) {
  int i, j, nx;
  double *x;
  FILE* f = hoc_obj_file_arg(1);
  nx = vector_instance_px(vv, &x);
  if (nx>scrsz) { 
    if (scrsz>0) { free(scr); scr=(unsigned int *)NULL; }
    scrsz=nx+10000;
    scr=(unsigned int *)ecalloc(scrsz, sizeof(int));
  }
  for (i=0;i<nx;i++) scr[i]=(int)x[i]; // copy into integer array 
  fwrite(&nx,sizeof(int),1,f);  // write out the size
  fwrite(scr,sizeof(int),nx,f);
  return (double)nx;
}
ENDVERBATIM

:* v.ird()
VERBATIM
static double ird(void* vv) {
  int i, j, nx, n;
  double *x;
  FILE* f = hoc_obj_file_arg(1);
  nx = vector_instance_px(vv, &x);
  fread(&n,sizeof(int),1,f);  // size
  if (n>scrsz) { 
    if (scrsz>0) { free(scr); scr=(unsigned int *)NULL; }
    scrsz=n+10000;
    scr=(unsigned int *)ecalloc(scrsz, sizeof(int));
  }
  if (n!=nx) { 
    nx=vector_buffer_size(vv);
    if (n<=nx) {
      vector_resize(vv, n); nx=n; 
    } else {
      printf("%d > %d :: ",n,nx);
      hoc_execerror("Vector max capacity too small for ird ", 0);
    }
  }
  fread(scr,sizeof(int),n,f);
  for (i=0;i<nx;i++) x[i]=(double)scr[i];
  return (double)n;
}
ENDVERBATIM

:* v.fread2()
VERBATIM
static double fread2(void* vv) {
  int i, j, nx, n, type, maxsz;
  double *x;
  FILE* fp;
  BYTEHEADER

  fp = hoc_obj_file_arg(1);
  nx = vector_instance_px(vv, &x);
  maxsz=vector_buffer_size(vv);
  n = (int)*getarg(2);
  type = (int)*getarg(3);
  if (n>maxsz) {
    printf("%d > %d :: ",n,maxsz);
    hoc_execerror("Vector max capacity too small for fread2 ", 0);
  } else {
    vector_resize(vv, n);
  }
  if (type==6 || type==16) {         // unsigned ints
    unsigned int *xs;
    if (n>scrsz) { 
      if (scrsz>0) { free(scr); scr=(unsigned int *)NULL; }
      scrsz=n+10000;
      scr=(unsigned int *)ecalloc(scrsz, sizeof(int));
    }
    xs=(unsigned int*)scr;
    fread(xs,sizeof(int),n,fp);
    if (type==16) BYTESWAP_FLAG=1;
    for (i=0;i<n;i++) {
      BYTESWAP(scr[i],int)
      x[i]=(double)scr[i];
    }
    return (double)n;
  } if (type==3 || type==13) { // straight float reads
    float *xf = (float *)malloc(n * (unsigned)sizeof(float));
    fread(xf,sizeof(float),n,fp);
    if (type==13) BYTESWAP_FLAG=1;
    for (i=0;i<n;i++) {
      BYTESWAP(xf[i],float)
      x[i]=(double)xf[i];
    }
    free((char *)xf);    
  } else hoc_execerror("Type unsupported in fread2 ", 0);
  return 0.0;
}
ENDVERBATIM

:* ind.insct(v1,v2)
: return v1 intersect v2
VERBATIM
static double insct (void* vv) {
	int i, j, k, nx, nv1, nv2, maxsz;
	double *x, *v1, *v2;
	nx = vector_instance_px(vv, &x);
        maxsz=vector_buffer_size(vv);
        vector_resize(vv, maxsz);
	nv1 = vector_arg_px(1, &v1);
	nv2 = vector_arg_px(2, &v2);
        for (i=0,k=0;i<nv1;i++) for (j=0;j<nv2;j++) if (v1[i]==v2[j]) {
          if (k<maxsz) { x[k++]=v1[i]; } else {k++;}}  // v1[i] found in both vectors 
        if (k>maxsz) { 
          printf("\tinsct WARNING: ran out of room: %d<%d\n",maxsz,k);
        } else { vector_resize(vv, k); }
	return (double)k;
}
ENDVERBATIM

:* vdest.vfill(vsrc)
: fill vdest with multiple instances of vsrc until reach size
VERBATIM
static double vfill (void* vv) {
	int i, nx, nv1;
	double *x, *v1;
	nx = vector_instance_px(vv, &x);
	nv1 = vector_arg_px(1, &v1);
        for (i=0;i<nx;i++) x[i]=v1[i%nv1];
  return 0.0;
}
ENDVERBATIM

:* vec.cull(src,key)
: remove numbers in vec that are found in the key
VERBATIM
static double cull (void* vv) {
	int i, j, k, nx, nv1, nv2, maxsz, flag;
	double *x, *v1, *v2;
	nx = vector_instance_px(vv, &x);
        maxsz=vector_buffer_size(vv);
        vector_resize(vv, maxsz);
	nv1 = vector_arg_px(1, &v1);
	nv2 = vector_arg_px(2, &v2);
        for (i=0,k=0;i<nv1;i++) {
          flag=1;
          for (j=0;j<nv2;j++) if (v1[i]==v2[j]) flag=0;
          if (flag) {if (k<maxsz) { x[k++]=v1[i]; } else { k++; }}
        }
        if (k>maxsz) { 
          printf("\tcull WARNING: ran out of room: %d<%d\n",maxsz,k);
        } else { vector_resize(vv, k); }
	return (double)k;
}
ENDVERBATIM

:* dest.redundout(src[,INDFLAG])
: flag redundant numbers; must sort src first
: with indflag set just returns indices of locations rather than values
VERBATIM
static double redundout (void* vv) {
  int i, j, nx, nv1, maxsz, ncntr, indflag, cntflag;
  double *x, *v1, *cntr, val;
  void *vc;
  if (ifarg(2)) indflag=(int)*getarg(2); else indflag=0; 
  if (ifarg(3)) { cntflag=1; 
    ncntr = vector_arg_px(3, &cntr);  vc=vector_arg(3); 
    ncntr = vector_buffer_size(vc);   vector_resize(vc, ncntr);
    for (i=0;i<ncntr;i++) cntr[i]=1.; // will be at least 1 of each #
  } else cntflag=0;
  nx = vector_instance_px(vv, &x);
  maxsz=vector_buffer_size(vv);  vector_resize(vv, maxsz);
  nv1 = vector_arg_px(1, &v1);
  val=v1[0]; x[0]=(indflag?0:val);
  if (cntflag) {
    for (j=1,i=1;i<nv1&&j<maxsz&&j<ncntr;i++) {
      if (v1[i]!=val) { val=v1[i]; x[j++]=(indflag?i:val); } else cntr[j-1]+=1;
    }
  } else {
    for (j=1,i=1;i<nv1&&j<maxsz;i++) if (v1[i]!=val) { val=v1[i]; x[j++]=(indflag?i:val); }
  }
  if (j>=maxsz) { 
    printf("\tredundout WARNING: ran out of room: %d<needed\n",maxsz);
  } else { vector_resize(vv, j); }
  if (cntflag) if (j>=ncntr) { 
    printf("\tredundout WARNING: cntr ran out of room: %d<needed\n",ncntr);
  } else { vector_resize(vc, j); }
  return (double)j;
}
ENDVERBATIM

:* ind.mredundout(veclistA[,INDFLAG,veclistB])
: check redundancy across multiple parallel vectors veclistA
: with indflag 1, returns index of matchs in ind but does not alter vecs in veclistA
: will also remove from veclistB in parallel 
: leaves the last of a series of matches 
: with indflag set just returns indices of locations rather than values
: NB using indices will to .remove match items will give results that differ from 
: direct use (last vs first of a series -- will be seen in the other columns -- ie veclistB)
VERBATIM
static double mredundout (void* vv) {
  int i, j, k, m, p, q, maxsz, ns, nx, av[VRRY], bv[VRRY], num, numb, indflag, match;
  Object *ob, *ob2;
  double *x, *avo[VRRY], *bvo[VRRY], val[VRRY];
  void *vva[VRRY],*vvb[VRRY];
  nx = vector_instance_px(vv, &x);
  ob = *hoc_objgetarg(1);
  if (ifarg(2)) indflag=(int)*getarg(2); else indflag=0; 
  if (ifarg(3)) { 
    ob2 = *hoc_objgetarg(3);
    numb = ivoc_list_count(ob2);
  } else numb=0;
  maxsz=vector_buffer_size(vv);
  if (indflag) vector_resize(vv, maxsz); // else vector is not used
  num = ivoc_list_count(ob);
  if (num>VRRY) hoc_execerror("mredundout ****ERRA****: can only handle VRRY vectors", 0);
  for (i=0;i<num;i++) { 
    av[i]=list_vector_px2(ob, i, &avo[i], &vva[i]);
    if (av[0]!=av[i]) { printf("mredundout ****ERRC**** %d %d %d\n",i,av[0],av[i]);
      hoc_execerror("Vectors must all be same size: ", 0); }
  }
  ns=av[0]; // size of source vecs
  for (i=0;i<numb;i++) { 
    bv[i]=list_vector_px2(ob2, i, &bvo[i], &vvb[i]);
    if (ns!=bv[i]) { printf("mredundout ****ERRC2**** %d %d %d\n",i,ns,bv[i]);
      hoc_execerror("Vectors must all be same size: ", 0); }
  }
  if (ns/4>scrsz) { 
    if (scrsz>0) { free(scr); scr=(unsigned int *)NULL; }
    scrsz=ns/4+10000;
    scr=(unsigned int *)ecalloc(scrsz, sizeof(int));
  }
  for (j=0;j<num;j++) val[j]=avo[j][0]; // initialize the val array
  for (i=1,k=0;i<ns;i++) { 
    for (j=0,match=1;j<num;j++) {
      if (val[j]!=avo[j][i]) { match=0; break; } // if no match say so
    }
    if (match) { // add this one to the list
      if (k>=scrsz){printf("mredundout****ERRD**** over scr size %d\n",k);hxe();}
      scr[k++]=i; // flag to get rid of this one
    } else for (j=0;j<num;j++) val[j]=avo[j][i]; // copy next set of vals
  }
  if (indflag) { // just fill ind with indices of the repeats
    if (k>maxsz){printf("mredundout****ERRE**** vec overflow %d>%d\n",k,maxsz);hxe();}
    for (i=0;i<k;i++) x[i]=(double)scr[i];
    vector_resize(vv, k);
  } else { // remove all the repeat rows
    if (k == 0) return (double)k;
    for (i=0,p=scr[0]; i<k-1; i++) { // iter thru the inds to remove
      for (m=scr[i],p--; m<scr[i+1]; m++,p++) { // move everything down till next ind
        for (j=0;j<num; j++) avo[j][p]=avo[j][m]; // go through all the A list vecs
        for (j=0;j<numb;j++) bvo[j][p]=bvo[j][m]; // go through B list vecs
      }
    }
    for (m=scr[i],p--; m<ns; m++,p++) { // finish up by moving down from last ind to end
      for (j=0;j<num; j++) avo[j][p]=avo[j][m]; 
      for (j=0;j<numb;j++) bvo[j][p]=bvo[j][m]; 
    }
    for (j=0;j<num; j++) vector_resize(vva[j], ns-k); // resize all the vectors
    for (j=0;j<numb;j++) vector_resize(vvb[j], ns-k);
  }
  return (double)k;
}
ENDVERBATIM


:* PIVOTA.join(PIVOTB,VECLISTA,VECLISTB)
:  NB must sort both pivots before use
VERBATIM
static double join (void* vv) {
  int i, j, k, m, p, q, maxsz, npiva, npivb, av[VRRY], bv[VRRY], num, numb, indflag, match;
  Object *ob, *ob2;
  double *piva, *pivb, *avo[VRRY], *bvo[VRRY], val[VRRY];
  void *vva[VRRY],*vvb[VRRY];
  npiva = vector_instance_px(vv, &piva);
  npivb = vector_arg_px(1, &pivb);
  ob = *hoc_objgetarg(2);
  ob2 = *hoc_objgetarg(3);
  num = ivoc_list_count(ob);
  numb = ivoc_list_count(ob2);
  if (num>VRRY) hoc_execerror("join ****ERRA****: can only handle VRRY vectors", 0);
  if (num!=numb) hoc_execerror("join ****ERRB****: different #of vecs in lists", 0);
  for (i=0;i<num;i++) { 
    av[i]=list_vector_px2(ob, i, &avo[i], &vva[i]);
    if (av[0]!=av[i]) { printf("join ****ERRC**** %d %d %d\n",i,av[0],av[i]);
      hoc_execerror("Vectors must all be same size: ", 0); }
  }
  for (i=0;i<numb;i++) { 
    bv[i]=list_vector_px2(ob2, i, &bvo[i], &vvb[i]);
    if (bv[0]!=bv[i]) { printf("join ****ERRC2**** %d %d %d\n",i,bv[0],bv[i]);
      hoc_execerror("Vectors must all be same size: ", 0); }
  }
  for (i=0,j=0;i<npiva;i++) {
    for (;piva[i]!=pivb[j] && j<npivb;j++); // move forward to find matching pivb[j]
    if (j==npivb) { printf("%g not found in PivotB\n",piva[i]);hxe();}
    for (k=0;k<num;k++) avo[k][i]=bvo[k][j];
  }
  return (double)k;
}
ENDVERBATIM

:* vscl() will scale a vector to -1,1
VERBATIM
static double vscl(double *x, double n) {
  int i; 
  double max,min,r,sf,b,a;
  max=-1e9; min=1e9; a=-1; b=1; 
  for (i=0;i<n;i++) { 
    if (x[i]>max) max=x[i];
    if (x[i]<min) min=x[i];
  }
  r=max-min;  // range
  sf = (b-a)/r; // scaling factor
  for (i=0;i<n;i++) x[i]=(x[i]-min)*sf+a;
  return 0.0;
}
ENDVERBATIM

:* vdest.scl(vsrc) nondestructive scaling
VERBATIM
static double scl(void* vv) {
  int i, j, k, nx, nsrc, nfilt, ntmp;
  double *x, *src, *filt, *tmp, sum, lpad, rpad;
  nx = vector_instance_px(vv, &x);
  nsrc = vector_arg_px(1, &src);
  if (nx!=nsrc) { hoc_execerror("scl:Vectors not same size: ", 0); }
  for (i=0;i<nx;i++) x[i]=src[i];
  vscl(x,nx);
  return 0.0;
}
ENDVERBATIM

:* vdest.sccvlv(vsrc,vfilt,tmp) // NOT CORRECT -- see scxing()
: scaled convolution -- scale section to -1,1 before multiplying by filter
VERBATIM
static double sccvlv(void* vv) {
  int i, j, k, nx, nsrc, nfilt, ntmp;
  double *x, *src, *filt, *tmp, sum, lpad, rpad;
  nx = vector_instance_px(vv, &x);
  nsrc = vector_arg_px(1, &src);
  nfilt = vector_arg_px(2, &filt);
  ntmp = vector_arg_px(3, &tmp);
  if (nx!=nsrc) { hoc_execerror("sccvlv:Vectors not same size: ", 0); }
  if (nfilt>nsrc) { hoc_execerror("sccvlv:Filter bigger than source ", 0); }
  if (nfilt!=ntmp){hoc_execerror("sccvlv:Filter (arg2) and tmp vector (arg3) diff size", 0);}
  for (i=0;i<nx;i++) {
    x[i]=0.0;
    for (j=0,k=i-(int)(nfilt/2);j<nfilt&&k>0&&k<nsrc;j++,k++) tmp[j]=src[k];
    vscl(tmp,j-1);
    for (k=0;k<j;k++) x[i]+=filt[k]*tmp[k];
  }
  return 0.0;
}
ENDVERBATIM

:* vdest.scxing(vsrc,tmp)
: scaled xing -- scale each section to -1,1 before checking # of 0-xing's
VERBATIM
static double scxing (void* vv) {
  int i, j, k, nx, nsrc, f, ntmp, maxsz;
  double *x, *src, *filt, *tmp, sum, maxsum, th;
  nx = vector_instance_px(vv, &x);
  nsrc = vector_arg_px(1, &src);
  ntmp = vector_arg_px(2, &tmp);
  if (nx!=nsrc) { hoc_execerror("scxing:Vectors not same size: ", 0); }
  th=0.0; 
  maxsum=-1e9;
  for (i=0;i<nx;i++) x[i]=0.; // clear
  for (i=ntmp/2+1;i<nx-ntmp/2-1;i++) {
    for (j=0,k=i-(int)(ntmp/2);j<ntmp;k++,j++) tmp[j]=src[k];
    vscl(tmp,j-1); // scale from -1 to 1
    for (k=0,f=0,sum=0.; k<nsrc; k++) {
      if (tmp[k]>th) { // ? passing thresh
        if (f==0) { 
          sum+=1; // just count the 0xing
          f=1; 
        }
      } else {       // below thresh 
        if (f==1) {
          f=0; // just passed going down 
          sum+=1;
        }
      }
    }
    if (sum>maxsum) maxsum=sum;
    x[i]=sum;
  }
  return (double)maxsum;
}
ENDVERBATIM

:* vdest.cvlv(vsrc,vfilt)
: convolution
VERBATIM
static double cvlv (void* vv) {
  int i, j, k, nx, nsrc, nfilt;
  double *x, *src, *filt, sum, lpad, rpad;
  nx = vector_instance_px(vv, &x);
  nsrc = vector_arg_px(1, &src);
  nfilt = vector_arg_px(2, &filt);
  if (nx!=nsrc) { hoc_execerror("Vectors not same size: ", 0); }
  if (nfilt>nsrc) { hoc_execerror("Filter bigger than source ", 0); }
  for (i=0;i<nx;i++) {
    x[i]=0.0;
    for (j=0,k=i-(int)(nfilt/2);j<nfilt;j++,k++) {
      if (k>0 && k<nsrc-1) x[i]+=filt[j]*src[k];
    }
  }
  return 0.0;
}
ENDVERBATIM

:* vdest.intrp(flag)
: interpolate numbers replacing numbers given as flag
VERBATIM
static double intrp(void* vv) {
  int i, la, lb, nx;
  double *x, fl, a, b;
  nx = vector_instance_px(vv, &x);
  fl = *getarg(1);
  i=0; a=x[0]; la=0;
  if (a==fl) a=0;
  while (i<nx-1) {
    for (i=la+1;x[i]==fl && i<nx-1; i++) ; // find the next one 
    b=x[i]; lb=i;
    for (i=la+1; i<lb; i++) x[i]= a + (b-a)/(lb-la)*(i-la);
    a=b; la=lb;
  }
  return (double)fl;
}
ENDVERBATIM

:* vec.sumabs()
: return sum of abs values
VERBATIM
static double sumabs(void* vv) {
  int i, nx;
  double *x, sum;
  nx = vector_instance_px(vv, &x);
  for (sum=0,i=0;i<nx; i++) sum+=fabs(x[i]);
  return sum;
}
ENDVERBATIM

:* vec.inv([NUMERATOR])
: element=NUMERATOR/element
VERBATIM
static double inv (void* vv) {
  int i, nx;
  double *x, nume;
  nx = vector_instance_px(vv, &x);
  if (ifarg(1)) nume=*getarg(1); else nume=1.; 
  for (i=0;i<nx; i++) x[i]=nume/x[i];
  return (double)i;
}
ENDVERBATIM

:* tmp.nind(ind,vec1,vec2[,vec3,vec4,...])
: picks out numbers not in ind from multiple vectors
: ind must be sorted
VERBATIM
static double nind(void* vv) {
	int i, j, k, m, nx, ni, nv[VRRY], num, c, last;
	double *x, *ind, *vvo[VRRY];
	nx = vector_instance_px(vv, &x);
        for (i=0;ifarg(i);i++);
	if (i>VRRY) hoc_execerror("ERR: nind can only handle VRRY vectors", 0);
	num = i-2; // number of vectors to be picked apart 
        for (i=0;i<num;i++) { 
          nv[i] = vector_arg_px(i+2, &vvo[i]);
          if (nx!=nv[i]) { printf("nind ERR %d %d %d\n",i,nx,nv[i]);
            hoc_execerror("Vectors must all be same size: ", 0); }
        }
        ni = vector_arg_px(1, &ind);
        c = nx-ni; // the elems indexed are to be eliminated 
        if (ni>scrsz) { 
          if (scrsz>0) { free(scr); scr=(unsigned int *)NULL; }
          scrsz=ni+10000;
          scr=(unsigned int *)ecalloc(scrsz, sizeof(int));
        }
        for (i=0,last=-1;i<ni;i++) { 
          scr[i]=(int)ind[i]; // copy into integer array 
          if (scr[i]<0 || scr[i]>=nx) hoc_execerror("nind(): Index out of bounds", 0);
          if (scr[i]<=last) hoc_execerror("nind(): indices should mono increase", 0);
          last=scr[i];
        }
        for (j=0;j<num;j++) { // each output vec 
          for (i=0,last=-1,m=0;i<ni;i++) { // the indices of ind 
            for (k=last+1;k<scr[i];k++) { x[m++]=vvo[j][k]; }
            last=scr[i];
          }
          for (k=last+1;k<nx;k++,m++) { x[m]=vvo[j][k]; }
          for (i=0;i<c;i++) vvo[j][i]=x[i];   
          vv=vector_arg(j+2); vector_resize(vv, c);
        }
	return c;
}
ENDVERBATIM

:* ind.keyind(key,vec1,vec2[,vec3,vec4,...])
: picks out indices of numbers in key from multiple vectors
VERBATIM
static double keyind(void* vv) {
  int i, j, k, ni, nk, nv[VRRY], num;
  double *ind, *key, *vvo[VRRY];
  ni = vector_instance_px(vv, &ind); // vv is ind
  for (i=0;ifarg(i);i++) {} i--; // drop back by one to get numarg()
  if (i>VRRY) hoc_execerror("ERR: keyind can only handle VRRY vectors", 0);
  num = i-1; // number of vectors to be picked apart 
  for (i=0;i<num;i++) { 
    nv[i] = vector_arg_px(i+2, &vvo[i]);
    if (ni!=nv[i]) { printf("keyind ERR %d %d %d\n",i,ni,nv[i]);
    hoc_execerror("Non-key vectors must be same size: ", 0); }
  }
  nk = vector_arg_px(1, &key);
  if (nk!=num) { printf("keyind ERR2 %d %d\n",nk,num); 
    hoc_execerror("Key length must be number of vecs",0); }
  k=0;
  for (j=0;j<ni;j++) { // j steps through elements of vectors
    for (i=0;i<nk;i++) { // i steps through the key
      if (key[i]==ALL) continue; // OK - do nothing
      if (key[i]==NOZ)        { if (vvo[i][j]==0.) break; else continue; 
      } else if (key[i]==POS) { if (vvo[i][j]<=0.) break; else continue; 
      } else if (key[i]==NEG) { if (vvo[i][j]>=0.) break; else continue; 
      } else if (key[i]!=vvo[i][j]) break; // default
    }
    if (i==nk) ind[k++]=j; // all equal
  }
  vector_resize(vv, k);
  return (double)k;
}
ENDVERBATIM
 
:* v1.thresh() threshold above and below thresh
VERBATIM
static double thresh(void* vv) {
  int i, nx, ny, cnt;
  double *x, *y, th;
  nx = vector_instance_px(vv, &x);
  cnt=0;
  if (hoc_is_object_arg(1)) { 
    ny = vector_arg_px(1, &y); th=0;
    if (nx!=ny) { hoc_execerror("Vector sizes don't match in thresh.", 0); }
    for (i=0; i<nx; i++) { if (x[i]>=y[i]) { x[i]= 1.; cnt++;} else { x[i]= BVBASE; } }
  } else { th = *getarg(1);
    for (i=0; i<nx; i++) { if (x[i] >= th) { x[i]= 1.; cnt++;} else { x[i]= BVBASE; } }
  }
  return cnt;
}
ENDVERBATIM
 
:* v1.triplet() return location of a triplet
VERBATIM
static double triplet(void* vv) {
  int i, nx;
  double *x, *y, a, b;
  nx = vector_instance_px(vv, &x);
  a = *getarg(1); b = *getarg(2);
  for (i=0; i<nx; i+=3) if (x[i]==a&&x[i+1]==b) break;
  if (i<nx) return (double)i; else return -1.;
}
ENDVERBATIM

:* state.onoff(volts,OBon,thresh,dur,refr)
:  looks at volts vector to decide if have reached threshold thresh
:  OBon takes account of burst dur and refractory period
VERBATIM
static double onoff(void* vv) {
  int i, j, n, nv, non, nt, nd, nr, num;
  double *st, *vol, *obon, *thr, *dur, *refr;
  n = vector_instance_px(vv, &st);
  nv   = vector_arg_px(1, &vol);
  non  = vector_arg_px(2, &obon);
  nt   = vector_arg_px(3, &thr);
  nd   = vector_arg_px(4, &dur);
  nr   = vector_arg_px(5, &refr);
  if (n!=nv||n!=non||n!=nt||n!=nd||n!=nr) {
    hoc_execerror("v.onoff: vectors not all same size", 0); }
  for (i=0,num=0;i<n;i++) {
    obon[i]--;
    if (obon[i]>0.) { st[i]=1.; continue; } // cell must fire 
    if (vol[i]>=thr[i] && obon[i]<= -refr[i]) { // past refractory period 
        st[i]=1.; obon[i]=dur[i]; num++;
    } else { st[i]= BVBASE; }
  }
  return (double)num;
}
ENDVERBATIM
 
:* vo.bpeval(outp,del)
:  service routine for back-prop: vo=outp*(1-outp)*del
VERBATIM
static double bpeval(void* vv) {
  int i, n, no, nd, flag=0;
  double add,div;
  double *vo, *outp, *del;
  n = vector_instance_px(vv, &vo);
  no   = vector_arg_px(1, &outp);
  nd   = vector_arg_px(2, &del);
  if (ifarg(3) && ifarg(4)) { add= *getarg(3); div= *getarg(4); flag=1;}
  if (n!=no||n!=nd) hoc_execerror("v.bpeval: vectors not all same size", 0);
  if (flag) {
    for (i=0;i<n;i++) vo[i]=((outp[i]+add)/div)*(1.-1.*((outp[i]+add)/div))*del[i];
  } else {
    for (i=0;i<n;i++) vo[i]=outp[i]*(1.-1.*outp[i])*del[i];
  }
  return 0.0;
}
ENDVERBATIM
 
:* v1.sedit() string edit
VERBATIM
static double sedit(void* vv) {
  int i, n, ni, f=0;
  double *x, *ind, th, val;
  Symbol* s; char *op;
  op = gargstr(1);
  n = vector_instance_px(vv, &x);
  sprintf(op,"hello world");
  return (double)n;
}
ENDVERBATIM

:* v1.w() a .where that sets elements in source vector
VERBATIM
static double w (void* vv) {
  int i, n, ni, f=0;
  double *x, *ind, th, val;
  Symbol* s; char *op;
  if (! ifarg(1)) { 
    printf("v1.w('op',thresh[,val,v2])\n"); 
    printf("  a .where that sets elements in v1 to val (default 0), if v2 => only look at these elements\n");
    printf("  'op'=='function name' is a .apply targeted by v2 called as func(x[i],thresh,val)\n");
    return -1.;
  }
  op = gargstr(1);
  n = vector_instance_px(vv, &x);
  th = *getarg(2);
  if (ifarg(3)) { val = *getarg(3); } else { val = 0.0; }
  if (ifarg(4)) {ni = vector_arg_px(4, &ind); f=1;} // just look at the spots indexed
  if (!strcmp(op,"==")) { 
    if (f==1) {for (i=0; i<ni;i++) {if (x[(int)ind[i]]==th) { x[(int)ind[i]] = val;}}
    } else {for (i=0; i<n; i++) {if (x[i]==th) { x[i] = val;}}}
  } else if (!strcmp(op,"!=")) {
    if (f==1) {for (i=0; i<ni;i++) {if (x[(int)ind[i]]!=th) { x[(int)ind[i]] = val;}}
    } else {for (i=0; i<n; i++) {if (x[i]!=th) { x[i] = val;}}}
  } else if (!strcmp(op,">")) {
    if (f==1) {for (i=0; i<ni;i++) {if (x[(int)ind[i]]>th) { x[(int)ind[i]] = val;}}
    } else {for (i=0; i<n; i++) {if (x[i]>th) { x[i] = val;}}}
  } else if (!strcmp(op,"<")) {
    if (f==1) {for (i=0; i<ni;i++) {if (x[(int)ind[i]]<th) { x[(int)ind[i]] = val;}}
    } else {for (i=0; i<n; i++) {if (x[i]<th) { x[i] = val;}}}
  } else if (!strcmp(op,">=")) {
    if (f==1) {for (i=0; i<ni;i++) {if (x[(int)ind[i]]>=th) { x[(int)ind[i]] = val;}}
    } else {for (i=0; i<n; i++) {if (x[i]>=th) { x[i] = val;}}}
  } else if (!strcmp(op,"<=")) {
    if (f==1) {for (i=0; i<ni;i++) {if (x[(int)ind[i]]<=th) { x[(int)ind[i]] = val;}}
    } else {for (i=0; i<n; i++) {if (x[i]<=th) { x[i] = val;}}}
  } else if ((s=hoc_lookup(op))) { // same as .apply but only does indexed ones
    if (f==1) {for (i=0; i<ni;i++) {
        hoc_pushx(x[(int)ind[i]]); hoc_pushx(th); hoc_pushx(val);
      x[(int)ind[i]]=hoc_call_func(s, 3);}
    } else {for (i=0; i<n;i++) {hoc_pushx(x[i]); hoc_pushx(th); hoc_pushx(val);
        x[i]=hoc_call_func(s, 3);}}
  }
  return (double)i;
}
ENDVERBATIM
 
:* ind.slone(SRC,VAL) select indices where ==VAL from sorted vector SRC
VERBATIM
static double slone (void* vv) {
  int i, j, n, ni, nsrc, maxsz;
  double *x, *src, val, max, min;
  n = vector_instance_px(vv, &x);
  maxsz=vector_buffer_size(vv);  vector_resize(vv, maxsz);
  nsrc = vector_arg_px(1, &src);
  val = *getarg(2);
  if (ifarg(3)) ni=(int)*getarg(3); else {
    min=src[0];  max=src[nsrc-1];
    ni=(int)(val-min)/(max-min)*(double)(nsrc-1); // where expect to find it
  }
  if (src[ni]<val) {
    for (i=ni;src[i]!=val&&i<nsrc;i++); 
  } else {
    for (i=ni;src[i]!=val&&i>=0;i--); 
    for (    ;src[i]==val&&i>=0;i--); 
    i++; // go back to the first one
  }
  for (j=0;src[i]==val && j<maxsz && i<nsrc;i++,j++) x[j]=i;
  if (j==maxsz) printf("vecst slone WARN: OOR %d %d\n",j,maxsz);
  vector_resize(vv, j);
  return (double)(i-1);
}
ENDVERBATIM
  
:* dest.xing(src,tvec,thresh) 
:  dest.xing(src,thresh)  -- returns indices, change into time by .mul(tstep)
:  dest.xing(src)  -- default thresh=0; returns indices
: a .where that looks for threshold crossings and then doesn't pick another till
: comes back down again; places values from tvec in dest; interpolates
VERBATIM
static double xing (void* vv) {
  int i, j, nsrc, ndest, ntvec, f, maxsz, tvf;
  double *src, *dest, *tvec, th;
  tvf=0;
  ndest = vector_instance_px(vv, &dest);
  nsrc = vector_arg_px(1, &src);
  if (ifarg(3)) {
    ntvec = vector_arg_px(2, &tvec);
    th = *getarg(3);
    tvf=1; // flag that tvec being used
  } else if (ifarg(2)) {
    th = *getarg(2);
  } else th=0.0; // default threshold
  maxsz=vector_buffer_size(vv);
  vector_resize(vv, maxsz);
  if (tvf && nsrc!=ntvec) hoc_execerror("v.xing: vectors not all same size", 0);
  for (i=0,f=0,j=0; i<nsrc; i++) {
    if (src[i]>th) { // ? passing thresh 
      if (f==0) { 
        if (j>=maxsz) {
          printf("(%d) :: ",maxsz);
          hoc_execerror("Dest vec too small in xing ", 0);
        }
        if (i>0) { // don't record if first value is above thresh 
          if (tvf) {
            dest[j++] = tvec[i-1] + (tvec[i]-tvec[i-1])*(th-src[i-1])/(src[i]-src[i-1]);
          } else {
            dest[j++] = (i-1) + (th-src[i-1])/(src[i]-src[i-1]);
          }
        }
        f=1; 
      }
    } else {       // below thresh 
      if (f==1) { f=0; } // just passed going down 
    }
  }
  vector_resize(vv, j);
  return (double)i;
}
ENDVERBATIM

:* src.updown(thresh,dlist)  -- returns indices, change into time by .mul(tstep)
:  dest.updown(src)  -- default thresh=0; returns indices
: a .where that looks for threshold crossings and then doesn't pick another till
: comes back down again; places values from tvec in dest; interpolates
VERBATIM
#define UDSL 50
#define UDNQ 10
static double updown (void* vv) {
  int i, k, m, n, nsrc, jj[UDSL], f[UDSL], lc, dsz[UDSL], nqsz[UDNQ], thsz, lc2, done, dbn;
  double *src, *tvec, *th, *dest[UDSL], *nq[UDNQ], *tmp, *dbx, lt, thdist;
  Object *ob, *ob2;
  void *vvd[UDSL], *vvth, *vnq[UDNQ];
  nsrc = vector_instance_px(vv, &src);
  thsz = vector_arg_px(1, &th);
  ob =  *hoc_objgetarg(2);
  ob2 = *hoc_objgetarg(3);
  tmp = (double *)ecalloc(nsrc, sizeof(double));  
  lc =  ivoc_list_count(ob);
  lc2 = ivoc_list_count(ob2);
  if (lc>UDSL) {printf("updown ERRF mismatch: max slice list:%d %d\n",UDSL,lc); hxe();}
  if (lc2>UDNQ) {printf("updown ERRG mismatch: max nq sz:%d %d\n",UDNQ,lc2); hxe();}
  // nq=new NQS("LOC","PEAK","WIDTH","BASE","HEIGHT","START","SLICES","SHARP","INDEX","FILE")
  if (lc2!=10)  {printf("updown ERRB mismatch: nq sz should be 10:%d\n",lc2); hxe();}
  if (nsrc<lc) {printf("updown ERRC mismatch: %d %d\n",lc,nsrc); hxe();}
  if (lc!=thsz) {printf("updown ERRA mismatch: %d %d\n",lc,thsz); hxe();}
  if (!ismono1(th,thsz,-1)) {printf("updown ERRD: not mono dec %g %d\n",th[0],thsz); hxe();}
  // thdist=(th[thsz-2]-th[thsz-1])/2; // NOT BEING USED: the smallest spike we will accept
  for (k=0;k<lc;k++) { 
    list_vector_px2(ob, k, &dest[k], &vvd[k]);
    dsz[k]=vector_buffer_size(vvd[k]);
    vector_resize(vvd[k], dsz[k]);
  }
  for (k=0;k<lc2;k++) { 
    list_vector_px2(ob2, k, &nq[k], &vnq[k]);
    nqsz[k]=vector_buffer_size(vnq[k]);
    vector_resize(vnq[k], nqsz[k]);
    if (nqsz[k]!=nqsz[0]) {printf("updown ERRE mismatch: %d %d %d\n",k,nqsz[0],nqsz[k]); hxe();}
  }

  // dest vectors will store crossing points and midpoints at each th[k] slice location
  // as triplets: up/max/down
  for (k=0; k<lc; k++) {     
    jj[k]=f[k]=0;
    for (i=0;i<nsrc && src[i]>th[k];i++) {} // make sure start below thresh
    for (; i<nsrc; i++) {
      if (src[i]>th[k]) { 
        if (f[k]==0) { // ? passing thresh 
          if (jj[k]>=dsz[k]){printf("(%d,%d,%d) :: ",k,jj[k],dsz[k]);
            hoc_execerror("Dest vec too small in updown ", 0); }
          dest[k][jj[k]++] = (i-1) + (th[k]-src[i-1])/(src[i]-src[i-1]);
          f[k]=1; 
          tmp[k]=-1e9; dest[k][jj[k]]=-1.;
        }
        if (f[k]==1 && src[i]>tmp[k]) { // use tmp[] temporarily
          tmp[k]=src[i]; // pick out max
          dest[k][jj[k]] = (double)i; // location of this peak
        }
      } else {       // below thresh 
        if (f[k]==1) {  // just passed going down 
          jj[k]++;
          dest[k][jj[k]++] = (i-1) + (src[i-1]-th[k])/(src[i-1]-src[i]);
          f[k]=0; 
        }
      }
    }
  }
  // truncate vectors to multiples of 3:
  for (k=0;k<lc;k++) vector_resize(vvd[k],(int)(floor((double)jj[k]/3.)*3.));
  for (i=0; i<nsrc; i++) tmp[i]=0.; // clear temp space
  // go through all the slices to find identical peaks and save widths and locations
  // tmp[] uses triplets centered around a location corresponding to a max loc in the
  // original vector; the widest flanks for each are then on either side of this loc
  for (k=0;k<lc;k++) { // need to go from top to bottom to widen flanks
    for (i=1;i<jj[k];i+=3) { // through centers
      m=(int)dest[k][i]; // hash: place center at location
      if (tmp[m-2]<0 || tmp[m-1]<0 || tmp[m+1]<0 || tmp[m+2]<0) continue; // too crowded
      tmp[m]--;  // how many slices have found this peak -- use negative
      tmp[m-1]=dest[k][i-1]; tmp[m+1]=dest[k][i+1]; // flanks
    }
  }
  // 1st of 2 loops through tmp[]
  // step through tmp[] looking for negatives which indicate the slice count and pick up 
  // flanks from these
  // nq=new NQS("LOC","PEAK","WIDTH","BASE","HEIGHT","START","SLICES","SHARP","INDEX","FILE")
  for (i=0,k=0; i<nsrc; i++) if (tmp[i]<0.) {
    if (k>=nqsz[0]) { printf("updown OOR: %d %d\n",k,nqsz[0]); hxe(); }
    nq[0][k]=(double)i; // approx location of the peak of the spike
    nq[2][k]=tmp[i+1]; // location of right side
    nq[5][k]=tmp[i-1]; // start of spike
    nq[6][k]=tmp[i]; // neg of # of slices
    k++;
  }
  if (DEBUG_VECST && ifarg(4)) { dbn = vector_arg_px(4, &dbx); // DEBUG
    if (dbn<nsrc) printf("Insufficient room in debug vec\n"); else {
      for (i=0;i<nsrc;i++) dbx[i]=tmp[i]; }}
  // use nq vecs to search through flanks comparing to neighboring centers
  // making sure that the flanks do not overlap the centers on LT or RT side
  // where they do overlap correct them by reference to original dest info
  for (i=0;i<k;i++) {
    if ((i-1)>0 && nq[5][i]<nq[0][i-1]) { // look at LT side
      if (DEBUG_VECST) printf("LT problem %d %g %g<%g\n",i,nq[0][i],nq[5][i],nq[0][i-1]);
      for (m=lc-1,done=0;m>=0 && !done;m--) { // go from bottom to top
        for (n=1;n<jj[m] && !done;n+=3) { // through centers
          if (floor(dest[m][n])==nq[0][i] && dest[m][n-1]>nq[0][i-1]) {
            // nq[6][i]=nq[5][i]; // overwrite #of slices with end of overlap
            nq[5][i]=dest[m][n-1]; nq[2][i]=dest[m][n+1]; done=1;
          }
        }
      }
    }
    if ((i+1)<k && nq[2][i]>nq[0][i+1]) { // look at RT side
      if (DEBUG_VECST) printf("RT problem %d %g %g>%g\n",i,nq[0][i],nq[2][i],nq[0][i+1]);
      for (m=lc-1,done=0;m>=0 && !done;m--) { // go from bottom to top
        for (n=1;n<jj[m] && !done;n+=3) { // through centers
          if (floor(dest[m][n])==nq[0][i] && dest[m][n+1]<nq[0][i+1]) {
            // nq[6][i]=nq[2][i]; // overwrite #of slices with end of overlap
            nq[5][i]=dest[m][n-1]; nq[2][i]=dest[m][n+1]; done=1;
          }
        }
      }        
    }
  }
  // 2nd loop through tmp[]
  // needed to split into 2 loops so that could check for overlaps and correct those
  // before filling in the rest of nq
  // nq=new NQS("LOC","PEAK","WIDTH","BASE","HEIGHT","START","SLICES","SHARP","INDEX","FILE")
  for (i=0,k=0; i<nsrc; i++) if (tmp[i]<0.) {
    // calculate base voltage lt as interpolated value on left side
    lt=src[(int)floor(nq[5][k])]+(nq[5][k]-floor(nq[5][k]))*\
      (src[(int)floor(nq[5][k]+1.)]-src[(int)floor(nq[5][k])]);
    nq[1][k]=src[i]; // peak voltage
    nq[2][k]-=nq[5][k]; // width
    nq[3][k]=lt; // base voltage
    nq[4][k]=src[i]-lt; // height
    nq[7][k]=(src[i]-src[i-(int)SHM_UPDOWN])-(src[i+(int)SHM_UPDOWN]-src[i]); // sharpness
    nq[8][k]=(double)k;  // index
    k++;
  }
  for (i=0;i<lc2;i++) vector_resize(vnq[i], k);
  free(tmp);
  return jj[0];
}
ENDVERBATIM

:* dest.snap(src,tvec,dt) 
: interpolate src with tvec to prior dt step, saves only highest value in each interval
: an .interpolate that doesn't loose spikes
VERBATIM
static double snap(void* vv) {
  int i, j, nsrc, ndest, ntvec, f, maxsz, size;
  double *src, *dest, *tvec, dtt, tstop, tt, val;
  ndest = vector_instance_px(vv, &dest);
  nsrc = vector_arg_px(1, &src);
  ntvec = vector_arg_px(2, &tvec);
  dtt = *getarg(3);
  maxsz=vector_buffer_size(vv);
  tstop = tvec[nsrc-1];
  size=(int)tstop/dtt;
  if (size>maxsz) { 
    printf("%g > %g\n",size,maxsz);
    hoc_execerror("v.snap: insufficient room in dest", 0); }
  vector_resize(vv, size);
  if (nsrc!=ntvec) hoc_execerror("v.snap: src and tvec not same size", 0);
  for (tt=0,i=0;i<size && tt<=tvec[0];i++,tt+=dtt) dest[i]=src[0];
  for (j=1, i--, tt-=dtt; i<size; i++, val=-1e9, tt+=dtt) {
    if (tvec[j]>tt) dest[i]=src[j-1]; else {
      for (;j<nsrc && tvec[j]<=tt;j++) if (src[j]>val) val=src[j];
      if (val==-1e9) printf("vecst:snap() internal ERROR\n");
      dest[i]=val;
    }
  }
  return (double)size;
}
ENDVERBATIM

:* v1.xzero() looks for zero crossings 
VERBATIM
static double xzero(void* vv) {
  int i, n, nv, up;
  double *x, *vc, th, cnt=0.;
  n = vector_instance_px(vv, &x);
  nv = vector_arg_px(1, &vc);
  if (ifarg(2)) { th = *getarg(2); } else { th=0.0;}
  if (vc[0]<th) up=0; else up=1;  // F or T 
  if (nv!=n) hoc_execerror("Vector size doesn't match.", 0);
  for (i=0; i<nv; i++) {
    x[i]=0.;
    if (up) { // look for passing down 
      if (vc[i]<th) { x[i]=-1; up=0; cnt++; }
    } else if (vc[i]>th) { up=x[i]=1; cnt++; }
  }
  return cnt;
}
ENDVERBATIM

:* v1.negwrap([FLAG]) wrap neg values to pos, FLAG==0 set them to 0, FLAG!=0 wrap
:  above FLAG
VERBATIM
static double negwrap(void* vv) {
  int i, n;
  double *x, cnt, sig;
  n = vector_instance_px(vv, &x);
  if (ifarg(1)) sig = (int)*getarg(1); else sig=1e9; // default: do wrap
  if (sig==0.) {
    for (i=0,cnt=0; i<n; i++) if (x[i]<0) { 
      x[i]=0.;
      cnt++; 
    }
  } else if (sig==1e9) {
    for (i=0,cnt=0; i<n; i++) if (x[i]<0) { 
      x[i]=-x[i];
      cnt++; 
    }
  } else {
    for (i=0,cnt=0; i<n; i++) if (x[i]<sig) { 
      x[i]=2*sig-x[i]; // sig+(sig-x[i]) wraps around sig
      cnt++; 
    }
  }
  return cnt;
}
ENDVERBATIM
  
:* v1.sw(FROM,TO) switchs all FROMs to TO
VERBATIM
static double sw(void* vv) {
  int i, n;
  double *x, fr, to;
  n = vector_instance_px(vv, &x);
  fr = *getarg(1);
  to = *getarg(2);
  for (i=0; i<n; i++) {
    if (x[i]==fr) { x[i] = to;}
  }
  return (double)n;
}
ENDVERBATIM
 
:* v.b2v(bytevec) copies from vector to bytevec
VERBATIM
static double b2v(void* vv) {
  int i, n, num;
  double *x; bvec* to; Object *ob;
  n = vector_instance_px(vv, &x);
  ob = *(hoc_objgetarg(1));
  to = (bvec*)ob->u.this_pointer; // doesn't check that this is actually a bvec
  if (to->size!=n) { hoc_execerror("Vector and bytevec sizes don't match.", 0); }
  for (i=0; i<n; i++) x[i] = (double)to->x[i];
  return (double)n;
}
ENDVERBATIM

:* v.v2d(&x) copies from vector to double area -- a seg error waiting to happen
VERBATIM
static double v2d(void* vv) {
  int i, n, num;
  double *x, *to;
  n = vector_instance_px(vv, &x);
  to = hoc_pgetarg(1);
  if (ifarg(2)) { num = *getarg(2); } else { num=-1;}
  if (num>-1 && num!=n) { hoc_execerror("Vector size doesn't match.", 0); }
  for (i=0; i<n; i++) {to[i] = x[i];}
  return (double)n;
}
ENDVERBATIM
 
:* v.d2v(&x) copies from double area to vector -- a seg error waiting to happen
VERBATIM
static double d2v(void* vv) {
  int i, n, num;
  double *x, *fr;
  n = vector_instance_px(vv, &x);
  fr = hoc_pgetarg(1);
  if (ifarg(2)) { num = *getarg(2); } else { num=-1;}
  if (num>-1 && num!=n) { hoc_execerror("Vector size doesn't match.", 0); }
  for (i=0; i<n; i++) {x[i] = fr[i];}
  return (double)n;
}
ENDVERBATIM
 
:* v.lcat(LIST)
VERBATIM
static double lcat(void* vv) {
  int i, j, k, n, lc, cap, maxsz;
  Object *ob1;
  double *x, *fr; 
  void *vw;
  n = vector_instance_px(vv, &x);
  vector_resize(vv,maxsz=vector_buffer_size(vv)); // open it up fully
  ob1 = *hoc_objgetarg(1);
  lc = ivoc_list_count(ob1);
  for (i=0,j=0;i<lc && j<maxsz;i++) {
    cap = list_vector_px2(ob1, i, &fr, &vw);
    for (k=0;k<cap && j<maxsz;k++,j++) x[j]=fr[k];
  }
  if (i<lc || k<cap) printf("vecst lcat WARN: not all vecs copied\n");
  vector_resize(vv,j);
  return (double)j;
}
ENDVERBATIM

:* v.mkcode(LIST,BITS) -- put together integer vectors from list by bit concatenating
VERBATIM
static double mkcode(void* vv) {
  int i, j, k, n, num, bits;
  Object *ob;
  double *x, *vvo[5];
  n = vector_instance_px(vv, &x);
  ob = *hoc_objgetarg(1);
  if (ifarg(2)) bits = *getarg(2); else bits=3;
  num = ivoc_list_count(ob);
  if (num!=5) hoc_execerror("mkcode ****ERRA****: can only handle 5 vectors", 0);
  for (i=0;i<num;i++) {  j=list_vector_px(ob, i, &vvo[i]);
    if (n!=j) { printf("mkcode ****ERRC**** %d %d %d\n",i,n,j);
      hoc_execerror("Vectors must all be same size: ", 0); }}
  for (i=0;i<n;i++) { // go through the vec length
    for (j=0,x[i]=0;j<5;j++) {
      if (vvo[j][i]<0. || vvo[j][i]>=sc[4] || floor(vvo[j][i]+0.5)!=vvo[j][i]) {
        printf("vec.mkcode OOB %g>%g in vec[%d].x[%d]\n",vvo[j][i],sc[4],j,i); hxe(); }
        x[i]+=vvo[j][i]*sc[j+1];
    }
  }
  return (double)i;
}
ENDVERBATIM

:* v.uncode(val) -- take apart val and place in vector
: v.uncode(VECLIST) -- take apart vector items and place in vectors in list (cf uncodf)
: v.uncode(vec,field) -- take apart v entries and place requested field in vec
: v.uncode(field,val) -- replace field in v with val (cf recodf)
: v.uncode(field,vec) -- replace field in v with values from vec
VERBATIM
static double uncode (void* vv) {
  int i, j, n, ny, num, field;
  Object *ob;
  double *x, *y, *vvo[5], val, old;
  n = vector_instance_px(vv, &x);
  field=0;
  if (!ifarg(1)) { // numarg()==0
    printf("\tv.uncode(val) -- take apart val and place in vector\n\tv.uncode(VECLIST) -- take apart vector items and place in vectors in list (cf uncodf)\n\tv.uncode(vec,field) -- take apart vector items and place requested field in vector\n\tv.uncode(field,val) -- replace field in v with val (cf recodf)\n\tv.uncode(field,vec) -- replace field in v with values from vec\n"); return 0.;
  } else if (!ifarg(2)) { // numarg()==1
    if (hoc_is_double_arg(1)) {
      val = *getarg(1);
      if (vector_buffer_size(vv)<5) {
        hoc_execerror("uncode ****ERRA****: vector too small to resize(5)", 0);}
      vector_resize(vv,5);
      for (i=1;i<=5;i++) UNCODE(val,i,x[i-1])
      return x[0];
    } else {
      ob = *hoc_objgetarg(1);
      num = ivoc_list_count(ob);
      if (num!=5) hoc_execerror("uncode ****ERRA****: can only handle 5 vectors", 0);
      for (i=0;i<num;i++) {  j=list_vector_px(ob, i, &vvo[i]);
        if (n!=j) { printf("uncode ****ERRC**** %d %d %d\n",i,n,j);
          hoc_execerror("Vectors must all be same size: ", 0); }
      }
      for (i=0;i<n;i++) for (j=1;j<=5;j++) UNCODE(x[i],j,vvo[j-1][i]);
      return (double)i;
    }
  } else { // numarg()==2
    if (hoc_is_double_arg(1)) { // replace values
      field = (int)chkarg(1,1.,5.);
      ny=-1;
      if (hoc_is_double_arg(2)) {
        val=chkarg(2,0.,sc[4]-1); 
        if (floor(val+0.5)!=val) hoc_execerror("uncode(vec) ****ERRG****: non-int val", 0);
      } else {
        ny=vector_arg_px(2, &y);
        if (ny!=n) hoc_execerror("uncode(vec) ****ERRH****: diff sized vecs", 0);
      }
      for (i=0;i<n;i++) {
        UNCODE(x[i],field,old)
        if (ny>0)  {
          if (y[i]<0.||y[i]>=sc[4]||floor(y[i]+0.5)!=y[i]) {
            printf("vec.uncode ERRJ OOB %g (%g max) at %d\n",y[i],sc[4],i);hxe();}
          x[i] += sc[field]*(y[i]-old);
        } else {
          x[i] += sc[field]*(val -old);
        }
      }
      return (double)i;
    } else {  // fill single vector with values
      ny = vector_arg_px(1, &y);
      field = (int)chkarg(2,1.,5.);
      if (ny!=n) hoc_execerror("uncode(vec) ****ERRI****: diff sized vecs", 0);
      for (i=0;i<n;i++) UNCODE(x[i],field,y[i])
      return (double)i;
    }
  }
}
ENDVERBATIM
 
VERBATIM
//* list_vector_px(LIST,ITEM#,DOUBLE PTR ADDRESS) 
// modeled on vector_arg_px() picks up a vec from a list
int list_vector_px(Object *ob, int i, double** px) {
  Object* obv;
  int sz;
  obv = ivoc_list_item(ob, i);
  sz = vector_capacity(obv->u.this_pointer);
  *px = vector_vec(obv->u.this_pointer);
  return sz;
}

//* list_vector_px2(LIST,ITEM#,DOUBLE PTR ADDRESS,VEC POINTER ADDRESS) 
//  returns the vector pointer as well as the double pointer
int list_vector_px2 (Object *ob, int i, double** px, void** vv) {
  Object* obv;
  int sz;
  obv = ivoc_list_item(ob, i);
  sz = vector_capacity(obv->u.this_pointer);
  *px = vector_vec(obv->u.this_pointer);
  *vv = (void*) obv->u.this_pointer;
  return sz;
}

//* list_vector_px3(LIST,ITEM#,DOUBLE PTR ADDRESS,VEC POINTER ADDRESS) 
//  same as px2 but returns max vec size instead of current vecsize
//  side effect -- increase vector size to maxsize
int list_vector_px3 (Object *ob, int i, double** px, void** vv) {
  Object* obv;
  int sz;
  obv = ivoc_list_item(ob, i);
  sz = vector_buffer_size(obv->u.this_pointer);
  *px = vector_vec(obv->u.this_pointer);
  *vv = (void*) obv->u.this_pointer;
  vector_resize(*vv,sz);
  return sz;
}

//* list_vector_resize(LIST,ITEM#,NEW SIZE)
int list_vector_resize (Object *ob, int i, int sz) {
  Object* obv;
  int maxsz;
  obv = ivoc_list_item(ob, i);
  maxsz = vector_buffer_size(obv->u.this_pointer);
  if (sz>maxsz) {
    printf("max:%d request:%d ",maxsz,sz);
    hoc_execerror("Can't grow vector in list_vector_resize ", 0);
    return -1;
  }
  vector_resize(obv->u.this_pointer,sz);
  return sz;
}
ENDVERBATIM

:* v1.ismono([arg]) asks whether is monotonically increasing, with arg==-1 - decreasing
:  with arg==0:all same; 2:no consec ==; 3: incrementing by 1
VERBATIM
static int ismono1 (double *x, int n, int flag) {
  int i; double last;
  last=x[0]; 
  if (flag==1) {
    for (i=1; i<n && x[i]>=last; i++) last=x[i];
  } else if (flag==-1) {
    for (i=1; i<n && x[i]<=last; i++) last=x[i];
  } else if (flag==0) {
    for (i=1; i<n && x[i]==last; i++) ;
  } else if (flag==2) {
    for (i=1; i<n && x[i]>last; i++) last=x[i];
  } else if (flag==-2) {
    for (i=1; i<n && x[i]<last; i++) last=x[i];
  } else  if (flag==3) {
    for (i=1; i<n && x[i]==last+1; i++) last=x[i];
  } else  if (flag==-3) {
    for (i=1; i<n && x[i]==last-1; i++) last=x[i];
  }
  if (i==n) return 1; else return 0;
}

static double ismono(void* vv) {
  int i, n, flag;
  double *x,last;
  n = vector_instance_px(vv, &x);
  if (ifarg(1)) { flag = (int)*getarg(1); } else { flag = 1; }
  return (double)ismono1(x,n,flag);
}
ENDVERBATIM
 
:* v1.count(num) returns number of instances of num
VERBATIM
static double count(void* vv) {
  int i, n, cnt=0;
  double *x,num;
  n = vector_instance_px(vv, &x);
  num = *getarg(1);
  for (i=0; i<n; i++) if (x[i]==num) cnt++;
  return cnt;
}
ENDVERBATIM

:* v1.rnd() rounds off to nearest integer
VERBATIM
static double rnd(void* vv) {
  int i, n;
  double *x;
  n = vector_instance_px(vv, &x);
  for (i=0; i<n; i++) x[i]=floor(x[i]+0.5);
  return (double)i;
}
ENDVERBATIM

:* v1.pop() removes last entry and shortens vector
VERBATIM
static double pop(void* vv) {
  int n;
  double *x;
  n = vector_instance_px(vv, &x);
  if (n==0) {printf("vec.pop ERR: empty vec\n");hxe();}
  vector_resize(vv,n-1);
  return x[n-1];
}
ENDVERBATIM

PROCEDURE Expo (x) {
  TABLE RES FROM -20 TO 0 WITH 5000
  RES = exp(x)
}

FUNCTION AAA (x) {
  Expo(x)
  AAA = RES
}

:* dest.smgs(src,low,high,step,var)
:  rewrite of v.sumgauss() in nrn5.3::ivoc/ivocvect.cpp:1078
VERBATIM
static double smgs (void* vv) {	
  int i, j, nx, xv, nsum, points, maxsz;
  double *x, *sum;
  double  low , high , step , var , svar , scale , arg;

  nsum = vector_instance_px(vv, &sum);
  nx = vector_arg_px(1,&x);
  low = *getarg(2);
  high = *getarg(3);
  step = *getarg(4);
  var = *getarg(5);

  points = (int)((high-low)/step+hoc_epsilon);
  if (nsum!=points) { 
    maxsz=vector_buffer_size(vv);
    if (points<=maxsz) {
      nsum=points;  vector_resize(vv, nsum); 
    } else {
      printf("%d > %d :: ",points,maxsz);
      hoc_execerror("Vector max capacity too small in smgs ", 0);
    }
  }

  svar = -2.*var*var/step/step;
  scale = 1./sqrt(2.*M_PI)/var;

  for (j=0; j<points;j++) sum[j] = 0.;
  for (i=0;i<nx;i++) {
    xv = (int)((x[i]-low)/step + 0.5);
    for (j=xv; j<points && (arg=(j-xv)*(j-xv)/svar)>-20;j++) {
      Expo(arg);
      sum[j] += RES;
    }
    for (j=xv-1; j>=0 && (arg=(j-xv)*(j-xv)/svar)>-20;j--) {
      Expo(arg);
      sum[j] += RES;
    }
  }
  for (j=0; j<points;j++) sum[j] *= scale;
  return svar;
}
ENDVERBATIM

:* dest.smsy(tvec,CVLV_VEC,tstop[,dt,del])
:  sum CVLV_VEC starting at each point given in tvec with optional delay
:  used for summing up syn potentials
VERBATIM
static double smsy (void* vv) {	
  int i, j, k, nx, nc, nsum, points, maxsz;
  double *x, *sum, *c;
  double del,tstop,dtt;

  if (! ifarg(1)) { printf("dest.smsy(tvec,CVLV_VEC,tstop[,dt,del])\n"); return -1.; }

  del=0.; dtt=0.2;
  nsum = vector_instance_px(vv, &sum);
  nx = vector_arg_px(1,&x);
  nc = vector_arg_px(2,&c);
  tstop = *getarg(3);
  if (ifarg(4)) dtt = *getarg(4);
  if (ifarg(5)) del = *getarg(5);

  points=(int)(tstop/dtt+hoc_epsilon);
  if (nsum!=points) { 
    maxsz=vector_buffer_size(vv);
    if (points<=maxsz) {
      vector_resize(vv, points); points=nsum; 
    } else {
      printf("%d > %d :: ",points,maxsz);
      hoc_execerror("Dest vector too small in smsy ", 0);
    }
  }

  // don't zero out dest vec
  for (i=0;i<nx;i++) for (j=0,k=(x[i]+del)/dtt;j<nc && k<nsum;j++,k++) sum[k] += c[j];
  return points;
}
ENDVERBATIM

:* int.vrdh(FILE,veclist,code)
: vector read header will read the headers from vecs saved with vread()
: needs to be generalized so reads code as well, also should do BYTESWAP
VERBATIM 
static double vrdh (void* vv) {	
  int code, i, num, n[2], maxsz;
  double *x;
  FILE* f;

  num = vector_instance_px(vv, &x);
  maxsz=vector_buffer_size(vv);
  f =     hoc_obj_file_arg(1);
  num = (int)*getarg(2); // number of vectors to look for

  if (maxsz<2*num){printf("vrdh ERR0 need %d room in vec\n",2*num);hxe();}
  vector_resize(vv, 2*num);

  for (i=0;i<num;i++) { 
    fread(&n,sizeof(int),2,f); // n[1] is type
    if (n[1]!=3){printf("vrdh ERRA code 3 only implemented %d:%d\n",i,n[1]);hxe();}
    x[2*i]=(double)n[0]; // size
    x[2*i+1]=(double)n[1];
    fseek(f,(long)n[1],SEEK_CUR);
  }
  return (double)num;
}
ENDVERBATIM

:* rdmany(FILE,{veclist or vec},code[,num])
VERBATIM
static double rdmany (void* vv) {	
  int code, i, j, ni, vsz, ny, nv, num, cnt, n[2], sz, hd, vflag, iflag, last;
  Object* ob;
  double *vvo[100], sf[2], *ind, *y;
  FILE* f;

  vflag=iflag=0;
  ni = vector_instance_px(vv, &ind);
  f =     hoc_obj_file_arg(1);
  ob =   *hoc_objgetarg(2);
  if (ifarg(3)) cnt=(int)*getarg(3); else { cnt=ni; iflag=1; }
  if (strncmp(hoc_object_name(ob),"Vector",6)==0) vflag=1;
  i=2*sizeof(int) + 2*sizeof(double); // size of header with scaling
  j=2*sizeof(int);  // size of header without scaling
  fread(&n,sizeof(int),2,f);
  vsz=n[0]; code=n[1];
  fseek(f,(long)-2*sizeof(int),SEEK_CUR);  // go back
  if (DEBUG_VECST) printf("rdmanyDBA: %ld %d %d\n",ftell(f),vsz,code);
  switch (code) {
    // case 1:sz=1; hd=i; break; // char
    case 2:sz=2; hd=i; break; // short
    case 3:sz=4; hd=j; break; // float
    case 4:sz=8; hd=j; break; // double
    // case 5:sz=4; hd=i; break; // int
    default: hoc_execerror("rdmany ERRE: code not recognized", 0);
  }
  if (vflag) {
    ny = vector_arg_px(2, &y);
    num= cnt;
    if (vsz*cnt!=ny) {
      printf("rdmany ERRD: wrong size vec: %d statt (%d*%d) %d\n",ny,vsz,cnt,vsz*cnt); hxe();}
  } else {
    num = ivoc_list_count(ob);
    if (num>100) hoc_execerror("rdmany ERRA: can only handle 100 vectors", 0);
    if (num!=cnt) {printf("rdmany ERRB: %d != %d",num,cnt); hxe();}
    for (i=0;i<num;i++) { 
      nv = list_vector_px(ob, i, &vvo[i]);
      if (vsz!=nv){printf("rdmany ERRC: Vectors must all be same size %d %d %d\n",i,vsz,nv);hxe();}
    }
  }
  if (vsz*sz>bufsz) { 
    if (scrsz>0) { free(scr); scr=(unsigned int *)NULL; }
    scrsz=vsz+10;
    scr=(unsigned int *)ecalloc(scrsz, sz);
    bufsz=scrsz*sz; // number of chars available
  }
  if (code==2) {
    unsigned short *xs;
    xs=(unsigned short *)scr;
    for (last=-1,i=0;i<num;i++) {
      if (iflag) {  // iflag // "i flag" not "if lag"
        fseek(f,(long)((int)ind[i]-last-1)*(hd+vsz*sizeof(short)),SEEK_CUR); 
        if (DEBUG_VECST) printf("rdmanyDBB %ld ",ftell(f));
        last=(int)ind[i]; 
      }
      fread(&n,sizeof(int),2,f);
      fread(&sf,sizeof(double),2,f);
      if (n[0]!=vsz){printf("rdmany ERRA vec(%d) %d vs %d\n",iflag?(int)ind[i]:i,vsz,n[0]);hxe();}
      if (n[1]!=code){printf("rdmany ERRB code mismatch %d %d\n",n[1],code);hxe();}
      fread(xs,sizeof(short),n[0],f);
      for (j=0;j<vsz;j++) if (vflag) {
            y[i*vsz+j]=(double)(xs[j]/sf[0] + sf[1]);
      } else vvo[i][j]=(double)(xs[j]/sf[0] + sf[1]);
    }
  } else if (code==3) {
    float *xs;
    xs=(float *)scr;
    for (last=-1,i=0;i<num;i++) {
      if (iflag) { 
        fseek(f,(long)((int)ind[i]-last-1)*(hd+vsz*sizeof(float)),SEEK_CUR); 
        last=(int)ind[i]; 
      }
      if (DEBUG_VECST) printf("rdmanyDBC:%ld ",ftell(f));
      fread(&n,sizeof(int),2,f);
      if (n[0]!=vsz){printf("rdmany ERRA vec(%d) %d vs %d\n",iflag?(int)ind[i]:i,vsz,n[0]);hxe();}
      if (n[1]!=code){printf("rdmany ERRB code mismatch %d %d\n",n[1],code);hxe();}
      fread(xs,sizeof(float),n[0],f);
      for (j=0;j<n[0];j++) if (vflag) {
            y[i*vsz+j]=(double)xs[j];
      } else vvo[i][j]=(double)xs[j];
    }
  } else if (code==4) {
    double *xs;
    xs=(double *)scr;
    for (last=-1,i=0;i<num;i++) {
      if (iflag) { 
        fseek(f,(long)((int)ind[i]-last-1)*(hd+vsz*sizeof(double)),SEEK_CUR); 
        last=(int)ind[i]; 
      }
      if (DEBUG_VECST) printf("rdmanyDBD %ld ",ftell(f));
      fread(&n,sizeof(int),2,f);
      if (n[0]!=vsz){printf("rdmany ERRA vec(%d) %d vs %d\n",iflag?(int)ind[i]:i,vsz,n[0]);hxe();}
      if (n[1]!=code){printf("rdmany ERRB code mismatch %d %d\n",n[1],code);hxe();}
      fread(xs,sizeof(double),n[0],f); // should just read directly into final array
      for (j=0;j<n[0];j++) if (vflag) y[i*vsz+j]=xs[j]; else vvo[i][j]=xs[j];
    }
  } else printf("rdmany() code %d not implemented\n",code);
  return (double)num;
}
ENDVERBATIM

:* PROCEDURE install_vecst()
PROCEDURE install_vecst () {
  if (VECST_INSTALLED==1) {
    printf("$Id: vecst.mod,v 1.326 2006/06/24 23:38:36 billl Exp $\n")
  } else {
  VECST_INSTALLED=1
  VERBATIM {
  int i,j;
  install_vector_method("indset", indset);
  install_vector_method("thresh", thresh);
  install_vector_method("triplet", triplet);
  install_vector_method("onoff", onoff);
  install_vector_method("bpeval", bpeval);
  install_vector_method("w", w);
  install_vector_method("sedit", sedit);
  install_vector_method("xing", xing);
  install_vector_method("updown", updown);
  install_vector_method("scxing", scxing);
  install_vector_method("cvlv", cvlv);
  install_vector_method("sccvlv", sccvlv);
  install_vector_method("scl", scl);
  install_vector_method("intrp", intrp);
  install_vector_method("xzero", xzero);
  install_vector_method("negwrap", negwrap);
  install_vector_method("sw", sw);
  install_vector_method("ismono", ismono);
  install_vector_method("count", count);
  install_vector_method("rnd", rnd);
  install_vector_method("fewind", fewind);
  install_vector_method("findx", findx);
  install_vector_method("sindx", sindx);
  install_vector_method("sindv", sindv);
  install_vector_method("nind", nind);
  install_vector_method("keyind", keyind);
  install_vector_method("slct", slct);
  install_vector_method("slor", slor);
  install_vector_method("insct", insct);
  install_vector_method("cull", cull);
  install_vector_method("redundout", redundout);
  install_vector_method("mredundout", mredundout);
  install_vector_method("d2v", d2v);
  install_vector_method("v2d", v2d);
  install_vector_method("b2v", b2v);
  install_vector_method("iwr", iwr);
  install_vector_method("ird", ird);
  install_vector_method("smgs", smgs);
  install_vector_method("smsy", smsy);
  install_vector_method("ident", ident);
  install_vector_method("lcat", lcat);
  install_vector_method("snap", snap);
  install_vector_method("fread2", fread2);
  install_vector_method("vfill", vfill);
  install_vector_method("vrdh", vrdh);
  install_vector_method("mkcode", mkcode);
  install_vector_method("uncode", uncode);
  install_vector_method("sumabs", sumabs);
  install_vector_method("inv", inv);
  install_vector_method("join", join);
  install_vector_method("slone", slone);
  install_vector_method("pop", pop);
  install_vector_method("rdmany", rdmany);
  for (i=0,j=5;i<=5;i++,j--) sc[i]=pow(2,10*j);
  }
  ENDVERBATIM
  }
}

: isojt(OB1,EXAMPLE_OBJ) return whether OB1 is an instance of EXAMPLE_OBJ
FUNCTION isojt () {
  VERBATIM {
  Object *ob1, *ob2;
  ob1 = *hoc_objgetarg(1); ob2 = *hoc_objgetarg(2);
  if (!ob1) if (!ob2) return 1; else return 0;
#define ctemplate template
#ifdef NRN_VERSION_GTEQ_8_2_0
#if NRN_VERSION_GTEQ(9, 0, 0)
#undef ctemplate
#define ctemplate ctemplate
#endif
#endif
  if (!ob2 || ob1->ctemplate != ob2->ctemplate) {
#undef ctemplate
    return 0;
  }
  return 1;
  }
  ENDVERBATIM
}

: isojn(OB1,NAME) return whether OB1 is an instance of EXAMPLE_OBJ
FUNCTION isojn () {
  VERBATIM {
  Object *ob1; char* name;
  ob1 = *hoc_objgetarg(1); name = gargstr(2);
  if (strncmp(hoc_object_name(ob1),name,3)==0) _lisojn=1.; else _lisojn=0.;
  }
  ENDVERBATIM
}

: ojtnum(OBJ) returns object number
: returns internal number of object, eg if vec[3] is Vector[432] returns 432
FUNCTION ojtnum () {
  VERBATIM {
  Object *ob1; char name[50]; int ii;
  ob1 = *hoc_objgetarg(1);
  if (!ob1) return -1;
  if (ifarg(2)) {
    strncpy(name, hoc_object_name(ob1),50);
    for (ii=strlen(name);ii>1;ii--) if (name[ii]==91) {name[ii]=0; break;} // 91 is [
    hoc_assign_str(hoc_pgargstr(2),name);
  }
  return (double)ob1->index;
  }
  ENDVERBATIM
}

: eqojt(OB1,OB2) return whether OB1 and OB2 point to same object
FUNCTION eqojt () {
  VERBATIM {
  Object *ob1, *ob2;
  ob1 = *hoc_objgetarg(1); ob2 = *hoc_objgetarg(2);
  if (ob1 && ob2 && ob1==ob2) {
    return 1;
  }
  return 0;
  }
  ENDVERBATIM
}

:* byteswap(FILE)
FUNCTION byteswap () {
  VERBATIM {
  int n[2];
  double ret;
  FILE* f;
  BYTEHEADER

  f =     hoc_obj_file_arg(1);
  fread(&n,sizeof(int),2,f);
  if (n[1] < 1 || n[1] > 5) {
    BYTESWAP_FLAG = 1;
    ret = 1.; 
  } else ret = 0.;
  BYTESWAP(n[1],int)
  if (n[1] < 1 || n[1] > 5) { 
    printf("byteswap: Something wrong with location sampled: %d\n",n[1]);
    ret = -1.;
  }
  fseek(f,-2*sizeof(int),SEEK_CUR); // go back to where we started
  return ret;
  }
  ENDVERBATIM
}


: mkcodf(val1,val2,val3,val4,val5) stuff 5 vals<=999 into a single double
FUNCTION mkcodf () {
  VERBATIM {
  int i;
  double x,a;
  for (x=0.,i=1;i<=5;i++) { 
    a=chkarg(i,0.,sc[4]-1);
    if (floor(a+0.5)!=a) {
      printf("mkcodf restricted to integers %g\n",a);hxe();}
    x+=a*sc[i];
  }
  return x;
  }
  ENDVERBATIM
}

: uncodf(code,i) returns field i (1-5) from code
FUNCTION uncodf () {
  VERBATIM {
  int i;
  double x,ret, *ptr;
  x=*getarg(1);
  if (hoc_is_double_arg(2)) {
    i=(int)chkarg(2,1.,5.); 
    UNCODE(x,i,ret);
    return ret;
  } else {
    for (i=2;i<=6;i++) if (ifarg(i)) {
      ptr = hoc_pgetarg(i);
      UNCODE(x,i-1,*ptr);
    } else break;
    return *ptr;
  }
  }
  ENDVERBATIM
}

: recodf(i,code,new) replaces field i (1-5) from code with new
FUNCTION recodf () {
  VERBATIM {
  int i;
  double x, y, old;
  i=(int)chkarg(1,1.,5.); x=*getarg(2); y=chkarg(3,0.,sc[4]-1);
  UNCODE(x,i,old);
  return x + sc[i]*(y-old);
  }
  ENDVERBATIM
}

: flor(val)
FUNCTION flor () {
  VERBATIM {
  return floor(*getarg(1));
  }
  ENDVERBATIM
}
