#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "routines.h"
#include "constants.h"
#define _USE_MATH_DEFINES
/*double dist=10;
int NInt=128;
int kNInt=128;
double kmax=1.8576485321061220;
double kratio=4.0;
double kdist = 25.0;*/
unsigned short mem1 =500;
void compile_fourier_func(double *** fourier_funcP, double *** fourier_funcN, double *** g_funcs, double *** f_funcs, double *** mbessels, unsigned short grandspinnum, unsigned short mkNmax, unsigned short kNmax, unsigned short pkNmax);
void dataforG(unsigned short grandspinnum, double **Nmat, double **Pmat, double *Nvals, double *Pvals, double *kn, double *mkn, double *pkn);
void momentanums(unsigned short grandspinnnum, unsigned short **momenta);
void g_function(double ** g1p, double **g2p, double **g1m, double **g2m, unsigned short grandspinnum, double ** Nmatrix, double ** Pmatrix, unsigned short mkNmax, unsigned short kNmax, unsigned short pkNmax, double *kn, double *pkn, double *mkn, double **jcc,double **jcp,double **jcm,double **jpc,double **jpp,double **jmc,double **jmm);
void f_function(double ** f1p, double **f2p, double **f1m, double **f2m, unsigned short grandspinnum, double ** Nmatrix, double ** Pmatrix, unsigned short mkNmax, unsigned short kNmax, unsigned short pkNmax, double *kn, double *pkn, double *mkn, double **jcc,double **jcp,double **jcm,double **jpc,double **jpp,double **jmc,double **jmm);
void bessels(unsigned short grandspinnum, unsigned short *momentanums, double *pks, double *ks, double *mks, double **jcc, double **jcp, double **jcm, double **jpc, double **jpp, double **jmc,double **jmm);
void compile_regvecs(double **** vectorsP, double **** vectorsN, double *** g_funcs, double *** f_funcs, unsigned short grandspinnum, unsigned short mkNmax, unsigned short kNmax, unsigned short pkNmax);
void Gspin_momentanums(unsigned short *granspinnum, unsigned short **pmomenta);
double k_inner_product(double **vec1, double **vec2, unsigned short grandspinnum);
void renormalise(double **vec, unsigned short grandspinnum);

void renormalise(double **vec, unsigned short grandspinnum) {
  double factor;
  unsigned short i, j;
  factor = k_inner_product(vec,vec,grandspinnum);
  factor = 1.0/sqrt(factor);
  if (grandspinnum==0) {
    for (i=0; i<2; i++) {
      for (j=0; j<=kNInt; j++) {
        vec[i][j] = factor*vec[i][j];
      }
    }
  } else {
    for (i=0; i<4; i++) {
      for (j=0; j<=kNInt; j++) {
        vec[i][j] = factor*vec[i][j];
      }
    }
  }
}
void Gspin_momentanums(unsigned short *grandspinnum, unsigned short **pmomenta) {
  FILE *momentanumsfile;
  unsigned short index, *momenta, i;
  unsigned short amomenta[200];
	char filename[200];

	sprintf(filename, "Data/expnum%d/momentanum_%d.b", expnum, expnum);
  momentanumsfile = fopen(filename, "rb");
  index =0;
  while (fread(&amomenta[index], sizeof(unsigned short), 1, momentanumsfile)) {
    index++;
  }
  momenta = (unsigned short *) malloc(sizeof(unsigned short)*(index+1));
  for (i=0; i<index; i++) {
    momenta[i] = amomenta[i];
  }
  momenta[index]=0;
  *pmomenta=momenta;
  *grandspinnum = index-2;
  fclose(momentanumsfile);
}
double k_inner_product(double **vec1, double **vec2, unsigned short grandspinnum) {
  unsigned long i, j, n;
  double val, tval, k;
  double intarr[kNInt+1];
  val=0;

  if (grandspinnum==0) {
    for (i=0; i<2; i++) {
      for (n=0; n<=kNInt; n++) {
        k=(double)n/(double)kNInt;
        intarr[n]=kdist*kdist*kdist*k*k*vec1[i][n]*vec2[i][n];
      }
      val+=booles(intarr, kNInt)/(double)kNInt;
    }
  } else {
    for (i=0; i<4; i++) {
      for (n=0; n<=kNInt; n++) {
        k=(double)n/(double)kNInt;
        intarr[n]=kdist*kdist*kdist*k*k*vec1[i][n]*vec2[i][n];
      }
      val+=booles(intarr, kNInt)/(double)kNInt;
    }
  }

  return val;
}
double r_inner_product(double **vec1, double **vec2, unsigned short grandspintot) {
   unsigned long i, j, n;
   double val, tval, r;
   double *intarr;
   intarr=(double*)malloc(sizeof(double)*(NInt+1));
   val=0;
   if (grandspintot==0) {
     for (i=0; i<2; i++) {
       for (n=0; n<=NInt; n++) {
         r=(double)n/(double)NInt;
         intarr[n]=r*r*vec1[i][n]*vec2[i][n];
       }
       val += booles(intarr,NInt)/(double)NInt;
     }
   } else {
     for (i=0; i<4; i++) {
       for (n=0; n<=NInt; n++) {
         r=(double)n/(double)NInt;
         intarr[n]=r*r*vec1[i][n]*vec2[i][n];
       }
       val += booles(intarr,NInt)/(double)NInt;
     }
   }

   free(intarr);
   return val;
}
void compile_regvecs(double **** pvectorsP, double **** pvectorsN, double *** g_funcs, double *** f_funcs, unsigned short grandspinnum, unsigned short mkNmax, unsigned short kNmax, unsigned short pkNmax){
//RE-ORGANISES THE ARRAYS TO A SENSICAL ORDER
  double ***vectorsP, ***vectorsN;
  unsigned long i, j, k;
  double *intarr, r, val;
  if (grandspinnum!=0) {
    vectorsP = (double***)malloc(sizeof(double**)*4*kNmax);
    for (i=0; i<4*kNmax; i++) {
      vectorsP[i]=(double**)malloc(sizeof(double*)*4);
      for (j=0; j<4; j++) {
        vectorsP[i][j]=(double*)malloc(sizeof(double)*(NInt+1));
      }
      for (k=0; k<=NInt; k++) {
        vectorsP[i][0][k] = g_funcs[0][i][k];
        vectorsP[i][1][k] = g_funcs[1][i][k];
        vectorsP[i][2][k] = f_funcs[0][i][k];
        vectorsP[i][3][k] = f_funcs[1][i][k];
      }
    }
    vectorsN = (double***)malloc(sizeof(double**)*(2*pkNmax+2*mkNmax));
    for (i=0; i<2*pkNmax+2*mkNmax; i++) {
      vectorsN[i]=(double**)malloc(sizeof(double*)*4);
      for (j=0; j<4; j++) {
        vectorsN[i][j]=(double*)malloc(sizeof(double)*(NInt+1));
      }
      for (k=0; k<=NInt; k++) {
        vectorsN[i][0][k] = g_funcs[2][i][k];
        vectorsN[i][1][k] = g_funcs[3][i][k];
        vectorsN[i][2][k] = f_funcs[2][i][k];
        vectorsN[i][3][k] = f_funcs[3][i][k];
      }
    }
  } else {
    vectorsP = (double***)malloc(sizeof(double**)*2*kNmax);
    for (i=0; i<2*kNmax; i++) {
      vectorsP[i]=(double**)malloc(sizeof(double*)*2);
      for (j=0; j<2; j++) {
        vectorsP[i][j]=(double*)malloc(sizeof(double)*(NInt+1));
      }
      for (k=0; k<=NInt; k++) {
        vectorsP[i][0][k] = g_funcs[0][i][k];
        vectorsP[i][1][k] = f_funcs[0][i][k];
      }
    }
    vectorsN = (double***)malloc(sizeof(double**)*(2*pkNmax));
    for (i=0; i<2*pkNmax; i++) {
      vectorsN[i]=(double**)malloc(sizeof(double*)*2);
      for (j=0; j<2; j++) {
        vectorsN[i][j]=(double*)malloc(sizeof(double)*(NInt+1));
      }
      for (k=0; k<=NInt; k++) {
        vectorsN[i][0][k] = g_funcs[2][i][k];
        vectorsN[i][1][k] = f_funcs[2][i][k];
      }
    }
  }
  *pvectorsP=vectorsP;
  *pvectorsN=vectorsN;
}
void compile_fourier_func(double *** fourier_funcP, double *** fourier_funcN, double *** g_funcs, double *** f_funcs, double ***mbessels, unsigned short grandspinnum, unsigned short mkNmax, unsigned short kNmax, unsigned short pkNmax) {
  double k, x;// *** fourier_funcP, *** fourier_funcN, k, x;
  unsigned short G;
  unsigned long i, j, n, temp;
  double intarrk[kNInt+1], intarr[NInt+1], val;
  G=grandspinnum;

  if (grandspinnum==0) {
    temp = 2*kNmax;
  } else {
    temp = 4*kNmax;
  }
  for (n=0; n<temp; n++) {
    for (i=0; i<=kNInt; i++) {
      k = (double)i/(double)kNInt;
      for (j=0; j<=NInt; j++) {
        x =(double)j/(double)NInt;
        intarr[j]= x*x*mbessels[j][i][G]*g_funcs[0][n][j];
      }
      fourier_funcP[n][0][i]=sqrt(2.0/M_PI)*dist*sqrt(dist)*booles(intarr, NInt)/(double)NInt;
    }
    if (grandspinnum==0) {
      for (i=0; i<=kNInt; i++) {
        k = (double)i/(double)kNInt;
        for (j=0; j<=NInt; j++) {
          x =(double)j/(double)NInt;
          intarr[j]= x*x*mbessels[j][i][G+1]*f_funcs[0][n][j];
        }
        fourier_funcP[n][1][i]=sqrt(2.0/M_PI)*dist*sqrt(dist)*booles(intarr, NInt)/(double)NInt;
      }
    } else {
      for (i=0; i<=kNInt; i++) {
        k = (double)i/(double)kNInt;
        for (j=0; j<=NInt; j++) {
          x =(double)j/(double)NInt;
          intarr[j]= x*x*mbessels[j][i][G]*g_funcs[1][n][j];
        }
        fourier_funcP[n][1][i]=sqrt(2.0/M_PI)*dist*sqrt(dist)*booles(intarr, NInt)/(double)NInt;
      }
      for (i=0; i<=kNInt; i++) {
        k = (double)i/(double)kNInt;
        for (j=0; j<=NInt; j++) {
          x =(double)j/(double)NInt;
          intarr[j]= x*x*mbessels[j][i][G+1]*f_funcs[0][n][j];
        }
        fourier_funcP[n][2][i]=sqrt(2.0/M_PI)*dist*sqrt(dist)*booles(intarr, NInt)/(double)NInt;
      }
      for (i=0; i<=kNInt; i++) {
        k = (double)i/(double)kNInt;
        for (j=0; j<=NInt; j++) {
          x =(double)j/(double)NInt;
          intarr[j]= x*x*mbessels[j][i][G-1]*f_funcs[1][n][j];
        }
        fourier_funcP[n][3][i]=sqrt(2.0/M_PI)*dist*sqrt(dist)*booles(intarr, NInt)/(double)NInt;
      }
    }
  }
  val=0;
  if (grandspinnum==0) {
    temp = pkNmax*2;
  } else {
    temp = pkNmax*2 + mkNmax*2;
  }
  for (n=0; n<temp; n++) {
    for (i=0; i<=kNInt; i++) {
      k = (double)i/(double)kNInt;
      for (j=0; j<=NInt; j++) {
        x =(double)j/(double)NInt;
        intarr[j]= x*x*mbessels[j][i][G+1]*g_funcs[2][n][j];
      }
      fourier_funcN[n][0][i]=sqrt(2.0/M_PI)*dist*sqrt(dist)*booles(intarr, NInt)/(double)NInt;
    }
    if (grandspinnum==0) {
      for (i=0; i<=kNInt; i++) {
        k = (double)i/(double)kNInt;
        for (j=0; j<=NInt; j++) {
          x =(double)j/(double)NInt;
          intarr[j]= x*x*mbessels[j][i][G]*f_funcs[2][n][j];
        }
        fourier_funcN[n][1][i]=sqrt(2.0/M_PI)*dist*sqrt(dist)*booles(intarr, NInt)/(double)NInt;
      }
    } else {
      for (i=0; i<=kNInt; i++) {
        k = (double)i/(double)kNInt;
        for (j=0; j<=NInt; j++) {
          x =(double)j/(double)NInt;
          intarr[j]= x*x*mbessels[j][i][G-1]*g_funcs[3][n][j];
        }
        fourier_funcN[n][1][i]=sqrt(2.0/M_PI)*dist*sqrt(dist)*booles(intarr, NInt)/(double)NInt;
      }
      for (i=0; i<=kNInt; i++) {
        k = (double)i/(double)kNInt;
        for (j=0; j<=NInt; j++) {
          x =(double)j/(double)NInt;
          intarr[j]= x*x*mbessels[j][i][G]*f_funcs[2][n][j];
        }
        fourier_funcN[n][2][i]=sqrt(2.0/M_PI)*dist*sqrt(dist)*booles(intarr, NInt)/(double)NInt;
      }
      for (i=0; i<=kNInt; i++) {
        k = (double)i/(double)kNInt;
        for (j=0; j<=NInt; j++) {
          x =(double)j/(double)NInt;
          intarr[j]= x*x*mbessels[j][i][G]*f_funcs[3][n][j];
        }
        fourier_funcN[n][3][i]=sqrt(2.0/M_PI)*dist*sqrt(dist)*booles(intarr, NInt)/(double)NInt;
      }
    }
  }
  val=0;
}
void bessels(unsigned short grandspinnum, unsigned short *momentanums, double *pks, double *ks, double *mks, double **jcc, double **jcp,double **jcm,double **jpc,double **jpp,double **jmc,double **jmm) {
  unsigned short pkNmax, kNmax, mkNmax, G;
  unsigned long i, j;
  double x;

  pkNmax = momentanums[grandspinnum+1];
  kNmax = momentanums[grandspinnum];
  mkNmax = 0;
  if (grandspinnum != 0) {
    mkNmax = momentanums[grandspinnum-1];
  }
  G = grandspinnum;

  for (i=0; i<kNmax; i++) {
    for (j=0; j<=NInt; j++) {
      x=(double)j/(double)NInt;
      jcc[i][j] = Jn(ks[i]*dist*x, G);
      jpc[i][j] = Jn(ks[i]*dist*x, G+1);
      jcc[2*kNmax-1-i][j]=jcc[i][j];
      jpc[2*kNmax-1-i][j]=jpc[i][j];
    }
  }
  for (i=0; i<pkNmax; i++) {
    for (j=0; j<=NInt; j++) {
      x = (double)j/(double)NInt;
      jcp[i][j] = Jn(pks[i]*dist*x, G);
      jpp[i][j] = Jn(pks[i]*dist*x, G+1);
      jcp[2*pkNmax -1-i][j] = jcp[i][j];
      jpp[2*pkNmax -1-i][j] = jpp[i][j];
    }
  }
  if (G != 0) {
    for (i=0; i<kNmax; i++) {
      for (j=0; j<=NInt; j++) {
        x = (double)j/(double)NInt;
        jmc[i][j] = Jn(ks[i]*dist*x, G-1);
        jmc[2*kNmax-1-i][j]=jmc[i][j];
      }
    }
    for (i=0; i<mkNmax; i++) {
      for (j=0; j<=NInt; j++) {
        x = (double)j/(double)NInt;
        jmm[i][j] = Jn(mks[i]*dist*x, G-1);
        jcm[i][j] = Jn(mks[i]*dist*x, G);
        jmm[2*mkNmax-1-i][j]=jmm[i][j];
        jcm[2*mkNmax-1-i][j]=jcm[i][j];
      }
    }
  }
}
void g_function(double **g1p, double **g2p, double **g1m, double **g2m, unsigned short grandspinnum, double ** Nmatrix, double ** Pmatrix, unsigned short mkNmax, unsigned short kNmax, unsigned short pkNmax, double *ks, double *pks, double *mks, double **jcc,double **jcp,double **jcm,double **jpc,double **jpp,double **jmc,double **jmm) {
  unsigned long i, j, k, n, temp1;
  double N, wp, something;//, **g1p, **g2p, **g1m, **g2m;
  double pEs[2*pkNmax], Es[2*kNmax], *mEs, x, wp1, wp2, N1, N2, *temp, *intarr, r;

  mEs = (double*)malloc(2*mkNmax*sizeof(double));
  for (i=0; i<pkNmax; i++) {
    pEs[i]=-sqrt(pks[i]*pks[i]+1);
    pEs[2*pkNmax-1-i]=-pEs[i];
  }
  for (i=0; i<kNmax; i++) {
    Es[i]=-sqrt(ks[i]*ks[i]+1);
    Es[2*kNmax-1-i]=-Es[i];
  }
  if (grandspinnum!=0) {
    for (i=0; i<mkNmax; i++) {
      mEs[i]=-sqrt(mks[i]*mks[i]+1);
      mEs[2*mkNmax-1-i]=-mEs[i];
    }
  }
  if (grandspinnum==0) {
    temp1 = 2*kNmax;
  } else {
    temp1 = 4*kNmax;
  }
  for (i=0; i<temp1; i++) {
    for (j=0; j<=NInt; j++) {
      g1p[i][j]=0;
      if (grandspinnum!=0) {
        g2p[i][j]=0;
      }
      for (n=0; n<kNmax*2; n++) {
        N=1/fabs(jpc[n][NInt]);
        wp = sqrt(1+1/Es[n]);
        g1p[i][j]+=Pmatrix[i][n]*N*wp*jcc[n][j];
        if (grandspinnum!=0) {
          g2p[i][j]+=Pmatrix[i][n+2*kNmax]*N*wp*jcc[n][j];
        }
      }
    }
  }
  if (grandspinnum==0) {
    temp1 = pkNmax*2;
  } else {
    temp1 = (2*pkNmax + 2*mkNmax);
  }
  for (i=0; i<temp1; i++) {
    for (j=0; j<=NInt; j++) {
      g1m[i][j]=0;
      if (grandspinnum!=0) {
        g2m[i][j]=0;
      }
      for (n=0; n<pkNmax*2; n++) {
        N=1/fabs(jcp[n][NInt]);
        wp = sqrt(1+1/pEs[n]);
        g1m[i][j]+=Nmatrix[i][n]*N*wp*jpp[n][j];
      }
      if (grandspinnum!=0) {
        for (n=0; n<mkNmax*2; n++) {
          N=1/fabs(jcm[n][NInt]);
          wp = sqrt(1+1/mEs[n]);
          g2m[i][j]+=Nmatrix[i][2*pkNmax+n]*N*wp*jmm[n][j];
        }
      }

    }
  }
  free(mEs);
}
void f_function(double **f1p, double **f2p, double **f1m, double **f2m, unsigned short grandspinnum, double ** Nmatrix, double ** Pmatrix, unsigned short mkNmax, unsigned short kNmax, unsigned short pkNmax, double *ks, double *pks, double *mks, double **jcc,double **jcp,double **jcm,double **jpc,double **jpp,double **jmc,double **jmm) {
  unsigned long i, j, k, n, temp1;
  double N, wm, N1, N2, wm1, wm2;//, **f1p, **f2p, **f1m, **f2m;
  double pEs[2*pkNmax], Es[2*kNmax], *mEs, *temp, x, something, *intarr, r;


  for (i=0; i<pkNmax; i++) {
    pEs[i]=-sqrt(pks[i]*pks[i]+1);
    pEs[2*pkNmax-1-i]=-pEs[i];
  }
  for (i=0; i<kNmax; i++) {
    Es[i]=-sqrt(ks[i]*ks[i]+1);
    Es[2*kNmax-1-i]=-Es[i];
  }
  if (grandspinnum!=0) {
    mEs = (double*)malloc(2*mkNmax*sizeof(double));
    for (i=0; i<mkNmax; i++) {
      mEs[i]=-sqrt(mks[i]*mks[i]+1);
      mEs[2*mkNmax-1-i]=-mEs[i];
    }
  }
  if (grandspinnum==0) {
    temp1 = kNmax*2;
  } else {
    temp1 = kNmax*4;
  }
  for (i=0; i<temp1; i++) {
    for (j=0; j<=NInt; j++) {
      f1p[i][j]=0;
      if (grandspinnum!=0){
        f2p[i][j]=0;
      }
      for (n=0; n<kNmax*2; n++) {
        N=1/fabs(jpc[n][NInt]);
        wm = (Es[n]/fabs(Es[n]))*sqrt(1-1/Es[n]);
        f1p[i][j]+=Pmatrix[i][n]*N*wm*jpc[n][j];
        if (grandspinnum!=0) {
          f2p[i][j]+=Pmatrix[i][n+2*kNmax]*N*wm*jmc[n][j];
        }
      }
    }
  }
  if (grandspinnum==0) {
    temp1 = 2*pkNmax;
  } else {
    temp1 = 2*pkNmax + 2*mkNmax;
  }
  for (i=0; i<temp1; i++) {
    for (j=0; j<=NInt; j++) {
      f1m[i][j]=0;
      if (grandspinnum!=0) {
        f2m[i][j]=0;
      }
      for (n=0; n<pkNmax*2; n++) {
        N=1/fabs(jcp[n][NInt]);
        wm = (pEs[n]/fabs(pEs[n]))*sqrt(1-1/pEs[n]);
        f1m[i][j]+=Nmatrix[i][n]*N*wm*jcp[n][j];
      }
      if (grandspinnum!=0) {
        for (n=0; n<mkNmax*2; n++) {
          N=1/fabs(jcm[n][NInt]);
          wm = (mEs[n]/fabs(mEs[n]))*sqrt(1-1/mEs[n]);
          f2m[i][j]+=Nmatrix[i][2*pkNmax+n]*N*wm*jcm[n][j];
        }
      }
    }
  }
  if (grandspinnum!=0) {
    free(mEs);
  }
}
void momentanums(unsigned short grandspinnum, unsigned short **momenta) {
  FILE * momentafile;
  unsigned short * momentanum, i;
  momentanum = (unsigned short *)malloc(sizeof(unsigned short)*(grandspinnum+2));
  momentafile = fopen("Data/expnum%d/momentanum_%d.b", "rb");

  for (i=0; i<=grandspinnum+1; i++) {
    fread(&momentanum[i], sizeof(unsigned short), 1, momentafile);
  }
  fclose(momentafile);
  *momenta = momentanum;
}
void dataforG(unsigned short grandspinnum, double **Nmat, double **Pmat, double *Nvals, double *Pvals, double *kn, double *mkn, double *pkn) {
  unsigned short * momentanum, i, j, lenP, lenN, plenk, lenk, mlenk;
  unsigned long seekto1, seekto2;
  double temp;
  FILE * momentafile, * eigvalfile, * eigvecfile, *knfile;
	char filename[100];

  momentanum = (unsigned short *)malloc(sizeof(unsigned short)*(grandspinnum+2));

	sprintf(filename, "Data/expnum%d/momentanum_%d.b", expnum, expnum);
  momentafile = fopen(filename, "rb");

	sprintf(filename, "Data/expnum%d/eigenergies_%d.b", expnum, expnum);
  eigvalfile = fopen(filename, "rb");

	sprintf(filename, "Data/expnum%d/eigvecs_%d.b", expnum, expnum);
  eigvecfile = fopen(filename, "rb");

	sprintf(filename, "Data/expnum%d/kn_%d.b", expnum, expnum);
  knfile = fopen(filename, "rb");

  for (i=0; i<=grandspinnum+1; i++) {
    fread(&momentanum[i], sizeof(unsigned short), 1, momentafile);
  }
  fclose(momentafile);
  lenk = momentanum[grandspinnum]*2;
  plenk = momentanum[grandspinnum+1]*2;
  mlenk = 0;
  if (grandspinnum!=0) {
    mlenk = momentanum[grandspinnum-1]*2;
  }
  if (grandspinnum==0) {
    lenP=momentanum[grandspinnum]*2;
    lenN=momentanum[grandspinnum+1]*2;
  } else {
    lenP= momentanum[grandspinnum]*4;
    lenN = momentanum[grandspinnum-1]*2 + momentanum[grandspinnum+1]*2;
  }
  seekto1 =0;
  for (i=0; i<=grandspinnum-2; i++) {
    seekto1 += momentanum[i]*2;
  }
  fseek(knfile, sizeof(double)*seekto1, SEEK_SET);
  if (grandspinnum!=0) {
    for (i=0; i<mlenk; i++) {
      fread(&mkn[i], sizeof(double), 1, knfile);
    }
  }
  for (i=0; i<lenk; i++) {
    fread(&kn[i], sizeof(double), 1, knfile);
  }
  for (i=0; i<plenk; i++) {
    fread(&pkn[i], sizeof(double), 1, knfile);
  }
  fclose(knfile);
  if (grandspinnum==0) {
    seekto1=0;
  } else {
    seekto1 =momentanum[0]*2+momentanum[1]*2;
  }
  for (i=1; i<=grandspinnum-1; i++) {
    seekto1 += momentanum[i]*4 + momentanum[i-1]*2 + momentanum[i+1]*2;
  }
  fseek(eigvalfile, sizeof(double)*seekto1, SEEK_SET);
  for (i=0; i<lenP; i++) {
    fread(&Pvals[i], sizeof(double), 1, eigvalfile);
  }
  for (i=0; i<lenN; i++) {
    fread(&Nvals[i], sizeof(double), 1, eigvalfile);
  }
  fclose(eigvalfile);
  if (grandspinnum==0) {
    seekto2= 0;
  } else {
    seekto2 = momentanum[0]*2*momentanum[0]*2 + momentanum[1]*2*momentanum[1]*2;
  }
  for (i=1; i<=grandspinnum-1; i++) {
    seekto2 += momentanum[i]*16*momentanum[i] + (momentanum[i-1]*2 + momentanum[i+1]*2)*(momentanum[i-1]*2 + momentanum[i+1]*2);
  }
  fseek(eigvecfile, sizeof(double)*seekto2, SEEK_SET);
  for (i=0; i<lenP; i++) {
    for (j=0; j<lenP; j++) {
      fread(&Pmat[i][j], sizeof(double), 1, eigvecfile);
    }
  }
  for (i=0; i<lenN; i++) {
    for (j=0; j<lenN; j++) {
      fread(&Nmat[i][j], sizeof(double), 1, eigvecfile);
    }
  }
  fclose(eigvecfile);
  free(momentanum);
}
int main(int argc, const char* argv[]) {
  FILE * fourierfile, *fouriertxt, *fourierfileF0;
  double **Nmatrix, **Pmatrix, *pks, *ks, *mks, *Pval, *Nval, ***fourier_funcP, ***fourier_funcN;
  double **jcc, **jpc, **jcp, **jpp, **jcm, **jmc, **jmm;
  unsigned short mkNmax, kNmax, pkNmax;
  double ***vectorsP, ***vectorsN;
  double **g1p, **g2p, **g1m, **g2m, **f1p, **f2p, **f1m, **f2m, ***g_funcs, ***f_funcs, ***mbessels;
  double val, x, r, c;
  double *temp, *intarr;
  unsigned short grandspinnum, *momenta, grandspintot, Gend, Gbegin;
  unsigned long i, j, k;
	char filename[100];
  grandspinnum=5;
  temp = (double*)malloc(sizeof(double)*(NInt+1));
  intarr = (double*)malloc(sizeof(double)*(NInt+1));

  fprintf(stderr, "kdist=%f\ndist=%f\nkratio=%f\nlambda=%f\ncqm=%f\nrep=%d\nexpnum=%d\n", kdist, dist, kratio, lambda, cqm, rep, expnum);
  Gspin_momentanums(&grandspintot, &momenta);
  fprintf(stderr, "grandspintot=%d\nkNInt=%lu\n", grandspintot, kNInt);
  Gbegin=0;
  Gend=1000;
  //return 0;
  g_funcs = (double***)malloc(sizeof(double**)*4);
  g_funcs[0]=(double**)malloc(mem1*sizeof(double*));
  g_funcs[1]=(double**)malloc(mem1*sizeof(double*));
  g_funcs[2]=(double**)malloc(mem1*sizeof(double*));
  g_funcs[3]=(double**)malloc(mem1*sizeof(double*));
  for (i=0; i<mem1; i++) {
    g_funcs[0][i]=(double*)malloc((NInt+1)*sizeof(double));
    g_funcs[1][i]=(double*)malloc((NInt+1)*sizeof(double));
    g_funcs[2][i]=(double*)malloc((NInt+1)*sizeof(double));
    g_funcs[3][i]=(double*)malloc((NInt+1)*sizeof(double));
  }
  f_funcs = (double***)malloc(sizeof(double**)*4);
  f_funcs[0]=(double**)malloc(mem1*sizeof(double*));
  f_funcs[1]=(double**)malloc(mem1*sizeof(double*));
  f_funcs[2]=(double**)malloc(mem1*sizeof(double*));
  f_funcs[3]=(double**)malloc(mem1*sizeof(double*));
  for (i=0; i<mem1; i++) {
    f_funcs[0][i]=(double*)malloc((NInt+1)*sizeof(double));
    f_funcs[1][i]=(double*)malloc((NInt+1)*sizeof(double));
    f_funcs[2][i]=(double*)malloc((NInt+1)*sizeof(double));
    f_funcs[3][i]=(double*)malloc((NInt+1)*sizeof(double));
  }

  fourier_funcP = (double***)malloc(sizeof(double**)*mem1);
  fourier_funcN = (double***)malloc(sizeof(double**)*mem1);
  for (j=0; j<mem1; j++) {
    fourier_funcP[j]=(double**)malloc(sizeof(double*)*4);
    fourier_funcN[j]=(double**)malloc(sizeof(double*)*4);
    for (i=0; i<4; i++) {
      fourier_funcP[j][i]=(double*)malloc(sizeof(double)*(kNInt+1));
      fourier_funcN[j][i]=(double*)malloc(sizeof(double)*(kNInt+1));
    }
  }

  Pmatrix=(double**)malloc(sizeof(double*)*mem1);
  Nmatrix=(double**)malloc(sizeof(double*)*mem1);
  for (i=0; i<mem1; i++) {
    Nmatrix[i]=(double*)malloc(sizeof(double)*mem1);
    Pmatrix[i]=(double*)malloc(sizeof(double)*mem1);
  }

  Nval = (double*)malloc(sizeof(double)*mem1);
  Pval = (double*)malloc(sizeof(double)*mem1);
  ks = (double*)malloc(sizeof(double)*mem1);
  mks= (double*)malloc(sizeof(double)*mem1);
  pks= (double*)malloc(sizeof(double)*mem1);

  jcc = (double**)malloc(250*sizeof(double*));
  jcp = (double**)malloc(250*sizeof(double*));
  jcm = (double**)malloc(250*sizeof(double*));
  jmc = (double**)malloc(250*sizeof(double*));
  jmm = (double**)malloc(250*sizeof(double*));
  jpc = (double**)malloc(250*sizeof(double*));
  jpp = (double**)malloc(250*sizeof(double*));
  for (i=0; i<250; i++) {
    jcc[i] = (double*)malloc((NInt+1)*sizeof(double));
    jcp[i] = (double*)malloc((NInt+1)*sizeof(double));
    jcm[i] = (double*)malloc((NInt+1)*sizeof(double));
    jmc[i] = (double*)malloc((NInt+1)*sizeof(double));
    jmm[i] = (double*)malloc((NInt+1)*sizeof(double));
    jpc[i] = (double*)malloc((NInt+1)*sizeof(double));
    jpp[i] = (double*)malloc((NInt+1)*sizeof(double));
  }
  mbessels = (double***)malloc(sizeof(double**)*(NInt+1));
  for (i=0; i<=NInt; i++) {
    x = i/(double)NInt * dist;
    mbessels[i] =(double**)malloc(sizeof(double*)*(kNInt+1));
    for (j=0; j<=kNInt; j++) {
      r = j/(double)kNInt * kdist;
      mbessels[i][j]=(double*)malloc(sizeof(double)*(grandspintot+2));
      Jnarr(mbessels[i][j], x*r, grandspintot+1);
    }
  }
  if (Gend>grandspintot) {
    Gend=grandspintot;
  }
  for (grandspinnum=Gbegin; grandspinnum<=Gend; grandspinnum++) {
		if (grandspinnum!=0) {
    	mkNmax = momenta[grandspinnum-1];
		}
		sprintf(filename, "Data/expnum%d/Fdata/fourierfileG%d_%d.b", expnum, grandspinnum, expnum);
  	fourierfile = fopen(filename, "wb");
		sprintf(filename, "Data/expnum%d/Fdata/fourierfileF0G%d_%d.b", expnum, grandspinnum, expnum);
  	fourierfileF0 = fopen(filename, "wb");

		kNmax = momenta[grandspinnum];
    pkNmax = momenta[grandspinnum+1];
    dataforG(grandspinnum, Nmatrix, Pmatrix, Nval, Pval, ks, mks, pks);
    for (i=0; i<2*momenta[0]; i++) {
    }
    bessels(grandspinnum, momenta, pks, ks, mks, jcc, jcp, jcm, jpc, jpp, jmc, jmm);
    g_function(g_funcs[0], g_funcs[1], g_funcs[2], g_funcs[3], grandspinnum, Nmatrix, Pmatrix, mkNmax, kNmax, pkNmax, ks, pks, mks, jcc, jcp, jcm, jpc, jpp, jmc, jmm);
    f_function(f_funcs[0], f_funcs[1], f_funcs[2], f_funcs[3], grandspinnum, Nmatrix, Pmatrix, mkNmax, kNmax, pkNmax, ks, pks, mks, jcc, jcp, jcm, jpc, jpp, jmc, jmm);
    compile_fourier_func(fourier_funcP, fourier_funcN, g_funcs, f_funcs, mbessels, grandspinnum, mkNmax, kNmax, pkNmax);
    if (grandspinnum == 0) {
      for (i=0; i<kNmax*2; i++) {
        if (renorm) {
          renormalise(fourier_funcP[i], grandspinnum);
        }
        val = k_inner_product(fourier_funcP[i], fourier_funcP[i], grandspinnum);
        fprintf(stderr, "inner product momentum P: %lu %f\n", i+1, val);
      }
      for (i=0; i<pkNmax*2; i++) {
        //renomalise(fourier_funcN[i])
        if (renorm) {
          renormalise(fourier_funcN[i], grandspinnum);
        }
        val = k_inner_product(fourier_funcN[i], fourier_funcN[i], grandspinnum);
        fprintf(stderr, "inner product momentum N: %lu %f\n", i+1, val);
      }
      for (i=0; i<kNmax*2; i++) {
        for (j=0; j<2; j++) {
          for (k=0; k<=kNInt; k++) {
            fwrite(&fourier_funcP[i][j][k], sizeof(double), 1, fourierfile);
          }
        }
      }
      for (i=0; i<pkNmax*2; i++) {//1.357492
        for (j=0; j<2; j++) {
          for (k=0; k<=kNInt; k++) {
            fwrite(&fourier_funcN[i][j][k], sizeof(double), 1, fourierfile);
          }
        }
      }
    } else {
      for (i=0; i<kNmax*4; i++) {
        if (renorm) {
          renormalise(fourier_funcP[i], grandspinnum);
        }
        val = k_inner_product(fourier_funcP[i], fourier_funcP[i], grandspinnum);
        fprintf(stderr, "inner product momentum P: %lu %f\n", i+1, val);
      }
      for (i=0; i<pkNmax*2+mkNmax*2; i++) {
        if (renorm) {
          renormalise(fourier_funcN[i], grandspinnum);
        }
        val = k_inner_product(fourier_funcN[i], fourier_funcN[i], grandspinnum);
        fprintf(stderr, "inner product momentum N: %lu %f\n", i+1, val);
      }
      for (i=0; i<kNmax*4; i++) {
        for (j=0; j<4; j++) {
          for (k=0; k<=kNInt; k++) {
            fwrite(&fourier_funcP[i][j][k], sizeof(double), 1, fourierfile);
          }
        }
      }
      for (i=0; i<2*mkNmax+2*pkNmax; i++) {
        for (j=0; j<4; j++) {
          for (k=0; k<=kNInt; k++) {
            fwrite(&fourier_funcN[i][j][k], sizeof(double), 1, fourierfile);
          }
        }
      }
    }
    fprintf(stderr, "Starting free case fourier transform\n");
    for (i=0; i<kNmax*4; i++) {
      for (j=0; j<kNmax*4;j++) {
        if (j==i) {
          Pmatrix[i][j]=1.0;
        } else {
          Pmatrix[i][j]=0.0;
        }
      }
    }
    if (grandspinnum==0) {
      for (i=0; i<2*pkNmax;i++) {
        for (j=0; j<2*pkNmax; j++) {
          if (j==i) {
            Nmatrix[i][j]=1.0;
          } else {
            Nmatrix[i][j]=0.0;
          }
        }
      }
    } else {
      for (i=0; i<2*pkNmax+2*mkNmax;i++) {
        for (j=0; j<2*pkNmax+2*mkNmax; j++) {
          if (j==i) {
            Nmatrix[i][j]=1.0;
          } else {
            Nmatrix[i][j]=0.0;
          }
        }
      }
    }
    g_function(g_funcs[0], g_funcs[1], g_funcs[2], g_funcs[3], grandspinnum, Nmatrix, Pmatrix, mkNmax, kNmax, pkNmax, ks, pks, mks, jcc, jcp, jcm, jpc, jpp, jmc, jmm);
    f_function(f_funcs[0], f_funcs[1], f_funcs[2], f_funcs[3], grandspinnum, Nmatrix, Pmatrix, mkNmax, kNmax, pkNmax, ks, pks, mks, jcc, jcp, jcm, jpc, jpp, jmc, jmm);
    compile_fourier_func(fourier_funcP, fourier_funcN, g_funcs, f_funcs, mbessels, grandspinnum, mkNmax, kNmax, pkNmax);
    if (grandspinnum == 0) {
      for (i=0; i<kNmax*2; i++) {
        if (renorm) {
          renormalise(fourier_funcP[i], grandspinnum);
        }
        val = k_inner_product(fourier_funcP[i], fourier_funcP[i], grandspinnum);
        fprintf(stderr, "inner product momentum P F0: %lu %f\n", i+1, val);
      }
      for (i=0; i<pkNmax*2; i++) {
        if (renorm) {
          renormalise(fourier_funcN[i], grandspinnum);
        }
        val = k_inner_product(fourier_funcN[i], fourier_funcN[i], grandspinnum);
        fprintf(stderr, "inner product momentum N F0: %lu %f\n", i+1, val);
      }
      for (i=0; i<kNmax*2; i++) {
        for (j=0; j<2; j++) {
          for (k=0; k<=kNInt; k++) {
            fwrite(&fourier_funcP[i][j][k], sizeof(double), 1, fourierfileF0);
          }
        }
      }
      for (i=0; i<pkNmax*2; i++) {
        for (j=0; j<2; j++) {
          for (k=0; k<=kNInt; k++) {
            fwrite(&fourier_funcN[i][j][k], sizeof(double), 1, fourierfileF0);
          }
        }
      }
    } else {
      for (i=0; i<kNmax*4; i++) {
        if (renorm) {
          renormalise(fourier_funcP[i], grandspinnum);
        }
        val = k_inner_product(fourier_funcP[i], fourier_funcP[i], grandspinnum);
        fprintf(stderr, "inner product momentum P F0: %lu %f\n", i+1, val);
      }
      for (i=0; i<pkNmax*2+mkNmax*2; i++) {
        if (renorm) {
          renormalise(fourier_funcN[i], grandspinnum);
        }
        val = k_inner_product(fourier_funcN[i], fourier_funcN[i], grandspinnum);
        fprintf(stderr, "inner product momentum N F0: %lu %f\n", i+1, val);
      }
      for (i=0; i<kNmax*4; i++) {
        for (j=0; j<4; j++) {
          for (k=0; k<=kNInt; k++) {
            fwrite(&fourier_funcP[i][j][k], sizeof(double), 1, fourierfileF0);
          }
        }
      }
      for (i=0; i<2*mkNmax+2*pkNmax; i++) {
        for (j=0; j<4; j++) {
          for (k=0; k<=kNInt; k++) {
            fwrite(&fourier_funcN[i][j][k], sizeof(double), 1, fourierfileF0);
          }
        }
      }

    }
    fprintf(stderr, "end grandspinnum = %d, grandspintot=%d\n", grandspinnum, grandspintot);
		fclose(fourierfile);
		fclose(fourierfileF0);
  }
  fprintf(stderr, "end of program\n");
  free(temp);
  return 0;

}
