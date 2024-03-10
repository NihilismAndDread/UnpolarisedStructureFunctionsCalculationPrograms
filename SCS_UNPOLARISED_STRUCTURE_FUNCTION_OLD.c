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
double lambda=1.8576485321061220;
double kratio=4.0;
double kdist = 25.0;
double M;
double xNInt;
double xdist=10;
double Nc=3;
double mpion = 135.0;
double cqm = 400.0;*/

void Gspin_momentanums(unsigned short *granspinnum, unsigned short **pmomenta);
void dataforGvals(unsigned short grandspinnum, double ** Nval, double ** Pval, double **mks, double ** ks, double **pks);
void dataforGfouriers(unsigned short grandspinnum, double ****pfourier_funcP, double **** pfourier_funcN, double **** pfourier_funcPF0, double **** pfourier_funcNF0);
void compute_USF_valence(double **USFvalP, double **USFvalN, unsigned short *momentanums, double mass);
double RF_to_IMF(double *F_IMF, double *F_RF);
void chiralangle_from_file(double **thetas);
double compute_soliton_energy(double * thetas, unsigned short grandspintot, unsigned short *momentanums);
void freeDataCollectedvals(double *Pval, double *Nval, double *PvalF0, double *NvalF0, double *mks, double *ks, double *pks, unsigned short G, unsigned short *momentanums);
void freeDataCollectedfouriers(double ***fourier_funcP, double ***fourier_funcN, double ***fourier_funcPF0, double ***fourier_funcNF0, unsigned short G, unsigned short *momentanums);
void createF0vals(double **pPvalF0, double **pNvalF0, unsigned short G, double *mks, double *ks, double *pks, unsigned short *momentanums);
double interpol(double y1, double y2, double Mx, unsigned short index);
void compute_USF_sea(double ** USFseaP, double ** USFseaN, double ** USFseaPF0, double ** USFseaNF0, unsigned short grandspintot, double mass, unsigned short *momentanums);

double interpol(double y1, double y2, double Mx, unsigned short index) {
  double p1, p2, result;

  p1=index/(double)kNInt *kdist;
  p2=(index+1)/(double)kNInt *kdist;
  result = y1+(y2-y1)/(p2-p1)*(fabs(Mx)-p1);
  return result;
}
void createF0vals(double **pPvalF0, double **pNvalF0, unsigned short G, double *mks, double *ks, double *pks, unsigned short *momentanums) {
  unsigned short i, j;
  double *PvalF0, *NvalF0, temp;
  if (G==0) {
    PvalF0=(double*)malloc(sizeof(double)*2*momentanums[G]);
    NvalF0=(double*)malloc(sizeof(double)*2*momentanums[G+1]);
    for (j=0; j<momentanums[G]; j++) {
      temp = sqrt(1.0+ks[j]*ks[j]);
      PvalF0[j]=-temp;
      PvalF0[momentanums[G]*2-1-j]=temp;
    }
    for (j=0; j<momentanums[G+1]; j++) {
      temp = sqrt(1.0+pks[j]*pks[j]);
      NvalF0[j]=-temp;
      NvalF0[2*momentanums[G+1]-1-j]=temp;
    }
  } else {
    PvalF0=(double*)malloc(sizeof(double)*4*momentanums[G]);
    NvalF0=(double*)malloc(sizeof(double)*(2*momentanums[G+1]+2*momentanums[G-1]));
    for (j=0; j<momentanums[G]; j++) {
      temp= sqrt(1.0+ks[j]*ks[j]);
      PvalF0[j]=-temp;
      PvalF0[momentanums[G]*2-1-j]=temp;
      PvalF0[momentanums[G]*2+j]=-temp;
      PvalF0[momentanums[G]*4-1-j]=temp;
    }
    for (j=0; j<momentanums[G+1]; j++) {
      temp= sqrt(1.0+pks[j]*pks[j]);
      NvalF0[j]=-temp;
      NvalF0[momentanums[G+1]*2-1-j]=temp;
    }
    for (j=0; j<momentanums[G-1]; j++) {
      temp = sqrt(1.0+mks[j]*mks[j]);
      NvalF0[2*momentanums[G+1]+j]=-temp;
      NvalF0[2*momentanums[G+1]+2*momentanums[G-1]-1-j]=temp;
    }
  }
  *pPvalF0=PvalF0;
  *pNvalF0=NvalF0;
}
void freeDataCollectedvals(double *Pval, double *Nval, double *PvalF0, double *NvalF0, double *mks, double *ks, double *pks, unsigned short G, unsigned short *momentanums) {
  unsigned short i, j;
  if (G==0) {
    free(ks);
    free(pks);
  } else {
    if (momentanums[G]!=0) {
      free(ks);
    }
    if (momentanums[G+1]!=0) {
      free(pks);
    }
    free(mks);
  }
  free(Pval);
  free(Nval);
  if (PvalF0!=NULL) {
    free(PvalF0);
  }
  if (NvalF0!=NULL) {
    free(NvalF0);
  }
}
void freeDataCollectedfouriers(double ***fourier_funcP, double ***fourier_funcN, double ***fourier_funcPF0, double ***fourier_funcNF0, unsigned short G, unsigned short *momentanums) {
  unsigned short i, j;
  if (G==0) {
    for (i=0; i<momentanums[0]*2; i++) {
      for (j=0; j<2; j++) {
        free(fourier_funcP[i][j]);
        free(fourier_funcPF0[i][j]);
      }
    }
    for (i=0; i<momentanums[1]*2; i++) {
      for (j=0; j<2; j++) {
        free(fourier_funcN[i][j]);
        free(fourier_funcNF0[i][j]);
      }
    }
  } else {
    for (i=0; i<momentanums[G]*4; i++) {
      for (j=0; j<4; j++) {
        free(fourier_funcP[i][j]);
        free(fourier_funcPF0[i][j]);
      }
      free(fourier_funcP[i]);
      free(fourier_funcPF0[i]);
    }

    for (i=0; i<momentanums[G+1]*2+momentanums[G-1]*2; i++) {
      for (j=0; j<4; j++) {
        free(fourier_funcN[i][j]);
        free(fourier_funcNF0[i][j]);
      }
      free(fourier_funcN[i]);
      free(fourier_funcNF0[i]);
    }
  }
  free(fourier_funcP);
  free(fourier_funcN);
  free(fourier_funcPF0);
  free(fourier_funcNF0);
}
void freeDataCollectedfourierIP(double ***fourierIPP, double ***fourierIPN, double ***fourierIPPF0, double ***fourierIPNF0, unsigned short G, unsigned short *momentanums) {
  unsigned short i, j;
  if (G==0) {
    for (i=0; i<momentanums[0]*2; i++) {
      for (j=0; j<2; j++) {
        free(fourierIPP[i][j]);
        free(fourierIPPF0[i][j]);
      }
    }
    for (i=0; i<momentanums[1]*2; i++) {
      for (j=0; j<2; j++) {
        free(fourierIPN[i][j]);
        free(fourierIPNF0[i][j]);
      }
    }
  } else {
    for (i=0; i<momentanums[G]*4; i++) {
      for (j=0; j<2; j++) {
        free(fourierIPP[i][j]);
        free(fourierIPPF0[i][j]);
      }
      free(fourierIPP[i]);
      free(fourierIPPF0[i]);
    }

    for (i=0; i<momentanums[G+1]*2+momentanums[G-1]*2; i++) {
      for (j=0; j<2; j++) {
        free(fourierIPN[i][j]);
        free(fourierIPNF0[i][j]);
      }
      free(fourierIPN[i]);
      free(fourierIPNF0[i]);
    }
  }
  free(fourierIPP);
  free(fourierIPN);
  free(fourierIPPF0);
  free(fourierIPNF0);
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
double compute_soliton_energy_val(double *Pvals, double momentanum) {
  unsigned short i;
  double Ev;
  for (i=0; i<2*momentanum; i++) {
    if (fabs(Pvals[i])<=1) {
      Ev= Pvals[i];
    }
  }
  return Ev;
}
double compute_soliton_energy(double * thetas, unsigned short grandspintot, unsigned short *momentanums) {
  double t1, t2, t3, Ev, t21, t22, *intarr, x, k, Nval, Pval, En, Esol;
  double *Pvals, *Nvals, *mks, *ks, *pks;
  unsigned long i, j;

  dataforGvals(0, &Nvals, &Pvals, &mks, &ks, &pks);
  intarr= (double*)malloc((NInt+1)*sizeof(double));
  for (i=0; i<2*momentanums[0]; i++) {
    if (fabs(Pvals[i])<=1) {
      Ev= Pvals[i];
    }
  }
  Esol = Ev;
  if (Ev<0) {
    Esol=0;
  }
  Esol=0;
  freeDataCollectedvals(Pvals, Nvals, NULL, NULL, mks, ks, pks, 0, momentanums);
  for (i=0; i<=grandspintot; i++) {
    t2=0;
    dataforGvals(i, &Nvals, &Pvals, &mks, &ks, &pks);
    if (i==0) {
      for (j=0; j<momentanums[i]*2; j++) {
        Pval = Pvals[j];
        k = ks[j];
        En = sqrt(k*k+1.0);
        t2 += 1.0/2.0*(fabs(Pval)-sqrt(Pval*Pval+lambda*lambda)+ 1.0/2.0*lambda*lambda/sqrt(Pval*Pval+lambda*lambda));
        t2 += -1.0/2.0*(fabs(En)- sqrt(En*En+lambda*lambda) + 1.0/2.0*lambda*lambda/sqrt(En*En+lambda*lambda));
      }
      for (j=0; j<momentanums[i+1]*2; j++) {
        Nval = Nvals[j];
        k = pks[j];
        En = sqrt(k*k+1);
        t2 += 1.0/2.0*(fabs(Nval)-sqrt(Nval*Nval+lambda*lambda)+ 1.0/2.0*lambda*lambda/sqrt(Nval*Nval+lambda*lambda));
        t2 += -1.0/2.0*(fabs(En)- sqrt(En*En+lambda*lambda) + 1.0/2.0*lambda*lambda/sqrt(En*En+lambda*lambda));
      }
    } else {
      for (j=0; j<momentanums[i]*4; j++) {
        Pval = Pvals[j];
        k = ks[j%(momentanums[i]*2)];
        En = sqrt(k*k+1);
        t2 += 1.0/2.0*(2*i+1)*(fabs(Pval)-sqrt(Pval*Pval+lambda*lambda)+ 1.0/2.0*lambda*lambda/sqrt(Pval*Pval+lambda*lambda));
        t2 += -1.0/2.0*(2*i+1)*(fabs(En)- sqrt(En*En+lambda*lambda) + 1.0/2.0*lambda*lambda/sqrt(En*En+lambda*lambda));
      }
      for (j=0; j<momentanums[i-1]*2+momentanums[i+1]*2; j++) {
        Nval = Nvals[j];
        t2 += 1.0/2.0*(2*i+1)*(fabs(Nval) - sqrt(Nval*Nval+lambda*lambda) + 1.0/2.0*lambda*lambda/sqrt(Nval*Nval+lambda*lambda));
      }
      for (j=0; j<momentanums[i-1]*2; j++) {
        k = mks[j];
        En = sqrt(k*k+1);
        t2 += -1.0/2.0*(2*i+1)*(fabs(En) -sqrt(En*En+lambda*lambda) + 1.0/2.0*lambda*lambda/sqrt(En*En+lambda*lambda));
      }
      for (j=0; j<momentanums[i+1]*2; j++) {
        k = pks[j];
        En = sqrt(k*k+1);
        t2 += -1.0/2.0*(2*i+1)*(fabs(En) -sqrt(En*En+lambda*lambda) + 1.0/2.0*lambda*lambda/sqrt(En*En+lambda*lambda));
      }
    }
    Esol += -t2;
    freeDataCollectedvals(Pvals, Nvals, NULL, NULL, mks, ks, pks, i, momentanums);
  }
  Esol *= Nc;
  for (i=0; i<=NInt; i++) {
    x= (double)i/(double)NInt;
    intarr[i]=x*x*(1.0-cos(thetas[i]));
  }
  t3 = 4*M_PI*mpion*mpion*93.0*93.0/(cqm*cqm)*dist*dist*dist*booles(intarr, NInt)/(double)NInt;
  //Esol += t3;
  return Esol;
}
double RF_to_IMF(double *F_IMF, double *F_RF) {
  double x_IMF, x_RF, x1, x2;
  unsigned short index, mark;
  unsigned long i, j;
  index =0;
  mark=0;
  for (i=0; i<=xNInt; i++) {
    x_IMF = (double)(i/xNInt);
    x_RF=-log(1-x_IMF);
    for (j=index; j<xNInt; j++) {
      x1 = (double)j/xNInt*xdist;
      x2 = (double)(j+1)/xNInt*xdist;
      if (x1<x_RF && x_RF<=x2) {
        F_IMF[i]= step_function(1-x_IMF)/(1-x_IMF) *F_RF[j+1];
        index = j;
        break;
      } else if (j+1==xNInt) {
        if (mark==0) {
          mark = i;
        }
        F_IMF[i]= step_function(1-x_IMF)/(1-x_IMF) *F_RF[j+1];
        break;
      }
    }
  }
  return mark;
}
void Gspin_momentanums(unsigned short *grandspinnum, unsigned short **pmomenta) {
  FILE *momentanumsfile;
  unsigned short index, *momenta, i;
  unsigned short amomenta[200];
	char filename[100];
	sprintf(filename, "Data/expnum%d/momentanum_%d.b", expnum, expnum);
  momentanumsfile = fopen(filename, "rb");
  index =0;
  while (fread(&amomenta[index], sizeof(unsigned short), 1, momentanumsfile)) {
    index++;
  }
  momenta = (unsigned short *) malloc(sizeof(unsigned short)*(index));
  for (i=0; i<index; i++) {
    momenta[i] = amomenta[i];
  }
  momenta[index-1]=0;
  *pmomenta=momenta;
  *grandspinnum = index-2;
  fclose(momentanumsfile);
}
void dataforGfouriers(unsigned short grandspinnum, double ****pfourier_funcP, double **** pfourier_funcN, double **** pfourier_funcPF0, double **** pfourier_funcNF0) {
  unsigned short *momentanum, lenP, lenN, plenk, lenk, mlenk;
  unsigned long seekto1, seekto2, i, j, k;
  double *** fourier_funcP, ***fourier_funcN, temp, ***fourier_funcPF0, ***fourier_funcNF0;
	char filename[100];

  FILE * momentafile, *wavemodesfile, *wavemodesfileF0;
  momentanum = (unsigned short *)malloc(sizeof(unsigned short)*(grandspinnum+2));
	sprintf(filename, "Data/expnum%d/momentanum_%d.b", expnum, expnum);
  momentafile = fopen(filename, "rb");

  sprintf(filename, "Data/expnum%d/Fdata/fourierfileF0G%d_%d.b", expnum, grandspinnum, expnum);
	wavemodesfileF0=fopen(filename, "rb");
	sprintf(filename, "Data/expnum%d/Fdata/fourierfileG%d_%d.b", expnum, grandspinnum, expnum);
	wavemodesfile  =fopen(filename, "rb");

  for (i=0; i<=grandspinnum+1; i++) {
    fread(&momentanum[i], sizeof(unsigned short), 1, momentafile);
  }
  fclose(momentafile);
  lenk = momentanum[grandspinnum]*2;
  plenk= momentanum[grandspinnum+1]*2;
  mlenk= 0;
  if (grandspinnum!=0) {
    mlenk=momentanum[grandspinnum-1]*2;
  }
  if (grandspinnum==0) {
    lenP=momentanum[grandspinnum]*2;
    lenN=momentanum[grandspinnum+1]*2;
  } else {
    lenP= momentanum[grandspinnum]*4;
    lenN = momentanum[grandspinnum-1]*2 + momentanum[grandspinnum+1]*2;
  }

  //ALLOCATING SPACE
  fourier_funcP   = (double***)malloc(sizeof(double**)*lenP);
  fourier_funcN   = (double***)malloc(sizeof(double**)*lenN);
  fourier_funcPF0 = (double***)malloc(sizeof(double**)*lenP);
  fourier_funcNF0 = (double***)malloc(sizeof(double**)*lenN);
  if (grandspinnum==0) {
    for (i=0; i<lenP; i++) {
      fourier_funcP[i]  =(double**)malloc(sizeof(double*)*2);
      fourier_funcPF0[i]=(double**)malloc(sizeof(double*)*2);
      for (j=0; j<2; j++) {
        fourier_funcP[i][j]  =(double*) malloc(sizeof(double)*(kNInt+1));
        fourier_funcPF0[i][j]=(double*) malloc(sizeof(double)*(kNInt+1));
      }
    }
    for (i=0; i<lenN; i++) {
      fourier_funcN[i]  =(double**)malloc(sizeof(double*)*2);
      fourier_funcNF0[i]=(double**)malloc(sizeof(double*)*2);
      for (j=0; j<2; j++) {
        fourier_funcN[i][j]  =(double*)malloc(sizeof(double)*(kNInt+1));
        fourier_funcNF0[i][j]=(double*)malloc(sizeof(double)*(kNInt+1));
      }
    }
  } else {
    for (i=0; i<lenP; i++) {
      fourier_funcP[i]  =(double**)malloc(sizeof(double*)*4);
      fourier_funcPF0[i]=(double**)malloc(sizeof(double*)*4);
      for (j=0; j<4; j++) {
        fourier_funcP[i][j]  =(double *) malloc(sizeof(double)*(kNInt+1));
        fourier_funcPF0[i][j]=(double *) malloc(sizeof(double)*(kNInt+1));
      }
    }
    for (i=0; i<lenN; i++) {
      fourier_funcN[i]  =(double**)malloc(sizeof(double*)*4);
      fourier_funcNF0[i]=(double**)malloc(sizeof(double*)*4);
      for (j=0; j<4; j++) {
        fourier_funcN[i][j]  =(double*)malloc(sizeof(double)*(kNInt+1));
        fourier_funcNF0[i][j]=(double*)malloc(sizeof(double)*(kNInt+1));
      }
    }

  }

  //SEEKING PLACE IN FILES
  /*if (grandspinnum==0) {
    seekto2= 0;
  } else {
    seekto2 = momentanum[0]*2*2*(kNInt+1)+ momentanum[1]*2*2*(kNInt+1);
  }
  for (i=1; i<=grandspinnum-1; i++) {
    seekto2 += momentanum[i]*4*4*(kNInt+1) + (momentanum[i-1]*2 + momentanum[i+1]*2)*4*(kNInt+1);
  }
  fseek(wavemodesfile, sizeof(double)*seekto2, SEEK_SET);
  fseek(wavemodesfileF0, sizeof(double)*seekto2, SEEK_SET);
  */

  //COLLECTING DATA FROM FILES
  if (grandspinnum==0) {
    for (i=0; i<lenP; i++) {
      for (j=0; j<2; j++) {
        for (k=0; k<=kNInt; k++) {
          fread(&temp, sizeof(double), 1, wavemodesfile);
          fourier_funcP[i][j][k]=temp;
          fread(&temp, sizeof(double), 1, wavemodesfileF0);
          fourier_funcPF0[i][j][k]=temp;
        }
      }
    }
    for (i=0; i<lenN; i++) {
      for (j=0; j<2; j++) {
        for (k=0; k<=kNInt; k++) {
          fread(&temp, sizeof(double), 1, wavemodesfile);
          fourier_funcN[i][j][k]=temp;
          fread(&temp, sizeof(double), 1, wavemodesfileF0);
          fourier_funcNF0[i][j][k]=temp;
        }
      }
    }
  } else {
    for (i=0; i<lenP; i++) {
      for (j=0; j<4; j++) {
        for (k=0; k<=kNInt; k++) {
          fread(&temp, sizeof(double), 1, wavemodesfile);
          fourier_funcP[i][j][k]=temp;
          fread(&temp, sizeof(double), 1, wavemodesfileF0);
          fourier_funcPF0[i][j][k]=temp;
        }
      }
    }
    for (i=0; i<lenN; i++) {
      for (j=0; j<4; j++) {
        for (k=0; k<=kNInt; k++) {
          fread(&temp, sizeof(double), 1, wavemodesfile);
          fourier_funcN[i][j][k]=temp;
          fread(&temp, sizeof(double), 1, wavemodesfileF0);
          fourier_funcNF0[i][j][k]=temp;
        }
      }
    }
  }
  *pfourier_funcP=fourier_funcP;
  *pfourier_funcN=fourier_funcN;
  *pfourier_funcPF0=fourier_funcPF0;
  *pfourier_funcNF0=fourier_funcNF0;

  fclose(wavemodesfile);
  fclose(wavemodesfileF0);

  free(momentanum);
}
void dataforGvals(unsigned short grandspinnum, double ** Nval, double ** Pval, double **mks, double ** ks, double **pks) {
  unsigned short * momentanum, i, j, k, lenP, lenN, plenk, lenk, mlenk;
  unsigned long seekto1, seekto2;
  double * Pvals, * Nvals, *kn, *pkn, *mkn, temp;
	char filename[100];
  FILE * momentafile, * eigvalfile, * eigvecfile, *knfile;

  momentanum = (unsigned short *)malloc(sizeof(unsigned short)*(grandspinnum+2));

	sprintf(filename, "Data/expnum%d/momentanum_%d.b", expnum, expnum);
  momentafile = fopen(filename, "rb");
	sprintf(filename, "Data/expnum%d/eigenergies_%d.b", expnum, expnum);
	eigvalfile = fopen(filename, "rb");
	sprintf(filename, "Data/expnum%d/kn_%d.b", expnum, expnum);
  knfile = fopen(filename, "rb");

  //fprintf(stderr, "dataforG 1\n");
  for (i=0; i<=grandspinnum+1; i++) {
    fread(&momentanum[i], sizeof(unsigned short), 1, momentafile);
  }
  fclose(momentafile);
  //fprintf(stderr, "dataforG 2\n");
  lenk = momentanum[grandspinnum]*2;
  plenk= momentanum[grandspinnum+1]*2;
  mlenk= 0;
  //fprintf(stderr, "dataforG 3\n");
  if (grandspinnum!=0) {
    mlenk=momentanum[grandspinnum-1]*2;
  }
  //fprintf(stderr, "dataforG 4\n");
  if (grandspinnum==0) {
    lenP=momentanum[grandspinnum]*2;
    lenN=momentanum[grandspinnum+1]*2;
  } else {
    lenP= momentanum[grandspinnum]*4;
    lenN = momentanum[grandspinnum-1]*2 + momentanum[grandspinnum+1]*2;
  }
  //fprintf(stderr, "dataforG 5\n");
  Pvals = (double*)malloc(sizeof(double)*lenP);
  Nvals = (double*)malloc(sizeof(double)*lenN);
  //fprintf(stderr, "dataforG 6\n");
  if (momentanum[grandspinnum+1]==0) {
    pkn = NULL;
  } else {
    pkn = (double*)malloc(sizeof(double)*plenk);
  }
  if (momentanum[grandspinnum]==0) {
    kn = NULL;
  } else {
    kn = (double*)malloc(sizeof(double)*lenk);
  }

  if (grandspinnum==0) {
    mkn = NULL;
  } else {
    mkn = (double*)malloc(sizeof(double)*mlenk);
  }
  //fprintf(stderr, "dataforG 7\n");

  //fprintf(stderr, "dataforG 8\n");
  seekto1 =0;
  for (i=0; i<=grandspinnum-2; i++) {
    seekto1 += momentanum[i]*2;
  }
  //fprintf(stderr, "dataforG 9\n");
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
  //fprintf(stderr, "dataforG 10\n");
  fseek(eigvalfile, sizeof(double)*seekto1, SEEK_SET);
  for (i=0; i<lenP; i++) {
    fread(&Pvals[i], sizeof(double), 1, eigvalfile);
  }
  for (i=0; i<lenN; i++) {
    fread(&Nvals[i], sizeof(double), 1, eigvalfile);
  }
  fclose(eigvalfile);
  free(momentanum);
  *Pval=Pvals;
  *Nval=Nvals;
  *ks=kn;
  *mks=mkn;
  *pks=pkn;
}
void chiralangle_from_file(double **pthetas) {
  FILE *thetas_file;
  double *thetas;
  unsigned short i;
	char filename[100];
  thetas = (double*)malloc((NInt+1)*sizeof(double));

	sprintf(filename, "Data/expnum%d/chiralangle_%d.b", expnum, expnum);
  thetas_file = fopen(filename, "rb");
  fseek(thetas_file, 0, SEEK_SET);
  for (i=0; i<=NInt; i++) {
    fread(&thetas[i], sizeof(double), 1, thetas_file);
  }
  fclose(thetas_file);
  *pthetas = thetas;

}
void compute_USF_valence(double **USFvalP, double **USFvalN, unsigned short *momentanums, double mass) {
  double t1, t2, t3, t4, omega0, *intarr1, *intarr2, *intarr3, *intarr4, x, Mx, p, val_e, **val_funcs, num;
	double ***fourier_funcP, ***fourier_funcN, ***fourier_funcPF0, ***fourier_funcNF0, *Nval, *Pval, *mks, *ks, *pks;
  unsigned long i, j, k, l, pnum, temp1, temp2, index;
	FILE *USFvalfile;
	char filename[100];
  fprintf(stderr, "starting USF valence\n");
  intarr1 = (double*)malloc(sizeof(double)*(kNInt+1));
  intarr2 = (double*)malloc(sizeof(double)*(kNInt+1));
	intarr3 = (double*)malloc(sizeof(double)*(xNInt+1));
	intarr4 = (double*)malloc(sizeof(double)*(xNInt+1));
	/*for (i=0; i<2*momentanum; i++) {
		num =k_inner_product(fourier_funcP[i], fourier_funcP[i], 0);
		fprintf(stderr, "%d %20.10f\n", i+1, num);
	}*/
	for (i=0; i<=xNInt; i++) {
    USFvalP[0][i]=0;
    USFvalN[0][i]=0;
		USFvalP[1][i]=0;
		USFvalN[1][i]=0;
  }
  dataforGfouriers(0, &fourier_funcP, &fourier_funcN, &fourier_funcPF0, &fourier_funcNF0);
  dataforGvals(0, &Nval, &Pval, &mks, &ks, &pks);
  num =k_inner_product(fourier_funcP[0], fourier_funcP[0], 0);
  fprintf(stderr, "%d %20.10f\n", 1, num);
  index =0;
  for (i=0; i<2*momentanums[0]; i++) {
		//fprintf(stderr, "%f\n", Pval[i]);
    if (fabs(Pval[i])<1) {
      index = i;
      break;
    } else if (i==2*momentanums[0]-1) {
      fprintf(stderr, "Valence quark not detected");
      exit(1);
    }
  }
  val_e = Pval[index];
  val_funcs = fourier_funcP[index];
  fprintf(stdout, "Eval=%f\n", val_e*3*cqm);
	/*fprintf(stderr, "index = %d, val=%f\n", index, val_e);
  for (i=0; i<20; i++) {
		//fprintf(stderr, "%d %20.10f\n", i, fourier_funcP[index][0][i]);
	}
  for (i=0; i<20; i++) {
		//fprintf(stderr, "%d %20.10f\n", i, fourier_funcP[index][1][i]);
	}*/

	for (l=0; l<=xNInt; l++) {
    x = (double) l/(double)xNInt*xdist;
    Mx = mass*x+val_e;
    if (kdist>fabs(Mx)) {
      for (pnum=0; pnum <=kNInt; pnum++) {
        p = (double)pnum/(double)kNInt *kdist;
        if (fabs(Mx) < p) {
          t1 = val_funcs[0][pnum]*val_funcs[0][pnum];
          t1 += val_funcs[1][pnum]*val_funcs[1][pnum];
          t2 = -2 *Mx*val_funcs[0][pnum]*val_funcs[1][pnum];
          intarr1[pnum]=-p*t1*kdist;
					intarr2[pnum]=t2*kdist;
        } else {
					intarr1[pnum]=0;
					intarr2[pnum]=0;
        }
      }
      USFvalN[0][l] = booles(intarr1, kNInt)/(double)kNInt;
    	USFvalN[1][l] = booles(intarr2,kNInt)/(double)kNInt;
		} else {
      USFvalN[0][l]=0;
			USFvalN[1][l]=0;
    }
    Mx=mass*x-val_e;
    if (kdist>fabs(Mx)) {
      for (pnum=0; pnum <=kNInt; pnum++) {
        p = (double)pnum/(double)kNInt *kdist;
        if (fabs(Mx)<p) {
          t1 = val_funcs[0][pnum]*val_funcs[0][pnum];
          t1 += val_funcs[1][pnum]*val_funcs[1][pnum];
          t2 = -2 *Mx*val_funcs[0][pnum]*val_funcs[1][pnum];
        	intarr1[pnum]=p*t1*kdist;
					intarr2[pnum]=t2*kdist;
				} else {
					intarr1[pnum]=0;
					intarr2[pnum]=0;
        }

      }
      USFvalP[0][l]=booles(intarr1, kNInt)/(double)(kNInt);
      USFvalP[1][l]=booles(intarr2, kNInt)/(double)(kNInt);
    } else {
			USFvalP[0][l]=0;
			USFvalP[1][l]=0;
    }

  }
	/*for (i=0; i<20; i++) {
		fprintf(stderr, "%d %20.10f\n", i, USFvalP[0][i]);
	}
	for (i=0; i<20; i++) {
		fprintf(stderr, "%d %20.10f\n", i, USFvalP[1][i]);
	}
	for (i=0; i<20; i++) {
		fprintf(stderr, "%d %20.10f\n", i, USFvalN[0][i]);
	}
	for (i=0; i<20; i++) {
		fprintf(stderr, "%d %20.10f\n", i, USFvalN[1][i]);
	}*/
	for (i=0; i<=xNInt; i++) {
		x=i*xdist/(double)xNInt;
		intarr3[i]=x*(USFvalP[0][i]+USFvalN[0][i]);
		intarr4[i]=x*(USFvalP[1][i]+USFvalN[1][i]);
	}
	for (i=0; i<=kNInt; i++) {
		p=i*kdist/(double)kNInt;
		intarr1[i]=-4.0/3.0/mass/mass*p*p*p*val_funcs[0][i]*val_funcs[1][i];
	}
	t1=booles(intarr3, xNInt)*xdist/(double)xNInt;
	t2=booles(intarr4, xNInt)*xdist/(double)xNInt;
	t3=2.0*val_e/mass/mass;
	t4=kdist/(double)kNInt*booles(intarr1, kNInt);
	fprintf(stderr, "%f %f %f %f\n", t3, t1, t4, t2);
	sprintf(filename, "Data/expnum%d/USFdata/USFval_%d.b", expnum, expnum);
	USFvalfile = fopen(filename, "wb");
	for (i=0; i<=xNInt; i++) {
		fwrite(&USFvalP[0][i], sizeof(double), 1, USFvalfile);
		//fprintf(stderr, "%f\n", USFvalP[0][i]);
	}
	for (i=0; i<=xNInt; i++) {
		fwrite(&USFvalP[1][i], sizeof(double), 1, USFvalfile);
	}
	for (i=0; i<=xNInt; i++) {
		fwrite(&USFvalN[0][i], sizeof(double), 1, USFvalfile);
	}
	for (i=0; i<=xNInt; i++) {
		fwrite(&USFvalN[1][i], sizeof(double), 1, USFvalfile);
	}
	fclose(USFvalfile);
	free(intarr1);
  free(intarr2);
	free(intarr3);
	free(intarr4);
  freeDataCollectedvals(Pval, Nval, NULL, NULL, mks, ks, pks, 0, momentanums);
  freeDataCollectedfouriers(fourier_funcP, fourier_funcN, fourier_funcPF0, fourier_funcNF0, 0, momentanums);
  fprintf(stderr, "complete USF valence\n");
}
void compute_USF_sea_Esep(double ** USFseaP, double ** USFseaN, double ** USFseaPF0, double ** USFseaNF0, unsigned short grandspintot, double mass, unsigned short *momentanums, unsigned short Gbegin, unsigned short Gend, unsigned short GELcut) {
  unsigned short G, tempN, tempP, pnum;
  unsigned long i, j, k, l, index1, index2, index3;
  double temp, p, c1, c2, c11, c12, c21, c22, c31, c3, c4, c5, Mx1, Mx2, Mx3, omega0, omega0F0, x, y1, y2, yr;
  double *intarr1x, *intarr2x, *intarr3x, *intarr4x, *intarr1, *intarr2, *intarrtemp1, *intarrtemp2, *intarrtemp3, *intarrtemp4, *intarrtemp5, *carr, *gpp, *gp0;
  double ***fourier_funcP, ***fourier_funcN, ***fourier_funcPF0, ***fourier_funcNF0, *Nval, *Pval, *NvalF0, *PvalF0, *mks, *ks, *pks;
	double ***USFseaPL, ***USFseaNL, ***USFseaPF0L, ***USFseaNF0L, **USFseaPEL, **USFseaNEL, **USFseaPF0EL, **USFseaNF0EL, **USFseaPELsep, **USFseaNELsep, **USFseaPF0ELsep, **USFseaNF0ELsep;
  double srule0, srule1, srule0l, srule1l, Esrule0, Esrule1, Esrule0m, Esrule0l, Esrule1l, EAsrule0, EAsrule0m, EAsrule1, EAsrule0l, EAsrule1l;
  double srule0x, srule1x, srule0lx, srule1lx, Esrule0x, Esrule1x, Esrule0mx, Esrule0lx, Esrule1lx, EAsrule0x, EAsrule0mx, EAsrule1x, EAsrule0lx, EAsrule1lx;
  double srule0a, srule1a, srule0la, srule1la, Esrule0a, Esrule1a, Esrule0ma, Esrule0la, Esrule1la, EAsrule0a, EAsrule0ma, EAsrule1a, EAsrule0la, EAsrule1la;
  double *sepsrule0l, *sepsrule0m, *sepsrule0, *sepsrule0lx, *sepsrule0mx, *sepsrule0x, *sepsrule0la, *sepsrule0ma, *sepsrule0a;
  double *sepEsrule0l, *sepEsrule0m, *sepEsrule0, *sepEsrule0lx, *sepEsrule0mx, *sepEsrule0x, *sepEsrule0la, *sepEsrule0ma, *sepEsrule0a;
  double EAsrule0lN, EAsrule0lP, EAsrule0lNx, EAsrule0lPx, EAsrule0lNa, EAsrule0lPa;
  double EAsrule0mN, EAsrule0mP, EAsrule0mNx, EAsrule0mPx, EAsrule0mNa, EAsrule0mPa;
  double Esrule0lN, Esrule0lP, Esrule0lNx, Esrule0lPx, Esrule0lNa, Esrule0lPa;
  double Esrule0mN, Esrule0mP, Esrule0mNx, Esrule0mPx, Esrule0mNa, Esrule0mPa;
  double *sepsrule0lP, *sepsrule0mP, *sepsrule0lPx, *sepsrule0mPx, *sepsrule0lPa, *sepsrule0mPa;
  double *sepsrule0lN, *sepsrule0mN, *sepsrule0lNx, *sepsrule0mNx, *sepsrule0lNa, *sepsrule0mNa;
  double *sepEsrule0lP, *sepEsrule0mP, *sepEsrule0lPx, *sepEsrule0mPx, *sepEsrule0lPa, *sepEsrule0mPa;
  double *sepEsrule0lN, *sepEsrule0mN, *sepEsrule0lNx, *sepEsrule0mNx, *sepEsrule0lNa, *sepEsrule0mNa;
  FILE *GsepUSFfile, *USFfile, *EsepSRfile, *EsepSRxfile, *EsepSRafile, *EsepSRsepfile, *EsepSRsepxfile, *EsepSRsepafile;
  FILE *EsepSRPNfile, *EsepSRPNxfile, *EsepSRPNafile, *EsepSRPNsepfile, *EsepSRPNsepxfile, *EsepSRPNsepafile;
	char filename[100];

	sprintf(filename, "Data/expnum%d/USFdata/USFfile_%d.b", expnum, expnum);
	USFfile = fopen(filename, "wb");
  //--------Total
  sepsrule0l = (double*)malloc(sizeof(double)*3);
  sepsrule0m = (double*)malloc(sizeof(double)*3);
  sepsrule0  = (double*)malloc(sizeof(double)*3);

  sepsrule0la = (double*)malloc(sizeof(double)*3);
  sepsrule0ma = (double*)malloc(sizeof(double)*3);
  sepsrule0a  = (double*)malloc(sizeof(double)*3);

  sepsrule0lx = (double*)malloc(sizeof(double)*3);
  sepsrule0mx = (double*)malloc(sizeof(double)*3);
  sepsrule0x  = (double*)malloc(sizeof(double)*3);

  sepEsrule0l = (double*)malloc(sizeof(double)*3);
  sepEsrule0m = (double*)malloc(sizeof(double)*3);
  sepEsrule0  = (double*)malloc(sizeof(double)*3);

  sepEsrule0la = (double*)malloc(sizeof(double)*3);
  sepEsrule0ma = (double*)malloc(sizeof(double)*3);
  sepEsrule0a  = (double*)malloc(sizeof(double)*3);

  sepEsrule0lx = (double*)malloc(sizeof(double)*3);
  sepEsrule0mx = (double*)malloc(sizeof(double)*3);
  sepEsrule0x  = (double*)malloc(sizeof(double)*3);
  //---------Quark
  sepsrule0lP = (double*)malloc(sizeof(double)*3);
  sepsrule0mP = (double*)malloc(sizeof(double)*3);

  sepsrule0lPa = (double*)malloc(sizeof(double)*3);
  sepsrule0mPa = (double*)malloc(sizeof(double)*3);

  sepsrule0lPx = (double*)malloc(sizeof(double)*3);
  sepsrule0mPx = (double*)malloc(sizeof(double)*3);

  sepEsrule0lP = (double*)malloc(sizeof(double)*3);
  sepEsrule0mP = (double*)malloc(sizeof(double)*3);

  sepEsrule0lPa = (double*)malloc(sizeof(double)*3);
  sepEsrule0mPa = (double*)malloc(sizeof(double)*3);

  sepEsrule0lPx = (double*)malloc(sizeof(double)*3);
  sepEsrule0mPx = (double*)malloc(sizeof(double)*3);
  //-----------Anti-Quark
  sepsrule0lN = (double*)malloc(sizeof(double)*3);
  sepsrule0mN = (double*)malloc(sizeof(double)*3);

  sepsrule0lNa = (double*)malloc(sizeof(double)*3);
  sepsrule0mNa = (double*)malloc(sizeof(double)*3);

  sepsrule0lNx = (double*)malloc(sizeof(double)*3);
  sepsrule0mNx = (double*)malloc(sizeof(double)*3);

  sepEsrule0lN = (double*)malloc(sizeof(double)*3);
  sepEsrule0mN = (double*)malloc(sizeof(double)*3);

  sepEsrule0lNa = (double*)malloc(sizeof(double)*3);
  sepEsrule0mNa = (double*)malloc(sizeof(double)*3);

  sepEsrule0lNx = (double*)malloc(sizeof(double)*3);
  sepEsrule0mNx = (double*)malloc(sizeof(double)*3);
  //---------
  //---------Allocate memory for arrays for integration
  intarr1=(double*)malloc(sizeof(double)*(kNInt+1));
  intarr2=(double*)malloc(sizeof(double)*(kNInt+1));
  intarrtemp3=(double*)malloc(sizeof(double)*(kNInt+1));
  intarrtemp2=(double*)malloc(sizeof(double)*(kNInt+1));
  intarrtemp1=(double*)malloc(sizeof(double)*(kNInt+1));
  intarrtemp4=(double*)malloc(sizeof(double)*(kNInt+1));
  intarrtemp5=(double*)malloc(sizeof(double)*(kNInt+1));
  intarr1x = (double*)malloc(sizeof(double)*(xNInt+1));
  intarr2x = (double*)malloc(sizeof(double)*(xNInt+1));
  intarr3x = (double*)malloc(sizeof(double)*(xNInt+1));
  intarr4x = (double*)malloc(sizeof(double)*(xNInt+1));
	carr=(double*)malloc(sizeof(double)*(kNInt+1));
	gpp = (double*)malloc(sizeof(double)*(kNInt+1));
	gp0 =(double*)malloc(sizeof(double)*(kNInt+1));
  //-----------Allocate memory for arrays for structure functions
	USFseaPL = (double***)malloc(2*sizeof(double**));
	USFseaNL = (double***)malloc(2*sizeof(double**));
	USFseaPF0L = (double***)malloc(2*sizeof(double**));
	USFseaNF0L = (double***)malloc(2*sizeof(double**));

	USFseaPL[0]=(double**)malloc(sizeof(double*)*(xNInt+1));
	USFseaPL[1]=(double**)malloc(sizeof(double*)*(xNInt+1));

	USFseaNL[0]=(double**)malloc(sizeof(double*)*(xNInt+1));
	USFseaNL[1]=(double**)malloc(sizeof(double*)*(xNInt+1));

	USFseaPF0L[0]=(double**)malloc(sizeof(double*)*(xNInt+1));
	USFseaPF0L[1]=(double**)malloc(sizeof(double*)*(xNInt+1));

	USFseaNF0L[0]=(double**)malloc(sizeof(double*)*(xNInt+1));
	USFseaNF0L[1]=(double**)malloc(sizeof(double*)*(xNInt+1));


  for (l=0; l<=xNInt; l++) {
    USFseaP[0][l]=0;
    USFseaN[0][l]=0;
    USFseaPF0[0][l]=0;
    USFseaNF0[0][l]=0;
    USFseaP[1][l]=0;
    USFseaN[1][l]=0;
    USFseaPF0[1][l]=0;
    USFseaNF0[1][l]=0;
  }

  srule0 =0;
	srule1 =0;
  srule0l=0;
	srule1l=0;

  EAsrule0 =0;
	EAsrule1 =0;
  EAsrule0l=0;
	EAsrule1l=0;
  EAsrule0m=0;

  srule0x =0;
	srule1x =0;
  srule0lx=0;
	srule1lx=0;

  EAsrule0x =0;
	EAsrule1x =0;
  EAsrule0lx=0;
	EAsrule1lx=0;
  EAsrule0mx=0;

  srule0a =0;
	srule1a =0;
  srule0la=0;
	srule1la=0;

  EAsrule0a =0;
	EAsrule1a =0;
  EAsrule0la=0;
	EAsrule1la=0;
  EAsrule0ma=0;

  EAsrule0lP =0;
  EAsrule0lN =0;
  EAsrule0lPx=0;
  EAsrule0lNx=0;
  EAsrule0lPa=0;
  EAsrule0lNa=0;

  EAsrule0mP =0;
  EAsrule0mN =0;
  EAsrule0mPx=0;
  EAsrule0mNx=0;
  EAsrule0mPa=0;
  EAsrule0mNa=0;

  for (i=0; i<3; i++) {
    sepsrule0l[i]=0;
    sepsrule0la[i]=0;
    sepsrule0lx[i]=0;

    sepsrule0m[i]=0;
    sepsrule0ma[i]=0;
    sepsrule0mx[i]=0;

    sepsrule0[i]=0;
    sepsrule0a[i]=0;
    sepsrule0x[i]=0;
    //--------Quark
    sepsrule0lP[i]=0;
    sepsrule0lPa[i]=0;
    sepsrule0lPx[i]=0;

    sepsrule0mP[i]=0;
    sepsrule0mPa[i]=0;
    sepsrule0mPx[i]=0;
    //--------Anti-Quark
    sepsrule0lN[i]=0;
    sepsrule0lNa[i]=0;
    sepsrule0lNx[i]=0;

    sepsrule0mN[i]=0;
    sepsrule0mNa[i]=0;
    sepsrule0mNx[i]=0;
  }

  if (Gend>grandspintot) {
    Gend = grandspintot;
  }
  //*************************************
  //*************************************
  //LOOP THROUGH G'S
  //*************************************
  //*************************************
  for (G=Gbegin; G<=Gend; G++) {
		//*************************************
		//*************************************
		//OPENING FILE
		sprintf(filename, "Data/expnum%d/USFdata/GsepUSFfileG%d_%d.b", expnum, G, expnum);
		GsepUSFfile = fopen(filename, "wb");
    //*************************************
    //*************************************
    //DATA COLLECTION
    dataforGfouriers(G, &fourier_funcP, &fourier_funcN, &fourier_funcPF0, &fourier_funcNF0);
    dataforGvals(G, &Nval, &Pval, &mks, &ks, &pks);
		//*************************************
    //*************************************
    //CREATING F0 VALS
    createF0vals(&PvalF0, &NvalF0, G, mks, ks, pks, momentanums);
    //*************************************

    if (G==0) {
      tempN=2*momentanums[1];
      tempP=2*momentanums[0];
    } else {
      tempN=2*momentanums[G-1]+2*momentanums[G+1];
      tempP=4*momentanums[G];
    }
		for (i=0; i<tempP; i++) {
			//c11 = k_inner_product(fourier_funcN[i], fourier_funcN[i], G);
			//fprintf(stderr, "%20.10f\n", c11);
		}
		//*************************************
		//*************************************
		//CLEANING VARIABLES AND ALLOCATING SPACE
		for (i=0; i<=xNInt; i++) {
			USFseaPL[0][i]=(double*)malloc(sizeof(double)*(tempP+tempN));
			USFseaNL[0][i]=(double*)malloc(sizeof(double)*(tempP+tempN));
			USFseaPF0L[0][i]=(double*)malloc(sizeof(double)*(tempP+tempN));
			USFseaNF0L[0][i]=(double*)malloc(sizeof(double)*(tempP+tempN));

      USFseaPL[1][i]=(double*)malloc(sizeof(double)*(tempP+tempN));
			USFseaNL[1][i]=(double*)malloc(sizeof(double)*(tempP+tempN));
			USFseaPF0L[1][i]=(double*)malloc(sizeof(double)*(tempP+tempN));
			USFseaNF0L[1][i]=(double*)malloc(sizeof(double)*(tempP+tempN));
		}
    if (G<GELcut) {
      USFseaPEL=(double**)malloc(sizeof(double*)*3);
      USFseaNEL=(double**)malloc(sizeof(double*)*3);
      USFseaPF0EL=(double**)malloc(sizeof(double*)*3);
      USFseaNF0EL=(double**)malloc(sizeof(double*)*3);

      USFseaPELsep=(double**)malloc(sizeof(double*)*3);
      USFseaNELsep=(double**)malloc(sizeof(double*)*3);
      USFseaPF0ELsep=(double**)malloc(sizeof(double*)*3);
      USFseaNF0ELsep=(double**)malloc(sizeof(double*)*3);
      for (i=0; i<3; i++) {
        USFseaPEL[i]=(double*)malloc(sizeof(double)*(xNInt+1));
        USFseaNEL[i]=(double*)malloc(sizeof(double)*(xNInt+1));
        USFseaPF0EL[i]=(double*)malloc(sizeof(double)*(xNInt+1));
        USFseaNF0EL[i]=(double*)malloc(sizeof(double)*(xNInt+1));

        USFseaPELsep[i]=(double*)malloc(sizeof(double)*(xNInt+1));
        USFseaNELsep[i]=(double*)malloc(sizeof(double)*(xNInt+1));
        USFseaPF0ELsep[i]=(double*)malloc(sizeof(double)*(xNInt+1));
        USFseaNF0ELsep[i]=(double*)malloc(sizeof(double)*(xNInt+1));
      }
    }
    //-*******************************
    //EVEN CHANNEL
    //-*******************************
    if (G<GELcut) {
      sprintf(filename, "Data/expnum%d/SRdata/GsepELSRfileEG%d_%d.txt", expnum, G, expnum);
      EsepSRfile = fopen(filename, "w+");
      sprintf(filename, "Data/expnum%d/SRdata/GsepELSRxfileEG%d_%d.txt", expnum, G, expnum);
      EsepSRxfile = fopen(filename, "w+");
      sprintf(filename, "Data/expnum%d/SRdata/GsepELSRafileEG%d_%d.txt", expnum, G, expnum);
      EsepSRafile = fopen(filename, "w+");

      sprintf(filename, "Data/expnum%d/SRdata/GsepELSRsepfileEG%d_%d.txt", expnum, G, expnum);
      EsepSRsepfile = fopen(filename, "w+");
      sprintf(filename, "Data/expnum%d/SRdata/GsepELSRsepafileEG%d_%d.txt", expnum, G, expnum);
      EsepSRsepafile = fopen(filename, "w+");
      sprintf(filename, "Data/expnum%d/SRdata/GsepELSRsepxfileEG%d_%d.txt", expnum, G, expnum);
      EsepSRsepxfile = fopen(filename, "w+");

      sprintf(filename, "Data/expnum%d/SRdata/GsepELSRPNfileEG%d_%d.txt", expnum, G, expnum);
      EsepSRPNfile = fopen(filename, "w+");
      sprintf(filename, "Data/expnum%d/SRdata/GsepELSRPNafileEG%d_%d.txt", expnum, G, expnum);
      EsepSRPNafile = fopen(filename, "w+");
      sprintf(filename, "Data/expnum%d/SRdata/GsepELSRPNxfileEG%d_%d.txt", expnum, G, expnum);
      EsepSRPNxfile = fopen(filename, "w+");

      sprintf(filename, "Data/expnum%d/SRdata/GsepELSRPNsepfileEG%d_%d.txt", expnum, G, expnum);
      EsepSRPNsepfile = fopen(filename, "w+");
      sprintf(filename, "Data/expnum%d/SRdata/GsepELSRPNsepafileEG%d_%d.txt", expnum, G, expnum);
      EsepSRPNsepafile = fopen(filename, "w+");
      sprintf(filename, "Data/expnum%d/SRdata/GsepELSRPNsepxfileEG%d_%d.txt", expnum, G, expnum);
      EsepSRPNsepxfile = fopen(filename, "w+");
    }
    for (i=0; i<tempP; i++) {
      Esrule1 =0;
      Esrule1x=0;
      Esrule1a=0;
      //LOOP THROUGH x's
      for (l=0; l<=xNInt; l++) {
        x=(double)l/(double)xNInt*xdist;
        //INTEGRALS--------------------
        for (pnum=0; pnum<=kNInt; pnum++) {
          p=(double)pnum/(double)kNInt*kdist;
          if (G==0) {
            intarr1[pnum] =fourier_funcP[i][0][pnum]*fourier_funcP[i][0][pnum];
            intarr1[pnum]+=fourier_funcP[i][1][pnum]*fourier_funcP[i][1][pnum];
            carr[pnum]=intarr1[pnum]*kdist;
  					intarr1[pnum]*=p*kdist;

            intarr2[pnum]=-2.0*kdist*fourier_funcP[i][0][pnum]*fourier_funcP[i][1][pnum];

          } else {
            intarr1[pnum] =fourier_funcP[i][0][pnum]*fourier_funcP[i][0][pnum];
            intarr1[pnum]+=fourier_funcP[i][1][pnum]*fourier_funcP[i][1][pnum];
            intarr1[pnum]+=fourier_funcP[i][2][pnum]*fourier_funcP[i][2][pnum];
            intarr1[pnum]+=fourier_funcP[i][3][pnum]*fourier_funcP[i][3][pnum];
            carr[pnum]=intarr1[pnum]*kdist*(G*2.0+1.0);
  					intarr1[pnum]*=p*kdist*(G*2.0+1.0);

            intarr2[pnum] =fourier_funcP[i][0][pnum]*fourier_funcP[i][2][pnum];
            intarr2[pnum]+=fourier_funcP[i][1][pnum]*fourier_funcP[i][3][pnum];
            intarr2[pnum]*=-2.0*kdist*(G*2.0+1.0);
          }
  				if (l==0) {
  					temp = sqrt(Pval[i]*Pval[i]+lambda*lambda);
  					gpp[pnum] = 2.0*Pval[i]/fabs(Pval[i])*(1.0-fabs(Pval[i])/temp-0.5*fabs(Pval[i])*lambda*lambda/temp/temp/temp);
  					gpp[pnum] *= -1.0/(3.0*mass*mass)*p*p*p*intarr2[pnum];
  				}
  			}
  			//fprintf(stderr, "vibe check4\n");
  			if (l==0) {
          c11=booles(gpp, kNInt)/(double)kNInt;
  				Esrule1 += c11;
          srule1 += c11;
          Esrule1x += c11;
          srule1x += c11;

  			}
        omega0 = sqrt(lambda*lambda+Pval[i]*Pval[i]);
        //-----------------------------
        //ANTI-QUARK CONTRIBUTION
        Mx1 = mass*x+fabs(Pval[i]);
        Mx2 = mass*x+omega0;
        Mx3 = mass*x;
        index1=0;
        index2=0;
        index3=0;
        for (pnum=0; pnum<=kNInt; pnum++) {
          p = (double)pnum/(double)kNInt*kdist;
          if (p<fabs(Mx1)) {
            index1 = pnum;
            intarrtemp1[pnum]=0;
            intarrtemp2[pnum]=0;
          } else {
            intarrtemp1[pnum]=intarr1[pnum];
            intarrtemp2[pnum]=intarr2[pnum];
          }

          if (p<fabs(Mx2)) {
            index2 = pnum;
            intarrtemp3[pnum]=0;
            intarrtemp4[pnum]=0;
          } else {
            intarrtemp3[pnum]=intarr1[pnum];
            intarrtemp4[pnum]=intarr2[pnum];
          }
          if (p<fabs(Mx3)) {
            index3=pnum;
            intarrtemp5[pnum]=0;
          } else {
            intarrtemp5[pnum]=intarr1[pnum];
          }
        }
  			c3 = lambda*lambda/(2.0*omega0)*Mx2*interpol(y1, y2, Mx2, index2);
        c11=0;
  			c12=0;
        c21=0;
  			c22=0;
        c31=0;
        c3=0;
        c4=0;
        c5=0;
        if (kNInt>=index1+1) {
          c11 = booles(intarrtemp1, kNInt)/(double)(kNInt);
          c12 = -Pval[i]/fabs(Pval[i])*Mx1*booles(intarrtemp2, kNInt)/(double)(kNInt);
        }
        if (kNInt>=index2+1) {
          c21 = booles(intarrtemp3, kNInt)/(double)(kNInt);
          c22 = -Pval[i]/omega0*Mx2*booles(intarrtemp4, kNInt)/(double)(kNInt);


          c5 = lambda*lambda*Pval[i]/omega0/omega0/omega0*mass*x/2.0 * booles(intarrtemp4, kNInt)/(double)(kNInt);

          y2 = carr[index2+1]/kdist;
          y1 = carr[index2]/kdist;
          c3 = lambda*lambda/(2.0*omega0)*Mx2*interpol(y1, y2, Mx2, index2);

  				y1 =intarr2[index2]/kdist;
          y2 =intarr2[index2+1]/kdist;
          c4 = lambda*lambda*Pval[i]/(2.0*omega0*omega0)*fabs(Mx2)*interpol(y1, y2, Mx2, index2);
        }
        if (kNInt>=index3+1) {
          c31 = booles(intarrtemp5, kNInt)/(double)kNInt;
        }
        if (G<GELcut) {
          USFseaNEL[0][l] = c11-c21-c3;
          USFseaNEL[1][l] = c12-c22+c5+c4;
          USFseaNEL[2][l] = c31;

          USFseaNELsep[0][l] = c11;
          USFseaNELsep[1][l] = -c21;
          USFseaNELsep[2][l] = -c3;
        }
        USFseaNL[0][l][i] = c11-c21-c3;
        USFseaNL[1][l][i] = c12-c22+c5+c4;
  			//-----------------------------
        //QUARK CONTRIBUTION
        Mx1 = mass*x-fabs(Pval[i]);
        Mx2 = mass*x-omega0;
        Mx3 = mass*x;
        index1=0;
        index2=0;
        index3=0;
        for (pnum=0; pnum<=kNInt; pnum++) {
          p = (double)pnum/(double)kNInt*kdist;
          if (p<fabs(Mx1)) {
            index1 = pnum;
            intarrtemp1[pnum]=0;
            intarrtemp2[pnum]=0;
          } else {
            intarrtemp1[pnum]=intarr1[pnum];
            intarrtemp2[pnum]=intarr2[pnum];
          }
          if (p<fabs(Mx2)) {
            index2 = pnum;
            intarrtemp3[pnum]=0;
            intarrtemp4[pnum]=0;
          } else {
            intarrtemp3[pnum]=intarr1[pnum];
            intarrtemp4[pnum]=intarr2[pnum];
          }
          if (p<fabs(Mx3)) {
            index3=0;
            intarrtemp5[pnum]=0;
          } else {
            intarrtemp5[pnum]=intarr1[pnum];
          }
        }
        c11=0;
  			c12=0;
        c21=0;
  			c22=0;
        c31=0;
        c3=0;
        c4=0;
        c5=0;
        if (kNInt>=index1+1) {
          c11 = -booles(intarrtemp1, kNInt)/(double)(kNInt);
          c12 = -Pval[i]/fabs(Pval[i])*Mx1*booles(intarrtemp2, kNInt)/(double)(kNInt);
        }
        if (kNInt>=index2+1) {
          c21 = -booles(intarrtemp3, kNInt)/(double)(kNInt);
          c22 = -Pval[i]/omega0*Mx2*booles(intarrtemp4, kNInt)/(double)(kNInt);

          p = (double)(index2+1.0)/(double)(kNInt)*kdist;
          c5 = lambda*lambda*Pval[i]/omega0/omega0/omega0*mass*x/2.0 * booles(intarrtemp4, kNInt)/(double)(kNInt);

          y2 = carr[index2+1]/kdist;
          y1 = carr[index2]/kdist;
          c3 = lambda*lambda/(2.0*omega0)*Mx2*interpol(y1, y2, Mx2, index2);

  				y1 = intarr2[index2]/kdist;
          y2 = intarr2[index2+1]/kdist;
          c4 = lambda*lambda*Pval[i]/(2.0*omega0*omega0)*fabs(Mx2)*interpol(y1, y2, Mx2, index2);
        }
        if (kNInt>=index3+1) {
          c31= -booles(intarrtemp5, kNInt)/(double)kNInt;
        }
        if (G<GELcut) {
          USFseaPEL[0][l] = c11-c21-c3;
          USFseaPEL[1][l] = c12-c22+c5-c4;
          USFseaPEL[2][l] = c31;

          USFseaPELsep[0][l] = c11;
          USFseaPELsep[1][l] = -c21;
          USFseaPELsep[2][l] = -c3;
        }
        USFseaPL[0][l][i] = c11-c21-c3;
        USFseaPL[1][l][i] = c12-c22+c5-c4;
  			//-----------------------------
      }
      //LOOP THROUGH x's F0
      for (l=0; l<=xNInt; l++) {
        x=(double)l/(double)xNInt*xdist;
        //INTEGRALS--------------------
        for (pnum=0; pnum<=kNInt; pnum++) {
          p=(double)pnum/(double)kNInt*kdist;
          if (G==0) {
            intarr1[pnum] =fourier_funcPF0[i][0][pnum]*fourier_funcPF0[i][0][pnum];
            intarr1[pnum]+=fourier_funcPF0[i][1][pnum]*fourier_funcPF0[i][1][pnum];
            carr[pnum]=intarr1[pnum]*kdist;
            intarr1[pnum]*=p*kdist;

            intarr2[pnum]=-2.0*kdist*fourier_funcPF0[i][0][pnum]*fourier_funcPF0[i][1][pnum];
          } else {
            intarr1[pnum] =fourier_funcPF0[i][0][pnum]*fourier_funcPF0[i][0][pnum];
            intarr1[pnum]+=fourier_funcPF0[i][1][pnum]*fourier_funcPF0[i][1][pnum];
            intarr1[pnum]+=fourier_funcPF0[i][2][pnum]*fourier_funcPF0[i][2][pnum];
            intarr1[pnum]+=fourier_funcPF0[i][3][pnum]*fourier_funcPF0[i][3][pnum];
            carr[pnum]=intarr1[pnum]*kdist*(G*2.0+1.0);
            intarr1[pnum]*=p*kdist*(G*2.0+1.0);

            intarr2[pnum] =fourier_funcPF0[i][0][pnum]*fourier_funcPF0[i][2][pnum];
            intarr2[pnum]+=fourier_funcPF0[i][1][pnum]*fourier_funcPF0[i][3][pnum];
            intarr2[pnum]*=-2.0*kdist*(G*2.0+1.0);
  				}
  				if (l==0) {
  					temp = sqrt(PvalF0[i]*PvalF0[i]+lambda*lambda);
  					gp0[pnum] = 2.0*PvalF0[i]/fabs(PvalF0[i])*(1.0-fabs(PvalF0[i])/temp-0.5*fabs(PvalF0[i])*lambda*lambda/temp/temp/temp);
  					gp0[pnum] *= -1.0/(3.0*mass*mass)*p*p*p*intarr2[pnum];
  				}
        }
  			if (l==0) {
          c11=-booles(gp0, kNInt)/(double)kNInt;
          Esrule1 += c11;
  				srule1  += c11;
          Esrule1a += c11;
          srule1a  += c11;
  			}
        omega0 = sqrt(lambda*lambda+PvalF0[i]*PvalF0[i]);
        //-----------------------------
        //ANTI-QUARK CONTRIBUTION
        Mx1 = mass*x+fabs(PvalF0[i]);
        Mx2 = mass*x+omega0;
        Mx3 = mass*x;
        index1=0;
        index2=0;
        index3=0;
        for (pnum=0; pnum<=kNInt; pnum++) {
          p = (double)pnum/(double)kNInt*kdist;
          if (p<fabs(Mx1)) {
            index1 = pnum;
            intarrtemp1[pnum]=0;
            intarrtemp2[pnum]=0;
          } else {
            intarrtemp1[pnum]=intarr1[pnum];
            intarrtemp2[pnum]=intarr2[pnum];
          }

          if (p<fabs(Mx2)) {
            index2 = pnum;
            intarrtemp3[pnum]=0;
            intarrtemp4[pnum]=0;
          } else {
            intarrtemp3[pnum]=intarr1[pnum];
            intarrtemp4[pnum]=intarr2[pnum];
          }
          if (p<fabs(Mx3)) {
            index3=0;
            intarrtemp5[pnum]=0;
          } else {
            intarrtemp5[pnum]=intarr1[pnum];
          }
        }
        c11=0;
        c21=0;
        c12=0;
        c22=0;
        c31=0;
        c3=0;
        c4=0;
        c5=0;
        if (kNInt>=index1+1) {
          c11 = booles(intarrtemp1, kNInt)/(double)(kNInt);
          c12 = -PvalF0[i]/fabs(PvalF0[i])*Mx1*booles(intarrtemp2, kNInt)/(double)(kNInt);
        }
        if (kNInt>=index2+1) {
          c21 = booles(intarrtemp3, kNInt)/(double)(kNInt);
          c22 = -PvalF0[i]/omega0*Mx2*booles(intarrtemp4, kNInt)/(double)(kNInt);

          c5 = lambda*lambda*PvalF0[i]/omega0/omega0/omega0*mass*x/2.0 * booles(intarrtemp4, kNInt)/(double)(kNInt);

          y2 = carr[index2+1]/kdist;
          y1 = carr[index2]/kdist;
  				c3 = lambda*lambda/(2.0*omega0)*Mx2*interpol(y1, y2, Mx2, index2);

  				y1 =intarr2[index2]/kdist;
          y2 =intarr2[index2+1]/kdist;
          c4 = lambda*lambda*PvalF0[i]/(2.0*omega0*omega0)*fabs(Mx2)*interpol(y1, y2, Mx2, index2);
        }
        if (kNInt>=index3+1) {
          c31=booles(intarrtemp5, kNInt)/(double)kNInt;
        }
        if (G<GELcut) {
          USFseaNF0EL[0][l] = c11-c21-c3;
          USFseaNF0EL[1][l] = c12-c22+c5+c4;
          USFseaNF0EL[2][l] = c31;

          USFseaNF0ELsep[0][l] = c11;
          USFseaNF0ELsep[1][l] = -c21;
          USFseaNF0ELsep[2][l] = -c3;
        }
        USFseaNF0L[0][l][i] = c11-c21-c3;
        USFseaNF0L[1][l][i] = c12-c22+c5+c4;
  			//-----------------------------
        //QUARK CONTRIBUTION
        Mx1 = mass*x-fabs(PvalF0[i]);
        Mx2 = mass*x-omega0;
        Mx3 = mass*x;
        index1=0;
        index2=0;
        index3=0;
        for (pnum=0; pnum<=kNInt; pnum++) {
          p = (double)pnum/(double)kNInt*kdist;
          if (p<fabs(Mx1)) {
            index1 = pnum;
            intarrtemp1[pnum]=0;
            intarrtemp2[pnum]=0;
          } else {
            intarrtemp1[pnum]=intarr1[pnum];
            intarrtemp2[pnum]=intarr2[pnum];
          }
          if (p<fabs(Mx2)) {
            index2 = pnum;
            intarrtemp3[pnum]=0;
            intarrtemp4[pnum]=0;
          } else {
            intarrtemp3[pnum]=intarr1[pnum];
            intarrtemp4[pnum]=intarr2[pnum];
          }
          if (p<fabs(Mx3)) {
            index3 = pnum;
            intarrtemp5[pnum]=0;
          } else {
            intarrtemp5[pnum]=intarr1[pnum];
          }
        }
        c11=0;
        c21=0;
        c12=0;
        c22=0;
        c31=0;
        c3=0;
        c4=0;
        c5=0;
        if (kNInt>=index1+1) {
          c11 = -booles(intarrtemp1, kNInt)/(double)(kNInt);
          c12 = -PvalF0[i]/fabs(PvalF0[i])*Mx1*booles(intarrtemp2, kNInt)/(double)(kNInt);
        }
        if (kNInt>=index2+1) {
          c21 = -booles(intarrtemp3, kNInt)/(double)(kNInt);
          c22 = -PvalF0[i]/omega0*Mx2*booles(intarrtemp4, kNInt)/(double)(kNInt);

          c5 = lambda*lambda*PvalF0[i]/omega0/omega0/omega0*mass*x/2.0 * booles(intarrtemp4, kNInt)/(double)(kNInt);
          y2 = carr[index2+1]/kdist;
          y1 = carr[index2]/kdist;
          c3 = lambda*lambda/(2.0*omega0)*Mx2*interpol(y1, y2, Mx2, index2);

  				y1 =intarr2[index2]/kdist;
          y2 =intarr2[index2+1]/kdist;
          c4 = lambda*lambda*PvalF0[i]/(2.0*omega0*omega0)*fabs(Mx2)*interpol(y1, y2, Mx2, index2);
        }
        if (kNInt>=index3+1) {
          c31=-booles(intarrtemp5, kNInt)/(double)kNInt;
        }
        if (G<GELcut) {
          USFseaPF0EL[0][l] = c11-c21-c3;
          USFseaPF0EL[1][l] = c12-c22+c5-c4;
          USFseaPF0EL[2][l] = c31;

          USFseaPF0ELsep[0][l] = c11;
          USFseaPF0ELsep[1][l] = -c21;
          USFseaPF0ELsep[2][l] = -c3;
        }
        USFseaPF0L[0][l][i]=c11-c21-c3;
        USFseaPF0L[1][l][i]=c12-c22+c5-c4;
  			//-----------------------------
      }

      //------------------------------------
      if (G<GELcut) {
        //fprintf(stderr, "here\n");
        omega0 = sqrt(lambda*lambda+Pval[i]*Pval[i]);
        omega0F0 = sqrt(lambda*lambda+PvalF0[i]*PvalF0[i]);

        Esrule0  = (2.0*G+1.0)*(fabs(Pval[i]) - sqrt(Pval[i]*Pval[i]+lambda*lambda) + 0.5*lambda*lambda/sqrt(Pval[i]*Pval[i]+lambda*lambda))/mass/mass;
        Esrule0x = Esrule0;
        c1 = -(2.0*G+1.0)*(fabs(PvalF0[i]) - sqrt(PvalF0[i]*PvalF0[i]+lambda*lambda) + 0.5*lambda*lambda/sqrt(PvalF0[i]*PvalF0[i]+lambda*lambda))/mass/mass;
        Esrule0 += c1;
        Esrule0a = c1;

        EAsrule0 += Esrule0;
        EAsrule1 += Esrule1;

        EAsrule0x+= Esrule0x;
        EAsrule1x+= Esrule1x;

        EAsrule0a+= Esrule0a;
        EAsrule1a+= Esrule1a;

        sepsrule0x[0] += (2.0*G+1.0)*fabs(Pval[i])/mass/mass;
        sepsrule0a[0] += -(2.0*G+1.0)*fabs(PvalF0[i])/mass/mass;
        sepsrule0[0] += (2.0*G+1.0)*fabs(Pval[i])/mass/mass-(2.0*G+1.0)*fabs(PvalF0[i])/mass/mass;

        sepsrule0x[1] += -(2.0*G+1.0)*sqrt(Pval[i]*Pval[i]+lambda*lambda)/mass/mass;
        sepsrule0a[1] += (2.0*G+1.0)*sqrt(PvalF0[i]*PvalF0[i]+lambda*lambda)/mass/mass;
        sepsrule0[1]  += -(2.0*G+1.0)*sqrt(Pval[i]*Pval[i]+lambda*lambda)/mass/mass+(2.0*G+1.0)*sqrt(PvalF0[i]*PvalF0[i]+lambda*lambda)/mass/mass;

        sepsrule0x[2] += (2.0*G+1.0)*0.5*lambda*lambda/sqrt(Pval[i]*Pval[i]+lambda*lambda)/mass/mass;
        sepsrule0a[2] += -(2.0*G+1.0)*0.5*lambda*lambda/sqrt(PvalF0[i]*PvalF0[i]+lambda*lambda)/mass/mass;
        sepsrule0[2]  += (2.0*G+1.0)*0.5*lambda*lambda/sqrt(Pval[i]*Pval[i]+lambda*lambda)/mass/mass-(2.0*G+1.0)*0.5*lambda*lambda/sqrt(PvalF0[i]*PvalF0[i]+lambda*lambda)/mass/mass;

        for (l=0; l<=xNInt; l++) {
    			x=(double)l/(double)xNInt*xdist;
    			intarr1x[l]=x*(USFseaNEL[0][l]-USFseaNF0EL[0][l]+USFseaPEL[0][l]-USFseaPF0EL[0][l]);
          intarr2x[l]=x*(USFseaNEL[0][l]-USFseaNF0EL[0][l]);
          intarr3x[l]=x*(USFseaPEL[0][l]-USFseaPF0EL[0][l]);

    			intarr4x[l]=x*(USFseaNEL[1][l]-USFseaNF0EL[1][l]+USFseaPEL[1][l]-USFseaPF0EL[1][l]);
        }
        Esrule0l = 0.5*booles(intarr1x, xNInt)*xdist/(double)xNInt;
        Esrule0lN = 0.5*booles(intarr2x, xNInt)*xdist/(double)xNInt;
        Esrule0lP = 0.5*booles(intarr3x, xNInt)*xdist/(double)xNInt;

    		Esrule1l = 0.5*booles(intarr4x, xNInt)*xdist/(double)xNInt;

        for (l=0; l<=xNInt; l++) {
    			x=(double)l/(double)xNInt*xdist;
    			intarr1x[l]=x*(USFseaNEL[0][l]+USFseaPEL[0][l]);
          intarr2x[l]=x*(USFseaNEL[0][l]);
          intarr3x[l]=x*(USFseaPEL[0][l]);

    			intarr4x[l]=x*(USFseaNEL[1][l]+USFseaPEL[1][l]);
    		}
        Esrule0lx = 0.5*booles(intarr1x, xNInt)*xdist/(double)xNInt;
        Esrule0lNx = 0.5*booles(intarr2x, xNInt)*xdist/(double)xNInt;
        Esrule0lPx = 0.5*booles(intarr3x, xNInt)*xdist/(double)xNInt;

    		Esrule1lx = 0.5*booles(intarr4x, xNInt)*xdist/(double)xNInt;

        for (l=0; l<=xNInt; l++) {
          x=(double)l/(double)xNInt*xdist;
          intarr1x[l]=x*(USFseaNF0EL[0][l]+USFseaPF0EL[0][l]);
          intarr2x[l]=x*(USFseaNF0EL[0][l]);
          intarr3x[l]=x*USFseaPF0EL[0][l];

          intarr4x[l]=x*(USFseaNF0EL[1][l]+USFseaPF0EL[1][l]);
        }
        Esrule0la = -0.5*booles(intarr1x, xNInt)*xdist/(double)xNInt;
        Esrule0lNa = -0.5*booles(intarr2x, xNInt)*xdist/(double)xNInt;
        Esrule0lPa = -0.5*booles(intarr3x, xNInt)*xdist/(double)xNInt;

        Esrule1la = -0.5*booles(intarr4x, xNInt)*xdist/(double)xNInt;

        EAsrule0l += Esrule0l;
        EAsrule1l += Esrule1l;
        EAsrule0lx+= Esrule0lx;
        EAsrule1lx+= Esrule1lx;
        EAsrule0la+= Esrule0la;
        EAsrule1la+= Esrule1la;

        EAsrule0lN+= Esrule0lN;
        EAsrule0lP+= Esrule0lP;
        EAsrule0lNx+=Esrule0lNx;
        EAsrule0lPx+=Esrule0lPx;
        EAsrule0lNa+=Esrule0lNa;
        EAsrule0lPa+=Esrule0lPa;

        for (j=0; j<3; j++) {
          for (l=0; l<=xNInt; l++) {
            x=(double)l/(double)xNInt*xdist;
            intarr1x[l]=x*(USFseaNELsep[j][l]+USFseaPELsep[j][l]);
            intarr2x[l]=x*(USFseaNELsep[j][l]);
            intarr3x[l]=x*(USFseaPELsep[j][l]);
          }
          sepsrule0lx[j]+=0.5*booles(intarr1x, xNInt)*xdist/(double)xNInt;
          sepsrule0lNx[j]+=0.5*booles(intarr2x, xNInt)*xdist/(double)xNInt;
          sepsrule0lPx[j]+=0.5*booles(intarr3x, xNInt)*xdist/(double)xNInt;
          sepEsrule0lx[j]=0.5*booles(intarr1x, xNInt)*xdist/(double)xNInt;
          sepEsrule0lNx[j]=0.5*booles(intarr2x, xNInt)*xdist/(double)xNInt;
          sepEsrule0lPx[j]=0.5*booles(intarr3x, xNInt)*xdist/(double)xNInt;
        }
        for (j=0; j<3; j++) {
          for (l=0; l<=xNInt; l++) {
            x=(double)l/(double)xNInt*xdist;
            intarr1x[l]=x*(USFseaNF0ELsep[j][l]+USFseaPF0ELsep[j][l]);
            intarr2x[l]=x*(USFseaNF0ELsep[j][l]);
            intarr3x[l]=x*(USFseaPF0ELsep[j][l]);
          }
          sepsrule0la[j]+=-0.5*booles(intarr1x, xNInt)*xdist/(double)xNInt;
          sepsrule0lNa[j]+=-0.5*booles(intarr2x, xNInt)*xdist/(double)xNInt;
          sepsrule0lPa[j]+=-0.5*booles(intarr3x, xNInt)*xdist/(double)xNInt;
          sepEsrule0la[j]=0.5*booles(intarr1x, xNInt)*xdist/(double)xNInt;
          sepEsrule0lNa[j]=0.5*booles(intarr2x, xNInt)*xdist/(double)xNInt;
          sepEsrule0lPa[j]=0.5*booles(intarr3x, xNInt)*xdist/(double)xNInt;
        }
        for (j=0; j<3; j++) {
          for (l=0; l<=xNInt; l++) {
            x=(double)l/(double)xNInt*xdist;
            intarr1x[l]=x*(USFseaNELsep[j][l]-USFseaNF0ELsep[j][l] +USFseaPELsep[j][l]-USFseaPF0ELsep[j][l]);
            intarr2x[l]=x*(USFseaNELsep[j][l]-USFseaNF0ELsep[j][l]);
            intarr3x[l]=x*(USFseaPELsep[j][l]-USFseaPF0ELsep[j][l]);
          }
          sepsrule0l[j]+=0.5*booles(intarr1x, xNInt)*xdist/(double)xNInt;
          sepsrule0lN[j]+=0.5*booles(intarr2x, xNInt)*xdist/(double)xNInt;
          sepsrule0lP[j]+=0.5*booles(intarr3x, xNInt)*xdist/(double)xNInt;
          sepEsrule0l[j]=0.5*booles(intarr1x, xNInt)*xdist/(double)xNInt;
          sepEsrule0lN[j]=0.5*booles(intarr2x, xNInt)*xdist/(double)xNInt;
          sepEsrule0lP[j]=0.5*booles(intarr3x, xNInt)*xdist/(double)xNInt;
        }
        //-********************
        for (l=0; l<=xNInt; l++) {
          x=(double)l/(double)xNInt*xdist;
          if (x<=fabs(Pval[i])/mass) {
            index1=l;
            intarr1x[l]=0;
          } else {
            intarr1x[l]=(x-fabs(Pval[i])/mass)*(USFseaNEL[2][l]);
          }
          if (x<=omega0/mass) {
            index2=l;
            intarr2x[l]=0;
            intarr3x[l]=0;
          } else {
            intarr2x[l]=(x-omega0/mass)*(USFseaNEL[2][l]);
            intarr3x[l]=0.5*lambda*lambda/(omega0*mass)*USFseaNEL[2][l];
          }
        }
        Esrule0m = xdist*(0.5*booles(intarr1x, xNInt)/(double)xNInt - 0.5*booles(intarr2x, xNInt)/(double)xNInt - 0.5*booles(intarr3x, xNInt)/(double)xNInt);
        Esrule0mx= Esrule0m;
        Esrule0mN =Esrule0m;
        Esrule0mNx=Esrule0m;
        sepEsrule0m[0]=  xdist*0.5*booles(intarr1x, xNInt)/(double)xNInt;
        sepEsrule0m[1]= -xdist*0.5*booles(intarr2x, xNInt)/(double)xNInt;
        sepEsrule0m[2]= -xdist*0.5*booles(intarr3x, xNInt)/(double)xNInt;
        sepEsrule0mx[0] = sepEsrule0m[0];
        sepEsrule0mx[1] = sepEsrule0m[1];
        sepEsrule0mx[2] = sepEsrule0m[2];
        sepEsrule0mN[0] =sepEsrule0m[0];
        sepEsrule0mN[1] =sepEsrule0m[1];
        sepEsrule0mN[2] =sepEsrule0m[2];
        sepEsrule0mNx[0] =sepEsrule0m[0];
        sepEsrule0mNx[1] =sepEsrule0m[1];
        sepEsrule0mNx[2] =sepEsrule0m[2];
        //-********************
        for (l=0; l<=xNInt; l++) {
          x=(double)l/(double)xNInt*xdist;
          if (x<=fabs(PvalF0[i])/mass) {
            index1=l;
            intarr1x[l]=0;
          } else {
            intarr1x[l]=(x-fabs(PvalF0[i])/mass)*(USFseaNF0EL[2][l]);
          }
          if (x<=omega0F0/mass) {
            index2=l;
            intarr2x[l]=0;
            intarr3x[l]=0;
          } else {
            intarr2x[l]=(x-omega0F0/mass)*(USFseaNF0EL[2][l]);
            intarr3x[l]=0.5*lambda*lambda/(omega0F0*mass)*(USFseaNF0EL[2][l]);
          }
        }
        c2 = -xdist*(0.5*booles(intarr1x, xNInt)/(double)xNInt);
        c3 = xdist*0.5*booles(intarr2x, xNInt)/(double)xNInt;
        c4 = xdist*0.5*booles(intarr3x, xNInt)/(double)xNInt;
        c1 = c2 + c3 + c4;
        Esrule0m += c1;
        Esrule0ma = c1;
        Esrule0mN+= c1;
        Esrule0mNa= c1;
        sepEsrule0m[0]+=c2;
        sepEsrule0m[1]+=c3;
        sepEsrule0m[2]+=c4;
        sepEsrule0ma[0]=c2;
        sepEsrule0ma[1]=c3;
        sepEsrule0ma[2]=c4;
        sepEsrule0mN[0]+=c2;
        sepEsrule0mN[1]+=c3;
        sepEsrule0mN[2]+=c4;
        sepEsrule0mNa[0]=c2;
        sepEsrule0mNa[1]=c3;
        sepEsrule0mNa[2]=c4;


        //-********************
        for (l=0; l<=xNInt; l++) {
          x=(double)l/(double)xNInt*xdist;
          intarr1x[l]=(x+fabs(Pval[i])/mass)*(USFseaPEL[2][l]);
          intarr2x[l]=(x+omega0/mass)*(USFseaPEL[2][l]);
          intarr3x[l]=-lambda*lambda*0.5/omega0/mass*(USFseaPEL[2][l]);
        }
        c1 = xdist*0.5*booles(intarr1x, xNInt)/(double)xNInt;
        c2 = -xdist*0.5*booles(intarr2x, xNInt)/(double)xNInt;
        c3 = -xdist*0.5*booles(intarr3x, xNInt)/(double)xNInt;
        //c1 = xdist*(0.5*booles(intarr1x, xNInt)/(double)xNInt - 0.5*booles(intarr2x, xNInt)/(double)xNInt - 0.5*booles(intarr3x, xNInt)/(double)xNInt);
        Esrule0m += c1+c2+c3;
        Esrule0mx+= c1+c2+c3;
        Esrule0mP = c1+c2+c3;
        Esrule0mPx= c1+c2+c3;


        sepEsrule0m[0]+=c1;
        sepEsrule0m[1]+=c2;
        sepEsrule0m[2]+=c3;
        sepEsrule0mx[0]+=c1;
        sepEsrule0mx[1]+=c2;
        sepEsrule0mx[2]+=c3;

        sepEsrule0mP[0] =c1;
        sepEsrule0mP[1] =c2;
        sepEsrule0mP[2] =c3;
        sepEsrule0mPx[0] =c1;
        sepEsrule0mPx[1] =c2;
        sepEsrule0mPx[2] =c3;
        for (l=0; l<=xNInt; l++) {
          x=(double)l/(double)xNInt*xdist;
          if (x<fabs(Pval[i])/mass) {
            index1=l;
            intarr1x[l]=(fabs(Pval[i])/mass-x)*(USFseaPEL[2][l]);
          } else {
            intarr1x[l]=0;
          }
          if (x<fabs(omega0/mass)) {
            index2=l;
            intarr2x[l]=(omega0/mass-x)*(USFseaPEL[2][l]);
            intarr3x[l]=-0.5*lambda*lambda/omega0/mass*USFseaPEL[2][l];
          } else {
            intarr2x[l]=0;
            intarr3x[l]=0;
          }
        }
        c1 = xdist*0.5*booles(intarr1x, xNInt)/(double)xNInt;
        c2 = -xdist*0.5*booles(intarr2x, xNInt)/(double)xNInt;
        c3 = -xdist*0.5*booles(intarr3x, xNInt)/(double)xNInt;
        //c1 = xdist*(0.5*booles(intarr1x, xNInt)/(double)xNInt - 0.5*booles(intarr2x, xNInt)/(double)xNInt);
        Esrule0m += c1+c2+c3;
        Esrule0mx+= c1+c2+c3;
        Esrule0mP += c1+c2+c3;
        Esrule0mPx+= c1+c2+c3;
        sepEsrule0m[0]+=c1;
        sepEsrule0m[1]+=c2;
        sepEsrule0m[2]+=c3;
        sepEsrule0mx[0]+=c1;
        sepEsrule0mx[1]+=c2;
        sepEsrule0mx[2]+=c3;
        sepEsrule0mP[0] +=c1;
        sepEsrule0mP[1] +=c2;
        sepEsrule0mP[2] +=c3;
        sepEsrule0mPx[0] +=c1;
        sepEsrule0mPx[1] +=c2;
        sepEsrule0mPx[2] +=c3;
        //-********************
        for (l=0; l<=xNInt; l++) {
          x=(double)l/(double)xNInt*xdist;
          intarr1x[l]=(x+fabs(PvalF0[i])/mass)*(USFseaPF0EL[2][l]);
          intarr2x[l]=(x+omega0F0/mass)*(USFseaPF0EL[2][l]);
          intarr3x[l]=-lambda*lambda*0.5/omega0F0/mass*USFseaPF0EL[2][l];
        }
        c1 = -xdist*0.5*booles(intarr1x, xNInt)/(double)xNInt;
        c2 = xdist*0.5*booles(intarr2x, xNInt)/(double)xNInt;
        c3 = xdist*0.5*booles(intarr3x, xNInt)/(double)xNInt;
        //c1 = -xdist*(0.5*booles(intarr1x, xNInt)/(double)xNInt - 0.5*booles(intarr2x, xNInt)/(double)xNInt);
        Esrule0m += c1+c2+c3;
        Esrule0ma+= c1+c2+c3;
        Esrule0mP+= c1+c2+c3;
        Esrule0mPa=c1+c2+c3;
        sepEsrule0m[0]+=c1;
        sepEsrule0m[1]+=c2;
        sepEsrule0m[2]+=c3;
        sepEsrule0ma[0]+=c1;
        sepEsrule0ma[1]+=c2;
        sepEsrule0ma[2]+=c3;

        sepEsrule0mP[0]+=c1;
        sepEsrule0mP[1]+=c2;
        sepEsrule0mP[2]+=c3;
        sepEsrule0mPa[0]=c1;
        sepEsrule0mPa[1]=c2;
        sepEsrule0mPa[2]=c3;
        for (l=0; l<=xNInt; l++) {
          x=(double)l/(double)xNInt*xdist;
          if (x<fabs(PvalF0[i])/mass) {
            index1=l;
            intarr1x[l]=(fabs(PvalF0[i])/mass-x)*(USFseaPF0EL[2][l]);
          } else {
            intarr1x[l]=0;
          }
          if (x<fabs(omega0F0/mass)) {
            index2=l;
            intarr2x[l]=(omega0F0/mass-x)*(USFseaPF0EL[2][l]);
            intarr3x[l]=-0.5*lambda*lambda/omega0F0/mass*USFseaPF0EL[2][l];
          } else {
            intarr2x[l]=0;
            intarr3x[l]=0;
          }
        }
        c1 =-xdist*0.5*booles(intarr1x, xNInt)/(double)xNInt;
        c2 = xdist*0.5*booles(intarr2x, xNInt)/(double)xNInt;
        c3 = xdist*0.5*booles(intarr3x, xNInt)/(double)xNInt;
        //c1 = -xdist*(0.5*booles(intarr1x, xNInt)/(double)xNInt - 0.5*booles(intarr2x, xNInt)/(double)xNInt);
        Esrule0m += c1+c2+c3;
        Esrule0ma+= c1+c2+c3;
        Esrule0mP+= c1+c2+c3;
        Esrule0mPa+=c1+c2+c3;

        sepEsrule0m[0]+=c1;
        sepEsrule0m[1]+=c2;
        sepEsrule0m[2]+=c3;
        sepEsrule0ma[0]+=c1;
        sepEsrule0ma[1]+=c2;
        sepEsrule0ma[2]+=c3;

        sepEsrule0mP[0]+=c1;
        sepEsrule0mP[1]+=c2;
        sepEsrule0mP[2]+=c3;
        sepEsrule0mPa[0]+=c1;
        sepEsrule0mPa[1]+=c2;
        sepEsrule0mPa[2]+=c3;
        //-********************

        EAsrule0m += Esrule0m;
        EAsrule0mx+= Esrule0mx;
        EAsrule0ma+= Esrule0ma;

        EAsrule0mP += Esrule0mP;
        EAsrule0mPx+= Esrule0mPx;
        EAsrule0mPa+= Esrule0mPa;

        EAsrule0mN += Esrule0mN;
        EAsrule0mNx+= Esrule0mNx;
        EAsrule0mNa+= Esrule0mNa;
        for (j=0; j<3; j++) {
          sepsrule0m[j]+=sepEsrule0m[j];
          sepsrule0ma[j]+=sepEsrule0ma[j];
          sepsrule0mx[j]+=sepEsrule0mx[j];

          sepsrule0mP[j]+=sepEsrule0mP[j];
          sepsrule0mPa[j]+=sepEsrule0mPa[j];
          sepsrule0mPx[j]+=sepEsrule0mPx[j];

          sepsrule0mN[j]+=sepEsrule0mN[j];
          sepsrule0mNa[j]+=sepEsrule0mNa[j];
          sepsrule0mNx[j]+=sepEsrule0mNx[j];
        }
        fprintf(EsepSRfile, "%20.16f %20.16f %20.16f %20.16f %20.16f\n", -Esrule0, Esrule0l, -Esrule0m, 0.5*Esrule1, Esrule1l);
        fprintf(EsepSRxfile, "%20.16f %20.16f %20.16f %20.16f %20.16f\n", -Esrule0x, Esrule0lx, -Esrule0mx, 0.5*Esrule1x, Esrule1lx);
        fprintf(EsepSRafile, "%20.16f %20.16f %20.16f %20.16f %20.16f\n", -Esrule0a, Esrule0la, -Esrule0ma, 0.5*Esrule1a, Esrule1la);
        fprintf(EsepSRsepfile, "%20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n", sepEsrule0[0], sepEsrule0[1], sepEsrule0[2], sepEsrule0l[0], sepEsrule0l[1], sepEsrule0l[2], sepEsrule0m[0], sepEsrule0m[1], sepEsrule0m[2]);
        fprintf(EsepSRsepxfile, "%20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n", sepEsrule0x[0], sepEsrule0x[1], sepEsrule0x[2], sepEsrule0lx[0], sepEsrule0lx[1], sepEsrule0lx[2], sepEsrule0mx[0], sepEsrule0mx[1], sepEsrule0mx[2]);
        fprintf(EsepSRsepafile, "%20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n", sepEsrule0a[0], sepEsrule0a[1], sepEsrule0a[2], sepEsrule0la[0], sepEsrule0la[1], sepEsrule0la[2], sepEsrule0ma[0], sepEsrule0ma[1], sepEsrule0ma[2]);


        fprintf(EsepSRPNfile, "%20.16f %20.16f %20.16f %20.16f\n", Esrule0lP, -Esrule0mP, Esrule0lN, -Esrule0mN);
        fprintf(EsepSRPNxfile, "%20.16f %20.16f %20.16f %20.16f\n", Esrule0lPx, -Esrule0mPx, Esrule0lNx, -Esrule0mNx);
        fprintf(EsepSRPNafile, "%20.16f %20.16f %20.16f %20.16f\n", Esrule0lPa, -Esrule0mPa, Esrule0lNa, -Esrule0mNa);
        fprintf(EsepSRPNsepfile, "%20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n", sepEsrule0lP[0], sepEsrule0lP[1], sepEsrule0lP[2], sepEsrule0mP[0], sepEsrule0mP[1], sepEsrule0mP[2], sepEsrule0lN[0], sepEsrule0lN[1], sepEsrule0lN[2], sepEsrule0mN[0], sepEsrule0mN[1], sepEsrule0mN[2]);
        fprintf(EsepSRPNsepxfile, "%20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n", sepEsrule0lPx[0], sepEsrule0lPx[1], sepEsrule0lPx[2], sepEsrule0mPx[0], sepEsrule0mPx[1], sepEsrule0mPx[2], sepEsrule0lNx[0], sepEsrule0lNx[1], sepEsrule0lNx[2], sepEsrule0mNx[0], sepEsrule0mNx[1], sepEsrule0mNx[2]);
        fprintf(EsepSRPNsepafile, "%20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n", sepEsrule0lPa[0], sepEsrule0lPa[1], sepEsrule0lPa[2], sepEsrule0mPa[0], sepEsrule0mPa[1], sepEsrule0mPa[2], sepEsrule0lNa[0], sepEsrule0lNa[1], sepEsrule0lNa[2], sepEsrule0mNa[0], sepEsrule0mNa[1], sepEsrule0mNa[2]);
      }
    }
    if (G<GELcut) {
      fclose(EsepSRfile);
      fclose(EsepSRxfile);
      fclose(EsepSRafile);
      fclose(EsepSRsepfile);
      fclose(EsepSRsepafile);
      fclose(EsepSRsepxfile);
      fclose(EsepSRPNfile);
      fclose(EsepSRPNafile);
      fclose(EsepSRPNxfile);
      fclose(EsepSRPNsepfile);
      fclose(EsepSRPNsepafile);
      fclose(EsepSRPNsepxfile);
    }

    //-*******************************
    //ODD CHANNEL
    //-*******************************
    if (G<GELcut) {
      sprintf(filename, "Data/expnum%d/SRdata/GsepELSRfileOG%d_%d.txt", expnum, G, expnum);
      EsepSRfile = fopen(filename, "w+");
      sprintf(filename, "Data/expnum%d/SRdata/GsepELSRxfileOG%d_%d.txt", expnum, G, expnum);
      EsepSRxfile = fopen(filename, "w+");
      sprintf(filename, "Data/expnum%d/SRdata/GsepELSRafileOG%d_%d.txt", expnum, G, expnum);
      EsepSRafile = fopen(filename, "w+");

      sprintf(filename, "Data/expnum%d/SRdata/GsepELSRsepfileOG%d_%d.txt", expnum, G, expnum);
      EsepSRsepfile = fopen(filename, "w+");
      sprintf(filename, "Data/expnum%d/SRdata/GsepELSRsepafileOG%d_%d.txt", expnum, G, expnum);
      EsepSRsepafile = fopen(filename, "w+");
      sprintf(filename, "Data/expnum%d/SRdata/GsepELSRsepxfileOG%d_%d.txt", expnum, G, expnum);
      EsepSRsepxfile = fopen(filename, "w+");

      sprintf(filename, "Data/expnum%d/SRdata/GsepELSRPNfileOG%d_%d.txt", expnum, G, expnum);
      EsepSRPNfile = fopen(filename, "w+");
      sprintf(filename, "Data/expnum%d/SRdata/GsepELSRPNafileOG%d_%d.txt", expnum, G, expnum);
      EsepSRPNafile = fopen(filename, "w+");
      sprintf(filename, "Data/expnum%d/SRdata/GsepELSRPNxfileOG%d_%d.txt", expnum, G, expnum);
      EsepSRPNxfile = fopen(filename, "w+");

      sprintf(filename, "Data/expnum%d/SRdata/GsepELSRPNsepfileOG%d_%d.txt", expnum, G, expnum);
      EsepSRPNsepfile = fopen(filename, "w+");
      sprintf(filename, "Data/expnum%d/SRdata/GsepELSRPNsepafileOG%d_%d.txt", expnum, G, expnum);
      EsepSRPNsepafile = fopen(filename, "w+");
      sprintf(filename, "Data/expnum%d/SRdata/GsepELSRPNsepxfileOG%d_%d.txt", expnum, G, expnum);
      EsepSRPNsepxfile = fopen(filename, "w+");
    }
    for (i=0; i<tempN; i++) {
      Esrule1=0;
      Esrule1x=0;
      Esrule1a=0;
      //LOOP THROUGH x's
      for (l=0; l<=xNInt; l++) {
        x=(double)l/(double)xNInt*xdist;
        //INTEGRALS--------------------
        for (pnum=0; pnum<=kNInt; pnum++) {
          p=(double)pnum/(double)kNInt*kdist;
          if (G==0) {
            intarr1[pnum] =fourier_funcN[i][0][pnum]*fourier_funcN[i][0][pnum];
            intarr1[pnum]+=fourier_funcN[i][1][pnum]*fourier_funcN[i][1][pnum];
            carr[pnum]=intarr1[pnum]*kdist;
            intarr1[pnum]*=p*kdist;

            intarr2[pnum] =-2.0*kdist*fourier_funcN[i][0][pnum]*fourier_funcN[i][1][pnum];
          } else {
            intarr1[pnum] =fourier_funcN[i][0][pnum]*fourier_funcN[i][0][pnum];
            intarr1[pnum]+=fourier_funcN[i][1][pnum]*fourier_funcN[i][1][pnum];
            intarr1[pnum]+=fourier_funcN[i][2][pnum]*fourier_funcN[i][2][pnum];
            intarr1[pnum]+=fourier_funcN[i][3][pnum]*fourier_funcN[i][3][pnum];
            carr[pnum]=intarr1[pnum]*kdist*(G*2.0+1.0);
            intarr1[pnum]*=p*kdist*(2.0*G+1.0);

            intarr2[pnum] =fourier_funcN[i][0][pnum]*fourier_funcN[i][2][pnum];
            intarr2[pnum]+=fourier_funcN[i][1][pnum]*fourier_funcN[i][3][pnum];
            intarr2[pnum]*=-2.0*kdist*(2.0*G+1.0);
          }
  				if (l==0) {
  					temp = sqrt(Nval[i]*Nval[i]+lambda*lambda);
  					gpp[pnum] = 2.0*Nval[i]/fabs(Nval[i])*(1.0-fabs(Nval[i])/temp-0.5*fabs(Nval[i])*lambda*lambda/temp/temp/temp);
  					gpp[pnum] *= -1.0/(3.0*mass*mass)*p*p*p*intarr2[pnum];
  				}
        }
  			if (l==0) {
          c11 = booles(gpp, kNInt)/(double)kNInt;
  				Esrule1 += c11;
          srule1 += c11;
          Esrule1x += c11;
          srule1x += c11;
        }
  			omega0=sqrt(lambda*lambda+Nval[i]*Nval[i]);
        //-----------------------------
        //ANTI-QUARK CONTRIBUTION
        Mx1 = mass*x+fabs(Nval[i]);
        Mx2 = mass*x+omega0;
        Mx3 = mass*x;
        index1=0;
        index2=0;
        index3=0;
        for (pnum=0; pnum<=kNInt; pnum++) {
          p = (double)pnum/(double)kNInt*kdist;
          if (p<fabs(Mx1)) {
            index1 = pnum;
            intarrtemp1[pnum]=0;
            intarrtemp2[pnum]=0;
          } else {
            intarrtemp1[pnum]=intarr1[pnum];
            intarrtemp2[pnum]=intarr2[pnum];
          }
          if (p<fabs(Mx2)) {
            index2 = pnum;
            intarrtemp3[pnum]=0;
            intarrtemp4[pnum]=0;
          } else {
            intarrtemp3[pnum]=intarr1[pnum];
            intarrtemp4[pnum]=intarr2[pnum];
          }
          if (p<fabs(Mx3)) {
            index3=pnum;
            intarrtemp5[pnum]=0;
          } else {
            intarrtemp5[pnum]=intarr1[pnum];
          }
        }
        c11=0;
        c12=0;
        c21=0;
        c22=0;
        c31=0;
        c3=0;
        c4=0;
        c5=0;
        if (kNInt>=index1+1) {
          c11 = booles(intarrtemp1, kNInt)/(double)(kNInt);
          c12 = -Nval[i]/fabs(Nval[i])*Mx1*booles(intarrtemp2, kNInt)/(double)(kNInt);
        }
        if (kNInt>=index2+1) {
          c21 = booles(intarrtemp3, kNInt)/(double)(kNInt);
          c22 = -Nval[i]/omega0*Mx2*booles(intarrtemp4, kNInt)/(double)(kNInt);

          c5 = lambda*lambda*Nval[i]/omega0/omega0/omega0*mass*x/2.0 * booles(intarrtemp4, kNInt)/(double)(kNInt);
          y2 = carr[index2+1]/kdist;
          y1 = carr[index2]/kdist;
          c3 = lambda*lambda/(2.0*omega0)*Mx2*interpol(y1, y2, Mx2, index2);
          y1 =intarr2[index2]/kdist;
          y2 =intarr2[index2+1]/kdist;
          c4 = lambda*lambda*Nval[i]/(2.0*omega0*omega0)*fabs(Mx2)*interpol(y1, y2, Mx2, index2);
        }
        if (kNInt>=index3+1) {
          c31=booles(intarrtemp5, kNInt)/(double)kNInt;
        }
        if (G<GELcut) {
          USFseaNEL[0][l]=c11-c21-c3;
          USFseaNEL[1][l]=c12-c22+c5+c4;
          USFseaNEL[2][l]=c31;

          USFseaNELsep[0][l]=c11;
          USFseaNELsep[1][l]=-c21;
          USFseaNELsep[2][l]=-c3;
        }
     		USFseaNL[0][l][tempP+i]= c11-c21-c3;
        USFseaNL[1][l][tempP+i]= c12-c22+c5+c4;
  			//-----------------------------
        //QUARK CONTRIBUTION
        Mx1 = mass*x-fabs(Nval[i]);
        Mx2 = mass*x-omega0;
        Mx3 = mass*x;
        index1=0;
        index2=0;
        index3=0;
        for (pnum=0; pnum<=kNInt; pnum++) {
          p = (double)pnum/(double)kNInt*kdist;
          if (p<fabs(Mx1)) {
            index1 = pnum;
            intarrtemp1[pnum]=0;
            intarrtemp2[pnum]=0;
          } else {
            intarrtemp1[pnum]=intarr1[pnum];
            intarrtemp2[pnum]=intarr2[pnum];
          }
          if (p<fabs(Mx2)) {
            index2 = pnum;
            intarrtemp3[pnum]=0;
            intarrtemp4[pnum]=0;
          } else {
            intarrtemp3[pnum]=intarr1[pnum];
            intarrtemp4[pnum]=intarr2[pnum];
          }
          if (p<fabs(Mx3)) {
            index3=pnum;
            intarrtemp5[pnum]=0;
          } else {
            intarrtemp5[pnum]=intarr1[pnum];
          }
        }
        c11=0;
        c21=0;
        c12=0;
        c22=0;
        c31=0;
        c3=0;
        c4=0;
        c5=0;
        if (kNInt>=index1+1) {
          c11 = -booles(intarrtemp1, kNInt)/(double)(kNInt);
          c12 = -Nval[i]/fabs(Nval[i])*Mx1*booles(intarrtemp2, kNInt)/(double)(kNInt);
        }
        if (kNInt>=index2+1) {
          c21 = -booles(intarrtemp3, kNInt)/(double)(kNInt);
          c22 = -Nval[i]/omega0*Mx2*booles(intarrtemp4, kNInt)/(double)(kNInt);

          c5 = lambda*lambda*Nval[i]/omega0/omega0/omega0*mass*x/2.0 * booles(intarrtemp4, kNInt)/(double)(kNInt);

          y2 = carr[index2+1]/kdist;
          y1 = carr[index2]/kdist;
          c3 = lambda*lambda/(2.0*omega0)*Mx2*interpol(y1, y2, Mx2, index2);
          y1 =intarr2[index2]/kdist;
          y2 =intarr2[index2+1]/kdist;
          c4 = lambda*lambda*Nval[i]/(2.0*omega0*omega0)*fabs(Mx2)*interpol(y1, y2, Mx2, index2);
        }
        if (kNInt>=index3+1) {
          c31=-booles(intarrtemp5, kNInt)/(double)kNInt;
        }
        if (G<GELcut) {
          USFseaPEL[0][l]=c11-c21-c3;
    			USFseaPEL[1][l]=c12-c22+c5-c4;
          USFseaPEL[2][l]=c31;

          USFseaPELsep[0][l]=c11;
          USFseaPELsep[1][l]=-c21;
          USFseaPELsep[2][l]=-c3;
        }
        USFseaPL[0][l][tempP+i]=c11-c21-c3;
  			USFseaPL[1][l][tempP+i]=c12-c22+c5-c4;
        //-----------------------------
      }
      //LOOP THROUGH x's F0
      for (l=0; l<=xNInt; l++) {
        x=(double)l/(double)xNInt*xdist;

        //INTEGRALS--------------------
        for (pnum=0; pnum<=kNInt; pnum++) {
          p=(double)pnum/(double)kNInt*kdist;
          if (G==0) {
            intarr1[pnum] =fourier_funcNF0[i][0][pnum]*fourier_funcNF0[i][0][pnum];
            intarr1[pnum]+=fourier_funcNF0[i][1][pnum]*fourier_funcNF0[i][1][pnum];
            carr[pnum]=intarr1[pnum]*kdist;
            intarr1[pnum]*=p*kdist;

            intarr2[pnum] =-2.0*kdist*fourier_funcNF0[i][0][pnum]*fourier_funcNF0[i][1][pnum];
          } else {
            intarr1[pnum] =fourier_funcNF0[i][0][pnum]*fourier_funcNF0[i][0][pnum];
            intarr1[pnum]+=fourier_funcNF0[i][1][pnum]*fourier_funcNF0[i][1][pnum];
            intarr1[pnum]+=fourier_funcNF0[i][2][pnum]*fourier_funcNF0[i][2][pnum];
            intarr1[pnum]+=fourier_funcNF0[i][3][pnum]*fourier_funcNF0[i][3][pnum];
            carr[pnum]=intarr1[pnum]*kdist*(G*2.0+1.0);
            intarr1[pnum]*=p*kdist*(2.0*G+1.0);

            intarr2[pnum] =fourier_funcNF0[i][0][pnum]*fourier_funcNF0[i][2][pnum];
            intarr2[pnum]+=fourier_funcNF0[i][1][pnum]*fourier_funcNF0[i][3][pnum];
            intarr2[pnum]*=-2.0*kdist*(2.0*G+1.0);
          }
  				if (l==0) {
  					temp = sqrt(NvalF0[i]*NvalF0[i]+lambda*lambda);
  					gp0[pnum] = 2.0*NvalF0[i]/fabs(NvalF0[i])*(1.0-fabs(NvalF0[i])/temp-0.5*fabs(NvalF0[i])*lambda*lambda/temp/temp/temp);
  					gp0[pnum] *= -1.0/(3.0*mass*mass)*p*p*p*intarr2[pnum];
  				}
        }
  			if (l==0) {
          c11 = -booles(gp0, kNInt)/(double)kNInt;
  				Esrule1 += c11;
          srule1 += c11;
          Esrule1a += c11;
          srule1a += c11;
  			}
        omega0=sqrt(lambda*lambda+NvalF0[i]*NvalF0[i]);
        //-----------------------------
        //ANTI-QUARK CONTRIBUTION
        Mx1 = mass*x+fabs(NvalF0[i]);
        Mx2 = mass*x+omega0;
        Mx3 = mass*x;
        index1=0;
        index2=0;
        index3=0;
        for (pnum=0; pnum<=kNInt; pnum++) {
          p = (double)pnum/(double)kNInt*kdist;
          if (p<fabs(Mx1)) {
            index1 = pnum;
            intarrtemp1[pnum]=0;
            intarrtemp2[pnum]=0;
          } else {
            intarrtemp1[pnum]=intarr1[pnum];
            intarrtemp2[pnum]=intarr2[pnum];
          }
          if (p<fabs(Mx2)) {
            index2 = pnum;
            intarrtemp3[pnum]=0;
            intarrtemp4[pnum]=0;
          } else {
            intarrtemp3[pnum]=intarr1[pnum];
            intarrtemp4[pnum]=intarr2[pnum];
          }
          if (p<fabs(Mx3)) {
            index3 = pnum;
            intarrtemp5[pnum]=0;
          } else {
            intarrtemp5[pnum]=intarr1[pnum];
          }
        }
        c11=0;
        c21=0;
        c12=0;
        c22=0;
        c31=0;
        c3=0;
        c4=0;
        c5=0;
        if (kNInt>=index1+1) {
          c11 = booles(intarrtemp1, kNInt)/(double)(kNInt);
          c12 = -NvalF0[i]/fabs(NvalF0[i])*Mx1*booles(intarrtemp2, kNInt)/(double)(kNInt);
        }
        if (kNInt>=index2+1) {
          c21 = booles(intarrtemp3, kNInt)/(double)(kNInt);
          c22 = -NvalF0[i]/omega0*Mx2*booles(intarrtemp4, kNInt)/(double)(kNInt);

          c5 = lambda*lambda*NvalF0[i]/omega0/omega0/omega0*mass*x/2.0 * booles(intarrtemp4, kNInt)/(double)(kNInt);
          y2 = carr[index2+1]/kdist;
          y1 = carr[index2]/kdist;
          c3 = lambda*lambda/(2.0*omega0)*Mx2*interpol(y1, y2, Mx2, index2);

  				y1 =intarr2[index2]/kdist;
          y2 =intarr2[index2+1]/kdist;
          c4 = lambda*lambda*NvalF0[i]/(2.0*omega0*omega0)*fabs(Mx2)*interpol(y1, y2, Mx2, index2);
        }
        if (kNInt>=index3+1) {
          c31=booles(intarrtemp5, kNInt)/(double)kNInt;
        }
        if (G<GELcut) {
          USFseaNF0EL[0][l]= c11-c21-c3;
    			USFseaNF0EL[1][l]= c12-c22+c5+c4;
          USFseaNF0EL[2][l]= c31;

          USFseaNF0ELsep[0][l]= c11;
          USFseaNF0ELsep[1][l]= -c21;
          USFseaNF0ELsep[2][l]= -c3;
        }
      	USFseaNF0L[0][l][tempP+i]= c11-c21-c3;
  			USFseaNF0L[1][l][tempP+i]= c12-c22+c5+c4;
        //-----------------------------
        //QUARK CONTRIBUTION
        Mx1 = mass*x-fabs(NvalF0[i]);
        Mx2 = mass*x-omega0;
        Mx3 = mass*x;
        index1=0;
        index2=0;
        index3=0;
        for (pnum=0; pnum<=kNInt; pnum++) {
          p = (double)pnum/(double)kNInt*kdist;
          if (p<fabs(Mx1)) {
            index1 = pnum;
            intarrtemp1[pnum]=0;
            intarrtemp2[pnum]=0;
          } else {
            intarrtemp1[pnum]=intarr1[pnum];
            intarrtemp2[pnum]=intarr2[pnum];
          }
          if (p<fabs(Mx2)) {
            index2 = pnum;
            intarrtemp3[pnum]=0;
            intarrtemp4[pnum]=0;
          } else {
            intarrtemp3[pnum]=intarr1[pnum];
            intarrtemp4[pnum]=intarr2[pnum];
          }
          if (p<fabs(Mx3)) {
            index3=0;
            intarrtemp5[pnum]=0;
          } else {
            intarrtemp5[pnum]=intarr1[pnum];
          }
        }
        c11=0;
        c21=0;
        c12=0;
        c22=0;
        c31=0;
        c3=0;
        c4=0;
        c5=0;
        if (kNInt>=index1+1) {
          c11 = -booles(intarrtemp1, kNInt)/(double)(kNInt);
          c12 = -NvalF0[i]/fabs(NvalF0[i])*Mx1*booles(intarrtemp2, kNInt)/(double)(kNInt);
        }
        if (kNInt>=index2+1) {
          c21 = -booles(intarrtemp3, kNInt)/(double)(kNInt);
          c22 = -NvalF0[i]/omega0*Mx2*booles(intarrtemp4, kNInt)/(double)(kNInt);

          c5 = lambda*lambda*NvalF0[i]/omega0/omega0/omega0*mass*x/2.0 * booles(intarrtemp4, kNInt)/(double)(kNInt);
          y2 = carr[index2+1]/kdist;
          y1 = carr[index2]/kdist;
          c3 = lambda*lambda/(2.0*omega0)*Mx2*interpol(y1, y2, Mx2, index2);

  				y1 =intarr2[index2]/kdist;
          y2 =intarr2[index2+1]/kdist;
          c4 = lambda*lambda*NvalF0[i]/(2.0*omega0*omega0)*fabs(Mx2)*interpol(y1, y2, Mx2, index2);
        }
        if (kNInt>=index3+1) {
          c31 = -booles(intarrtemp5, kNInt)/(double)kNInt;
        }
        if (G<GELcut) {
          USFseaPF0EL[0][l]=c11-c21-c3;
    			USFseaPF0EL[1][l]=c12-c22+c5-c4;
          USFseaPF0EL[2][l]=c31;

          USFseaPF0ELsep[0][l]=c11;
          USFseaPF0ELsep[1][l]=-c21;
          USFseaPF0ELsep[2][l]=-c3;
        }
        USFseaPF0L[0][l][tempP+i]=c11-c21-c3;
  			USFseaPF0L[1][l][tempP+i]=c12-c22+c5-c4;
        //-----------------------------
      }
      //------------------------------------
      if (G<GELcut) {
        //fprintf(stderr, "here\n");
        omega0 = sqrt(lambda*lambda+Nval[i]*Nval[i]);
        omega0F0 = sqrt(lambda*lambda+NvalF0[i]*NvalF0[i]);

        Esrule0  = (2.0*G+1.0)*(fabs(Nval[i]) - sqrt(Nval[i]*Nval[i]+lambda*lambda) + 0.5*lambda*lambda/sqrt(Nval[i]*Nval[i]+lambda*lambda))/mass/mass;
        Esrule0x = Esrule0;
        c1 = -(2.0*G+1.0)*(fabs(NvalF0[i]) - sqrt(NvalF0[i]*NvalF0[i]+lambda*lambda) + 0.5*lambda*lambda/sqrt(NvalF0[i]*NvalF0[i]+lambda*lambda))/mass/mass;
        Esrule0 += c1;
        Esrule0a = c1;

        EAsrule0 += Esrule0;
        EAsrule1 += Esrule1;

        EAsrule0x+= Esrule0x;
        EAsrule1x+= Esrule1x;

        EAsrule0a+= Esrule0a;
        EAsrule1a+= Esrule1a;

        sepsrule0x[0] += (2.0*G+1.0)*fabs(Nval[i])/mass/mass;
        sepsrule0a[0] += -(2.0*G+1.0)*fabs(NvalF0[i])/mass/mass;
        sepsrule0[0] += (2.0*G+1.0)*fabs(Nval[i])/mass/mass-(2.0*G+1.0)*fabs(NvalF0[i])/mass/mass;

        sepsrule0x[1] += -(2.0*G+1.0)*sqrt(Nval[i]*Nval[i]+lambda*lambda)/mass/mass;
        sepsrule0a[1] += (2.0*G+1.0)*sqrt(NvalF0[i]*NvalF0[i]+lambda*lambda)/mass/mass;
        sepsrule0[1]  += -(2.0*G+1.0)*sqrt(Nval[i]*Nval[i]+lambda*lambda)/mass/mass+(2.0*G+1.0)*sqrt(NvalF0[i]*NvalF0[i]+lambda*lambda)/mass/mass;

        sepsrule0x[2] += (2.0*G+1.0)*0.5*lambda*lambda/sqrt(Nval[i]*Nval[i]+lambda*lambda)/mass/mass;
        sepsrule0a[2] += -(2.0*G+1.0)*0.5*lambda*lambda/sqrt(NvalF0[i]*NvalF0[i]+lambda*lambda)/mass/mass;
        sepsrule0[2]  += (2.0*G+1.0)*0.5*lambda*lambda/sqrt(Nval[i]*Nval[i]+lambda*lambda)/mass/mass-(2.0*G+1.0)*0.5*lambda*lambda/sqrt(NvalF0[i]*NvalF0[i]+lambda*lambda)/mass/mass;
        for (l=0; l<=xNInt; l++) {
          x=(double)l/(double)xNInt*xdist;
          intarr1x[l]=x*(USFseaNEL[0][l]-USFseaNF0EL[0][l]+USFseaPEL[0][l]-USFseaPF0EL[0][l]);
          intarr2x[l]=x*(USFseaNEL[0][l]-USFseaNF0EL[0][l]);
          intarr3x[l]=x*(USFseaPEL[0][l]-USFseaPF0EL[0][l]);

          intarr4x[l]=x*(USFseaNEL[1][l]-USFseaNF0EL[1][l]+USFseaPEL[1][l]-USFseaPF0EL[1][l]);
        }
        Esrule0l = 0.5*booles(intarr1x, xNInt)*xdist/(double)xNInt;
        Esrule0lN = 0.5*booles(intarr2x, xNInt)*xdist/(double)xNInt;
        Esrule0lP = 0.5*booles(intarr3x, xNInt)*xdist/(double)xNInt;

        Esrule1l = 0.5*booles(intarr4x, xNInt)*xdist/(double)xNInt;

        for (l=0; l<=xNInt; l++) {
          x=(double)l/(double)xNInt*xdist;
          intarr1x[l]=x*(USFseaNEL[0][l]+USFseaPEL[0][l]);
          intarr2x[l]=x*USFseaNEL[0][l];
          intarr3x[l]=x*USFseaPEL[0][l];
          intarr4x[l]=x*(USFseaNEL[1][l]+USFseaPEL[1][l]);
        }
        Esrule0lx = 0.5*booles(intarr1x, xNInt)*xdist/(double)xNInt;
        Esrule0lNx = 0.5*booles(intarr2x, xNInt)*xdist/(double)xNInt;
        Esrule0lPx = 0.5*booles(intarr3x, xNInt)*xdist/(double)xNInt;

        Esrule1lx = 0.5*booles(intarr4x, xNInt)*xdist/(double)xNInt;

        for (l=0; l<=xNInt; l++) {
          x=(double)l/(double)xNInt*xdist;
          intarr1x[l]=x*(USFseaNF0EL[0][l]+USFseaPF0EL[0][l]);
          intarr2x[l]=x*USFseaNF0EL[0][l];
          intarr3x[l]=x*USFseaPF0EL[0][l];

          intarr4x[l]=x*(USFseaNF0EL[1][l]+USFseaPF0EL[1][l]);
        }
        Esrule0la = -0.5*booles(intarr1x, xNInt)*xdist/(double)xNInt;
        Esrule0lNa = -0.5*booles(intarr2x, xNInt)*xdist/(double)xNInt;
        Esrule0lPa = -0.5*booles(intarr3x, xNInt)*xdist/(double)xNInt;

        Esrule1la = -0.5*booles(intarr4x, xNInt)*xdist/(double)xNInt;

        EAsrule0l += Esrule0l;
        EAsrule1l += Esrule1l;
        EAsrule0lx+= Esrule0lx;
        EAsrule1lx+= Esrule1lx;
        EAsrule0la+= Esrule0la;
        EAsrule1la+= Esrule1la;

        EAsrule0lN+= Esrule0lN;
        EAsrule0lP+= Esrule0lP;
        EAsrule0lNx+=Esrule0lNx;
        EAsrule0lPx+=Esrule0lPx;
        EAsrule0lNa+=Esrule0lNa;
        EAsrule0lPa+=Esrule0lPa;

        for (j=0; j<3; j++) {
          for (l=0; l<=xNInt; l++) {
            x=(double)l/(double)xNInt*xdist;
            intarr1x[l]=x*(USFseaNELsep[j][l]+USFseaPELsep[j][l]);
            intarr2x[l]=x*USFseaNELsep[j][l];
            intarr3x[l]=x*USFseaPELsep[j][l];
          }
          sepsrule0lx[j]+=0.5*booles(intarr1x, xNInt)*xdist/(double)xNInt;
          sepsrule0lNx[j]+=0.5*booles(intarr2x, xNInt)*xdist/(double)xNInt;
          sepsrule0lPx[j]+=0.5*booles(intarr3x, xNInt)*xdist/(double)xNInt;
          sepEsrule0lx[j]=0.5*booles(intarr1x, xNInt)*xdist/(double)xNInt;
          sepEsrule0lNx[j]=0.5*booles(intarr2x, xNInt)*xdist/(double)xNInt;
          sepEsrule0lPx[j]=0.5*booles(intarr3x, xNInt)*xdist/(double)xNInt;
        }
        for (j=0; j<3; j++) {
          for (l=0; l<=xNInt; l++) {
            x=(double)l/(double)xNInt*xdist;
            intarr1x[l]=x*(USFseaNF0ELsep[j][l]+USFseaPF0ELsep[j][l]);
            intarr2x[l]=x*USFseaNF0ELsep[j][l];
            intarr3x[l]=x*USFseaPF0ELsep[j][l];
          }
          sepsrule0la[j]+=-0.5*booles(intarr1x, xNInt)*xdist/(double)xNInt;
          sepsrule0lNa[j]+=-0.5*booles(intarr2x, xNInt)*xdist/(double)xNInt;
          sepsrule0lPa[j]+=-0.5*booles(intarr3x, xNInt)*xdist/(double)xNInt;
          sepEsrule0la[j]=-0.5*booles(intarr1x, xNInt)*xdist/(double)xNInt;
          sepEsrule0lNa[j]=-0.5*booles(intarr2x, xNInt)*xdist/(double)xNInt;
          sepEsrule0lPa[j]=-0.5*booles(intarr3x, xNInt)*xdist/(double)xNInt;
        }
        for (j=0; j<3; j++) {
          for (l=0; l<=xNInt; l++) {
            x=(double)l/(double)xNInt*xdist;
            intarr1x[l]=x*(USFseaNELsep[j][l]-USFseaNF0ELsep[j][l] +USFseaPELsep[j][l]-USFseaPF0ELsep[j][l]);
            intarr2x[l]=x*(USFseaNELsep[j][l]-USFseaNF0ELsep[j][l]);
            intarr3x[l]=x*(USFseaPELsep[j][l]-USFseaPF0ELsep[j][l]);
          }
          sepsrule0l[j]+=0.5*booles(intarr1x, xNInt)*xdist/(double)xNInt;
          sepsrule0lN[j]+=0.5*booles(intarr2x, xNInt)*xdist/(double)xNInt;
          sepsrule0lP[j]+=0.5*booles(intarr3x, xNInt)*xdist/(double)xNInt;
          sepEsrule0l[j]=0.5*booles(intarr1x, xNInt)*xdist/(double)xNInt;
          sepEsrule0lN[j]=0.5*booles(intarr2x, xNInt)*xdist/(double)xNInt;
          sepEsrule0lP[j]=0.5*booles(intarr3x, xNInt)*xdist/(double)xNInt;
        }
        //********************
        for (l=0; l<=xNInt; l++) {
          x=(double)l/(double)xNInt*xdist;
          if (x<=fabs(Nval[i])/mass) {
            index1=l;
            intarr1x[l]=0;
          } else {
            intarr1x[l]=(x-fabs(Nval[i])/mass)*(USFseaNEL[2][l]);
          }
          if (x<=omega0/mass) {
            index2=l;
            intarr2x[l]=0;
            intarr3x[l]=0;
          } else {
            intarr2x[l]=(x-omega0/mass)*(USFseaNEL[2][l]);
            intarr3x[l]=0.5*lambda*lambda/(omega0*mass)*USFseaNEL[2][l];
          }
        }
        c1 = xdist*0.5*booles(intarr1x, xNInt)/(double)xNInt;
        c2 = -xdist*0.5*booles(intarr2x, xNInt)/(double)xNInt;
        c3 = -xdist*0.5*booles(intarr3x, xNInt)/(double)xNInt;
        Esrule0m=c1+c2+c3;
        Esrule0mx= Esrule0m;
        Esrule0mN=c1+c2+c3;
        Esrule0mNx=Esrule0m;

        sepEsrule0m[0]=c1;
        sepEsrule0m[1]=c2;
        sepEsrule0m[2]=c3;
        sepEsrule0mx[0]=c1;
        sepEsrule0mx[1]=c2;
        sepEsrule0mx[2]=c3;

        sepEsrule0mN[0] =c1;
        sepEsrule0mN[1] =c2;
        sepEsrule0mN[2] =c3;
        sepEsrule0mNx[0]=c1;
        sepEsrule0mNx[1]=c2;
        sepEsrule0mNx[2]=c3;

        //********************
        for (l=0; l<=xNInt; l++) {
          x=(double)l/(double)xNInt*xdist;
          if (x<=fabs(NvalF0[i])/mass) {
            index1=l;
            intarr1x[l]=0;
          } else {
            intarr1x[l]=(x-fabs(NvalF0[i])/mass)*(USFseaNF0EL[2][l]);
          }
          if (x<=omega0F0/mass) {
            index2=l;
            intarr2x[l]=0;
            intarr3x[l]=0;
          } else {
            intarr2x[l]=(x-omega0F0/mass)*(USFseaNF0EL[2][l]);
            intarr3x[l]=0.5*lambda*lambda/(omega0F0*mass)*USFseaNF0EL[2][l];
          }
        }
        c1 = -xdist*0.5*booles(intarr1x, xNInt)/(double)xNInt;
        c2 = xdist*0.5*booles(intarr2x, xNInt)/(double)xNInt;
        c3 = xdist*0.5*booles(intarr3x, xNInt)/(double)xNInt;
        //c1 = -xdist*(0.5*booles(intarr1x, xNInt)/(double)xNInt - 0.5*booles(intarr2x, xNInt)/(double)xNInt);
        Esrule0m += c1+c2+c3;
        Esrule0ma = c1+c2+c3;
        Esrule0mN+=c1+c2+c3;
        Esrule0mNa=c1+c2+c3;


        sepEsrule0m[0]+=c1;
        sepEsrule0m[1]+=c2;
        sepEsrule0m[2]+=c3;
        sepEsrule0ma[0]=c1;
        sepEsrule0ma[1]=c2;
        sepEsrule0ma[2]=c3;

        sepEsrule0mN[0]+=c1;
        sepEsrule0mN[1]+=c2;
        sepEsrule0mN[2]+=c3;
        sepEsrule0mNa[0]=c1;
        sepEsrule0mNa[1]=c2;
        sepEsrule0mNa[2]=c3;
        //********************
        for (l=0; l<=xNInt; l++) {
          x=(double)l/(double)xNInt*xdist;
          intarr1x[l]=(x+fabs(Nval[i])/mass)*(USFseaPEL[2][l]);
          intarr2x[l]=(x+omega0/mass)*(USFseaPEL[2][l]);
          intarr3x[l]=-lambda*lambda*0.5/omega0/mass*USFseaPEL[2][l];
        }
        c1= xdist*0.5*booles(intarr1x, xNInt)/(double)xNInt;
        c2= -xdist*0.5*booles(intarr2x, xNInt)/(double)xNInt;
        c3= -xdist*0.5*booles(intarr3x, xNInt)/(double)xNInt;

        //c1 = xdist*(0.5*booles(intarr1x, xNInt)/(double)xNInt - 0.5*booles(intarr2x, xNInt)/(double)xNInt);
        Esrule0m +=c1+c2+c3;
        Esrule0mx +=c1+c2+c3;
        Esrule0mP =c1+c2+c3;
        Esrule0mPx=c1+c2+c3;

        sepEsrule0m[0]+=c1;
        sepEsrule0m[1]+=c2;
        sepEsrule0m[2]+=c3;
        sepEsrule0mx[0]+=c1;
        sepEsrule0mx[1]+=c2;
        sepEsrule0mx[2]+=c3;

        sepEsrule0mP[0]=c1;
        sepEsrule0mP[1]=c2;
        sepEsrule0mP[2]=c3;
        sepEsrule0mPx[0]=c1;
        sepEsrule0mPx[1]=c2;
        sepEsrule0mPx[2]=c3;
        for (l=0; l<=xNInt; l++) {
          x=(double)l/(double)xNInt*xdist;
          if (x<fabs(Nval[i])/mass) {
            index1=l;
            intarr1x[l]=(fabs(Nval[i])/mass-x)*(USFseaPEL[2][l]);
          } else {
            intarr1x[l]=0;
          }
          if (x<fabs(omega0/mass)) {
            index2=l;
            intarr2x[l]=(omega0/mass-x)*(USFseaPEL[2][l]);
            intarr3x[l]=(-0.5*lambda*lambda/omega0/mass)*USFseaPEL[2][l];
          } else {
            intarr2x[l]=0;
            intarr3x[l]=0;
          }
        }
        c1= xdist*0.5*booles(intarr1x, xNInt)/(double)xNInt;
        c2= -xdist*0.5*booles(intarr2x, xNInt)/(double)xNInt;
        c3= -xdist*0.5*booles(intarr3x, xNInt)/(double)xNInt;

        //c1 = xdist*(0.5*booles(intarr1x, xNInt)/(double)xNInt - 0.5*booles(intarr2x, xNInt)/(double)xNInt);
        Esrule0m += c1+c2+c3;
        Esrule0mx+= c1+c2+c3;
        Esrule0mP += c1+c2+c3;
        Esrule0mPx+= c1+c2+c3;

        sepEsrule0m[0]+=c1;
        sepEsrule0m[1]+=c2;
        sepEsrule0m[2]+=c3;
        sepEsrule0mx[0]+=c1;
        sepEsrule0mx[1]+=c2;
        sepEsrule0mx[2]+=c3;

        sepEsrule0mP[0]+=c1;
        sepEsrule0mP[1]+=c2;
        sepEsrule0mP[2]+=c3;
        sepEsrule0mPx[0]+=c1;
        sepEsrule0mPx[1]+=c2;
        sepEsrule0mPx[2]+=c3;
        //********************
        for (l=0; l<=xNInt; l++) {
          x=(double)l/(double)xNInt*xdist;
          intarr1x[l]=(x+fabs(NvalF0[i])/mass)*(USFseaPF0EL[2][l]);
          intarr2x[l]=(x+omega0F0/mass)*(USFseaPF0EL[2][l]);
          intarr3x[l]=-lambda*lambda*0.5/omega0F0/mass*(USFseaPF0EL[2][l]);
        }
        c1= -xdist*0.5*booles(intarr1x, xNInt)/(double)xNInt;
        c2= xdist*0.5*booles(intarr2x, xNInt)/(double)xNInt;
        c3= xdist*0.5*booles(intarr3x, xNInt)/(double)xNInt;
        //c1 = -xdist*(0.5*booles(intarr1x, xNInt)/(double)xNInt - 0.5*booles(intarr2x, xNInt)/(double)xNInt);
        Esrule0m += c1+c2+c3;
        Esrule0ma+= c1+c2+c3;
        Esrule0mP+=c1+c2+c3;
        Esrule0mPa=c1+c2+c3;


        sepEsrule0m[0]+=c1;
        sepEsrule0m[1]+=c2;
        sepEsrule0m[2]+=c3;
        sepEsrule0ma[0]+=c1;
        sepEsrule0ma[1]+=c2;
        sepEsrule0ma[2]+=c3;

        sepEsrule0mP[0]+=c1;
        sepEsrule0mP[1]+=c2;
        sepEsrule0mP[2]+=c3;
        sepEsrule0mPa[0]=c1;
        sepEsrule0mPa[1]=c2;
        sepEsrule0mPa[2]=c3;
        for (l=0; l<=xNInt; l++) {
          x=(double)l/(double)xNInt*xdist;
          if (x<fabs(NvalF0[i])/mass) {
            index1=l;
            intarr1x[l]=(fabs(NvalF0[i])/mass-x)*(USFseaPF0EL[2][l]);
          } else {
            intarr1x[l]=0;
          }
          if (x<fabs(omega0F0/mass)) {
            index2=l;
            intarr2x[l]=(omega0F0/mass-x)*(USFseaPF0EL[2][l]);
            intarr3x[l]= -0.5*lambda*lambda/omega0F0/mass*USFseaPF0EL[2][l];
          } else {
            intarr2x[l]=0;
            intarr3x[l]=0;
          }
        }
        c1=-xdist*0.5*booles(intarr1x, xNInt)/(double)xNInt;
        c2= xdist*0.5*booles(intarr2x, xNInt)/(double)xNInt;
        c3= xdist*0.5*booles(intarr3x, xNInt)/(double)xNInt;
        //c1 = -xdist*(0.5*booles(intarr1x, xNInt)/(double)xNInt - 0.5*booles(intarr2x, xNInt)/(double)xNInt);
        Esrule0m += c1+c2+c3;
        Esrule0ma+= c1+c2+c3;
        Esrule0mP+=c1+c2+c3;
        Esrule0mPa+=c1+c2+c3;

        sepEsrule0m[0]+=c1;
        sepEsrule0m[1]+=c2;
        sepEsrule0m[2]+=c3;
        sepEsrule0ma[0]+=c1;
        sepEsrule0ma[1]+=c2;
        sepEsrule0ma[2]+=c3;

        sepEsrule0mP[0]+=c1;
        sepEsrule0mP[1]+=c2;
        sepEsrule0mP[2]+=c3;
        sepEsrule0mPx[0]+=c1;
        sepEsrule0mPx[1]+=c2;
        sepEsrule0mPx[2]+=c3;
        //********************

        EAsrule0m += Esrule0m;
        EAsrule0mx+= Esrule0mx;
        EAsrule0ma+= Esrule0ma;

        EAsrule0mP += Esrule0mP;
        EAsrule0mPx+= Esrule0mPx;
        EAsrule0mPa+= Esrule0mPa;

        EAsrule0mN += Esrule0mN;
        EAsrule0mNx+= Esrule0mNx;
        EAsrule0mNa+= Esrule0mNa;

        for (j=0; j<3; j++) {
          sepsrule0m[j] +=sepEsrule0m[j];
          sepsrule0ma[j]+=sepEsrule0ma[j];
          sepsrule0mx[j]+=sepEsrule0mx[j];

          sepsrule0mN[j] +=sepEsrule0mN[j];
          sepsrule0mNa[j]+=sepEsrule0mNa[j];
          sepsrule0mNx[j]+=sepEsrule0mNx[j];

          sepsrule0mP[j] +=sepEsrule0mP[j];
          sepsrule0mPa[j]+=sepEsrule0mPa[j];
          sepsrule0mPx[j]+=sepEsrule0mPx[j];
        }

        fprintf(EsepSRfile, "%20.16f %20.16f %20.16f %20.16f %20.16f\n", -Esrule0, Esrule0l, -Esrule0m, 0.5*Esrule1, Esrule1l);
        fprintf(EsepSRxfile, "%20.16f %20.16f %20.16f %20.16f %20.16f\n", -Esrule0x, Esrule0lx, -Esrule0mx, 0.5*Esrule1x, Esrule1lx);
        fprintf(EsepSRafile, "%20.16f %20.16f %20.16f %20.16f %20.16f\n", -Esrule0a, Esrule0la, -Esrule0ma, 0.5*Esrule1a, Esrule1la);
        fprintf(EsepSRsepfile, "%20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n", sepEsrule0[0], sepEsrule0[1], sepEsrule0[2], sepEsrule0l[0], sepEsrule0l[1], sepEsrule0l[2], sepEsrule0m[0], sepEsrule0m[1], sepEsrule0m[2]);
        fprintf(EsepSRsepxfile, "%20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n", sepEsrule0x[0], sepEsrule0x[1], sepEsrule0x[2], sepEsrule0lx[0], sepEsrule0lx[1], sepEsrule0lx[2], sepEsrule0mx[0], sepEsrule0mx[1], sepEsrule0mx[2]);
        fprintf(EsepSRsepafile, "%20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n", sepEsrule0a[0], sepEsrule0a[1], sepEsrule0a[2], sepEsrule0la[0], sepEsrule0la[1], sepEsrule0la[2], sepEsrule0ma[0], sepEsrule0ma[1], sepEsrule0ma[2]);


        fprintf(EsepSRPNfile, "%20.16f %20.16f %20.16f %20.16f\n", Esrule0lP, -Esrule0mP, Esrule0lN, -Esrule0mN);
        fprintf(EsepSRPNxfile, "%20.16f %20.16f %20.16f %20.16f\n", Esrule0lPx, -Esrule0mPx, Esrule0lN, -Esrule0mN);
        fprintf(EsepSRPNafile, "%20.16f %20.16f %20.16f %20.16f\n", Esrule0lPa, -Esrule0mPa, Esrule0lN, -Esrule0mN);
        fprintf(EsepSRPNsepfile, "%20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n", sepEsrule0lP[0], sepEsrule0lP[1], sepEsrule0lP[2], sepEsrule0mP[0], sepEsrule0mP[1], sepEsrule0mP[2], sepEsrule0lN[0], sepEsrule0lN[1], sepEsrule0lN[2], sepEsrule0mN[0], sepEsrule0mN[1], sepEsrule0mN[2]);
        fprintf(EsepSRPNsepxfile, "%20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n", sepEsrule0lPx[0], sepEsrule0lPx[1], sepEsrule0lPx[2], sepEsrule0mPx[0], sepEsrule0mPx[1], sepEsrule0mPx[2], sepEsrule0lNx[0], sepEsrule0lNx[1], sepEsrule0lNx[2], sepEsrule0mNx[0], sepEsrule0mNx[1], sepEsrule0mNx[2]);
        fprintf(EsepSRPNsepafile, "%20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n", sepEsrule0lPa[0], sepEsrule0lPa[1], sepEsrule0lPa[2], sepEsrule0mPa[0], sepEsrule0mPa[1], sepEsrule0mPa[2], sepEsrule0lNa[0], sepEsrule0lNa[1], sepEsrule0lNa[2], sepEsrule0mNa[0], sepEsrule0mNa[1], sepEsrule0mNa[2]);
      }
    }
    if (G<GELcut) {
      fclose(EsepSRfile);
      fclose(EsepSRxfile);
      fclose(EsepSRafile);
      fclose(EsepSRsepfile);
      fclose(EsepSRsepafile);
      fclose(EsepSRsepxfile);
      fclose(EsepSRPNfile);
      fclose(EsepSRPNafile);
      fclose(EsepSRPNxfile);
      fclose(EsepSRPNsepfile);
      fclose(EsepSRPNsepafile);
      fclose(EsepSRPNsepxfile);
    }

    //*******************************

		//srule0 =0;
		for (i=0; i<tempP; i++) {
			srule0 += (2.0*G+1.0)*(fabs(Pval[i]) - sqrt(Pval[i]*Pval[i]+lambda*lambda) + 0.5*lambda*lambda/sqrt(Pval[i]*Pval[i]+lambda*lambda))/mass/mass;
			srule0 += -(2.0*G+1.0)*(fabs(PvalF0[i]) - sqrt(PvalF0[i]*PvalF0[i]+lambda*lambda) + 0.5*lambda*lambda/sqrt(PvalF0[i]*PvalF0[i]+lambda*lambda))/mass/mass;
      srule0x += (2.0*G+1.0)*(fabs(Pval[i]) - sqrt(Pval[i]*Pval[i]+lambda*lambda) + 0.5*lambda*lambda/sqrt(Pval[i]*Pval[i]+lambda*lambda))/mass/mass;
      srule0a += -(2.0*G+1.0)*(fabs(PvalF0[i]) - sqrt(PvalF0[i]*PvalF0[i]+lambda*lambda) + 0.5*lambda*lambda/sqrt(PvalF0[i]*PvalF0[i]+lambda*lambda))/mass/mass;
    }
		for (i=0; i<tempN; i++) {
			srule0 +=  (2.0*G+1.0)*(fabs(Nval[i]) - sqrt(Nval[i]*Nval[i]+lambda*lambda) + 0.5*lambda*lambda/sqrt(Nval[i]*Nval[i]+lambda*lambda))/mass/mass;
			srule0 += -(2.0*G+1.0)*(fabs(NvalF0[i]) - sqrt(NvalF0[i]*NvalF0[i]+lambda*lambda) + 0.5*lambda*lambda/sqrt(NvalF0[i]*NvalF0[i]+lambda*lambda))/mass/mass;
      srule0x += (2.0*G+1.0)*(fabs(Nval[i]) - sqrt(Nval[i]*Nval[i]+lambda*lambda) + 0.5*lambda*lambda/sqrt(Nval[i]*Nval[i]+lambda*lambda))/mass/mass;
      srule0a +=-(2.0*G+1.0)*(fabs(NvalF0[i]) - sqrt(NvalF0[i]*NvalF0[i]+lambda*lambda) + 0.5*lambda*lambda/sqrt(NvalF0[i]*NvalF0[i]+lambda*lambda))/mass/mass;
    }
		//srule0 *= -1/mass/mass;
		for (i=0; i<=xNInt; i++) {
			for (j=0;j<tempP+tempN;j++) {
				USFseaP[0][i] += USFseaPL[0][i][j];
				USFseaN[0][i] += USFseaNL[0][i][j];
				USFseaP[1][i] += USFseaPL[1][i][j];
				USFseaN[1][i] += USFseaNL[1][i][j];

				USFseaPF0[0][i] += USFseaPF0L[0][i][j];
				USFseaNF0[0][i] += USFseaNF0L[0][i][j];
				USFseaPF0[1][i] += USFseaPF0L[1][i][j];
				USFseaNF0[1][i] += USFseaNF0L[1][i][j];
			}
		}
		for (i=0; i<=xNInt; i++) {
			intarr1x[i]=0;
			intarr2x[i]=0;
			x=(double)i/(double)xNInt*xdist;
			for (j=0;j<tempP+tempN;j++) {
				intarr1x[i]+=x*(USFseaNL[0][i][j]-USFseaNF0L[0][i][j]+USFseaPL[0][i][j]-USFseaPF0L[0][i][j]);
        intarr2x[i]+=x*(USFseaNL[1][i][j]-USFseaNF0L[1][i][j]+USFseaPL[1][i][j]-USFseaPF0L[1][i][j]);
			}
		}
		srule0l += 0.5*booles(intarr1x, xNInt)*xdist/(double)xNInt;
		srule1l += 0.5*booles(intarr2x, xNInt)*xdist/(double)xNInt;

    for (i=0; i<=xNInt; i++) {
			intarr1x[i]=0;
			intarr2x[i]=0;
			x=(double)i/(double)xNInt*xdist;
			for (j=0;j<tempP+tempN;j++) {
				intarr1x[i]+=x*(USFseaNL[0][i][j]+USFseaPL[0][i][j]);
				intarr2x[i]+=x*(USFseaNL[1][i][j]+USFseaPL[1][i][j]);
			}
		}
		srule0lx += 0.5*booles(intarr1x, xNInt)*xdist/(double)xNInt;
		srule1lx += 0.5*booles(intarr2x, xNInt)*xdist/(double)xNInt;

    for (i=0; i<=xNInt; i++) {
			intarr1x[i]=0;
			intarr2x[i]=0;
			x=(double)i/(double)xNInt*xdist;
			for (j=0;j<tempP+tempN;j++) {
				intarr1x[i]+=x*(USFseaNF0L[0][i][j]+USFseaPF0L[0][i][j]);
				intarr2x[i]+=x*(USFseaNF0L[1][i][j]+USFseaPF0L[1][i][j]);
			}
		}
		srule0la += -0.5*booles(intarr1x, xNInt)*xdist/(double)xNInt;
		srule1la += -0.5*booles(intarr2x, xNInt)*xdist/(double)xNInt;

    fprintf(stderr, "---------------------------------\n");
		fprintf(stderr, "%d %20.10f %20.10f %20.10f %20.10f (G sr0 sr0l sr1 sr1l)\n", G, -srule0, srule0l, 0.5*srule1, srule1l);
    fprintf(stderr, "%d %20.10f %20.10f %20.10f %20.10f (G sr0x sr0lx sr1x sr1lx)\n", G, -srule0x, srule0lx, 0.5*srule1x, srule1lx);
    fprintf(stderr, "%d %20.10f %20.10f %20.10f %20.10f (G sr0a sr0la sr1a sr1la)\n\n", G, -srule0a, srule0la, 0.5*srule1a, srule1la);
    if (G<GELcut) {
      fprintf(stderr, "%12.7f %12.7f %12.7f %12.7f %12.7f (easr0 easr0l easr0m easr1 easr1l)\n", -EAsrule0, EAsrule0l, EAsrule0m, 0.5*EAsrule1, EAsrule1l);
      fprintf(stderr, "%12.7f %12.7f %12.7f %12.7f %12.7f (easr0x easr0lx easr0mx easr1x easr1lx)\n", -EAsrule0x, EAsrule0lx, EAsrule0mx, 0.5*EAsrule1x, EAsrule1lx);
      fprintf(stderr, "%12.7f %12.7f %12.7f %12.7f %12.7f (easr0a easr0la easr0ma easr1a easr1la)\n\n", -EAsrule0a, EAsrule0la, EAsrule0ma, 0.5*EAsrule1a, EAsrule1la);
      fprintf(stderr, "%12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f (ssr0 ssr1 ssr2 ssrx0 ssrx1 ssrx2 ssra0 ssra1 ssra2)\n", -sepsrule0[0], -sepsrule0[1], -sepsrule0[2], sepsrule0x[0], sepsrule0x[1], sepsrule0x[2], sepsrule0a[0], sepsrule0a[1], sepsrule0a[2]);
      fprintf(stderr, "%12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f (ssrl0 ssrl1 ssrl2 ssrlx0 ssrlx1 ssrlx2 ssrla0 ssrla1 ssrla2)\n", sepsrule0l[0], sepsrule0l[1], sepsrule0l[2], sepsrule0lx[0], sepsrule0lx[1], sepsrule0lx[2], sepsrule0la[0], sepsrule0la[1], sepsrule0la[2]);
      fprintf(stderr, "%12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f (ssrm0 ssrm1 ssrm2 ssrmx0 ssrmx1 ssrmx2 ssrma0 ssrma1 ssrma2)\n\n", sepsrule0m[0], sepsrule0m[1], sepsrule0m[2], sepsrule0mx[0], sepsrule0mx[1], sepsrule0mx[2], sepsrule0ma[0], sepsrule0ma[1], sepsrule0ma[2]);

      fprintf(stderr, "%12.7f %12.7f %12.7f %12.7f (easr0lN easr0mN easr0lP ears0mP)\n", EAsrule0lN, EAsrule0mN, EAsrule0lP, EAsrule0mP);
      fprintf(stderr, "%12.7f %12.7f %12.7f %12.7f (easr0lNx easr0mNx easr0lPx easr0mPx)\n", EAsrule0lNx, EAsrule0mNx, EAsrule0lPx, EAsrule0mPx);
      fprintf(stderr, "%12.7f %12.7f %12.7f %12.7f (easr0lNa easr0mNa easr0lPa easr0mPa)\n\n", EAsrule0lNa, EAsrule0mNa, EAsrule0lPa, EAsrule0mPa);

      fprintf(stderr, "%12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f (ssrlP0 ssrlP1 ssrlP2 ssrlPx0 ssrlPx1 ssrlPx2 ssrlPa0 ssrlPa1 ssrlPa2)\n", sepsrule0lP[0], sepsrule0lP[1], sepsrule0lP[2], sepsrule0lPx[0], sepsrule0lPx[1], sepsrule0lPx[2], sepsrule0lPa[0], sepsrule0lPa[1], sepsrule0lPa[2]);
      fprintf(stderr, "%12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f (ssrmP0 ssrmP1 ssrmP2 ssrmPx0 ssrmPx1 ssrmPx2 ssrmPa0 ssrmPa1 ssrmPa2)\n\n", sepsrule0mP[0], sepsrule0mP[1], sepsrule0mP[2], sepsrule0mPx[0], sepsrule0mPx[1], sepsrule0mPx[2], sepsrule0mPa[0], sepsrule0mPa[1], sepsrule0mPa[2]);

      fprintf(stderr, "%12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f (ssrlN0 ssrlN1 ssrlN2 ssrlNx0 ssrlNx1 ssrlNx2 ssrlNa0 ssrlNa1 ssrlNa2)\n", sepsrule0lN[0], sepsrule0lN[1], sepsrule0lN[2], sepsrule0lNx[0], sepsrule0lNx[1], sepsrule0lNx[2], sepsrule0lNa[0], sepsrule0lNa[1], sepsrule0lNa[2]);
      fprintf(stderr, "%12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f (ssrmN0 ssrmN1 ssrmN2 ssrmNx0 ssrmNx1 ssrmNx2 ssrmNa0 ssrmNa1 ssrmNa2)\n\n", sepsrule0mN[0], sepsrule0mN[1], sepsrule0mN[2], sepsrule0mNx[0], sepsrule0mNx[1], sepsrule0mNx[2], sepsrule0mNa[0], sepsrule0mNa[1], sepsrule0mNa[2]);
      fprintf(stderr, "---------------------------------\n");
    } else {
      fprintf(stderr, "---------------------------------\n");
    }
    for (l=0; l<tempP+tempN; l++) {
			for (i=0; i<=xNInt; i++) {
				fwrite(&USFseaPL[0][i][l], sizeof(double), 1, GsepUSFfile);
			}
			for (i=0; i<=xNInt; i++) {
				fwrite(&USFseaNL[0][i][l], sizeof(double), 1, GsepUSFfile);
			}
			for (i=0; i<=xNInt; i++) {
				fwrite(&USFseaPF0L[0][i][l], sizeof(double), 1, GsepUSFfile);
			}
			for (i=0; i<=xNInt; i++) {
				fwrite(&USFseaNF0L[0][i][l], sizeof(double), 1, GsepUSFfile);
			}
			for (i=0; i<=xNInt; i++) {
				fwrite(&USFseaPL[1][i][l], sizeof(double), 1, GsepUSFfile);
			}
			for (i=0; i<=xNInt; i++) {
				fwrite(&USFseaNL[1][i][l], sizeof(double), 1, GsepUSFfile);
			}
			for (i=0; i<=xNInt; i++) {
				fwrite(&USFseaPF0L[1][i][l], sizeof(double), 1, GsepUSFfile);
			}
			for (i=0; i<=xNInt; i++) {
				fwrite(&USFseaNF0L[1][i][l], sizeof(double), 1, GsepUSFfile);
			}
		}


    //*************************************
    //*************************************
    //FREEING MEMORY
		fclose(GsepUSFfile);
    freeDataCollectedvals(Pval, Nval, PvalF0, NvalF0, mks, ks, pks, G, momentanums);
    freeDataCollectedfouriers(fourier_funcP, fourier_funcN, fourier_funcPF0, fourier_funcNF0, G, momentanums);
    for (i=0;i<=xNInt;i++) {
			free(USFseaPL[0][i]);
			free(USFseaNL[0][i]);
			free(USFseaPL[1][i]);
			free(USFseaNL[1][i]);

			free(USFseaPF0L[0][i]);
			free(USFseaNF0L[0][i]);
			free(USFseaPF0L[1][i]);
			free(USFseaNF0L[1][i]);
		}
    if (G<GELcut) {
      for (i=0; i<3; i++) {
        free(USFseaPEL[i]);
        free(USFseaNEL[i]);
        free(USFseaNF0EL[i]);
        free(USFseaPF0EL[i]);

        free(USFseaNELsep[i]);
        free(USFseaPELsep[i]);
        free(USFseaNF0ELsep[i]);
        free(USFseaPF0ELsep[i]);
      }
      free(USFseaPEL);
      free(USFseaNEL);
      free(USFseaNF0EL);
      free(USFseaPF0EL);
    }

		//*************************************
  }
  for (l=0; l<=xNInt; l++) {
    //USFseaP[l]*=5*mass*Nc/144.0;
    //USFseaN[l]*=5*mass*Nc/144.0;
    //USFseaPF0[l]*=5*mass*Nc/144.0;
    //USFseaNF0[l]*=5*mass*Nc/144.0;
  }
	for (i=0; i<=xNInt; i++) {
		fwrite(&USFseaP[0][i], sizeof(double), 1, USFfile);
	}
	for (i=0; i<=xNInt; i++) {
		fwrite(&USFseaN[0][i], sizeof(double), 1, USFfile);
	}
	for (i=0; i<=xNInt; i++) {
		fwrite(&USFseaPF0[0][i], sizeof(double), 1, USFfile);
	}
	for (i=0; i<=xNInt; i++) {
		fwrite(&USFseaNF0[0][i], sizeof(double), 1, USFfile);
	}

	for (i=0; i<=xNInt; i++) {
		fwrite(&USFseaP[1][i], sizeof(double), 1, USFfile);
	}
	for (i=0; i<=xNInt; i++) {
		fwrite(&USFseaN[1][i], sizeof(double), 1, USFfile);
	}
	for (i=0; i<=xNInt; i++) {
		fwrite(&USFseaPF0[1][i], sizeof(double), 1, USFfile);
	}
	for (i=0; i<=xNInt; i++) {
		fwrite(&USFseaNF0[1][i], sizeof(double), 1, USFfile);
	}
	fclose(USFfile);
  free(intarr1);
  free(intarr2);
  free(intarrtemp3);
  free(intarrtemp1);
  free(intarrtemp2);
  free(intarrtemp4);
}



int main(int argc, const char* argv[]) {
  unsigned short grandspintot, grandspinnum, *momentanums;
  unsigned long i, j, k, l;
  double x, Esol, Ev, mass, temp, c1, c2;
  double *mks, *ks, *pks, *Pval, *Nval, *PvalF0, *NvalF0, ***fourier_funcP, ***fourier_funcN, **USFvalN, **USFseaN, **USFvalP, **USFseaP, **USFseaPF0, **USFseaNF0, *USF_RF, *USF_IMF;
  double *thetas, ***fourier_funcPF0, ***fourier_funcNF0, *USFseaPIMF, *USFvalPIMF, *USFseaNIMF, *USFvalNIMF, *arrtemp1, *arrtemp2;
  FILE *ks_check, *Evals_check, *theta_check, *fourier_check, *fourier_checktxt, *USFval_file, *USFsea_file, *USFseaF0_file, *USFvalIMF_file, *USFseaIMF_file;
  //fprintf(stderr, "main 1\n");
fprintf(stderr,"start Gspin_momentanum\n");
  Gspin_momentanums(&grandspintot, &momentanums);
fprintf(stderr,"start valence\n");
  //fprintf(stdout, "grandspintot=%d\n", grandspintot);
  //return 0;

  thetas = (double*)malloc(sizeof(double)*(NInt+1));
  USFvalN = (double**)malloc(sizeof(double*)*2);
	USFvalP = (double**)malloc(sizeof(double*)*2);
	USFseaN = (double**)malloc(sizeof(double*)*2);
	USFseaP = (double**)malloc(sizeof(double*)*2);
	USFseaNF0 = (double**)malloc(sizeof(double*)*2);
	USFseaPF0 = (double**)malloc(sizeof(double*)*2);

	USFvalP[0] = (double*)malloc(sizeof(double)*(xNInt+1));
	USFvalP[1] = (double*)malloc(sizeof(double)*(xNInt+1));
	USFvalN[0] = (double*)malloc(sizeof(double)*(xNInt+1));
	USFvalN[1] = (double*)malloc(sizeof(double)*(xNInt+1));

	USFseaN[0] = (double*)malloc(sizeof(double)*(xNInt+1));
  USFseaP[0] = (double*)malloc(sizeof(double)*(xNInt+1));
  USFseaNF0[0] = (double*)malloc(sizeof(double)*(xNInt+1));
  USFseaPF0[0] = (double*)malloc(sizeof(double)*(xNInt+1));

	USFseaN[1] = (double*)malloc(sizeof(double)*(xNInt+1));
  USFseaP[1] = (double*)malloc(sizeof(double)*(xNInt+1));
  USFseaNF0[1] = (double*)malloc(sizeof(double)*(xNInt+1));
  USFseaPF0[1] = (double*)malloc(sizeof(double)*(xNInt+1));

	arrtemp1 = (double*)malloc(sizeof(double)*(xNInt+1));
  arrtemp2 = (double*)malloc(sizeof(double)*(xNInt+1));
  chiralangle_from_file(&thetas);
  mpion = mpion/cqm;
  fprintf(stderr, "kdist=%f\ndist=%f\nkratio=%f\nlambda=%f\ncqm=%f\nrep=%d\nexpnum=%d\nxdist=%f\nxNInt=%f\n", kdist, dist, kratio, lambda, cqm, rep, expnum, xdist, xNInt);
	compute_USF_valence(USFvalP, USFvalN, momentanums, 940.0/cqm);
  compute_USF_sea_Esep(USFseaP, USFseaN, USFseaPF0, USFseaNF0, grandspintot, 940.0/cqm, momentanums, 0, 1000, 1000);
  //Esol=compute_soliton_energy(thetas, grandspintot, momentanums);
  fprintf(stderr, "%f\n", 400.0*Esol);
  free(arrtemp1);
	free(arrtemp2);
}
