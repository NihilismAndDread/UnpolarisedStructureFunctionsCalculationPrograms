#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "constants.h"
#include "routines.h"

void chiralangle_from_file(double **pthetas);
void freeDataCollectedvals(double *Pval, double *Nval, double *PvalF0, double *NvalF0, double *mks, double *ks, double *pks, unsigned short G, unsigned short *momentanums);
void freeDataCollectedfouriers(double ***fourier_funcP, double ***fourier_funcN, double ***fourier_funcPF0, double ***fourier_funcNF0, unsigned short G, unsigned short *momentanums);
void dataforGfouriers(unsigned short grandspinnum, double ****pfourier_funcP, double **** pfourier_funcN, double **** pfourier_funcPF0, double **** pfourier_funcNF0);
void dataforGvals(unsigned short grandspinnum, double ** Pval, double ** Nval, double **mks, double ** ks, double **pks);
void chiralangle_from_file(double **pthetas);
void Gspin_momentanums(unsigned short *grandspinnum, unsigned short **pmomenta);
void USFseaAMdata(double ***pUSFsea);
void sepUSFseaAMdata(unsigned short grandspinnum, double ****pUSFseaL);
void sepUSFseadata(unsigned short grandspinnum, double ****pUSFseaL, double ****pUSFseaF0L);
void freeUSFLs(unsigned short grandspinnum, unsigned short *momentanum, double ***USFseaL, double ***USFseaF0L);
void SumRuleFitter(double *USFseaF0alt, double *USFseaF0);
void flatten_function(double **pansfunc, double *func, unsigned short len, double dist);
//void alt_function(double **pansfunc, double *func1, double *func2, unsigned short len, double dist, unsigned short altnum);
void USFvaldata(double ***pUSFval);
void freeUSFvaldata(double **USFval);
void createF0vals(double **pPvalF0, double **pNvalF0, unsigned short G, double *mks, double *ks, double *pks, unsigned short *momentanums);
void altMassFormat();
void altAverageFormat(char printIMF);
double compute_sea_energy(double * thetas, unsigned short grandspintot, unsigned short *momentanums);
void chiral_to_file(unsigned short texpnum, unsigned short labelmax);
void chiralangle_from_labeledfile(double **pthetas, unsigned short texpnum, unsigned short labelnum);
double RF_to_IMF(double **pF_IMF, double *F_RF);
//void altFormat(unsigned short altnum);

double RF_to_IMF(double **pF_IMF, double *F_RF) {
  double x_IMF, x_RF, x1, x2, *F_IMF;
  unsigned short index, mark;
  unsigned long i, j;
  index =0;
  mark=0;
  F_IMF = (double*)malloc(sizeof(double)*(xNInt+1));
  for (i=0; i<xNInt; i++) {
    x_IMF = (double)i/(double)xNInt;
    x_RF=-log(1.0-x_IMF);
    for (j=0; j<xNInt; j++) {
      x1 = (double)j/xNInt*xdist;
      x2 = (double)(j+1)/xNInt*xdist;
      if (x1<x_RF && x_RF<=x2) {
        F_IMF[i]= step_function(1-x_IMF)/(1-x_IMF) *F_RF[j+1];
        break;
      } else if (j+1==xNInt) {
        if (mark==0) {
          mark = i;
        }
        F_IMF[i]= step_function(1-x_IMF)/(1-x_IMF) *F_RF[j+1];
        break;
      }
    }

    //fprintf(stderr, "x_RF=%f\n", x_RF);
  }
  *pF_IMF=F_IMF;
  return mark;
}
void chiral_to_file(unsigned short texpnum, unsigned short labelmax) {
  double *chiral;
  FILE *chiral_file;
  char sfile[100];
  unsigned short i, j;
  for (j=0; j<labelmax; j++) {
    chiralangle_from_labeledfile(&chiral, texpnum, j);
    sprintf(sfile, "formatted_Data/expnum%d/chiral_%d.txt", texpnum, j);
    chiral_file = fopen(sfile, "w+");
    fprintf(stderr, "chiral angle for %d: %f", j, chiral[0]);
    fprintf(chiral_file, "%20.16f", chiral[0]);
    for (i=1; i<= NInt; i++) {
      fprintf(stderr, ", %f", chiral[i]);
      fprintf(chiral_file, "\n%20.16f", chiral[i]);
    }
    fprintf(stderr, "\n");
    fclose(chiral_file);
  }
}
double compute_sea_energy(double * thetas, unsigned short grandspintot, unsigned short *momentanums) {
  double t1, t2, t3, Ev, t21, t22, *intarr, x, k, Nval, Pval, En, Esol, smpion;
  double *Pvals, *Nvals, *mks, *ks, *pks;
  unsigned short i, j;
  smpion=mpion/cqm;
  Esol=0;
  for (i=0; i<=grandspintot; i++) {
    t2=0;
    dataforGvals(i, &Pvals, &Nvals, &mks, &ks, &pks);
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
  fprintf(stderr, "Completed eigen values contributions\n");
  intarr=(double*)malloc(sizeof(double)*(NInt+1));
  //Esol *= Nc;
  //fprintf(stderr, "E0=%f\n", cqm*Nc*Esol/940.0);
  for (i=0; i<=NInt; i++) {
    x= (double)i/(double)NInt;
    intarr[i]=x*x*(1.0-cos(thetas[i]));
  }
  t3 = 4*M_PI*smpion*smpion*93.0*93.0/(cqm*cqm)*dist*dist*dist;
  t2 = simps(intarr, NInt)/(double)NInt;
  //fprintf(stderr, "E_mesonic=%f\n", cqm*t3*t2);
  //Esol += t3;
  free(intarr);
  return Esol;
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
void dataforGfouriers(unsigned short grandspinnum, double ****pfourier_funcP, double **** pfourier_funcN, double **** pfourier_funcPF0, double **** pfourier_funcNF0) {
  unsigned short *momentanum, i, j, k, lenP, lenN, plenk, lenk, mlenk;
  unsigned int seekto1, seekto2;
  double *** fourier_funcP, ***fourier_funcN, temp, ***fourier_funcPF0, ***fourier_funcNF0;
	char filename[100];
  FILE * momentafile, *wavemodesfile, *wavemodesfileF0;

  momentanum = (unsigned short *)malloc(sizeof(unsigned short)*(grandspinnum+2));
	sprintf(filename, "Data/expnum%d/momentanum_%d.b", expnum, expnum);
  momentafile = fopen(filename, "rb");

  sprintf(filename, "Data/expnum%d/fourierfileF0G%d_%d.b", expnum, grandspinnum, expnum);
	wavemodesfileF0=fopen(filename, "rb");
	sprintf(filename, "Data/expnum%d/fourierfileG%d_%d.b", expnum, grandspinnum, expnum);
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
void dataforGvals(unsigned short grandspinnum, double ** Pval, double ** Nval, double **mks, double ** ks, double **pks) {
  unsigned short * momentanum, i, j, k, lenP, lenN, plenk, lenk, mlenk, index2;
  unsigned int seekto1, seekto2;
  double * Pvals, * Nvals, *kn, *pkn, *mkn, temp;
	char filename[100];
  FILE * momentafile, * eigvalfile, * eigvecfile, *knfile;

  momentanum = (unsigned short *)malloc(sizeof(unsigned short)*(grandspinnum+2));
  //fprintf(stderr, "vibe check Gvals %d\n", index2++);
	sprintf(filename, "Data/expnum%d/eigenergies_%d.b", expnum, expnum);
	eigvalfile = fopen(filename, "rb");
	sprintf(filename, "Data/expnum%d/kn_%d.b", expnum, expnum);
  knfile = fopen(filename, "rb");
  sprintf(filename,"Data/expnum%d/momentanum_%d.b", expnum, expnum);
  momentafile = fopen(filename, "rb");
  //fprintf(stderr, "vibe check Gvals %d\n", index2++);
  //fprintf(stderr, "dataforG 1\n");
  for (i=0; i<=grandspinnum+1; i++) {
    fread(&momentanum[i], sizeof(unsigned short), 1, momentafile);
  }
  fclose(momentafile);
  //fprintf(stderr, "vibe check Gvals %d\n", index2++);
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
  //fprintf(stderr, "vibe check Gvals %d\n", index2++);
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
  //fprintf(stderr, "vibe check Gvals %d\n", index2++);
  if (grandspinnum==0) {
    seekto1=0;
  } else {
    seekto1 =momentanum[0]*2+momentanum[1]*2;
  }
  for (i=1; i<=grandspinnum-1; i++) {
    seekto1 += momentanum[i]*4 + momentanum[i-1]*2 + momentanum[i+1]*2;
  }
  //fprintf(stderr, "vibe check Gvals %d %d\n", index2++, seekto1);
  //fprintf(stderr, "dataforG 10\n");
  fseek(eigvalfile, sizeof(double)*seekto1, SEEK_SET);
  for (i=0; i<lenP; i++) {
    fread(&Pvals[i], sizeof(double), 1, eigvalfile);
  }
  for (i=0; i<lenN; i++) {
    fread(&Nvals[i], sizeof(double), 1, eigvalfile);
  }
  //fprintf(stderr, "vibe check Gvals %d\n", index2++);
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
  //fprintf(stderr, "file opened\n");
  fseek(thetas_file, 0, SEEK_SET);
  for (i=0; i<=NInt; i++) {
    fread(&thetas[i], sizeof(double), 1, thetas_file);
  }
  fclose(thetas_file);
  *pthetas = thetas;
}
void chiralangle_from_labeledfile(double **pthetas, unsigned short texpnum, unsigned short labelnum) {
  FILE *thetas_file;
  double *thetas;
  unsigned short i;
	char filename[100];
  thetas = (double*)malloc((NInt+1)*sizeof(double));

	sprintf(filename, "Data/expnum%d/chiralangle_%d.b", texpnum, labelnum);
  thetas_file = fopen(filename, "rb");
  //fprintf(stderr, "file opened\n");
  fseek(thetas_file, 0, SEEK_SET);
  for (i=0; i<=NInt; i++) {
    fread(&thetas[i], sizeof(double), 1, thetas_file);
  }
  fclose(thetas_file);
  *pthetas = thetas;
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
  momenta[index]=0;
  *pmomenta=momenta;
  *grandspinnum = index-2;
  fclose(momentanumsfile);
}
void USFseaAMdata(double ***pUSFsea) {
  unsigned short i, j, k;
  double **USFsea, **USFseaF0;
  FILE *USFseabfile;
  char filename[100];

  sprintf(filename, "Data/expnum%d/USFfile_%d_%d.b", expnum, expnum, massnum);
  USFseabfile = fopen(filename, "rb");

  USFsea = (double**)malloc(sizeof(double*)*4);
  for (i=0; i<4; i++) {
    USFsea[i] = (double*)malloc(sizeof(double)*(kNInt+1));
  }
  for (i=0; i<=xNInt; i++) {
    fread(&USFsea[0][i], sizeof(double), 1, USFseabfile);
  }
  for (i=0; i<=xNInt; i++) {
    fread(&USFsea[2][i], sizeof(double), 1, USFseabfile);
  }

  for (i=0; i<=xNInt; i++) {
    fread(&USFsea[1][i], sizeof(double), 1, USFseabfile);
  }
  for (i=0; i<=xNInt; i++) {
    fread(&USFsea[3][i], sizeof(double), 1, USFseabfile);
  }
  fclose(USFseabfile);
  *pUSFsea=USFsea;
}
void USFseadata(double ***pUSFsea, double ***pUSFseaF0) {
  unsigned short i, j, k;
  double **USFsea, **USFseaF0;
  FILE *USFseabfile;
  char filename[100];

  sprintf(filename, "Data/expnum%d/USFdata/USFfile_%d.b", expnum, expnum);
  USFseabfile = fopen(filename, "rb");

  USFsea = (double**)malloc(sizeof(double*)*4);
  USFseaF0 = (double**)malloc(sizeof(double*)*4);
  for (i=0; i<4; i++) {
    USFsea[i] = (double*)malloc(sizeof(double)*(kNInt+1));
    USFseaF0[i] = (double*)malloc(sizeof(double)*(kNInt+1));
  }
  for (i=0; i<=xNInt; i++) {
    fread(&USFsea[0][i], sizeof(double), 1, USFseabfile);
  }
  for (i=0; i<=xNInt; i++) {
    fread(&USFsea[2][i], sizeof(double), 1, USFseabfile);
  }
  for (i=0; i<=xNInt; i++) {
    fread(&USFseaF0[0][i], sizeof(double), 1, USFseabfile);
  }
  for (i=0; i<=xNInt; i++) {
    fread(&USFseaF0[2][i], sizeof(double), 1, USFseabfile);
  }

  for (i=0; i<=xNInt; i++) {
    fread(&USFsea[1][i], sizeof(double), 1, USFseabfile);
  }
  for (i=0; i<=xNInt; i++) {
    fread(&USFsea[3][i], sizeof(double), 1, USFseabfile);
  }
  for (i=0; i<=xNInt; i++) {
    fread(&USFseaF0[1][i], sizeof(double), 1, USFseabfile);
  }
  for (i=0; i<=xNInt; i++) {
    fread(&USFseaF0[3][i], sizeof(double), 1, USFseabfile);
  }
  fclose(USFseabfile);
  *pUSFsea=USFsea;
  *pUSFseaF0=USFseaF0;
}
void sepUSFseaAMdata(unsigned short grandspinnum, double ****pUSFseaL) {
  unsigned short *momentanum, i, j, k, lenP, lenN, plenk, lenk, mlenk;
  unsigned int seekto1, seekto2;
  double temp, ***USFseaL, ***USFseaF0L;
  char filename[100];
  FILE * momentafile, * GsepUSFfile;

  momentanum = (unsigned short *)malloc(sizeof(unsigned short)*(grandspinnum+2));
  sprintf(filename, "Data/expnum%d/momentanum_%d.b", expnum, expnum);
  momentafile = fopen(filename, "rb");
  for (i=0; i<=grandspinnum+1; i++) {
    fread(&momentanum[i], sizeof(unsigned short), 1, momentafile);
  }
  fclose(momentafile);

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

  USFseaL = (double***)malloc(sizeof(double**)*4);
  for (i=0; i<4; i++) {
    USFseaL[i]=(double**)malloc(sizeof(double*)*(lenP+lenN));
    for (j=0; j<lenN+lenP; j++) {
      USFseaL[i][j]=(double*)malloc(sizeof(double)*(xNInt+1));
    }
  }
  sprintf(filename, "Data/expnum%d/GsepUSFfileG%d_%d_%d.b", expnum, grandspinnum, expnum, massnum);
  GsepUSFfile = fopen(filename, "rb");
  for (i=0; i<lenP+lenN; i++) {
    for (j=0; j<=xNInt; j++) {
      fread(&USFseaL[0][i][j], sizeof(double), 1, GsepUSFfile);
    }
    for (j=0; j<=xNInt; j++) {
      fread(&USFseaL[2][i][j], sizeof(double), 1, GsepUSFfile);
    }
    for (j=0; j<=xNInt; j++) {
      fread(&USFseaL[1][i][j], sizeof(double), 1, GsepUSFfile);
    }
    for (j=0; j<=xNInt; j++) {
      fread(&USFseaL[3][i][j], sizeof(double), 1, GsepUSFfile);
    }
  }
  fclose(GsepUSFfile);
  free(momentanum);
  *pUSFseaL=USFseaL;
}
void sepUSFseadata(unsigned short grandspinnum, double ****pUSFseaL, double ****pUSFseaF0L) {
  unsigned short * momentanum, i, j, k, lenP, lenN, plenk, lenk, mlenk;
  unsigned int seekto1, seekto2;
  double temp, ***USFseaL, ***USFseaF0L;
  char filename[100];
  FILE * momentafile, * GsepUSFfile;

  momentanum = (unsigned short *)malloc(sizeof(unsigned short)*(grandspinnum+2));
  sprintf(filename, "Data/expnum%d/momentanum_%d.b", expnum, expnum);
  momentafile = fopen(filename, "rb");
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
  USFseaL = (double***)malloc(sizeof(double**)*4);
  USFseaF0L = (double***)malloc(sizeof(double**)*4);
  for (i=0; i<4; i++) {
    USFseaL[i]=(double**)malloc(sizeof(double*)*(lenP+lenN));
    USFseaF0L[i]=(double**)malloc(sizeof(double*)*(lenP+lenN));
    for (j=0; j<lenN+lenP; j++) {
      USFseaL[i][j]=(double*)malloc(sizeof(double)*(xNInt+1));
      USFseaF0L[i][j]=(double*)malloc(sizeof(double)*(xNInt+1));
    }
  }
  sprintf(filename, "Data/expnum%d/USFdata/GsepUSFfileG%d_%d.b", expnum, grandspinnum, expnum);
  GsepUSFfile = fopen(filename, "rb");
  for (i=0; i<lenP+lenN; i++) {
    for (j=0; j<=xNInt; j++) {
      fread(&USFseaL[0][i][j], sizeof(double), 1, GsepUSFfile);
    }
    for (j=0; j<=xNInt; j++) {
      fread(&USFseaL[2][i][j], sizeof(double), 1, GsepUSFfile);
    }
    for (j=0; j<=xNInt; j++) {
      fread(&USFseaF0L[0][i][j], sizeof(double), 1, GsepUSFfile);
    }
    for (j=0; j<=xNInt; j++) {
      fread(&USFseaF0L[2][i][j], sizeof(double), 1, GsepUSFfile);
    }
    for (j=0; j<=xNInt; j++) {
      fread(&USFseaL[1][i][j], sizeof(double), 1, GsepUSFfile);
    }
    for (j=0; j<=xNInt; j++) {
      fread(&USFseaL[3][i][j], sizeof(double), 1, GsepUSFfile);
    }
    for (j=0; j<=xNInt; j++) {
      fread(&USFseaF0L[1][i][j], sizeof(double), 1, GsepUSFfile);
    }
    for (j=0; j<=xNInt; j++) {
      fread(&USFseaF0L[3][i][j], sizeof(double), 1, GsepUSFfile);
    }
  }
  fclose(GsepUSFfile);
  free(momentanum);
  *pUSFseaL=USFseaL;
  *pUSFseaF0L=USFseaF0L;
}
void freeUSFLs(unsigned short grandspinnum, unsigned short *momentanum, double ***USFseaL, double ***USFseaF0L) {
  unsigned short i, j, k, l;
  unsigned short lenk, plenk, mlenk, lenP, lenN;
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

  for (i=0; i<4; i++) {
    for (j=0; j<lenP+lenN; j++) {
      free(USFseaL[i][j]);
      free(USFseaF0L[i][j]);
    }
    free(USFseaL[i]);
    free(USFseaF0L[i]);
  }
  free(USFseaL);
  free(USFseaF0L);
}
void SumRuleFitter(double *USFseaF0alt, double *USFseaF0) {
  unsigned short i;
  double *intarr1, *intarr2, x , c1, c2;

  intarr1 = (double*)malloc(sizeof(double)*(xNInt+1));
  intarr2 = (double*)malloc(sizeof(double)*(xNInt+1));

  for (i=0; i<=xNInt; i++) {
    x=(double)i/(double)xNInt *xdist;
    intarr1[i]=x*USFseaF0[i];
    intarr2[i]=x*USFseaF0alt[i];
  }
  c1 = booles(intarr1, xNInt)*xdist/(double)xNInt;
  c2 = booles(intarr2, xNInt)*xdist/(double)xNInt;
  for (i=0; i<=xNInt; i++) {
    USFseaF0alt[i]=c1/c2*USFseaF0alt[i];
  }

}
void flatten_function(double **pansfunc, double *func, unsigned short len, double dist) {
  unsigned short i;
  double average, *ansfunc, *intarr, x, c;
  ansfunc = (double*)malloc(sizeof(double)*(len+1));
  intarr = (double*)malloc(sizeof(double)*(len+1));
  
  for (i=0; i<=len; i++) {
    x=(double)i/(double)len*dist;
    intarr[i]=x*func[i];
  }
  c = 2.0*booles(intarr, len)/(double)len/dist;
  for (i=0; i<=len; i++) {
    ansfunc[i]=c;
  }
  *pansfunc=ansfunc;
}
/*void alt_function(double **pansfunc, double *func1, double *func2, unsigned short len, double dist, unsigned short altnum) {
  double c, *ansfunc, *intarr, s, x;
  unsigned short i;
  ansfunc = (double*)malloc(sizeof(double)*(len+1));
  intarr = (double*)malloc(sizeof(double)*(len+1));
  for (i=0; i<=len; i++) {
    x=(double)i/(double)len*dist;
    intarr[i]=x*func2[i];
  }
  c = booles(intarr, len)/(double)len;
  switch(altnum) {
    case 1:

    default:
      fprintf(stderr, "altnum %d not valid, printing back default F", altnum);
      *pansfunc = func2;
  }

}*/
void USFvaldata(double ***pUSFval) {
  unsigned short i, j ,k;
  char filename[100];
  FILE *USFvalfile;
  double **USFval;

  USFval = (double**)malloc(sizeof(double*)*4);
  for (i=0; i<4; i++) {
    USFval[i]=(double*)malloc(sizeof(double)*(xNInt+1));
  }
  sprintf(filename, "Data/expnum%d/USFdata/USFval_%d.b", expnum, expnum);
  USFvalfile = fopen(filename, "rb");
  for (i=0; i<=xNInt; i++) {
    fread(&USFval[0][i], sizeof(double), 1, USFvalfile);
  }
  for (i=0; i<=xNInt; i++) {
    fread(&USFval[1][i], sizeof(double), 1, USFvalfile);
  }
  for (i=0; i<=xNInt; i++) {
    fread(&USFval[2][i], sizeof(double), 1, USFvalfile);
  }
  for (i=0; i<=xNInt; i++) {
    fread(&USFval[3][i], sizeof(double), 1, USFvalfile);
    //fprintf(stderr, "%f\n", USFval[3][i]);
  }
  //fprintf(stderr, "cobwebs in my hair\n");
  fclose(USFvalfile);
  *pUSFval=USFval;
}
void freeUSFvaldata(double **USFval) {
  unsigned short i, j;
  for (i=0; i<4; i++) {
    free(USFval[i]);
  }
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
void EigVals() {
  double *Pvals, *Nvals, *ks, *pks, *mks;
  unsigned short i, j, k, G, grandspintot, *momentanums;
  char filename[50];
  FILE *eigvalfile;
  fprintf(stderr, "Starting EigVal printing\n");
  Gspin_momentanums(&grandspintot, &momentanums);
  for (G=0; G<=grandspintot; G++) {
    //fprintf(stderr, "Beginning for G=%d\n", G);
    sprintf(filename, "formatted_Data/expnum%d/EigvalsPG%d_%d.txt", expnum, G, expnum);
    eigvalfile = fopen(filename, "w+");
    dataforGvals(G, &Pvals, &Nvals, &mks, &ks, &pks);
    if (G!=0) {
      if (momentanums[G]!=0) {
        fprintf(eigvalfile, "%20.16f", Pvals[0]);
      }
      for (i=1; i<4*momentanums[G]; i++) {
        fprintf(eigvalfile, "\n%20.16f", Pvals[i]);
      }
    } else {
      fprintf(eigvalfile, "%20.16f", Pvals[0]);
      for (i=1; i<2*momentanums[G]; i++) {
        fprintf(eigvalfile, "\n%20.16f", Pvals[i]);
        if (fabs(Pvals[i])<1) {
          fprintf(stderr, "valence number i=%d, %f\n", i, Pvals[i]);
        }
      }
    }
    if (G<6) {
      j=4*momentanums[G];
      if (G==0) {
        j=2*momentanums[G];
      }
      fprintf(stderr, "%d\n", j);
    }
    fclose(eigvalfile);
    sprintf(filename, "formatted_Data/expnum%d/EigvalsNG%d_%d.txt", expnum, G, expnum);
    eigvalfile = fopen(filename, "w+");
    if (G!=0){
      if (2*momentanums[G+1]+2*momentanums[G-1]!=0) {
        fprintf(eigvalfile, "%20.16f", Nvals[0]);
      }
      for (i=1; i<2*momentanums[G+1]+2*momentanums[G-1]; i++) {
        fprintf(eigvalfile, "\n%20.16f", Nvals[i]);
      }
    } else {
      fprintf(eigvalfile, "%20.16f", Nvals[0]);
      for (i=1; i<2*momentanums[G+1]; i++) {
        fprintf(eigvalfile, "\n%20.16f", Nvals[i]);
      }
    }
    fclose(eigvalfile);
    //fprintf(stderr, "Completed G=%d\n", G);
  }
  fprintf(stderr, "Done!\n");
}
void printmomentsnums() {
  FILE *momentanumfile;
  unsigned short *momentanums, grandspintot, grandspinnum;
  char sfile[100];
  Gspin_momentanums(&grandspintot, &momentanums);
  sprintf(sfile, "formatted_Data/expnum%d/momentanum_%d.txt", expnum, expnum);
  momentanumfile = fopen(sfile, "w+");
  for (grandspinnum=0; grandspinnum<=grandspintot; grandspinnum++) {
    fprintf(momentanumfile, "%d\n", momentanums[grandspinnum]);
  }
  fclose(momentanumfile);
}
void altMassFormat() {
  unsigned short grandspinmax, *momentanums, i, j, k, l, G, lenP, lenN, lenk, plenk, mlenk, index, pnum, index2, xNIntp;
  double ***USFseaL, ***USFseaF0L, ***USFseaF0altL, **USFsea, **USFseaF0, **USFseaF0alt, **USFval, *thetas;
  double ***fourier_funcP, ***fourier_funcN, ***fourier_funcPF0, ***fourier_funcNF0, **val_funcs;
  double *Pvals, *Nvals, *mks, *ks, *pks, *PvalsF0, *NvalsF0;
  double *intarr1, *intarr2, *intarr3, *intarr4; //ARRAYS FOR INTEGRALS
  double val_e, x, p, mass, t1, t2, t3, t4, temp, srule0, srule1, srule0l, srule1l, srule0lalt, srule1lalt;
  FILE *USFseafile, *USFseaF0file, *USFseaF0altfile, *USFseaLfile, *USFseaF0Lfile, *USFseaF0altLfile;
  char filename[100];
  index2=0;
  mass=940.0/cqm;
  xNIntp=floor(graphcut/xdist*xNInt);
  //fprintf(stderr, "vibe check %d\n", index2++);
  Gspin_momentanums(&grandspinmax, &momentanums);
  //fprintf(stderr, "vibe check %d\n", index2++);
  USFseadata(&USFsea, &USFseaF0);
  USFseaF0alt = (double**)malloc(sizeof(double*)*4);
  for (i=0; i<4; i++) {
    USFseaF0alt[i] = (double*)malloc(sizeof(double)*(xNInt+1));
    for (j=0; j<=xNInt;j++) {
      USFseaF0alt[i][j]=0;
    }
  }
  mass = 940.0/cqm;
  //fprintf(stderr, "vibe check %d\n", index2++);
  USFvaldata(&USFval);
  for (i=0; i<=xNInt; i++) {
    t1=USFval[3][i];
    //fprintf(stderr, "%d %e\n", i, t1);
  }

  //fprintf(stderr, "vibe check %d\n", index2++);
  dataforGfouriers(0, &fourier_funcP,  &fourier_funcN, &fourier_funcPF0, &fourier_funcNF0);
  //fprintf(stderr, "vibe check here %d\n", index2++);
  dataforGvals(0, &Pvals, &Nvals, &mks, &ks, &pks);
  //fprintf(stderr, "vibe check %d\n", index2++);
  index =0;
  for (i=0; i<2*momentanums[0]; i++) {
    //fprintf(stderr, "Pvals=%f\n", Pvals[i]);
    if (fabs(Pvals[i])<1) {
      index = i;
      break;
    } else if (i==2*momentanums[0]-1) {
      fprintf(stderr, "Valence quark not detected");
      exit(1);
    }
  }
  val_e = Pvals[index];
  //fprintf(stderr, "vibe check %d\n", index2++);
  val_funcs = fourier_funcP[index];
  //fprintf(stderr, "vibe check %d\n", index2++);
  intarr3 = (double*)malloc(sizeof(double)*(xNInt+1));
  intarr4 = (double*)malloc(sizeof(double)*(xNInt+1));
  intarr1 = (double*)malloc(sizeof(double)*(kNInt+1));

  for (i=0; i<=xNInt; i++) {
    x=i*xdist/(double)xNInt;
    intarr3[i]=x*(USFval[0][i]+USFval[2][i]);
    intarr4[i]=x*(USFval[1][i]+USFval[3][i]);
  }
  //fprintf(stderr, "vibe check %d\n", index2++);
  for (i=0; i<=kNInt; i++) {
    p=i*kdist/(double)kNInt;
    intarr1[i]=-4.0/3.0/mass/mass*p*p*p*val_funcs[0][i]*val_funcs[1][i];
  }
  //fprintf(stderr, "vibe check %d\n", index2++);
  t1=booles(intarr3, xNInt)*xdist/(double)xNInt;
  t2=booles(intarr4, xNInt)*xdist/(double)xNInt;
  t3=2.0*val_e/mass/mass;
  t4=kdist/(double)kNInt*booles(intarr1, kNInt);
  fprintf(stderr, "%f %f %f %f\n", t3, t1, t4, t2);

  t1 = 3.0*5.0*mass/72.0;
  for (i=0; i<4; i++) {
    for (j=0; j<=xNInt; j++) {
      USFval[i][j]*=t1;
    }
  }
  for (i=0; i<=xNInt; i++) {
    x=i*xdist/(double)xNInt;
    intarr3[i]=x*(USFval[0][i]+USFval[2][i]);
  }
  t1=36.0/5.0*booles(intarr3, xNInt)*xdist/(double)xNInt;
  t2=val_e/mass*3.0;
  fprintf(stderr, "E_v/M=%f, [M^0_G]_v=%f\n", t2 ,t1);

  freeDataCollectedfouriers(fourier_funcP, fourier_funcN, fourier_funcPF0, fourier_funcNF0, 0, momentanums);
  freeDataCollectedvals(Pvals, Nvals, NULL, NULL, mks, ks, pks, 0, momentanums);
  free(intarr3);
  free(intarr4);
  free(intarr1);
  //fprintf(stderr, "E_v/M=%f, [M^0_G]_v=%f\n", t2 ,t1);
  //exit(0);
  srule0=0;
  srule1=0;
  srule0l=0;
  srule1l=0;
  srule0lalt=0;
  srule1lalt=0;
  for (G=0; G<grandspinmax; G++) {
    lenk = momentanums[G]*2;
    plenk= momentanums[G+1]*2;
    mlenk= 0;
    //fprintf(stderr, "dataforG 3\n");
    if (G!=0) {
      mlenk=momentanums[G-1]*2;
    }
    //fprintf(stderr, "dataforG 4\n");
    if (G==0) {
      lenP=momentanums[G]*2;
      lenN=momentanums[G+1]*2;
    } else {
      lenP= momentanums[G]*4;
      lenN = momentanums[G-1]*2 + momentanums[G+1]*2;
    }
    sepUSFseadata(G, &USFseaL, &USFseaF0L);
    sepUSFseaAMdata(G, &USFseaF0altL);
    dataforGfouriers(G, &fourier_funcP,  &fourier_funcN, &fourier_funcPF0, &fourier_funcNF0);
    dataforGvals(G, &Pvals, &Nvals, &mks, &ks, &pks);
    createF0vals(&PvalsF0, &NvalsF0, G, mks, ks, pks, momentanums);
    for (i=0; i<4; i++) {
      for (j=0; j<lenP+lenN; j++) {
        SumRuleFitter(USFseaF0altL[i][j], USFseaF0L[i][j]);
      }
    }
    for (i=0; i<4; i++) {
      for (j=0; j<=xNInt; j++) {
        USFsea[i][j]  =0;
        USFseaF0[i][j]=0;
        for (k=0; k<lenP+lenN; k++) {
          USFseaF0alt[i][j]+=USFseaF0altL[i][k][j];
          USFsea[i][j]+=USFseaL[i][k][j];
          USFseaF0[i][j]+=USFseaF0L[i][k][j];
        }
      }
    }
    intarr2=(double*)malloc(sizeof(double)*(kNInt+1));
    intarr1=(double*)malloc(sizeof(double)*(kNInt+1));
    for (i=0; i<lenP; i++) {
      for (pnum=0;pnum<=kNInt;pnum++) {
        p=(double)pnum/(double)kNInt*kdist;
        if (G==0) {
          intarr2[pnum]=-2.0*kdist*fourier_funcP[i][0][pnum]*fourier_funcP[i][1][pnum];
          intarr1[pnum]=-2.0*kdist*fourier_funcPF0[i][0][pnum]*fourier_funcPF0[i][1][pnum];
        } else {
          intarr2[pnum] =fourier_funcP[i][0][pnum]*fourier_funcP[i][2][pnum];
          intarr2[pnum]+=fourier_funcP[i][1][pnum]*fourier_funcP[i][3][pnum];
          intarr2[pnum]*=-2.0*kdist*(G*2.0+1.0);

          intarr1[pnum] =fourier_funcPF0[i][0][pnum]*fourier_funcPF0[i][2][pnum];
          intarr1[pnum]+=fourier_funcPF0[i][1][pnum]*fourier_funcPF0[i][3][pnum];
          intarr1[pnum]*=-2.0*kdist*(G*2.0+1.0);
        }
        temp = sqrt(Pvals[i]*Pvals[i]+lambda*lambda);
        intarr2[pnum] *= 2.0*Pvals[i]/fabs(Pvals[i])*(1.0-fabs(Pvals[i])/temp-0.5*fabs(Pvals[i])*lambda*lambda/temp/temp/temp);
        intarr2[pnum] *= -1.0/(3.0*mass*mass)*p*p*p;

        temp = sqrt(PvalsF0[i]*PvalsF0[i]+lambda*lambda);
        intarr1[pnum] *= 2.0*PvalsF0[i]/fabs(PvalsF0[i])*(1.0-fabs(PvalsF0[i])/temp-0.5*fabs(PvalsF0[i])*lambda*lambda/temp/temp/temp);
        intarr1[pnum] *= -1.0/(3.0*mass*mass)*p*p*p;
      }
      srule1 += booles(intarr2, kNInt)/(double)kNInt - booles(intarr1,kNInt)/(double)kNInt;

    }
    for (i=0; i<lenN; i++) {
      for (pnum=0;pnum<=kNInt;pnum++) {
        p=(double)pnum*kdist/(double)kNInt;
        if (G==0) {
          intarr2[pnum] =-2.0*kdist*fourier_funcN[i][0][pnum]*fourier_funcN[i][1][pnum];
          intarr1[pnum] =-2.0*kdist*fourier_funcNF0[i][0][pnum]*fourier_funcNF0[i][1][pnum];
        } else {
          intarr2[pnum] =fourier_funcN[i][0][pnum]*fourier_funcN[i][2][pnum];
          intarr2[pnum]+=fourier_funcN[i][1][pnum]*fourier_funcN[i][3][pnum];
          intarr2[pnum]*=-2.0*kdist*(2.0*G+1.0);

          intarr1[pnum] =fourier_funcNF0[i][0][pnum]*fourier_funcNF0[i][2][pnum];
          intarr1[pnum]+=fourier_funcNF0[i][1][pnum]*fourier_funcNF0[i][3][pnum];
          intarr1[pnum]*=-2.0*kdist*(2.0*G+1.0);
        }
        temp = sqrt(Nvals[i]*Nvals[i]+lambda*lambda);
        intarr2[pnum] *= 2.0*Nvals[i]/fabs(Nvals[i])*(1.0-fabs(Nvals[i])/temp-0.5*fabs(Nvals[i])*lambda*lambda/temp/temp/temp);
        intarr2[pnum] *= -1.0/(3.0*mass*mass)*p*p*p;

        temp = sqrt(NvalsF0[i]*NvalsF0[i]+lambda*lambda);
        intarr1[pnum] *= 2.0*NvalsF0[i]/fabs(NvalsF0[i])*(1.0-fabs(NvalsF0[i])/temp-0.5*fabs(NvalsF0[i])*lambda*lambda/temp/temp/temp);
        intarr1[pnum] *= -1.0/(3.0*mass*mass)*p*p*p;
      }
      srule1 += booles(intarr2, kNInt)/(double)kNInt - booles(intarr1,kNInt)/(double)kNInt;
    }
    free(intarr1);
    free(intarr2);
    intarr1=(double*)malloc(sizeof(double)*(xNInt+1));
    intarr2=(double*)malloc(sizeof(double)*(xNInt+1));
    for (i=0; i<=xNInt; i++) {
      intarr1[i]=0;
      intarr2[i]=0;
      x=(double)i/(double)xNInt*xdist;
      for (j=0; j<lenP+lenN; j++) {
        intarr1[i]+=x*(USFseaL[0][j][i]-USFseaF0L[0][j][i] + USFseaL[2][j][i]-USFseaF0L[2][j][i]);
        intarr2[i]+=x*(USFseaL[1][j][i]-USFseaF0L[1][j][i] + USFseaL[3][j][i]-USFseaF0L[3][j][i]);
      }
    }
    for (i=0; i<lenP; i++) {
      srule0 += (2.0*G+1.0)*(fabs(Pvals[i]) - sqrt(Pvals[i]*Pvals[i]+lambda*lambda) + 0.5*lambda*lambda/sqrt(Pvals[i]*Pvals[i]+lambda*lambda))/mass/mass;
      srule0 += -(2.0*G+1.0)*(fabs(PvalsF0[i]) - sqrt(PvalsF0[i]*PvalsF0[i]+lambda*lambda) + 0.5*lambda*lambda/sqrt(PvalsF0[i]*PvalsF0[i]+lambda*lambda))/mass/mass;
    }
    for (i=0; i<lenN; i++) {
      srule0 +=(2.0*G+1.0)*(fabs(Nvals[i]) - sqrt(Nvals[i]*Nvals[i]+lambda*lambda) + 0.5*lambda*lambda/sqrt(Nvals[i]*Nvals[i]+lambda*lambda))/mass/mass;
      srule0 += -(2.0*G+1.0)*(fabs(NvalsF0[i]) - sqrt(NvalsF0[i]*NvalsF0[i]+lambda*lambda) + 0.5*lambda*lambda/sqrt(NvalsF0[i]*NvalsF0[i]+lambda*lambda))/mass/mass;
    }
    srule0l+=0.5*booles(intarr1, xNInt)*xdist/(double)xNInt;
    srule1l+=0.5*booles(intarr2, xNInt)*xdist/(double)xNInt;
    //fprintf(stderr, "%20.10f %20.10f %20.10f %20.10f\n", -srule0, srule0l, 0.5*srule1, srule1l);
    for (i=0; i<=xNInt; i++) {
      intarr1[i]=0;
      intarr2[i]=0;
      x=(double)i/(double)xNInt*xdist;
      for (j=0; j<lenP+lenN; j++) {
        intarr1[i]+=x*(USFseaL[0][j][i]-USFseaF0altL[0][j][i] + USFseaL[2][j][i]-USFseaF0altL[2][j][i]);
        intarr2[i]+=x*(USFseaL[1][j][i]-USFseaF0altL[1][j][i] + USFseaL[3][j][i]-USFseaF0altL[3][j][i]);
      }
    }
    srule0lalt+=0.5*booles(intarr1, xNInt)*xdist/(double)xNInt;
    srule1lalt+=0.5*booles(intarr2, xNInt)*xdist/(double)xNInt;
    fprintf(stderr, "(%20.10f %20.10f %20.10f) (%20.10f %20.10f %20.10f)\n", -srule0, srule0l, srule0lalt, 0.5*srule1, srule1l, srule1lalt);

    t1=(3.0*mass*5.0/144.0);
    for (i=0; i<lenP+lenN; i++) {
      sprintf(filename, "formatted_data/USFseaLG%dL%d_%d.txt", G, i, expnum);
      USFseaLfile = fopen(filename, "w+");
      sprintf(filename, "formatted_data/USFseaF0LG%dL%d_%d.txt", G, i, expnum);
      USFseaF0Lfile = fopen(filename, "w+");
      sprintf(filename, "formatted_data/USFseaF0altMLG%dL%d_%d_%d.txt", G, i, expnum, massnum);
      USFseaF0altLfile = fopen(filename, "w+");

      fprintf(USFseaLfile, "%20.16f", t1*USFseaL[0][i][0]);
      fprintf(USFseaF0Lfile, "%20.16f", t1*USFseaF0L[0][i][0]);
      fprintf(USFseaF0altLfile, "%20.16f", t1*USFseaF0altL[0][i][0]);
      for (j=1; j<=xNIntp; j++) {
        fprintf(USFseaLfile, ",%20.16f", t1*USFseaL[0][i][j]);
        fprintf(USFseaF0Lfile, ",%20.16f", t1*USFseaF0L[0][i][j]);
        fprintf(USFseaF0altLfile, ",%20.16f", t1*USFseaF0altL[0][i][j]);
      }
      for (j=1; j<4; j++) {
        fprintf(USFseaLfile, "\n%20.16f", t1*USFseaL[j][i][0]);
        fprintf(USFseaF0Lfile, "\n%20.16f", t1*USFseaF0L[j][i][0]);
        fprintf(USFseaF0altLfile, "\n%20.16f", t1*USFseaF0altL[j][i][0]);
        for (k=1; k<=xNIntp; k++) {
          fprintf(USFseaLfile, ",%20.16f", t1*USFseaL[j][i][k]);
          fprintf(USFseaF0Lfile, ",%20.16f", t1*USFseaF0L[j][i][k]);
          fprintf(USFseaF0altLfile, ",%20.16f", t1*USFseaF0altL[j][i][k]);
        }
      }
      fclose(USFseaLfile);
      fclose(USFseaF0Lfile);
      fclose(USFseaF0altLfile);
    }

    for (i=0; i<4; i++) {
      for (j=0; j<lenP+lenN; j++) {
        free(USFseaL[i][j]);
        free(USFseaF0L[i][j]);
        free(USFseaF0altL[i][j]);
      }
      free(USFseaL[i]);
      free(USFseaF0L[i]);
      free(USFseaF0altL[i]);
    }
    free(USFseaL);
    free(USFseaF0L);
    free(USFseaF0altL);

    freeDataCollectedfouriers(fourier_funcP, fourier_funcN, fourier_funcPF0, fourier_funcNF0, G, momentanums);
    freeDataCollectedvals(Pvals, Nvals, PvalsF0, NvalsF0, mks, ks, pks, G, momentanums);
    free(intarr1);
    free(intarr2);
    //fprintf(stderr, "here\n");
  }

  intarr1=(double*)malloc(sizeof(double)*(xNInt+1));
  intarr2=(double*)malloc(sizeof(double)*(xNInt+1));
  for (i=0; i<=xNInt; i++) {
    x=(double)i/(double)xNInt*xdist;
    intarr1[i]=x*(USFsea[0][i]-USFseaF0[0][i] + USFsea[2][i]-USFseaF0[2][i]);
    intarr2[i]=x*(USFsea[1][i]-USFseaF0[1][i] + USFsea[3][i]-USFseaF0[3][i]);
  }
  t1=0.5*booles(intarr1, xNInt)*xdist/(double)xNInt;
  t2=0.5*booles(intarr2, xNInt)*xdist/(double)xNInt;
  for (i=0; i<=xNInt; i++) {
    x=(double)i/(double)xNInt*xdist;
    intarr1[i]=x*(USFsea[0][i]-USFseaF0alt[0][i] + USFsea[2][i]-USFseaF0alt[2][i]);
    intarr2[i]=x*(USFsea[1][i]-USFseaF0alt[1][i] + USFsea[3][i]-USFseaF0alt[3][i]);
  }
  t3=0.5*booles(intarr1, xNInt)*xdist/(double)xNInt;
  t4=0.5*booles(intarr2, xNInt)*xdist/(double)xNInt;
  fprintf(stderr, "%f %f %f %f\n", t1, t3, t2, t4);


  free(intarr1);
  free(intarr2);
  t1=(3.0*mass*5.0/144.0);
  for (i=0; i<4; i++) {
    for (j=0; j<=xNInt; j++) {
      USFsea[i][j]*=t1;
      USFseaF0[i][j]*=t1;
      USFseaF0alt[i][j]*=t1;
    }
  }
  intarr1=(double*)malloc(sizeof(double)*(xNInt+1));
  intarr2=(double*)malloc(sizeof(double)*(xNInt+1));
  for (i=0; i<=xNInt; i++) {
    x=(double)i/(double)xNInt*xdist;
    intarr1[i]=x*(USFsea[0][i]+USFsea[2][i]-USFseaF0[0][i]-USFseaF0[2][i]);
    intarr2[i]=x*(USFsea[0][i]+USFsea[2][i]-USFseaF0alt[0][i]-USFseaF0alt[2][i]);
  }
  t1=36.0/5.0*booles(intarr1, xNInt)*xdist/(double)xNInt;
  t2=36.0/5.0*booles(intarr2, xNInt)*xdist/(double)xNInt;
  chiralangle_from_file(&thetas);

  t3 = compute_sea_energy(thetas, grandspinmax, momentanums);
  fprintf(stderr, "%f %f %f\n", t3*cqm*Nc/mass, t1, t2);
  free(thetas);
  sprintf(filename, "formatted_data/USFsea_%d.txt", expnum);
  USFseafile = fopen(filename, "w+");
  sprintf(filename, "formatted_data/USFseaF0_%d.txt", expnum);
  USFseaF0file = fopen(filename, "w+");
  sprintf(filename, "formatted_data/USFseaF0altM_%d_%d.txt", expnum, massnum);
  USFseaF0altfile = fopen(filename, "w+");
  fprintf(USFseafile, "%20.16f", USFsea[0][0]);
  fprintf(USFseaF0file, "%20.16f", USFseaF0[0][0]);
  fprintf(USFseaF0altfile, "%20.16f", USFseaF0alt[0][0]);
  for (i=1; i<=xNIntp; i++) {
    fprintf(USFseafile, ",%20.16f", USFsea[0][i]);
    fprintf(USFseaF0file, ",%20.16f", USFseaF0[0][i]);
    fprintf(USFseaF0altfile, ",%20.16f", USFseaF0alt[0][i]);
  }
  for (i=1; i<4; i++) {
    fprintf(USFseafile, "\n%20.16f", USFsea[i][0]);
    fprintf(USFseaF0file, "\n%20.16f", USFseaF0[i][0]);
    fprintf(USFseaF0altfile, "\n%20.16f", USFseaF0alt[i][0]);
    for (j=1; j<=xNIntp; j++) {
      fprintf(USFseafile, ",%20.16f", USFsea[i][j]);
      fprintf(USFseaF0file, ",%20.16f", USFseaF0[i][j]);
      fprintf(USFseaF0altfile, ",%20.16f", USFseaF0alt[i][j]);
    }
  }

  for (i=0; i<4; i++) {
    free(USFsea[i]);
    free(USFseaF0[i]);
    free(USFseaF0alt[i]);
  }
  free(USFsea);
  free(USFseaF0);
  free(USFseaF0alt);

  fclose(USFseafile);
  fclose(USFseaF0file);
  fclose(USFseaF0altfile);
}
void altAverageFormat(char printIMF) {
  unsigned short grandspinmax, *momentanums, i, j, k, l, G, lenP, lenN, lenk, plenk, mlenk, index, pnum, index2, xNIntp;
  double ***USFseaL, ***USFseaF0L, ***USFseaF0Lalt, **USFsea, **USFseaF0, **USFseaF0alt, **USFval, **USFseaG, **USFseaF0G, **USFseaF0Galt;
  double ***fourier_funcP, ***fourier_funcN, ***fourier_funcPF0, ***fourier_funcNF0, **val_funcs;
  double *Pvals, *Nvals, *mks, *ks, *pks, *PvalsF0, *NvalsF0, *pholder;
  double *intarr1, *intarr2, *intarr3, *intarr4; //ARRAYS FOR INTEGRALS
  double val_e, x, p, mass, t1, t2, t3, t4, temp, srule0, srule1, srule0l, srule1l, srule0lalt, srule1lalt, srule0ltot, srule1ltot, srule0lalttot, srule1lalttot;
  FILE *USFseafile, *USFseaF0file, *USFseaF0altfile;
  FILE *USFseaLfile, *USFseaF0Lfile, *USFseaF0altLfile;
  FILE *USFseaGfile, *USFseaF0Gfile, *USFseaF0altGfile;
  char filename[100];
  index2=0;
  mass=940.0/cqm;
  xNIntp=floor(graphcut/xdist*xNInt);
  if (printIMF) {
    xNIntp=floor(IMFgraphcut*xNInt);
  }
  //fprintf(stderr, "vibe check %d\n", index2++);
  Gspin_momentanums(&grandspinmax, &momentanums);
  //fprintf(stderr, "vibe check %d\n", index2++);
  USFseaF0alt = (double**)malloc(sizeof(double*)*4);
  for (i=0; i<4; i++) {
    USFseaF0alt[i] = (double*)malloc(sizeof(double)*(xNInt+1));
    for (j=0; j<=xNInt;j++) {
      USFseaF0alt[i][j]=0;
    }
  }
  //fprintf(stderr, "vibe check %d\n", index2++);
  USFseadata(&USFsea, &USFseaF0);
  
  mass = 940.0/cqm;
  //fprintf(stderr, "vibe check %d\n", index2++);
  USFvaldata(&USFval);
  dataforGvals(0, &Pvals, &Nvals, &mks, &ks, &pks);
  index =0;
  for (i=0; i<2*momentanums[0]; i++) {
    //fprintf(stderr, "Pvals=%f\n", Pvals[i]);
    if (fabs(Pvals[i])<1) {
      index = i;
      break;
    } else if (i==2*momentanums[0]-1) {
      fprintf(stderr, "Valence quark not detected");
      exit(1);
    }
  }
  val_e = Pvals[index];
  //fprintf(stderr, "vibe check %d\n", index2++);
  intarr3 = (double*)malloc(sizeof(double)*(xNInt+1));
  intarr4 = (double*)malloc(sizeof(double)*(xNInt+1));
  intarr1 = (double*)malloc(sizeof(double)*(kNInt+1));
  for (i=0; i<=xNInt; i++) {
    x=i*xdist/(double)xNInt;
    intarr3[i]=x*(USFval[0][i]+USFval[2][i]);
    intarr4[i]=x*(USFval[1][i]+USFval[3][i]);
  }
  t1=booles(intarr3, xNInt)*xdist/(double)xNInt;
  t2=booles(intarr4, xNInt)*xdist/(double)xNInt;
  t3=2.0*val_e/mass/mass;
  fprintf(stderr, "%f %f\n", t1, t2);

  t1 = 3.0*5.0*mass/72.0;
  for (i=0; i<4; i++) {
    for (j=0; j<=xNInt; j++) {
      USFval[i][j]*=t1;
    }
  }
  for (i=0; i<=xNInt; i++) {
    x=i*xdist/(double)xNInt;
    intarr3[i]=x*(USFval[0][i]+USFval[2][i]);
  }
  t1=36.0/5.0*booles(intarr3, xNInt)*xdist/(double)xNInt;
  t2=val_e/mass*3.0;
  fprintf(stderr, "E_v/M=%f, [M^0_G]_v=%f\n", t2 ,t1);
  freeDataCollectedvals(Pvals, Nvals, NULL, NULL, mks, ks, pks, 0, momentanums);
  free(intarr3);
  free(intarr4);
  free(intarr1);
  srule0=0;
  srule1=0;
  srule0l=0;
  srule1l=0;
  srule0lalt=0;
  srule1lalt=0;
  srule0ltot=0;
  srule1ltot=0;
  srule0lalttot=0;
  srule1lalttot=0;
  USFseaG = (double**)malloc(sizeof(double*)*4);
  USFseaF0G = (double**)malloc(sizeof(double*)*4);
  USFseaF0Galt = (double**)malloc(sizeof(double*)*4);
  for (i=0; i<4; i++) {
    USFseaG[i]=(double*)malloc(sizeof(double)*(xNInt+1));
    USFseaF0G[i]=(double*)malloc(sizeof(double)*(xNInt+1));
    USFseaF0Galt[i]=(double*)malloc(sizeof(double)*(xNInt+1));
  }
  for (G=0; G<=10; G++) {
    lenk = momentanums[G]*2;
    plenk= momentanums[G+1]*2;
    mlenk= 0;
    if (G!=0) {
      mlenk=momentanums[G-1]*2;
    }
    if (G==0) {
      lenP=momentanums[G]*2;
      lenN=momentanums[G+1]*2;
    } else {
      lenP= momentanums[G]*4;
      lenN = momentanums[G-1]*2 + momentanums[G+1]*2;
    }
    fprintf(stderr, "Reading in files for G=%d.\n", G);
    sepUSFseadata(G, &USFseaL, &USFseaF0L);
    dataforGvals(G, &Pvals, &Nvals, &mks, &ks, &pks);
    createF0vals(&PvalsF0, &NvalsF0, G, mks, ks, pks, momentanums);
    fprintf(stderr, "Done reading in files for G=%d.\n", G);
    USFseaF0Lalt = (double***)malloc(sizeof(double**)*4);
    for (i=0; i<4; i++) {
      USFseaF0Lalt[i]=(double**)malloc(sizeof(double*)*(lenP+lenN));
      for (j=0; j<lenP+lenN; j++) {
        flatten_function(&USFseaF0Lalt[i][j], USFseaF0L[i][j], xNInt, xdist);
      }
    }
    fprintf(stderr, "done Flatten function\n");
    for (i=0; i<4; i++) {
      for (j=0; j<=xNInt; j++) {
        for (k=0; k<lenP+lenN; k++) {
          USFseaF0alt[i][j]+=USFseaF0Lalt[i][k][j];
        }
      }
    }
    intarr1=(double*)malloc(sizeof(double)*(xNInt+1));
    intarr2=(double*)malloc(sizeof(double)*(xNInt+1));
    intarr3=(double*)malloc(sizeof(double)*(xNInt+1));
    intarr4=(double*)malloc(sizeof(double)*(xNInt+1));
    for (i=0; i<=xNInt; i++) {
      intarr1[i]=0;
      intarr2[i]=0;
      intarr3[i]=0;
      intarr4[i]=0;
      x=(double)i/(double)xNInt*xdist;
      for (j=0; j<lenP+lenN; j++) {
        intarr1[i]+=x*(USFseaL[0][j][i]-USFseaF0L[0][j][i] + USFseaL[2][j][i]-USFseaF0L[2][j][i]);
        intarr2[i]+=x*(USFseaL[1][j][i]-USFseaF0L[1][j][i] + USFseaL[3][j][i]-USFseaF0L[3][j][i]);
        intarr3[i]+=x*(USFseaL[0][j][i]-USFseaF0Lalt[0][j][i] + USFseaL[2][j][i]-USFseaF0Lalt[2][j][i]);
        intarr4[i]+=x*(USFseaL[1][j][i]-USFseaF0Lalt[1][j][i] + USFseaL[3][j][i]-USFseaF0Lalt[3][j][i]);
      }
    }
    for (i=0; i<lenP; i++) {
      srule0 += (2.0*G+1.0)*(fabs(Pvals[i]) - sqrt(Pvals[i]*Pvals[i]+lambda*lambda) + 0.5*lambda*lambda/sqrt(Pvals[i]*Pvals[i]+lambda*lambda))/mass/mass;
      srule0 += -(2.0*G+1.0)*(fabs(PvalsF0[i]) - sqrt(PvalsF0[i]*PvalsF0[i]+lambda*lambda) + 0.5*lambda*lambda/sqrt(PvalsF0[i]*PvalsF0[i]+lambda*lambda))/mass/mass;
    }
    for (i=0; i<lenN; i++) {
      srule0 +=(2.0*G+1.0)*(fabs(Nvals[i]) - sqrt(Nvals[i]*Nvals[i]+lambda*lambda) + 0.5*lambda*lambda/sqrt(Nvals[i]*Nvals[i]+lambda*lambda))/mass/mass;
      srule0 += -(2.0*G+1.0)*(fabs(NvalsF0[i]) - sqrt(NvalsF0[i]*NvalsF0[i]+lambda*lambda) + 0.5*lambda*lambda/sqrt(NvalsF0[i]*NvalsF0[i]+lambda*lambda))/mass/mass;
    }
    srule0l+=0.5*booles(intarr1, xNInt)*xdist/(double)xNInt;
    srule1l+=0.5*booles(intarr2, xNInt)*xdist/(double)xNInt;
    srule0lalt+=0.5*booles(intarr3, xNInt)*xdist/(double)xNInt;
    srule1lalt+=0.5*booles(intarr4, xNInt)*xdist/(double)xNInt;
    fprintf(stderr, "%20.10f %20.10f %20.16f %20.16f\n", srule0l, srule0lalt, srule1l, srule1lalt);
    t1=(3.0*mass*5.0/144.0);
    for (i=0; i<4; i++) {
      for (j=0; j<=xNInt; j++) {
        USFseaG[i][j]=0;
        USFseaF0G[i][j]=0;
        USFseaF0Galt[i][j]=0;
      }
    }

    for (i=0; i<lenP+lenN; i++) {
      for (j=0; j<4; j++) {
        for (k=0; k<=xNInt; k++) {
          USFseaG[j][k]+=USFseaL[j][i][k];
          USFseaF0G[j][k]+=USFseaF0L[j][i][k];
          USFseaF0Galt[j][k]+=USFseaF0Lalt[j][i][k];
        }
      }
      if (!printIMF) {
        sprintf(filename, "formatted_data/expnum%d/USFseaLG%dL%d_%d.txt", expnum, G, i, expnum);
        USFseaLfile = fopen(filename, "w+");
        sprintf(filename, "formatted_data/expnum%d/USFseaF0LG%dL%d_%d.txt", expnum, G, i, expnum);
        USFseaF0Lfile = fopen(filename, "w+");
        sprintf(filename, "formatted_data/expnum%d/USFseaF0altLG%dL%d_%d.txt", expnum, G, i, expnum);
        USFseaF0altLfile = fopen(filename, "w+");
      
        fprintf(USFseaLfile, "%20.16f", t1*USFseaL[0][i][0]);
        fprintf(USFseaF0Lfile, "%20.16f", t1*USFseaF0L[0][i][0]);
        fprintf(USFseaF0altLfile, "%20.16f", t1*USFseaF0Lalt[0][i][0]);
        for (j=1; j<=xNIntp; j++) {
          fprintf(USFseaLfile, ",%20.16f", t1*USFseaL[0][i][j]);
          fprintf(USFseaF0Lfile, ",%20.16f", t1*USFseaF0L[0][i][j]);
          fprintf(USFseaF0altLfile, ",%20.16f", t1*USFseaF0Lalt[0][i][j]);
        }
        for (j=1; j<4; j++) {
          fprintf(USFseaLfile, "\n%20.16f", t1*USFseaL[j][i][0]);
          fprintf(USFseaF0Lfile, "\n%20.16f", t1*USFseaF0L[j][i][0]);
          fprintf(USFseaF0altLfile, "\n%20.16f", t1*USFseaF0Lalt[j][i][0]);
          for (k=1; k<=xNIntp; k++) {
            fprintf(USFseaLfile, ",%20.16f", t1*USFseaL[j][i][k]);
            fprintf(USFseaF0Lfile, ",%20.16f", t1*USFseaF0L[j][i][k]);
            fprintf(USFseaF0altLfile, ",%20.16f", t1*USFseaF0Lalt[j][i][k]);
          }
        }
        fclose(USFseaLfile);
        fclose(USFseaF0Lfile);
        fclose(USFseaF0altLfile);
      }
    }

    for (i=0; i<=xNInt; i++) {
      x=(double)i/(double)xNInt*xdist;
      intarr1[i]=x*(USFseaG[0][i]-USFseaF0G[0][i] + USFseaG[2][i]-USFseaF0G[2][i]);
      intarr2[i]=x*(USFseaG[1][i]-USFseaF0G[1][i] + USFseaG[3][i]-USFseaF0G[3][i]);
      intarr3[i]=x*(USFseaG[0][i]-USFseaF0Galt[0][i] + USFseaG[2][i]-USFseaF0Galt[2][i]);
      intarr4[i]=x*(USFseaG[1][i]-USFseaF0Galt[1][i] + USFseaG[3][i]-USFseaF0Galt[3][i]);
    }
    srule0ltot+=0.5*booles(intarr1, xNInt)*xdist/(double)xNInt;
    srule1ltot+=0.5*booles(intarr2, xNInt)*xdist/(double)xNInt;
    srule0lalttot+=0.5*booles(intarr3, xNInt)*xdist/(double)xNInt;
    srule1lalttot+=0.5*booles(intarr4, xNInt)*xdist/(double)xNInt;
    fprintf(stderr, "%20.10f %20.10f %20.16f %20.16f, %d\n\n", srule0ltot, srule0lalttot, srule1ltot, srule1lalttot, xNIntp);
    
    if (printIMF) {
      sprintf(filename, "formatted_data/expnum%d/USFseaIMFG%d_%d.txt", expnum, G, expnum);
      USFseaGfile=fopen(filename, "w+");
      sprintf(filename, "formatted_data/expnum%d/USFseaIMFF0G%d_%d.txt", expnum, G, expnum);
      USFseaF0Gfile=fopen(filename, "w+");
      sprintf(filename, "formatted_data/expnum%d/USFseaIMFF0altG%d_%d.txt", expnum, G, expnum);
      USFseaF0altGfile=fopen(filename, "w+");
      fprintf(stderr, "converting to IMF for G=%d\n", G);
      for (i=0; i<4; i++) {
        RF_to_IMF(&pholder, USFseaG[i]);
        free(USFseaG[i]);
        USFseaG[i]=pholder;
        RF_to_IMF(&pholder, USFseaF0G[i]);
        free(USFseaF0G[i]);
        USFseaF0G[i]=pholder;
        RF_to_IMF(&pholder, USFseaF0Galt[i]);
        free(USFseaF0Galt[i]);
        USFseaF0Galt[i]=pholder;
        fprintf(stderr, "%d converted\n", i+1);
      }
      fprintf(stderr, "done converting to IMF\n");
    } else {
      sprintf(filename, "formatted_data/expnum%d/USFseaG%d_%d.txt", expnum, G, expnum);
      USFseaGfile=fopen(filename, "w+");
      sprintf(filename, "formatted_data/expnum%d/USFseaF0G%d_%d.txt", expnum, G, expnum);
      USFseaF0Gfile=fopen(filename, "w+");
      sprintf(filename, "formatted_data/expnum%d/USFseaF0altG%d_%d.txt", expnum, G, expnum);
      USFseaF0altGfile=fopen(filename, "w+");
    }
    //t1=1.0;
    
    fprintf(USFseaGfile, "%20.16f", t1*USFseaG[0][0]);
    fprintf(USFseaF0Gfile, "%20.16f", t1*USFseaF0G[0][0]);
    fprintf(USFseaF0altGfile, "%20.16f", t1*USFseaF0Galt[0][0]);
    for (k=1; k<=xNIntp; k++) {
      fprintf(USFseaGfile, ",%20.16f", t1*USFseaG[0][k]);
      fprintf(USFseaF0Gfile, ",%20.16f", t1*USFseaF0G[0][k]);
      fprintf(USFseaF0altGfile, ",%20.16f", t1*USFseaF0Galt[0][k]);
    }
    for (j=1; j<4; j++) {
      fprintf(USFseaGfile, "\n%20.16f", t1*USFseaG[j][0]);
      fprintf(USFseaF0Gfile, "\n%20.16f", t1*USFseaF0G[j][0]);
      fprintf(USFseaF0altGfile, "\n%20.16f", t1*USFseaF0Galt[j][0]);
      for (k=1; k<=xNIntp; k++) {
        fprintf(USFseaGfile, ",%20.16f", t1*USFseaG[j][k]);
        fprintf(USFseaF0Gfile, ",%20.16f", t1*USFseaF0G[j][k]);
        fprintf(USFseaF0altGfile, ",%20.16f", t1*USFseaF0Galt[j][k]);
      }
    }
    fclose(USFseaGfile);
    fclose(USFseaF0Gfile);
    fclose(USFseaF0altGfile);
    //freeDataCollectedfouriers(fourier_funcP, fourier_funcN, fourier_funcPF0, fourier_funcNF0, G, momentanums);
    freeDataCollectedvals(Pvals, Nvals, PvalsF0, NvalsF0, mks, ks, pks, G, momentanums);
    free(intarr1);
    free(intarr2);
    free(intarr3);
    free(intarr4);
    //fprintf(stderr, "here\n");
  }
  fprintf(stderr, "Completed grandspins\n");
  for (i=0; i<4; i++) {
    free(USFseaG[i]);
    free(USFseaF0G[i]);
    free(USFseaF0Galt[i]);
  }
  free(USFseaG);
  free(USFseaF0G);
  free(USFseaF0Galt);
  intarr1=(double*)malloc(sizeof(double)*(xNInt+1));
  intarr2=(double*)malloc(sizeof(double)*(xNInt+1));
  for (i=0; i<=xNInt; i++) {
    x=(double)i/(double)xNInt*xdist;
    intarr1[i]=x*(USFsea[0][i]-USFseaF0[0][i] + USFsea[2][i]-USFseaF0[2][i]);
    intarr2[i]=x*(USFsea[1][i]-USFseaF0[1][i] + USFsea[3][i]-USFseaF0[3][i]);
  }
  t1=0.5*booles(intarr1, xNInt)*xdist/(double)xNInt;
  t2=0.5*booles(intarr2, xNInt)*xdist/(double)xNInt;
  for (i=0; i<=xNInt; i++) {
    x=(double)i/(double)xNInt*xdist;
    intarr1[i]=x*(USFsea[0][i]-USFseaF0alt[0][i] + USFsea[2][i]-USFseaF0alt[2][i]);
    intarr2[i]=x*(USFsea[1][i]-USFseaF0alt[1][i] + USFsea[3][i]-USFseaF0alt[3][i]);
  }
  t3=0.5*booles(intarr1, xNInt)*xdist/(double)xNInt;
  t4=0.5*booles(intarr2, xNInt)*xdist/(double)xNInt;
  fprintf(stderr, "%f %f %f %f\n", t1, t3, t2, t4);


  free(intarr1);
  free(intarr2);
  t1=(3.0*mass*5.0/144.0);
  for (i=0; i<4; i++) {
    for (j=0; j<=xNInt; j++) {
      USFsea[i][j]*=t1;
      USFseaF0[i][j]*=t1;
      USFseaF0alt[i][j]*=t1;
    }
  }
  if (printIMF) {
    for (i=0; i<4; i++) {
      RF_to_IMF(&pholder, USFsea[i]);
      free(USFsea[i]);
      USFsea[i]=pholder;

      RF_to_IMF(&pholder, USFseaF0[i]);
      free(USFseaF0[i]);
      USFseaF0[i]=pholder;

      RF_to_IMF(&pholder, USFseaF0alt[i]);
      free(USFseaF0alt[i]);
      USFseaF0alt[i]=pholder;
    }
    sprintf(filename, "formatted_data/expnum%d/USFseaIMF_%d.txt", expnum, expnum);
    USFseafile = fopen(filename, "w+");
    sprintf(filename, "formatted_data/expnum%d/USFseaIMFF0_%d.txt", expnum, expnum);
    USFseaF0file = fopen(filename, "w+");
    sprintf(filename, "formatted_data/expnum%d/USFseaIMFF0alt_%d.txt", expnum, expnum);
  } else {
    sprintf(filename, "formatted_data/expnum%d/USFsea_%d.txt", expnum, expnum);
    USFseafile = fopen(filename, "w+");
    sprintf(filename, "formatted_data/expnum%d/USFseaF0_%d.txt", expnum, expnum);
    USFseaF0file = fopen(filename, "w+");
    sprintf(filename, "formatted_data/expnum%d/USFseaF0alt_%d.txt", expnum, expnum);
  }
  
  USFseaF0altfile = fopen(filename, "w+");
  fprintf(USFseafile, "%20.16f", USFsea[0][0]);
  fprintf(USFseaF0file, "%20.16f", USFseaF0[0][0]);
  fprintf(USFseaF0altfile, "%20.16f", USFseaF0alt[0][0]);
  for (i=1; i<=xNIntp; i++) {
    fprintf(USFseafile, ",%20.16f", USFsea[0][i]);
    fprintf(USFseaF0file, ",%20.16f", USFseaF0[0][i]);
    fprintf(USFseaF0altfile, ",%20.16f", USFseaF0alt[0][i]);
  }
  for (i=1; i<4; i++) {
    fprintf(USFseafile, "\n%20.16f", USFsea[i][0]);
    fprintf(USFseaF0file, "\n%20.16f", USFseaF0[i][0]);
    fprintf(USFseaF0altfile, "\n%20.16f", USFseaF0alt[i][0]);
    for (j=1; j<=xNIntp; j++) {
      fprintf(USFseafile, ",%20.16f", USFsea[i][j]);
      fprintf(USFseaF0file, ",%20.16f", USFseaF0[i][j]);
      fprintf(USFseaF0altfile, ",%20.16f", USFseaF0alt[i][j]);
    }
  }

  for (i=0; i<4; i++) {
    free(USFsea[i]);
    free(USFseaF0[i]);
    free(USFseaF0alt[i]);
  }
  free(USFsea);
  free(USFseaF0);
  free(USFseaF0alt);

  fclose(USFseafile);
  fclose(USFseaF0file);
  fclose(USFseaF0altfile);
}
void DefaultFormat() {
  unsigned short grandspinmax, *momentanums, i, j, k, l, G, lenP, lenN, lenk, plenk, mlenk, index, pnum, index2, xNIntp;
  double ***USFseaL, ***USFseaF0L, **USFsea, **USFseaF0, **USFval, **USFseaG, **USFseaF0G;
  double ***fourier_funcP, ***fourier_funcN, ***fourier_funcPF0, ***fourier_funcNF0, **val_funcs;
  double *Pvals, *Nvals, *mks, *ks, *pks, *PvalsF0, *NvalsF0;
  double *intarr1, *intarr2, *intarr3, *intarr4; //ARRAYS FOR INTEGRALS
  double val_e, x, p, mass, t1, t2, t3, t4, temp, srule0, srule1, srule0l, srule1l;
  FILE *USFseafile, *USFseaF0file;
  FILE *USFseaLfile, *USFseaF0Lfile;
  FILE *USFseaGfile, *USFseaF0Gfile;
  char filename[100];
  index2=0;
  mass=940.0/cqm;
  xNIntp=floor(graphcut/xdist*xNInt);
  //fprintf(stderr, "vibe check %d\n", index2++);
  Gspin_momentanums(&grandspinmax, &momentanums);
  //fprintf(stderr, "vibe check %d\n", index2++);
  USFseadata(&USFsea, &USFseaF0);
  mass = 940.0/cqm;
  //fprintf(stderr, "vibe check %d\n", index2++);
  USFvaldata(&USFval);
  for (i=0; i<=xNInt; i++) {
    t1=USFval[3][i];
    //fprintf(stderr, "%d %e\n", i, t1);
  }
  dataforGvals(0, &Pvals, &Nvals, &mks, &ks, &pks);
  index =0;
  for (i=0; i<2*momentanums[0]; i++) {
    //fprintf(stderr, "Pvals=%f\n", Pvals[i]);
    if (fabs(Pvals[i])<1) {
      index = i;
      break;
    } else if (i==2*momentanums[0]-1) {
      fprintf(stderr, "Valence quark not detected");
      exit(1);
    }
  }
  val_e = Pvals[index];
  //fprintf(stderr, "vibe check %d\n", index2++);
  intarr3 = (double*)malloc(sizeof(double)*(xNInt+1));
  intarr4 = (double*)malloc(sizeof(double)*(xNInt+1));
  intarr1 = (double*)malloc(sizeof(double)*(kNInt+1));
  for (i=0; i<=xNInt; i++) {
    x=i*xdist/(double)xNInt;
    intarr3[i]=x*(USFval[0][i]+USFval[2][i]);
    intarr4[i]=x*(USFval[1][i]+USFval[3][i]);
  }
  t1=booles(intarr3, xNInt)*xdist/(double)xNInt;
  t2=booles(intarr4, xNInt)*xdist/(double)xNInt;
  t3=2.0*val_e/mass/mass;
  fprintf(stderr, "%f %f\n", t1, t2);

  t1 = 3.0*5.0*mass/72.0;
  for (i=0; i<4; i++) {
    for (j=0; j<=xNInt; j++) {
      USFval[i][j]*=t1;
    }
  }
  for (i=0; i<=xNInt; i++) {
    x=i*xdist/(double)xNInt;
    intarr3[i]=x*(USFval[0][i]+USFval[2][i]);
  }
  t1=36.0/5.0*booles(intarr3, xNInt)*xdist/(double)xNInt;
  t2=val_e/mass*3.0;
  fprintf(stderr, "E_v/M=%f, [M^0_G]_v=%f\n", t2 ,t1);
  freeDataCollectedvals(Pvals, Nvals, NULL, NULL, mks, ks, pks, 0, momentanums);
  free(intarr3);
  free(intarr4);
  free(intarr1);
  srule0=0;
  srule1=0;
  srule0l=0;
  srule1l=0;
  USFseaG = (double**)malloc(sizeof(double*)*4);
  USFseaF0G = (double**)malloc(sizeof(double*)*4);
  for (i=0; i<4; i++) {
    USFseaG[i]=(double*)malloc(sizeof(double)*(xNIntp+1));
    USFseaF0G[i]=(double*)malloc(sizeof(double)*(xNIntp+1));
  }
  for (G=0; G<=-1; G++) {
    lenk = momentanums[G]*2;
    plenk= momentanums[G+1]*2;
    mlenk= 0;
    if (G!=0) {
      mlenk=momentanums[G-1]*2;
    }
    if (G==0) {
      lenP=momentanums[G]*2;
      lenN=momentanums[G+1]*2;
    } else {
      lenP= momentanums[G]*4;
      lenN = momentanums[G-1]*2 + momentanums[G+1]*2;
    }
    sepUSFseadata(G, &USFseaL, &USFseaF0L);
    dataforGvals(G, &Pvals, &Nvals, &mks, &ks, &pks);
    createF0vals(&PvalsF0, &NvalsF0, G, mks, ks, pks, momentanums);
    intarr1=(double*)malloc(sizeof(double)*(xNInt+1));
    intarr2=(double*)malloc(sizeof(double)*(xNInt+1));
    for (i=0; i<=xNInt; i++) {
      intarr1[i]=0;
      intarr2[i]=0;
      x=(double)i/(double)xNInt*xdist;
      for (j=0; j<lenP+lenN; j++) {
        intarr1[i]+=x*(USFseaL[0][j][i]-USFseaF0L[0][j][i] + USFseaL[2][j][i]-USFseaF0L[2][j][i]);
        intarr2[i]+=x*(USFseaL[1][j][i]-USFseaF0L[1][j][i] + USFseaL[3][j][i]-USFseaF0L[3][j][i]);
      }
    }
    for (i=0; i<lenP; i++) {
      srule0 += (2.0*G+1.0)*(fabs(Pvals[i]) - sqrt(Pvals[i]*Pvals[i]+lambda*lambda) + 0.5*lambda*lambda/sqrt(Pvals[i]*Pvals[i]+lambda*lambda))/mass/mass;
      srule0 += -(2.0*G+1.0)*(fabs(PvalsF0[i]) - sqrt(PvalsF0[i]*PvalsF0[i]+lambda*lambda) + 0.5*lambda*lambda/sqrt(PvalsF0[i]*PvalsF0[i]+lambda*lambda))/mass/mass;
    }
    for (i=0; i<lenN; i++) {
      srule0 +=(2.0*G+1.0)*(fabs(Nvals[i]) - sqrt(Nvals[i]*Nvals[i]+lambda*lambda) + 0.5*lambda*lambda/sqrt(Nvals[i]*Nvals[i]+lambda*lambda))/mass/mass;
      srule0 += -(2.0*G+1.0)*(fabs(NvalsF0[i]) - sqrt(NvalsF0[i]*NvalsF0[i]+lambda*lambda) + 0.5*lambda*lambda/sqrt(NvalsF0[i]*NvalsF0[i]+lambda*lambda))/mass/mass;
    }
    srule0l+=0.5*booles(intarr1, xNInt)*xdist/(double)xNInt;
    srule1l+=0.5*booles(intarr2, xNInt)*xdist/(double)xNInt;
    fprintf(stderr, "%d %20.10f %20.10f\n", G, srule0l, srule1l);
    t1=(3.0*mass*5.0/144.0);
    for (i=0; i<4; i++) {
      for (j=0; j<=xNIntp; j++) {
        USFseaG[i][j]=0;
        USFseaF0G[i][j]=0;
      }
    }
    for (i=0; i<lenP+lenN; i++) {
      for (j=0; j<4; j++) {
        for (k=0; k<=xNIntp; k++) {
          USFseaG[j][k]+=USFseaL[j][i][k];
          USFseaF0G[j][k]+=USFseaF0L[j][i][k];
        }
      }
      sprintf(filename, "formatted_data/expnum%d/USFseaLG%dL%d_%d.txt", expnum, G, i, expnum);
      USFseaLfile = fopen(filename, "w+");
      sprintf(filename, "formatted_data/expnum%d/USFseaF0LG%dL%d_%d.txt", expnum, G, i, expnum);
      USFseaF0Lfile = fopen(filename, "w+");
      fprintf(USFseaLfile, "%20.16f", t1*USFseaL[0][i][0]);
      fprintf(USFseaF0Lfile, "%20.16f", t1*USFseaF0L[0][i][0]);
      for (j=1; j<=xNIntp; j++) {
        fprintf(USFseaLfile, ",%20.16f", t1*USFseaL[0][i][j]);
        fprintf(USFseaF0Lfile, ",%20.16f", t1*USFseaF0L[0][i][j]);
      }
      for (j=1; j<4; j++) {
        fprintf(USFseaLfile, "\n%20.16f", t1*USFseaL[j][i][0]);
        fprintf(USFseaF0Lfile, "\n%20.16f", t1*USFseaF0L[j][i][0]);
        for (k=1; k<=xNIntp; k++) {
          fprintf(USFseaLfile, ",%20.16f", t1*USFseaL[j][i][k]);
          fprintf(USFseaF0Lfile, ",%20.16f", t1*USFseaF0L[j][i][k]);
        }
      }
      fclose(USFseaLfile);
      fclose(USFseaF0Lfile);
    }
    sprintf(filename, "formatted_data/expnum%d/USFseaG%d_%d.txt", expnum, G, expnum);
    USFseaGfile=fopen(filename, "w+");
    sprintf(filename, "formatted_data/expnum%d/USFseaF0G%d_%d.txt", expnum, G, expnum);
    USFseaF0Gfile=fopen(filename, "w+");
    fprintf(USFseaGfile, "%20.16f", t1*USFseaG[0][0]);
    fprintf(USFseaF0Gfile, "%20.16f", t1*USFseaF0G[0][0]);
    for (k=1; k<=xNIntp; k++) {
      fprintf(USFseaGfile, ",%20.16f", t1*USFseaG[0][k]);
      fprintf(USFseaF0Gfile, ",%20.16f", t1*USFseaF0G[0][k]);
    }
    for (j=1; j<4; j++) {
      fprintf(USFseaGfile, "\n%20.16f", t1*USFseaG[j][0]);
      fprintf(USFseaF0Gfile, "\n%20.16f", t1*USFseaF0G[j][0]);
      for (k=1; k<=xNIntp; k++) {
        fprintf(USFseaGfile, ",%20.16f", t1*USFseaG[j][k]);
        fprintf(USFseaF0Gfile, ",%20.16f", t1*USFseaF0G[j][k]);
      }
    }
    fclose(USFseaGfile);
    fclose(USFseaF0Gfile);
    //freeDataCollectedfouriers(fourier_funcP, fourier_funcN, fourier_funcPF0, fourier_funcNF0, G, momentanums);
    freeDataCollectedvals(Pvals, Nvals, PvalsF0, NvalsF0, mks, ks, pks, G, momentanums);
    free(intarr1);
    free(intarr2);
    //fprintf(stderr, "here\n");
  }
  for (i=0; i<4; i++) {
    free(USFseaG[i]);
    free(USFseaF0G[i]);
  }
  free(USFseaG);
  free(USFseaF0G);
  intarr1=(double*)malloc(sizeof(double)*(xNInt+1));
  intarr2=(double*)malloc(sizeof(double)*(xNInt+1));
  for (i=0; i<=xNInt; i++) {
    x=(double)i/(double)xNInt*xdist;
    intarr1[i]=x*(USFsea[0][i]-USFseaF0[0][i] + USFsea[2][i]-USFseaF0[2][i]);
    intarr2[i]=x*(USFsea[1][i]-USFseaF0[1][i] + USFsea[3][i]-USFseaF0[3][i]);
  }
  t1=0.5*booles(intarr1, xNInt)*xdist/(double)xNInt;
  t2=0.5*booles(intarr2, xNInt)*xdist/(double)xNInt;
  fprintf(stderr, "%f %f\n", t1, t2);


  free(intarr1);
  free(intarr2);
  t1=(3.0*mass*5.0/144.0);
  for (i=0; i<4; i++) {
    for (j=0; j<=xNInt; j++) {
      USFsea[i][j]*=t1;
      USFseaF0[i][j]*=t1;
    }
  }

  sprintf(filename, "formatted_data/expnum%d/USFsea_%d.txt", expnum, expnum);
  USFseafile = fopen(filename, "w+");
  sprintf(filename, "formatted_data/expnum%d/USFseaF0_%d.txt", expnum, expnum);
  USFseaF0file = fopen(filename, "w+");
  fprintf(USFseafile, "%20.16f", USFsea[0][0]);
  fprintf(USFseaF0file, "%20.16f", USFseaF0[0][0]);
  for (i=1; i<=xNIntp; i++) {
    fprintf(USFseafile, ",%20.16f", USFsea[0][i]);
    fprintf(USFseaF0file, ",%20.16f", USFseaF0[0][i]);
  }
  for (i=1; i<4; i++) {
    fprintf(USFseafile, "\n%20.16f", USFsea[i][0]);
    fprintf(USFseaF0file, "\n%20.16f", USFseaF0[i][0]);
    for (j=1; j<=xNIntp; j++) {
      fprintf(USFseafile, ",%20.16f", USFsea[i][j]);
      fprintf(USFseaF0file, ",%20.16f", USFseaF0[i][j]);
    }
  }

  for (i=0; i<4; i++) {
    free(USFsea[i]);
    free(USFseaF0[i]);
  }
  free(USFsea);
  free(USFseaF0);

  fclose(USFseafile);
  fclose(USFseaF0file);
}
/*
void altFormat(unsigned short altnum) {
  unsigned short grandspinmax, *momentanums, i, j, k, l, G, lenP, lenN, lenk, plenk, mlenk, index, pnum, index2, xNIntp;
  double ***USFseaL, ***USFseaF0L, ***USFseaF0Lalt, **USFsea, **USFseaF0, **USFseaF0alt, **USFval, **USFseaG, **USFseaF0G, **USFseaF0Galt;
  double ***fourier_funcP, ***fourier_funcN, ***fourier_funcPF0, ***fourier_funcNF0, **val_funcs;
  double *Pvals, *Nvals, *mks, *ks, *pks, *PvalsF0, *NvalsF0;
  double *intarr1, *intarr2, *intarr3, *intarr4; //ARRAYS FOR INTEGRALS
  double val_e, x, p, mass, t1, t2, t3, t4, temp, srule0, srule1, srule0l, srule1l, srule0lalt, srule1lalt, srule0ltot, srule1ltot, srule0lalttot, srule1lalttot;
  FILE *USFseafile, *USFseaF0file, *USFseaF0altfile;
  FILE *USFseaLfile, *USFseaF0Lfile, *USFseaF0altLfile;
  FILE *USFseaGfile, *USFseaF0Gfile, *USFseaF0altGfile;
  char filename[100];
  index2=0;
  mass=940.0/cqm;
  xNIntp=floor(graphcut/xdist*xNInt);
  //fprintf(stderr, "vibe check %d\n", index2++);
  Gspin_momentanums(&grandspinmax, &momentanums);
  //fprintf(stderr, "vibe check %d\n", index2++);
  USFseaF0alt = (double**)malloc(sizeof(double*)*4);
  for (i=0; i<4; i++) {
    USFseaF0alt[i] = (double*)malloc(sizeof(double)*(xNInt+1));
    for (j=0; j<=xNInt;j++) {
      USFseaF0alt[i][j]=0;
    }
  }
  //fprintf(stderr, "vibe check %d\n", index2++);
  USFseadata(&USFsea, &USFseaF0);
  mass = 940.0/cqm;
  //fprintf(stderr, "vibe check %d\n", index2++);
  USFvaldata(&USFval);
  for (i=0; i<=xNInt; i++) {
    t1=USFval[3][i];
    //fprintf(stderr, "%d %e\n", i, t1);
  }
  dataforGvals(0, &Pvals, &Nvals, &mks, &ks, &pks);
  index =0;
  for (i=0; i<2*momentanums[0]; i++) {
    //fprintf(stderr, "Pvals=%f\n", Pvals[i]);
    if (fabs(Pvals[i])<1) {
      index = i;
      break;
    } else if (i==2*momentanums[0]-1) {
      fprintf(stderr, "Valence quark not detected");
      exit(1);
    }
  }
  val_e = Pvals[index];
  //fprintf(stderr, "vibe check %d\n", index2++);
  intarr3 = (double*)malloc(sizeof(double)*(xNInt+1));
  intarr4 = (double*)malloc(sizeof(double)*(xNInt+1));
  intarr1 = (double*)malloc(sizeof(double)*(kNInt+1));
  for (i=0; i<=xNInt; i++) {
    x=i*xdist/(double)xNInt;
    intarr3[i]=x*(USFval[0][i]+USFval[2][i]);
    intarr4[i]=x*(USFval[1][i]+USFval[3][i]);
  }
  t1=booles(intarr3, xNInt)*xdist/(double)xNInt;
  t2=booles(intarr4, xNInt)*xdist/(double)xNInt;
  t3=2.0*val_e/mass/mass;
  fprintf(stderr, "%f %f\n", t1, t2);

  t1 = 3.0*5.0*mass/72.0;
  for (i=0; i<4; i++) {
    for (j=0; j<=xNInt; j++) {
      USFval[i][j]*=t1;
    }
  }
  for (i=0; i<=xNInt; i++) {
    x=i*xdist/(double)xNInt;
    intarr3[i]=x*(USFval[0][i]+USFval[2][i]);
  }
  t1=36.0/5.0*booles(intarr3, xNInt)*xdist/(double)xNInt;
  t2=val_e/mass*3.0;
  fprintf(stderr, "E_v/M=%f, [M^0_G]_v=%f\n", t2 ,t1);
  freeDataCollectedvals(Pvals, Nvals, NULL, NULL, mks, ks, pks, 0, momentanums);
  free(intarr3);
  free(intarr4);
  free(intarr1);
  srule0=0;
  srule1=0;
  srule0l=0;
  srule1l=0;
  srule0lalt=0;
  srule1lalt=0;
  srule0ltot=0;
  srule1ltot=0;
  srule0lalttot=0;
  srule1lalttot=0;
  USFseaG = (double**)malloc(sizeof(double*)*4);
  USFseaF0G = (double**)malloc(sizeof(double*)*4);
  USFseaF0Galt = (double**)malloc(sizeof(double*)*4);
  for (i=0; i<4; i++) {
    USFseaG[i]=(double*)malloc(sizeof(double)*(xNInt+1));
    USFseaF0G[i]=(double*)malloc(sizeof(double)*(xNInt+1));
    USFseaF0Galt[i]=(double*)malloc(sizeof(double)*(xNInt+1));
  }
  for (G=0; G<=grandspinmax; G++) {
    lenk = momentanums[G]*2;
    plenk= momentanums[G+1]*2;
    mlenk= 0;
    if (G!=0) {
      mlenk=momentanums[G-1]*2;
    }
    if (G==0) {
      lenP=momentanums[G]*2;
      lenN=momentanums[G+1]*2;
    } else {
      lenP= momentanums[G]*4;
      lenN = momentanums[G-1]*2 + momentanums[G+1]*2;
    }
    sepUSFseadata(G, &USFseaL, &USFseaF0L);
    dataforGvals(G, &Pvals, &Nvals, &mks, &ks, &pks);
    createF0vals(&PvalsF0, &NvalsF0, G, mks, ks, pks, momentanums);
    USFseaF0Lalt = (double***)malloc(sizeof(double**)*4);
    for (i=0; i<4; i++) {
      USFseaF0Lalt[i]=(double**)malloc(sizeof(double*)*(lenP+lenN));
      for (j=0; j<lenP+lenN; j++) {
        USFseaF0Lalt[i][j]=(double*)malloc(sizeof(double)*(xNInt+1));
        alt_function(&USFseaF0Lalt[i][j], USFseaL[i][j], USFseaF0L[i][j], xNInt, xdist, altnum);
      }
    }
    for (i=0; i<4; i++) {
      for (j=0; j<=xNInt; j++) {
        for (k=0; k<lenP+lenN; k++) {
          USFseaF0alt[i][j]+=USFseaF0Lalt[i][k][j];
        }
      }
    }
    intarr1=(double*)malloc(sizeof(double)*(xNInt+1));
    intarr2=(double*)malloc(sizeof(double)*(xNInt+1));
    intarr3=(double*)malloc(sizeof(double)*(xNInt+1));
    intarr4=(double*)malloc(sizeof(double)*(xNInt+1));
    for (i=0; i<=xNInt; i++) {
      intarr1[i]=0;
      intarr2[i]=0;
      intarr3[i]=0;
      intarr4[i]=0;
      x=(double)i/(double)xNInt*xdist;
      for (j=0; j<lenP+lenN; j++) {
        intarr1[i]+=x*(USFseaL[0][j][i]-USFseaF0L[0][j][i] + USFseaL[2][j][i]-USFseaF0L[2][j][i]);
        intarr2[i]+=x*(USFseaL[1][j][i]-USFseaF0L[1][j][i] + USFseaL[3][j][i]-USFseaF0L[3][j][i]);
        intarr3[i]+=x*(USFseaL[0][j][i]-USFseaF0Lalt[0][j][i] + USFseaL[2][j][i]-USFseaF0Lalt[2][j][i]);
        intarr4[i]+=x*(USFseaL[1][j][i]-USFseaF0Lalt[1][j][i] + USFseaL[3][j][i]-USFseaF0Lalt[3][j][i]);
      }
    }
    for (i=0; i<lenP; i++) {
      srule0 += (2.0*G+1.0)*(fabs(Pvals[i]) - sqrt(Pvals[i]*Pvals[i]+lambda*lambda) + 0.5*lambda*lambda/sqrt(Pvals[i]*Pvals[i]+lambda*lambda))/mass/mass;
      srule0 += -(2.0*G+1.0)*(fabs(PvalsF0[i]) - sqrt(PvalsF0[i]*PvalsF0[i]+lambda*lambda) + 0.5*lambda*lambda/sqrt(PvalsF0[i]*PvalsF0[i]+lambda*lambda))/mass/mass;
    }
    for (i=0; i<lenN; i++) {
      srule0 +=(2.0*G+1.0)*(fabs(Nvals[i]) - sqrt(Nvals[i]*Nvals[i]+lambda*lambda) + 0.5*lambda*lambda/sqrt(Nvals[i]*Nvals[i]+lambda*lambda))/mass/mass;
      srule0 += -(2.0*G+1.0)*(fabs(NvalsF0[i]) - sqrt(NvalsF0[i]*NvalsF0[i]+lambda*lambda) + 0.5*lambda*lambda/sqrt(NvalsF0[i]*NvalsF0[i]+lambda*lambda))/mass/mass;
    }
    srule0l+=0.5*booles(intarr1, xNInt)*xdist/(double)xNInt;
    srule1l+=0.5*booles(intarr2, xNInt)*xdist/(double)xNInt;
    srule0lalt+=0.5*booles(intarr3, xNInt)*xdist/(double)xNInt;
    srule1lalt+=0.5*booles(intarr4, xNInt)*xdist/(double)xNInt;
    fprintf(stderr, "%20.10f %20.10f %20.16f %20.16f\n", srule0l, srule0lalt, srule1l, srule1lalt);
    t1=(3.0*mass*5.0/144.0);
    for (i=0; i<4; i++) {
      for (j=0; j<=xNInt; j++) {
        USFseaG[i][j]=0;
        USFseaF0G[i][j]=0;
        USFseaF0Galt[i][j]=0;
      }
    }

    for (i=0; i<lenP+lenN; i++) {
      for (j=0; j<4; j++) {
        for (k=0; k<=xNInt; k++) {
          USFseaG[j][k]+=USFseaL[j][i][k];
          USFseaF0G[j][k]+=USFseaF0L[j][i][k];
          USFseaF0Galt[j][k]+=USFseaF0Lalt[j][i][k];
        }
      }
      sprintf(filename, "formatted_data/expnum%d/USFseaLG%dL%d_%d.txt", expnum, G, i, expnum);
      USFseaLfile = fopen(filename, "w+");
      sprintf(filename, "formatted_data/expnum%d/USFseaF0LG%dL%d_%d.txt", expnum, G, i, expnum);
      USFseaF0Lfile = fopen(filename, "w+");
      sprintf(filename, "formatted_data/expnum%d/USFseaF0altLG%dL%d_%d.txt", expnum, G, i, expnum);
      USFseaF0altLfile = fopen(filename, "w+");
      fprintf(USFseaLfile, "%20.16f", t1*USFseaL[0][i][0]);
      fprintf(USFseaF0Lfile, "%20.16f", t1*USFseaF0L[0][i][0]);
      fprintf(USFseaF0altLfile, "%20.16f", t1*USFseaF0Lalt[0][i][0]);
      for (j=1; j<=xNIntp; j++) {
        fprintf(USFseaLfile, ",%20.16f", t1*USFseaL[0][i][j]);
        fprintf(USFseaF0Lfile, ",%20.16f", t1*USFseaF0L[0][i][j]);
        fprintf(USFseaF0altLfile, ",%20.16f", t1*USFseaF0Lalt[0][i][j]);
      }
      for (j=1; j<4; j++) {
        fprintf(USFseaLfile, "\n%20.16f", t1*USFseaL[j][i][0]);
        fprintf(USFseaF0Lfile, "\n%20.16f", t1*USFseaF0L[j][i][0]);
        fprintf(USFseaF0altLfile, "\n%20.16f", t1*USFseaF0Lalt[j][i][0]);
        for (k=1; k<=xNIntp; k++) {
          fprintf(USFseaLfile, ",%20.16f", t1*USFseaL[j][i][k]);
          fprintf(USFseaF0Lfile, ",%20.16f", t1*USFseaF0L[j][i][k]);
          fprintf(USFseaF0altLfile, ",%20.16f", t1*USFseaF0Lalt[j][i][k]);
        }
      }
      fclose(USFseaLfile);
      fclose(USFseaF0Lfile);
      fclose(USFseaF0altLfile);
    }
    for (i=0; i<=xNInt; i++) {
      x=(double)i/(double)xNInt*xdist;
      intarr1[i]=x*(USFseaG[0][i]-USFseaF0G[0][i] + USFseaG[2][i]-USFseaF0G[2][i]);
      intarr2[i]=x*(USFseaG[1][i]-USFseaF0G[1][i] + USFseaG[3][i]-USFseaF0G[3][i]);
      intarr3[i]=x*(USFseaG[0][i]-USFseaF0Galt[0][i] + USFseaG[2][i]-USFseaF0Galt[2][i]);
      intarr4[i]=x*(USFseaG[1][i]-USFseaF0Galt[1][i] + USFseaG[3][i]-USFseaF0Galt[3][i]);
    }
    srule0ltot+=0.5*booles(intarr1, xNInt)*xdist/(double)xNInt;
    srule1ltot+=0.5*booles(intarr2, xNInt)*xdist/(double)xNInt;
    srule0lalttot+=0.5*booles(intarr3, xNInt)*xdist/(double)xNInt;
    srule1lalttot+=0.5*booles(intarr4, xNInt)*xdist/(double)xNInt;
    fprintf(stderr, "%20.10f %20.10f %20.16f %20.16f, %d\n\n", srule0ltot, srule0lalttot, srule1ltot, srule1lalttot, xNIntp);
    sprintf(filename, "formatted_data/expnum%d/USFseaG%d_%d.txt", expnum, G, expnum);
    USFseaGfile=fopen(filename, "w+");
    sprintf(filename, "formatted_data/expnum%d/USFseaF0G%d_%d.txt", expnum, G, expnum);
    USFseaF0Gfile=fopen(filename, "w+");
    sprintf(filename, "formatted_data/expnum%d/USFseaF0altG%d_%d.txt", expnum, G, expnum);
    USFseaF0altGfile=fopen(filename, "w+");
    //t1=1.0;
    
    fprintf(USFseaGfile, "%20.16f", t1*USFseaG[0][0]);
    fprintf(USFseaF0Gfile, "%20.16f", t1*USFseaF0G[0][0]);
    fprintf(USFseaF0altGfile, "%20.16f", t1*USFseaF0Galt[0][0]);
    for (k=1; k<=xNIntp; k++) {
      fprintf(USFseaGfile, ",%20.16f", t1*USFseaG[0][k]);
      fprintf(USFseaF0Gfile, ",%20.16f", t1*USFseaF0G[0][k]);
      fprintf(USFseaF0altGfile, ",%20.16f", t1*USFseaF0Galt[0][k]);
    }
    for (j=1; j<4; j++) {
      fprintf(USFseaGfile, "\n%20.16f", t1*USFseaG[j][0]);
      fprintf(USFseaF0Gfile, "\n%20.16f", t1*USFseaF0G[j][0]);
      fprintf(USFseaF0altGfile, "\n%20.16f", t1*USFseaF0Galt[j][0]);
      for (k=1; k<=xNIntp; k++) {
        fprintf(USFseaGfile, ",%20.16f", t1*USFseaG[j][k]);
        fprintf(USFseaF0Gfile, ",%20.16f", t1*USFseaF0G[j][k]);
        fprintf(USFseaF0altGfile, ",%20.16f", t1*USFseaF0Galt[j][k]);
      }
    }
    fclose(USFseaGfile);
    fclose(USFseaF0Gfile);
    fclose(USFseaF0altGfile);
    //freeDataCollectedfouriers(fourier_funcP, fourier_funcN, fourier_funcPF0, fourier_funcNF0, G, momentanums);
    freeDataCollectedvals(Pvals, Nvals, PvalsF0, NvalsF0, mks, ks, pks, G, momentanums);
    free(intarr1);
    free(intarr2);
    free(intarr3);
    free(intarr4);
    //fprintf(stderr, "here\n");
  }
  for (i=0; i<4; i++) {
    free(USFseaG[i]);
    free(USFseaF0G[i]);
    free(USFseaF0Galt[i]);
  }
  free(USFseaG);
  free(USFseaF0G);
  free(USFseaF0Galt);
  intarr1=(double*)malloc(sizeof(double)*(xNInt+1));
  intarr2=(double*)malloc(sizeof(double)*(xNInt+1));
  for (i=0; i<=xNInt; i++) {
    x=(double)i/(double)xNInt*xdist;
    intarr1[i]=x*(USFsea[0][i]-USFseaF0[0][i] + USFsea[2][i]-USFseaF0[2][i]);
    intarr2[i]=x*(USFsea[1][i]-USFseaF0[1][i] + USFsea[3][i]-USFseaF0[3][i]);
  }
  t1=0.5*booles(intarr1, xNInt)*xdist/(double)xNInt;
  t2=0.5*booles(intarr2, xNInt)*xdist/(double)xNInt;
  for (i=0; i<=xNInt; i++) {
    x=(double)i/(double)xNInt*xdist;
    intarr1[i]=x*(USFsea[0][i]-USFseaF0alt[0][i] + USFsea[2][i]-USFseaF0alt[2][i]);
    intarr2[i]=x*(USFsea[1][i]-USFseaF0alt[1][i] + USFsea[3][i]-USFseaF0alt[3][i]);
  }
  t3=0.5*booles(intarr1, xNInt)*xdist/(double)xNInt;
  t4=0.5*booles(intarr2, xNInt)*xdist/(double)xNInt;
  fprintf(stderr, "%f %f %f %f\n", t1, t3, t2, t4);


  free(intarr1);
  free(intarr2);
  t1=(3.0*mass*5.0/144.0);
  for (i=0; i<4; i++) {
    for (j=0; j<=xNInt; j++) {
      USFsea[i][j]*=t1;
      USFseaF0[i][j]*=t1;
      USFseaF0alt[i][j]*=t1;
    }
  }

  sprintf(filename, "formatted_data/expnum%d/USFsea_%d.txt", expnum, expnum);
  USFseafile = fopen(filename, "w+");
  sprintf(filename, "formatted_data/expnum%d/USFseaF0_%d.txt", expnum, expnum);
  USFseaF0file = fopen(filename, "w+");
  sprintf(filename, "formatted_data/expnum%d/USFseaF0alt_%d.txt", expnum, expnum);
  USFseaF0altfile = fopen(filename, "w+");
  fprintf(USFseafile, "%20.16f", USFsea[0][0]);
  fprintf(USFseaF0file, "%20.16f", USFseaF0[0][0]);
  fprintf(USFseaF0altfile, "%20.16f", USFseaF0alt[0][0]);
  for (i=1; i<=xNIntp; i++) {
    fprintf(USFseafile, ",%20.16f", USFsea[0][i]);
    fprintf(USFseaF0file, ",%20.16f", USFseaF0[0][i]);
    fprintf(USFseaF0altfile, ",%20.16f", USFseaF0alt[0][i]);
  }
  for (i=1; i<4; i++) {
    fprintf(USFseafile, "\n%20.16f", USFsea[i][0]);
    fprintf(USFseaF0file, "\n%20.16f", USFseaF0[i][0]);
    fprintf(USFseaF0altfile, "\n%20.16f", USFseaF0alt[i][0]);
    for (j=1; j<=xNIntp; j++) {
      fprintf(USFseafile, ",%20.16f", USFsea[i][j]);
      fprintf(USFseaF0file, ",%20.16f", USFseaF0[i][j]);
      fprintf(USFseaF0altfile, ",%20.16f", USFseaF0alt[i][j]);
    }
  }

  for (i=0; i<4; i++) {
    free(USFsea[i]);
    free(USFseaF0[i]);
    free(USFseaF0alt[i]);
  }
  free(USFsea);
  free(USFseaF0);
  free(USFseaF0alt);

  fclose(USFseafile);
  fclose(USFseaF0file);
  fclose(USFseaF0altfile);

}*/

void ValenceFormat(char printIMF) {
  FILE * valFile;
  char sfile[100];
  double **USFval, *pholder;
  unsigned short i, j, xNIntp;

  sprintf(sfile, "formatted_Data/expnum%d/USFval_%d.txt", expnum, expnum);
  if (printIMF) {
    sprintf(sfile, "formatted_Data/expnum%d/USFvalIMF_%d.txt", expnum, expnum);
  }
  valFile = fopen(sfile, "w");
  xNIntp=floor(graphcut/xdist*xNInt);
  if (printIMF) {
    xNIntp=floor(IMFgraphcut*xNInt);
  }
  USFvaldata(&USFval);
  if (printIMF) {
    for (i=0; i<4; i++) {
      RF_to_IMF(&pholder, USFval[i]);
      free(USFval[i]);
      USFval[i]=pholder;
    }
  }
  for (i=0; i<4; i++) {
    fprintf(valFile, "%20.16f", USFval[i][0]);
    for (j=1; j<=xNIntp; j++) {
      fprintf(valFile, ",%20.16f", USFval[i][j]);
    }
    fprintf(valFile, "\n");
  }
  fclose(valFile);
}
int main(int argc, const char *argv[]) {
  double Esol, *thetas;
  unsigned short grandspinmax, *momentanums, i;

  //altMassFormat();
  //altAverageFormat(0);
  //altAverageFormat(1);
  //ValenceFormat(0);
  //ValenceFormat(1);
  DefaultFormat();
  //EigVals();
  //printmomentsnums();
  //chiral_to_file(200, 3);
}
