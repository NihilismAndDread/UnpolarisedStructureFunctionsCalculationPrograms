#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "routines.h"
#include "constants.h"
#define _USE_MATH_DEFINES
/*int NInt=128;
double chiral_angle=0.5;
double lambda=1.8576485321061220;
double kratio=4.0;
double mpion = 135.0;
double cqm = 400.0;
int rep = 100;
unsigned short Nc=3;*/
unsigned short printtofile=1;
double theta(double x);
double dereg(double x);

double dereg(double x) {
  double temp;
  temp = sqrt(1+1/(x*x));
  return 1-(1+1.5/(x*x))/(temp*temp*temp);
}
int theta_from_file(char sfile[], double thetas[]) {
  FILE *file;
  int check, i;
  if (!(file=fopen(sfile, "r"))) {
    return 0;
  }
  i=0;
  while (check!=EOF) {
    check = fscanf(file, "%lf", &thetas[i]);
    i++;
  }
  fclose(file);
  return 1;
}
double theta(double x) {
  double temp;
  temp = -M_PI*exp(-x/chiral_angle);
  return temp;
}
int main(int argc, const char* argv[]) {
  unsigned long i, j, k, n, m;
  unsigned short kNmax, mkNmax, pkNmax, G, u, index;
  double c1, c2, c3, c4, x1, x2, rn, rm, kcut, x, tval, Nconst;
  double ktemp[250], eigvalsP[500], eigvalsN[500], s[500], p[500], st[500], pt[500], EsolGsep1[500], EsolGsep2[500];
  double intarr[500], thetas[500], newthetas[500];
  double f1p[500], g1p[500], f2p[500], g2p[500];
  double f1m[500], g1m[500], f2m[500], g2m[500];
  double **Hmat1, **Hmat, **eigvecsP, **eigvecsN;
  double *ks, *mks, *pks, *Es, *mEs, *pEs, *kEtemp;
  double **jcc, **jcp, **jcm, **jpc, **jpp, **jmc, **jmm, **jtemp;
  double tracet, trace, Eval, Esol;
  FILE * eigvecsfile, * eigenergiesfile, * momentanumberfile, * chiralanglefile, * knfile, * matcheckm;
  unsigned short jerr=0;
	char filename[100];
  trace=0;
  mpion = mpion/cqm;
  c1 = mpion*mpion*93.0*93.0/(cqm*cqm);
  //READING THETAS FROM FILE
if (usedefault) {
  if (NInt==256) {
    if (theta_from_file("input2.in", thetas) == 0) {
      fprintf(stderr, "No such input file for chiral angle\n");
      return 0;
    }
  }  else if (NInt==128) {
    if (theta_from_file("input.in", thetas) == 0) {
      fprintf(stderr, "No such input file for chiral angle\n");
      return 0;
    }
  } else {
    fprintf(stderr, "No fitting input file for chiral angle, aborting program\n");
    exit(0);
  }
} else {
  for (i=0; i<=NInt; i++) {
    thetas[i]=-M_PI*exp(-10*(double)i/(double)NInt);
  }
}




  Hmat1 = (double**)malloc(500*sizeof(double*));
  Hmat = (double**)malloc(500*sizeof(double*));
  eigvecsP = (double**)malloc(500*sizeof(double*));
  eigvecsN = (double**)malloc(500*sizeof(double*));
  for (i=0; i<500; i++) {
    Hmat[i]=(double*)malloc(500*sizeof(double));
    eigvecsP[i] = (double*)malloc(500*sizeof(double));
    eigvecsN[i] = (double*)malloc(500*sizeof(double));
    Hmat1[i]= (double*)malloc(500*sizeof(double));
  }
  ks = (double*)malloc(250*sizeof(double));
  mks = (double*)malloc(250*sizeof(double));
  pks = (double*)malloc(250*sizeof(double));
  Es = (double*)malloc(250*sizeof(double));
  mEs = (double*)malloc(250*sizeof(double));
  pEs = (double*)malloc(250*sizeof(double));
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
  fprintf(stderr, "Confirmation of Constants\n");
  fprintf(stderr, "kdist=%f\ndist=%f\nkratio=%f\nlambda=%f\ncqm=%f\nrep=%d\nexpnum=%d\n", kdist, dist, kratio, lambda, cqm, rep, expnum);
  for (k=0; k<rep; k++) {
    fprintf(stdout, "Iteration %lu of %d\n", k+1, rep);
		//OPENING BINARY FILES
    if (printtofile==1) {
			sprintf(filename, "Data/expnum%d/eigvecs_%d.b", expnum, expnum);
      eigvecsfile = fopen(filename, "wb");

			sprintf(filename, "Data/expnum%d/eigenergies_%d.b", expnum, expnum);
			eigenergiesfile = fopen(filename, "wb");

			sprintf(filename, "Data/expnum%d/momentanum_%d.b", expnum, expnum);
      momentanumberfile = fopen(filename, "wb");

			sprintf(filename, "Data/expnum%d/chiralangle_%d.b", expnum, expnum);
			chiralanglefile = fopen(filename, "wb");

			sprintf(filename, "Data/expnum%d/kn_%d.b", expnum, expnum);
      knfile = fopen(filename, "wb");
    }
    //PREPARING FOR SOLITON ENERGY CHECK
    if (k!=0) {
      for (i=0; i<=G+1; i++) {
        EsolGsep2[i] = EsolGsep1[i];
      }
    }
  //MOMENTA AND ENERGY FOR G=0 POSITIVE PARITY
    kNmax=0;
    kcut = M_PI/dist;
    while (kcut <= lambda*kratio) {
      kNmax+=1;
      kcut += M_PI/dist;
    }
    for (i=0; i<kNmax; i++) {
      ks[i] = ((double)(kNmax-i))*M_PI/dist;
      ks[2*kNmax-1-i]=ks[i];
      Es[i] = -sqrt(1.0+ks[i]*ks[i]);
      Es[2*kNmax-1-i]=-Es[i];
    }
		fprintf(stderr, "kNmax=%d\n", kNmax);
    if (printtofile==1) {
      fwrite(&kNmax, sizeof(unsigned short), 1, momentanumberfile);
    }

    //CREATING BESSEL FUNCTIONS FOR G=0 POSITIVE PARITY
    for (i=0; i<kNmax; i++) {
      for (j=0; j<=NInt; j++) {
        x=(double)j/(double)NInt;
        jcc[i][j] = Jn(ks[i]*dist*x, 0);
        jpc[i][j] = Jn(ks[i]*dist*x, 1);
        jcc[2*kNmax-1-i][j]=jcc[i][j];
        jpc[2*kNmax-1-i][j]=jpc[i][j];
      }
    }
    //CREATING MATRIX ELEMENTS FOR G=0 POSITIVE PARITY
    for (n=0; n<kNmax*2; n++) {
      rn=ks[n]/(1.0+Es[n]);
      for (m=0; m<=n; m++) {
        rm=ks[m]/(1.0+Es[m]);
        intarr[0]=0.0;
        for (i=1; i<=NInt; i++) {
          x=(double)i/(double)NInt;
          tval=thetas[i];

          c1 = jcc[n][i]*jcc[m][i]-jpc[n][i]*jpc[m][i]*rn*rm;
          c2 = jpc[n][i]*jcc[m][i]*rn+jcc[n][i]*jpc[m][i]*rm;
          intarr[i]=x*x*((cos(tval)-1)*c1+sin(tval)*c2);
        }
        c3 = (Es[n]+1)*(Es[m]+1)/(Es[n]*Es[m]);
        c4 = (ks[n]*dist*ks[m]*dist)*sqrt(c3);


        Hmat[n][m]=c4*simps(intarr, NInt)/(double)NInt;
        //fprintf(stderr, "%f", Hmat[n][m]);
        if (n!=m) {
          Hmat[m][n]=Hmat[n][m];
        } else if (n==m){
          Hmat[m][n]=Hmat[n][m]+Es[n];
        }
      }
      //fprintf(stderr, "\n");

    }
    tracet=0;
    for (n=0; n<kNmax*2; n++) {
      tracet+=Hmat[n][n];
    }
    trace+=tracet;
    //DIAGONALIZING MATRIX FOR G=0 POSITIVE PARITY
    tracet=0;
    jacobi(Hmat, eigvalsP, eigvecsP, kNmax*2, jerr);
		if (printtofile==1) {
      for (i=0; i<kNmax*2; i++) {
        fwrite(&eigvalsP[i], sizeof(double), 1, eigenergiesfile);
        fwrite(&ks[i], sizeof(double), 1, knfile);
        //fprintf(stderr, "eigvalsP[%d]=%f\n", i, eigvalsP[i]);
      }
      for (i=0; i<kNmax*2; i++) {
        for (j=0; j<kNmax*2; j++) {
          fwrite(&eigvecsP[j][i], sizeof(double), 1, eigvecsfile);
        }
      }
    }
  //INITIALIZE DENSITIES
    for (i=0; i<=NInt; i++) {
      s[i]=0;
      p[i]=0;
    }

  //DENSITY CALCULATIONS FOR POSITIVE PARITY G=0

    for (u=0; u<kNmax*2; u++) {
      c1 = -0.5 * (eigvalsP[u]/fabs(eigvalsP[u])) * dereg(eigvalsP[u]/lambda);
      for (i=0; i<=NInt; i++) {
        x=(double)i/(double)NInt;
        g1p[i] =0;
        f1p[i] =0;
        for (n=0; n<kNmax*2; n++) {
          rn = ks[n]/(1.0+Es[n]);
          Nconst = sqrt(1+1/Es[n])/fabs(jpc[n][NInt]);
          g1p[i] += eigvecsP[n][u]*Nconst*jcc[n][i];
          f1p[i] += eigvecsP[n][u]*Nconst*rn*jpc[n][i];
        }
        s[i] += c1*(g1p[i]*g1p[i]-f1p[i]*f1p[i]);
        p[i] += c1*(2*g1p[i]*f1p[i]);
      }
    }
  //FIND VALENCE LEVEL AND COMPUTE THE DENSITY CONTRIBUTION FOR IT

    Eval=Es[2*kNmax-1];
    j=kNmax-1;

    for (i=0; i<2*kNmax; i++) {
      if (fabs(eigvalsP[i])<=fabs(Eval)) {
        Eval = eigvalsP[i];
        j=i;
      }
    }
    if (fabs(Eval)>=1.0) {
      fprintf(stderr, "Valence level not detected\n");
    }
    for (i=0; i<=NInt; i++) {
      x=(double)i/(double)NInt;
      g1p[i] =0;
      f1p[i] =0;
      for (n=0; n<kNmax*2; n++) {
        rn = ks[n]/(1+Es[n]);
        Nconst = sqrt(1+1/Es[n])/fabs(jpc[n][NInt]);
        g1p[i] += eigvecsP[n][j]*Nconst*jcc[n][i];
        f1p[i] += eigvecsP[n][j]*Nconst*rn*jpc[n][i];
      }
      s[i] += (g1p[i]*g1p[i]-f1p[i]*f1p[i]);
      p[i] += (2*g1p[i]*f1p[i]);
    }
  //MOMENTA AND ENERGY FOR NEGATIVE PARITY G=0
    for (i=0; i<kNmax; i++) {
      ktemp[i] = ks[kNmax + i];
    }
    ktemp[kNmax+1]=0.0;
    for (i=0; i<kNmax; i++) {
      x1 = ktemp[i]*dist;
      x2 = max(x1+M_PI, ktemp[i+1]*dist);
      ktemp[i] = findJnZero(1, x1, x2, 1E-16)/dist;
    }
    pkNmax =0;
    for (i=0; i<kNmax; i++) {
      if (ktemp[i] < kratio*lambda) {
        pkNmax++;
      }
    }
    for (i=0; i<pkNmax; i++) {
      pks[i]=ktemp[pkNmax-1-i];
      pks[2*pkNmax-1-i]=ktemp[pkNmax-1-i];
      pEs[i]=-sqrt(pks[i]*pks[i]+1);
      pEs[2*pkNmax-1-i]=-pEs[i];
    }
    if  (printtofile==1) {
      fwrite(&pkNmax, sizeof(unsigned short), 1, momentanumberfile);
    }

    //CREATING BESSEL FUNCTIONS FOR G=0 NEGATIVE PARITY
    for (i=0; i<pkNmax; i++) {
      for (j=0; j<=NInt; j++) {
        x = (double)j/(double)NInt;
        jcp[i][j] = Jn(pks[i]*dist*x, 0);
        jpp[i][j] = Jn(pks[i]*dist*x, 1);
        jcp[2*pkNmax -1-i][j] = jcp[i][j];
        jpp[2*pkNmax -1-i][j] = jpp[i][j];
      }
    }
  //CREATING MATRIX ELEMENTS FOR G=0 NEGATIVE PARITY
    for (n=0; n<2*pkNmax; n++) {
      rn = pks[n]/(1+pEs[n]);
      for (m=0; m<=n; m++) {
        rm = pks[m]/(1+pEs[m]);
        intarr[0]=0.0;
        for (i=1; i<=NInt; i++) {
          x = (double)i/(double)NInt;
          tval = thetas[i];
          c1 = jpp[n][i]*jpp[m][i]-jcp[n][i]*jcp[m][i]*rn*rm;
          c2 = jpp[n][i]*jcp[m][i]*rm+jpp[m][i]*jcp[n][i]*rn;
          intarr[i]=x*x*((cos(tval)-1)*c1-sin(tval)*c2);
        }
        c3 = ((pEs[n]+1)*(pEs[m]+1))/(pEs[n]*pEs[m]);
        c4 = sqrt(c3)/fabs(jcp[n][NInt]*jcp[m][NInt]);
        Hmat[n][m]=c4*simps(intarr, NInt)/(double)NInt;
        if (n==m) {
          Hmat[n][m]+= pEs[n];
        } else {
          Hmat[m][n]=Hmat[n][m];
        }
      }
    }

    tracet=0;
    for (n=0; n<pkNmax*2; n++) {
      tracet+=Hmat[n][n];
    }
    trace+=tracet;

  //DIAGONALIZING MATRIX FOR G=0 NEGATIVE PARITY
    jacobi(Hmat, eigvalsN, eigvecsN, pkNmax*2, jerr);
    if (printtofile==1) {
      for (i=0; i<pkNmax*2; i++) {
        fwrite(&eigvalsN[i], sizeof(double), 1, eigenergiesfile);
        fwrite(&pks[i], sizeof(double), 1, knfile);
      }
      for (i=0; i<pkNmax*2; i++) {
        for (j=0; j<pkNmax*2; j++) {
          fwrite(&eigvecsN[j][i], sizeof(double), 1, eigvecsfile);
        }
      }
    }

  //CONTRIBUTION TO DENSITIES G=0 NEGATIVE PARITY
    for (u=0; u<pkNmax*2; u++) {
      c1 = -0.5 * (eigvalsN[u]/fabs(eigvalsN[u])) * dereg(eigvalsN[u]/lambda);
      for (i=0; i<=NInt; i++) {
        x=(double)i/(double)NInt;
        g1m[i] =0;
        f1m[i] =0;
        for (n=0; n<pkNmax*2; n++) {
          rn = pks[n]/(1+pEs[n]);
          Nconst = sqrt(1+1/pEs[n])/fabs(jcp[n][NInt]);
          g1m[i] += eigvecsN[n][u]*Nconst*jpp[n][i];
          f1m[i] += eigvecsN[n][u]*Nconst*rn*jcp[n][i];
        }
        s[i] += c1*(g1m[i]*g1m[i]-f1m[i]*f1m[i]);
        p[i] += -c1*(2*g1m[i]*f1m[i]);
      }
    }

  //SOLITON ENERGY FOR G=0
    Esol =0;
    //if (Eval > 0) {
    //  Esol = Eval;
    //}
    c1=0;
    c2=0;
    for (i=0; i<kNmax*2; i++) {
      c1 += -1.0/2.0*(fabs(eigvalsP[i])-sqrt(eigvalsP[i]*eigvalsP[i]+lambda*lambda)+ 1.0/2.0*lambda*lambda/sqrt(eigvalsP[i]*eigvalsP[i]+lambda*lambda));
      c1 += 1.0/2.0*(fabs(Es[i])-sqrt(Es[i]*Es[i]+lambda*lambda)+ 1.0/2.0*lambda*lambda/sqrt(Es[i]*Es[i]+lambda*lambda));
    }

    for (i=0; i<pkNmax*2; i++) {
      c1 += -1.0/2.0*(fabs(eigvalsN[i])-sqrt(eigvalsN[i]*eigvalsN[i]+lambda*lambda)+ 1.0/2.0*lambda*lambda/sqrt(eigvalsN[i]*eigvalsN[i]+lambda*lambda));
      c1 += 1.0/2.0*(fabs(pEs[i])-sqrt(pEs[i]*pEs[i]+lambda*lambda)+ 1.0/2.0*lambda*lambda/sqrt(pEs[i]*pEs[i]+lambda*lambda));
    }
    fprintf(stderr, "Esol[0]=%f\n", c1);
    EsolGsep1[0] = c1;
    Esol += c1;
  //GENERAL G
    G=0;

    while (kNmax!=0) {
      G++;
      //CALCULATING MOMENTA FOR GENERAL G
      kEtemp = mks;
      mks = ks;
      ks = pks;
      pks = kEtemp;

      kEtemp =mEs;
      mEs = Es;
      Es = pEs;
      pEs = kEtemp;

      mkNmax=kNmax;
      kNmax=pkNmax;
      for (i=0; i<kNmax; i++) {
        ktemp[i] = ks[kNmax + i];
      }
      ktemp[kNmax]=0.0;
      for (i=0; i<kNmax; i++) {
        x1 = ktemp[i]*dist;
        x2 = max(x1+M_PI, ktemp[i+1]*dist);
        while (Jn(x1, G+1)*Jn(x2, G+1)>=0) {
          x2 +=-0.000001;
          x1 += 0.000001;
          if (x1>x2) {
            fprintf(stderr, "x1 greater than x2\n");
          }
        }
        ktemp[i] = findJnZero(G+1, x1, x2, 1E-16)/dist;
      }
      pkNmax =0;
      for (i=0; i<kNmax; i++) {
        if (ktemp[i] < lambda*kratio) {
          pkNmax++;
        } else {
          break;
        }
      }
      for (i=0; i<pkNmax; i++) {
        pks[i]=ktemp[pkNmax-1-i];
        pks[2*pkNmax-1-i]=ktemp[pkNmax-1-i];
        pEs[i]=-sqrt(pks[i]*pks[i]+1);
        pEs[2*pkNmax-1-i]=-pEs[i];
      }
      if (printtofile==1) {
        fwrite(&pkNmax, sizeof(unsigned short), 1, momentanumberfile);
        for (i=0; i<2*pkNmax; i++) {
          fwrite(&pks[i], sizeof(double), 1, knfile);
        }
      }
      //CALCULATING BESSEL FUNCTIONS FOR GENERAL G
      jtemp = jmm;
      jmm = jcc;
      jcc = jpp;
      jpp = jtemp;

      jtemp = jcm;
      jcm = jpc;
      jpc = jtemp;

      jtemp = jmc;
      jmc = jcp;
      jcp = jtemp;
      jcp[0][NInt]=1;
      for (i=0; i<pkNmax; i++) {
        for (j=0; j<=NInt; j++) {
          x = (double)j/(double)NInt;
          jcp[i][j] = Jn(pks[i]*dist*x, G);
          jpp[i][j] = Jn(pks[i]*dist*x, G+1);
          jcp[2*pkNmax -1-i][j] = jcp[i][j];
          jpp[2*pkNmax -1-i][j] = jpp[i][j];
        }
      }
      for (i=0; i<kNmax; i++) {
        for (j=0; j<=NInt; j++) {
          x = (double)j/(double)NInt;
          jpc[i][j] = Jn(ks[i]*dist*x, G+1);

          jpc[2*kNmax -1-i][j] = jpc[i][j];
        }
      }
      //CALCULATING MATRIX FOR GENERAL G, POSITIVE PARITY
     for (n=0; n<kNmax*2; n++) {
        rn=ks[n]/(1.0+Es[n]);
        for (m=0; m<=n; m++) {
          rm=ks[m]/(1.0+Es[m]);
          intarr[0]=0.0;
          for (i=1; i<=NInt; i++) {
            x=(double)i/(double)NInt;
            tval = thetas[i];
            c1 = jcc[n][i]*jcc[m][i]-jpc[n][i]*jpc[m][i]*rn*rm;
            c2 = jpc[n][i]*jcc[m][i]*rn+jcc[n][i]*jpc[m][i]*rm;
            intarr[i]=x*x*((cos(tval)-1)*c1+sin(tval)*c2/(double)(2.0*G+1.0));
          }
          c3 = (Es[n]+1)*(Es[m]+1)/(Es[n]*Es[m]);
          c4 = sqrt(c3)/fabs(jmc[n][NInt]*jmc[m][NInt]);
          Hmat[n][m]=c4*simps(intarr, NInt)/(double)NInt;
          Hmat[m][n]=Hmat[n][m];

          for (i=1; i<=NInt; i++) {
            x=(double)i/(double)NInt;
            tval = thetas[i];
            c1 = jcc[n][i]*jcc[m][i]-jmc[n][i]*jmc[m][i]*rn*rm;
            c2 = jmc[n][i]*jcc[m][i]*rn+jcc[n][i]*jmc[m][i]*rm;
            intarr[i]=x*x*((cos(tval)-1.0)*c1+sin(tval)*c2/(double)(2.0*G+1.0));
          }
          Hmat[kNmax*2+m][kNmax*2+n]=c4*simps(intarr, NInt)/(double)NInt;
          Hmat[kNmax*2+n][kNmax*2+m]=Hmat[kNmax*2+m][kNmax*2+n];
        }
      }
      for (n=0; n<kNmax*2; n++) {
        Hmat[n][n]+=Es[n];
        Hmat[kNmax*2+n][kNmax*2+n]+=Es[n];
      }
      for (n=0; n<kNmax*2; n++) {
        rn=ks[n]/(1.0+Es[n]);
        for (m=0; m<kNmax*2; m++) {
          rm=ks[m]/(1.0+Es[m]);
          intarr[0]=0.0;
          for (i=1; i<=NInt; i++) {
            x=(double)i/(double)NInt;
            tval = thetas[i];
            c2 = rm*jmc[m][i]*jcc[n][i]-rn*jcc[m][i]*jpc[n][i];
            intarr[i]=x*x*sin(tval)*c2*2.0*sqrt((double)G*(G+1))/(double)(2*G+1);
          }
          c3 = (Es[n]+1)*(Es[m]+1)/(Es[n]*Es[m]);
          c4 = sqrt(c3)/fabs(jpc[n][NInt]*jpc[m][NInt]);
          Hmat[kNmax*2+m][n]=c4*simps(intarr, NInt)/(double)NInt;
          Hmat[n][kNmax*2+m]=Hmat[kNmax*2+m][n];
        }
      }
      tracet=0;
      for (i=0; i<4*kNmax; i++) {
        tracet+=Hmat[i][i];
      }
      trace += tracet;
      //DIAGONALISE MATRIX FOR GENERAL G POSITIVE PARITY

			jacobi(Hmat, eigvalsP, eigvecsP, kNmax*4, jerr);
      for (i=0; i<kNmax*4; i++) {
        tracet += eigvalsP[i];
      }
      if (printtofile==1) {
        for (i=0; i<kNmax*4; i++) {
          fwrite(&eigvalsP[i], sizeof(double), 1, eigenergiesfile);
				}
        for (i=0; i<kNmax*4; i++) {
          for (j=0; j<kNmax*4; j++) {
            fwrite(&eigvecsP[j][i], sizeof(double), 1, eigvecsfile);
          }
        }
      }
      //CALCULATE DENSITIES FOR GENERAL G, POSITIVE PARITY
      for (u=0; u<kNmax*4; u++) {
        c1 = -0.5 * (eigvalsP[u]/fabs(eigvalsP[u])) * dereg(eigvalsP[u]/lambda);
        for (i=0; i<=NInt; i++) {
          x=(double)i/(double)NInt;
          g1p[i] =0;
          f1p[i] =0;
          g2p[i] =0;
          f2p[i] =0;
          for (n=0; n<kNmax*2; n++) {
            rn = ks[n]/(1+Es[n]);
            Nconst = sqrt(1+1/Es[n])/fabs(jpc[n][NInt]);
            g1p[i] += eigvecsP[n][u]*Nconst*jcc[n][i];
            f1p[i] += eigvecsP[n][u]*Nconst*rn*jpc[n][i];
            g2p[i] += eigvecsP[n+kNmax*2][u]*Nconst*jcc[n][i];
            f2p[i] += eigvecsP[n+kNmax*2][u]*rn*Nconst*jmc[n][i];
          }
          s[i] += c1*(2*G+1)*(g1p[i]*g1p[i]-f1p[i]*f1p[i] + g2p[i]*g2p[i] - f2p[i]*f2p[i]);
          p[i] += c1*(2*(g1p[i]*f1p[i] + g2p[i]*f2p[i]) + 4*sqrt(G*(G+1))*(g1p[i]*f2p[i]-g2p[i]*f1p[i]));
        }
      }
      //compute matrix elements for general G negative parity

      for (n=0; n<2*mkNmax; n++) {
        rn=mks[n]/(1.0+mEs[n]);
        for (m=0; m<=n; m++) {
          rm=mks[m]/(1.0+mEs[m]);
          intarr[0]=0;
          for (i=1; i<=NInt; i++) {
            x=(double)i/(double)NInt;
            tval = thetas[i];
            c1=jmm[n][i]*jmm[m][i]-rn*rm*jcm[n][i]*jcm[m][i];
            c2=jmm[n][i]*jcm[m][i]*rm+jcm[n][i]*jmm[m][i]*rn;
            intarr[i]=x*x*((cos(tval)-1)*c1-sin(tval)*c2/(double)(2*G+1));
          }
          c3=(mEs[n]+1)*(mEs[m]+1)/(mEs[n]*mEs[m]);
          c4=sqrt(c3)/fabs(jcm[m][NInt]*jcm[n][NInt]);
          Hmat[m+2*pkNmax][n+2*pkNmax]=c4*simps(intarr, NInt)/(double)NInt;
          if (n==m) {
            Hmat[m+2*pkNmax][n+2*pkNmax]+=mEs[n];
          } else {
            Hmat[n+2*pkNmax][m+2*pkNmax]=Hmat[m+2*pkNmax][n+2*pkNmax];
          }
        }
      }

      for (n=0; n<2*pkNmax; n++) {
        rn = pks[n]/(1.0+pEs[n]);
        for (m=0; m<=n; m++) {
          rm = pks[m]/(1.0+pEs[m]);
          for (i=1;i<=NInt; i++) {
            x=(double)i/(double)NInt;
            tval = thetas[i];
            c1=jpp[n][i]*jpp[m][i]-rn*rm*jcp[n][i]*jcp[m][i];
            c2=jpp[n][i]*jcp[m][i]*rm+jcp[n][i]*jpp[m][i]*rn;
            intarr[i]=x*x*((cos(thetas[i])-1)*c1-sin(thetas[i])*c2/(double)(2*G+1));
          }
          c3=(pEs[n]+1)*(pEs[m]+1)/(pEs[n]*pEs[m]);
          c4=sqrt(c3)/fabs(jcp[m][NInt]*jcp[n][NInt]);
          Hmat[m][n]=c4*simps(intarr, NInt)/(double)NInt;
          if (n==m) {
            Hmat[m][n]+=pEs[n];
          } else {
            Hmat[n][m]=Hmat[m][n];
          }
        }
      }
      for (n=0; n<2*mkNmax;n++) {
        rn = mks[n]/(1.0+mEs[n]);
        for (m=0; m<2*pkNmax; m++) {
          rm = pks[m]/(1.0+pEs[m]);
          for (i=1; i<=NInt; i++) {
            x=(double)i/(double)(NInt);
            tval = thetas[i];
            c2=rm*jmm[n][i]*jcp[m][i]-rn*jcm[n][i]*jpp[m][i];
            intarr[i]=x*x*(sin(tval)*c2)*2.0*sqrt((double)G*(G+1))/(double)(2*G+1);
          }
          c3=(1.0+pEs[m])*(1.0+mEs[n])/(mEs[n]*pEs[m]);
          c4=sqrt(c3)/fabs(jcm[n][NInt]*jcp[m][NInt]);
          Hmat[m][n+2*pkNmax]=c4*simps(intarr, NInt)/(double)NInt;
          Hmat[n+2*pkNmax][m]=Hmat[m][n+2*pkNmax];
        }
      }
      tracet=0;
      for (i=0; i< 2*mkNmax +2*pkNmax; i++) {
        tracet+=Hmat[i][i];
      }
      trace+=tracet;
      jacobi(Hmat, eigvalsN, eigvecsN, 2*mkNmax+2*pkNmax, jerr);
      if (printtofile==1) {
        for (i=0; i<2*mkNmax +2*pkNmax; i++) {
          fwrite(&eigvalsN[i], sizeof(double), 1, eigenergiesfile);
        }
        for (i=0; i<2*mkNmax +2*pkNmax; i++) {
          for (j=0; j<2*mkNmax +2*pkNmax; j++) {
            fwrite(&eigvecsN[j][i], sizeof(double), 1, eigvecsfile);
          }
        }
      }
      //CALCULATE DENSITIES FOR GENERAL G, NEGATIVE PARITY
      for (u=0; u<(pkNmax*2+mkNmax*2); u++) {
        c1 = -0.5 * (eigvalsN[u]/fabs(eigvalsN[u])) * dereg(eigvalsN[u]/lambda);
        for (i=0; i<=NInt; i++) {
          x=(double)i/(double)NInt;
          g1m[i] =0;
          f1m[i] =0;
          g2m[i] =0;
          f2m[i] =0;
          for (n=0; n<pkNmax*2; n++) {
            rn = pks[n]/(1+pEs[n]);
            Nconst = sqrt(1+1/pEs[n])/fabs(jcp[n][NInt]);
            g1m[i] += eigvecsN[n][u]*Nconst*jpp[n][i];
            f1m[i] += eigvecsN[n][u]*Nconst*rn*jcp[n][i];
          }
          for (n=0; n<mkNmax*2; n++) {
            rn = mks[n]/(1+mEs[n]);
            Nconst = sqrt(1+1/mEs[n])/fabs(jcm[n][NInt]);
            g2m[i] += eigvecsN[pkNmax*2+n][u]*Nconst*jmm[n][i];
            f2m[i] += eigvecsN[pkNmax*2+n][u]*Nconst*rn*jcm[n][i];
          }
          s[i] += c1*(2*G+1)*(g1m[i]*g1m[i]-f1m[i]*f1m[i] + g2m[i]*g2m[i] - f2m[i]*f2m[i]);
          p[i] += c1*(-2*(g1m[i]*f1m[i] + g2m[i]*f2m[i]) + 4*sqrt(G*(G+1))*(g2m[i]*f1m[i]-g1m[i]*f2m[i]));
        }
      }
      //CALCLATING SOLITON MASS CONTRIBUTION FROM GENERAL G
      c1=0;
      for (i=0; i<kNmax*4; i++) {
        c1 += -1.0/2.0*(G*2+1)*(fabs(eigvalsP[i])-sqrt(eigvalsP[i]*eigvalsP[i]+lambda*lambda)+ 1.0/2.0*lambda*lambda/sqrt(eigvalsP[i]*eigvalsP[i]+lambda*lambda));
        c1 += 1.0/2.0*(G*2+1)*(fabs(Es[i%(2*kNmax)])-sqrt(Es[i%(2*kNmax)]*Es[i%(2*kNmax)]+lambda*lambda)+ 1.0/2.0*lambda*lambda/sqrt(Es[i%(2*kNmax)]*Es[i%(2*kNmax)]+lambda*lambda));
      }
      for (i=0; i<pkNmax*2+mkNmax*2; i++) {
        c1 += -1.0/2.0*(G*2+1)*(fabs(eigvalsN[i])-sqrt(eigvalsN[i]*eigvalsN[i]+lambda*lambda)+ 1.0/2.0*lambda*lambda/sqrt(eigvalsN[i]*eigvalsN[i]+lambda*lambda));
      }
      for (i=0; i<pkNmax*2; i++) {
        c1 += 1.0/2.0*(G*2+1)*(fabs(pEs[i])-sqrt(pEs[i]*pEs[i]+lambda*lambda)+ 1.0/2.0*lambda*lambda/sqrt(pEs[i]*pEs[i]+lambda*lambda));
      }
      for (i=0; i<mkNmax*2; i++) {
        c1 += 1.0/2.0*(G*2+1)*(fabs(mEs[i])-sqrt(mEs[i]*mEs[i]+lambda*lambda)+ 1.0/2.0*lambda*lambda/sqrt(mEs[i]*mEs[i]+lambda*lambda));
      }
      Esol += c1;
      EsolGsep1[G] = c1;
    }
    //Esol = (double)Nc*Esol;
    printf("highest grand spin = %d, total trace: %f\n", G, trace);

    //CALCULATING NEW THETA
    c1 = mpion*mpion*93.0*93.0/(cqm*cqm);
    c2 =0;
    for (i=0; i<=NInt; i++) {
      s[i] = s[i]/(dist*dist*dist);
      s[i] = s[i] - (4.0*M_PI/3.0)*c1;
      p[i] = p[i]/(dist*dist*dist);
      c3 =thetas[i];
      if (p[i]<0) {
        thetas[i] = acos(-s[i]/sqrt(s[i]*s[i]+p[i]*p[i]));
      } else if (p[i]>=0) {
        thetas[i] = -acos(-s[i]/sqrt(s[i]*s[i]+p[i]*p[i]));
      }
      c3 = c3-thetas[i];
      if (c2<fabs(c3)) {
        c2 = c3;
      }

    }
    trace=0;
    //CALCULATING SOLITON ENERGY
    for (i=0; i<=NInt; i++) {
      x= i/(double)NInt;
      intarr[i]=x*x*(1.0-cos(thetas[i]));
    }
    //Esol *= cqm;
    fprintf(stdout, "E0 = %f\n", Esol*cqm*Nc);
    fprintf(stdout, "EV+E0 = %f\n", Esol*3*cqm+Eval*3*cqm);
    c1 = 4*M_PI*mpion*mpion*93.0*93.0/(cqm*cqm)*dist*dist*dist;
    c2 = simps(intarr, NInt)/(double)NInt;
    fprintf(stdout, "EMESONIC = %f\n", cqm*c2*c1);
    //Esol += c1*c2;
    fprintf(stdout, "Total energy: %f\n", Esol*3*cqm+Eval*3*cqm+cqm*c2*c1);
    fprintf(stdout, "Eval=%f\n", Eval*3*cqm);

    fprintf(stdout, "New thetas:");
    for (i=0; i<=NInt; i++) {
      fprintf(stdout, " %f", thetas[i]);
      if (printtofile==1) {
        fwrite(&thetas[i], sizeof(double), 1, chiralanglefile);
      }
		}
    fprintf(stdout, "\n");
    if (printtofile==1) {
      fclose(eigvecsfile);
      fclose(eigenergiesfile);
      fclose(momentanumberfile);
      fclose(chiralanglefile);
      fclose(knfile);
    }
    if (k!=0) {
      c1=0;
      index=0;
      for (i=0; i<=G; i++) {
        //fprintf(stderr, "EsolGsep1[%d]=%f\n", i, EsolGsep1[i]);
        if (c1<=fabs((EsolGsep1[i]-EsolGsep2[i])/EsolGsep1[i])) {
          c1=fabs((EsolGsep1[i]-EsolGsep2[i])/EsolGsep1[i]);
          c2=i;
        }
      }
      //fprintf(stdout, "Esol Max diff on G spin %d, %f\n", index, c1);
    }
    if (0) {
      fprintf(stdout, "Desired accuracy achieved.\n");
      for (i=0; i<250; i++) {
        free(Hmat[i]);
        free(eigvecsP[i]);
        free(eigvecsN[i]);
      }
      free(Hmat);
      free(eigvecsP);
      free(eigvecsN);

      for (i=0; i<125; i++) {
        free(jcc[i]);
        free(jcp[i]);
        free(jcm[i]);
        free(jmc[i]);
        free(jmm[i]);
        free(jpc[i]);
        free(jpp[i]);
      }
      free(ks);
      free(mks);
      free(pks);
      free(Es);
      free(mEs);
      free(pEs);
      free(jcc);
      free(jcp);
      free(jcm);
      free(jmc);
      free(jmm);
      free(jpc);
      free(jpp);
      return 0;
    }
  }
  fprintf(stdout, "Maximum number of iterations reached, may not have the desired accuracy.\n");
  for (i=0; i<250; i++) {
    free(Hmat[i]);
    free(eigvecsP[i]);
    free(eigvecsN[i]);
  }
  free(Hmat);
  free(eigvecsP);
  free(eigvecsN);

  for (i=0; i<125; i++) {
    free(jcc[i]);
    free(jcp[i]);
    free(jcm[i]);
    free(jmc[i]);
    free(jmm[i]);
    free(jpc[i]);
    free(jpp[i]);
  }
  free(ks);
  free(mks);
  free(pks);
  free(Es);
  free(mEs);
  free(pEs);
  free(jcc);
  free(jcp);
  free(jcm);
  free(jmc);
  free(jmm);
  free(jpc);
  free(jpp);
}
