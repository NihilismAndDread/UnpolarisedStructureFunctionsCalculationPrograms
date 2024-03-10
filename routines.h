#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#define _USE_MATH_DEFINES
#define BIG_NO 1e10
#define BIG_NI 1e-10
#define THRESH 1E-10

short unsigned MAXIT=1000;
double dJn(double x, unsigned short n);
double uJn(double x, unsigned short n);
double Jn(double x, unsigned short n);
void dJnarr(double *jns, double x, unsigned short index, unsigned short n);
void uJnarr(double *jns, double x, unsigned short n);
void Jnarr(double *jns, double x, unsigned short n);
unsigned short JnZeros(unsigned short n, double xcut, double* prezeros, unsigned short prenum, double** zeros);
double simpsons(double *fs, unsigned short N);
double findJnZero(unsigned short n, double a, double b, double xcut);
double max(double a, double b);
double min(double a, double b);
void updated(unsigned short k, double t, double c, double* val);
int maxind(double** mat, unsigned short k, unsigned short mlen);
void rotate(double** mat, double c, double s, unsigned short k, unsigned short l, unsigned short i, unsigned short j);
void jacobi(double** mat, double* val, double** vec, unsigned short len, unsigned short showerr);
short sgn(double num);
double step_function(double x);
void jrotate1(double**mat, double s, double tau, unsigned short i, unsigned short j, unsigned short k, unsigned short l);
void jacobi1(double** mat, double* val, double** vec, unsigned short len, unsigned short showerr);


double step_function(double x) {
  if (x<=0) {
    return 0;
  }
  return 1;
}
short sgn(double num) {
  if (num>=0) {
    return 1;
  } else {
    return -1;
  }
  return 0;
}
void jacobi(double** mat, double* val, double**vec,unsigned short n, unsigned short showerr) {
  double *b, *z;
	double theta, tresh, tau, s, c, h, t, sum, g;
	unsigned short i, j, k, l, o, p, nrot, index1, index2;
	b = (double*)malloc(sizeof(double)*n);
	z = (double*)malloc(sizeof(double)*n);
	for (i=0; i<n; i++) {
		b[i] = mat[i][i];
		val[i] = b[i];
		z[i]=0.0;
		for (j=0; j<n; j++) {
			if(i==j) {
				vec[i][j] = 1.0;
			} else {
				vec[i][j] = 0.0;
			}
		}
	}
	for (i=0; i<n; i++) {
		for (j=0; j<n; j++) {
			//fprintf(stderr, "%20.6f ", mat[i][j]);
		}
		//fprintf(stderr, "\n");
	}
	nrot=0;
	for (i=0; i<50; i++) {
		sum =0;
		for (j=0; j<n-1; j++) {
			for (k=j+1; k<n; k++) {
				sum += fabs(mat[j][k]);
			}
		}
		//fprintf(stderr, "sum=%f\n", sum);
		if (sum == 0.0) {return;}
		//fprintf(stderr, "here\n");
		if (i<3) {
			tresh=0.2*sum/((double)512*512);
		} else {
			tresh = 0.0;
		}
		//fprintf(stderr, "n=%d tresh=%20.5f\n", n, tresh);
		index1=0;
		index2=0;
		for (j=0; j<n-1; j++) {
			for (k=j+1; k<n; k++) {
				if (fabs(mat[j][k])!=0.0) {
					//fprintf(stderr, "elseif\n");
					index2++;
					/*for (o=0; o<n; o++) {
						for (p=0; p<n; p++) {
							fprintf(stderr, "%2d %2d %20.5f\n", o+1, p+1, mat[o][p]);
						}
					}
					exit(0);*/
					h=val[k]-val[j];
					theta = 0.5*h/mat[j][k];
					t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
					if (theta < 0.0) {t=-t;}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*mat[j][k];
					//fprintf(stderr, "t=%20.5f\n", t);
					z[j]+=-h;
					z[k]+=h;
					val[j]+=-h;
					val[k]+=h;
					mat[j][k]=0.0;
					for (l=0; l<=j-1; l++) {
						jrotate1(mat, s, tau, l, j, l, k);
					}
					for (l=j+1; l<=k-1; l++) {
						jrotate1(mat, s, tau, j, l, l, k);
					}
					for (l=k+1; l<n; l++) {
						jrotate1(mat, s, tau, j, l, k, l);
					}

					for (l=0; l<n; l++) {
						jrotate1(vec, s, tau, l, j, l, k);
					}
					nrot++;
					//exit(0);
				}
			}
		}
		/*for (o=0; o<n; o++) {
			for (p=0; p<n; p++) {
				fprintf(stderr, "%2d %2d %20.5f\n", o+1, p+1, vec[o][p]);
			}
		}
		exit(0);*/
		//fprintf(stderr, "index1=%d index2=%d\n", index1, index2);
		/*sum=0;
		for (j=0; j<n-1; j++) {
			for (k=j+1; k<n; k++) {
				sum += fabs(mat[j][k]);
			}
		}
		fprintf(stderr, "sum=%f\n", sum);*/
		//return;
		for (j=0; j<n; j++) {
			b[j]+=z[j];
			val[j]=b[j];
		}
		for (j=0; j<n; j++) {
			z[j]=0;
		}

	}
	fprintf(stderr, "Max iteration exceeded in diagonalisation program");
}
void jrotate1(double**mat, double s, double tau, unsigned short i, unsigned short j, unsigned short k, unsigned short l) {
	double temp;
	temp = mat[i][j];
	mat[i][j] += -s*(mat[k][l]+mat[i][j]*tau);
	mat[k][l] +=  s*(temp-mat[k][l]*tau);

}
void jacobi1(double** mat, double* val, double** vec, unsigned short n, unsigned short showerr) {
  unsigned short i, j, k, l, m, state, zs;
  double sum, s, c, t, p, y, d, r, Eik, Eil;
  unsigned short ind[n];
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      if (i==j) {
        vec[i][j]=1.0;
      } else {
        vec[i][j]=0.0;
        vec[i][j]=0.0;
      }
    }
  }
  state = n;
  for (k=0; k<n; k++) {
    ind[k]=maxind(mat, k, n);
    val[k]=mat[k][k];
  }
  while (state!=0) {
    m = 0;
    for (k=1; k<n-1; k++) {
      if (fabs(mat[m][ind[m]])<fabs(mat[k][ind[k]])) {
        m=k;
      }
    }
    k=m;
    l=ind[m];
    if (fabs(mat[k][l])==0) {
      state=0;
      sum=0;
      for (i=0; i<n; i++) {
        for (j=i+1; j<n; j++) {
          sum += mat[i][j];
        }
      }
      if (sum!=0) {
        fprintf(stderr, "wtf\n");
      }
      break;
    }
    p=mat[k][l];
    //fprintf(stderr, "p=%f\n", p);
    y=(val[l]-val[k])/2.0;
    d=fabs(y)+sqrt(p*p+y*y);
    r=sqrt(p*p+d*d);
    c=d/r;
    s=p/r;
    t=p*p/d;
    if (showerr) {
      fprintf(stderr, "c=%f, s=%f, t=%f\n", c, s, t);
    }
    if (y<0) {
      s=-s;
      t=-t;
    }
    mat[k][l]=0.0;
    mat[l][k]=0.0;
    updated(k, -t, c, val);
    updated(l, t, c, val);
    for (i=1; i<=k-1; i++) {
      rotate(mat, c, s, i, k, i, l);
    }
    for (i=k+1; i<=l-1; i++) {
      rotate(mat, c, s, k, i, i, l);
    }
    for (i=l+1; i<n; i++) {
      rotate(mat, c, s, k, i, l, i);
    }
    for (i=0; i<n; i++) {
      Eik=vec[i][k];
      Eil=vec[i][l];
      vec[i][k]=c*Eik-s*Eil;
      vec[i][l]=s*Eik+c*Eil;
    }
    //ind[k]=maxind(mat, k, n);
    //ind[l]=maxind(mat, l, n);
    for (i=0; i<n-1; i++) {
      ind[i]=maxind(mat, i, n);
    }
  }

}
void updated(unsigned short k, double t, double c, double* val) {
  double y;
  y=val[k];
  val[k]=y+t;
}
int maxind(double** mat, unsigned short k, unsigned short matlength) {
  unsigned short m, i;
  m=k+1;
  for (i=k+2; i<matlength; i++) {
    if (fabs(mat[k][m])<=fabs(mat[k][i])) {
      m = i;
    }
  }
  return m;
}
void rotate(double** mat, double c, double s, unsigned short k, unsigned short l, unsigned short i, unsigned short j) {
  double Skl, Sij;
  Skl=mat[k][l];
  Sij=mat[i][j];
  mat[k][l]=c*Skl-s*Sij;
  mat[i][j]=s*Skl+c*Sij;
}
double max(double a, double b){
  if (a>b) {
    return a;
  } else {
    return b;
  }
  return b;
}
double min(double a, double b) {
  if (a<b) {
    return a;
  } else {
    return b;
  }
  return b;
}
double findJnZero(unsigned short n, double a, double b, double xcut) {
  double dx, midf, midx, f, zero;
  unsigned short i;
  f=Jn(a, n);
  midf=Jn(b, n);
  if (midf*f >=0) {
    fprintf(stderr, "findJnZero: must not have the same sign\n");
  }
  if (f < 0.0) {
    dx = b-a;
    zero = a;
  } else {
    dx = a-b;
    zero = b;
  }
  for (i=0; i<MAXIT; i++) {
    dx = dx*0.5;
    midx = zero + dx;
    midf = Jn(midx, n);
    if (midf <= 0.0) {
      zero = midx;
    }
    if (fabs(dx) < xcut || midf == 0.0) {
      return zero;
    }
  }
  fprintf(stderr, "findJnZero: max iterations reached\n");
  return 0;
}
void dJnarr(double *jns, double x, unsigned short index, unsigned short n) {
  double a, b, c, j, sum;
  int i, k;
  unsigned short flag, num=40, start;
  j=0;
  b=0;
  c=1;
  start = 2*((n+(unsigned short)sqrt((double)num*n))/2);
  //fprintf(stderr, "x=%f", x);
	if (x==0.0) {
    if (n==0) {
      jns[0]=1.0;
      return;
    } else {
      //jns[0]=1.0;
      for (i=1; i<=n; i++) {
        jns[i]=0.0;
      }
			if (index==0) {
				jns[0]=1.0;
			}
      return;
    }
  }
  for (i=start; i>=1; i=i-1) {
    //fprintf(stderr, "M=%d\n", M);
    a=((double)2.0*i+1.0)*b/x - c;
    c=b;
    b=a;
    if (fabs(b)>BIG_NO) {
      b=b*BIG_NI;
      c=c*BIG_NI;
      for (k=n; k>=index; k--) {
        if (k<i) {
          break;
        }
        jns[k-index]=jns[k-index]*BIG_NI;
      }
    }
    if ((i-1)<=n && (i-1)>=index) {jns[i-1-index]=b;}
  }
  for (i=0; i<=n-index; i++) {
    jns[i] = jns[i]*sin(x)/(x*b);

  }
}
void uJnarr(double *jns, double x, unsigned short n) {
  double temp, temp1, temp2;

  int i;
  if (x==0.0) {
    jns[0] = 1.0;
    for (i=1; i<=n; i++) {
      jns[i]=0.0;
    }
  } else {
    jns[0] = sin(x)/x;
    if (n>=1) {
      jns[1] = sin(x)/(x*x) - cos(x)/x;
    }
    temp1 = sin(x)/x;
    temp = (temp1-cos(x))/x;
    for (i=2; i<=n; i++) {
      temp2=(2*(double)i-1)*temp/x - temp1;
      temp1 = temp;
      temp = temp2;
      jns[i] = temp;
    }
  }
}
void Jnarr(double *jns, double x, unsigned short n) {
  double *jns1, *jns2;
  unsigned short i, index;
  /*jns[0] = uJn(x, 0);
  if (1<=n) {
    jns[1]=uJn(x, 1);
  }*/
  index=0;
  for (i=0; i<=n; i++) {
    if (x > (double)i) {
      index=i;
    }
  }
  jns1 = (double*)malloc(sizeof(double)*(index+1));
  jns2 = (double*)malloc(sizeof(double)*(n-index));
  //fprintf(stderr, "index=%d\n", index);
	uJnarr(jns1, x, index);
  dJnarr(jns2, x, index+1, n);
  for (i=0; i<=n; i++) {
    if (i<=index) {
      jns[i]=jns1[i];
    } else {
      jns[i]=jns2[i-index-1];
    }
  }
  free(jns1);
  free(jns2);
}
double dJn(double x, unsigned short n) {
  double a, b, c, j, sum;
  int i;
  unsigned short iacc=40, start;
  j=0;
  b=0;
  c=1;
  start = 2*((n+(unsigned short)sqrt((double)iacc*n))/2);
  if (x==0.0) {
    if (n==0) {
      return 1;
    } else {
      return 0;
    }
  }

  for (i=start; i>=1; i=i-1) {
    //fprintf(stderr, "M=%d\n", M);
    a=((double)2.0*(double)i+1.0)*b/x - c;
    c=b;
    b=a;
    if (fabs(b)>BIG_NO) {
      b=b*BIG_NI;
      c=c*BIG_NI;
      j=j*BIG_NI;
    }
    if ((i-1)==n) {j=b;}
  }
  return j*sin(x)/(x*b);
}
double uJn(double x, unsigned short n) {
  double temp, temp1, temp2;

  int i;
  if (n==0) {
    if (x==0.0) {
      return 1.0;
    }
    temp = sin(x)/x;
    //fprintf(stderr, "n=%d, temp=%.15Lf\n", n, temp);
    return temp;
  } else {
    if (x==0.0) {
      return 0.0;
    } else {
      temp1 = sin(x)/x;
      temp = (temp1-cos(x))/x;
      for (i=2; i<=n; i++) {
        temp2=(2*(double)i-1)*temp/x - temp1;
        temp1 = temp;
        temp = temp2;
      }
      return temp;

    }
  }

  return -1.0;
}
double Jn(double x, unsigned short n) {
  double x1, x2;
  if (n==0 || n==1) {
    return uJn(x, n);
  }
  /*x1=uJn(x, n);
  x2=dJn(x, n);
  if (fabs(x1)<=fabs(x2)) {
    return x1;
  } else {
    return x2;
  }
  return -1;*/
  if (x>(double)n) {
    return uJn(x, n);
  } else {
    return dJn(x, n);
  }
}
double simps(double *fs, unsigned long N) {
  double ans;
  unsigned long i;
  ans = fs[0] + fs[N];
  for (i=1; i<=N/2; i++) {
    ans += 4.0*fs[2*i-1];
  }
  for (i=1; i<= N/2 -1; i++) {
    ans += 2.0*fs[2*i];
  }
  ans = ans/3.0;
  return ans;
}
double booles(double *fs, unsigned long N) {
  double ans;
  unsigned long i;
  ans = 7.0*fs[0]+7.0*fs[N];
  for (i=1; i<=N-1; i+=2) {
    ans+=32.0*fs[i];
  }
  for (i=2;i<=N-2; i+=4) {
    ans+=12.0*fs[i];
  }
  for (i=4; i<=N-4; i+=4) {
    ans+=14.0*fs[i];
  }
  ans *= 2.0/45.0;
  return ans;
}
//int main(int argc, const char *argv[]) {
//  printf("zero: %f", Jn(0.0, 0));
//}
