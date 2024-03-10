unsigned short NInt=128;
double chiral_angle=0.5;
//double lambda=1.8576485321061220;//old lambda
double lambda=1.8477410815682611;//new lambda
//double lambda=1.6176198091949803;//lambda m=450
//double lambda=1.9;
//double lambda=1.0755495809106856;//m=700
char usedefault=0;
char renorm=1;
double kratio=8.5;
double mpion = 135.0;
double cqm = 400.0;
double alpha = 12.0;
unsigned short rep = 25;
unsigned short Nc=3;
double dist=10.0;
unsigned long kNInt=33000;
double IMFgraphcut=0.99;
//double kdist=1.5*4.0*1.8576485321061220; //old lambda
double kdist=8.5*1.5*1.8477410815682611;// new lambda
//double kdist=5.0*1.5*4.0*1.6176198091949803;//lambda m=450
//double kdist=4.0*1.5*1.9;
//double kdist=2.5*4.0*1.0755495809106856;
unsigned short xNInt=3500;
double xdist=4.0;
unsigned short expnum=118;
double graphcut=2.0;
unsigned short massnum=1;
//Expnum =4 NInt=128, dist=10, new lambda, cqm=400, rep=20 (1 kNInt=33000 ratio=3.0)
//Expnum =5 NInt=256, dist=12, new lambda, cqm=400, rep=30 (1 kNInt=33000 rep=25)
//Expnum =6 NInt=256, dist=15, new lambda, cqm=400, rep=30 (1 kNInt=33000 rep=25 2.5)
//Expnum =7 NInt=128, dist=10, new lambda, cqm=400, rep=30 (1 kNInt=33000 ratio=3.0)

//Expnum =8 NInt=128, dist=10, new lambda, cqm=400, rep=30, kratio=4
//Expnum =9 NInt=128, dist=10, new lambda, cqm=400, rep=30, kratio=6
//Expnum =10 NInt=128, dist=10, old lambda, cqm=400, rep=20, kratio=4 (1 kNInt=33000 ratio=3.0)

//Expnum =11 NInt=128, dist=10, old lambda, cqm=450, rep=20, kratio=4
//Expnum =12 NInt=128, dist=10, new lambda, cqm=450, rep=50, kratio=6
//Expnum =13 NInt=256, dist=12, old lambda, cqm=400, rep=50, kratio=4
//Expnum =14 NInt=256, dist=12, new lambda, cqm=450, rep=50, kratio=6
//Expnum =15 NInt=256, dist=12, old lambda, cqm=450, rep=50, kratio=4

//Expnum =16 NInt=128, dist=10, new lambda, cqm=500, rep=50, kratio=5
//Expnum =17 NInt=256, dist=12, new lambda, cqm=500, rep=50, kratio=5
//Expnum =18 NInt=128, dist=10, new lambda, cqm=350, rep=50, kratio=5
//Expnum =19 NInt 256, dist=12, new lambda, cqm=350, rep=50, kratio=5

//Expnum =20 NInt=128, dist=10, new lambda, cqm=400, rep=20, kratio=6
//Expnum =21 NInt=256, dist=12, new lambda, cqm=400, rep=50, kratio=6

//Expnum =22 NInt=128, dist=8, new lambda, cqm=500, rep=50, kratio=5 

//Expnum =23 NInt=128, dist=10, new lambda, cqm=400, rep=35, kratio=6
//Expnum =24 NInt=128, dist=10, new lambda, cqm=400, rep=30, kratio=6
//Expnum =25 NInt=128, dist=10, new lambda, cqm=400, rep=40, kratio=6
//Expnum =26 NInt=128, dist=10, new lambda, cqm=400, rep=45, kratio=6
//Expnum =27 NInt=128, dist=10, new lambda, cqm=400, rep=25, kratio=6
//Expnum =28 NInt=128, dist=10, new lambda, cqm=400, rep=55, kratio=6
//Expnum =29 NInt=128, dist=10, new lambda, cqm=400, rep=60, kratio=6

//Expnum =30 NInt=128, dist=10, new lambda, cqm=400, rep=25, kratio=6, (kNInt=33000)
//Expnum =31 NInt=128, dist=8, new lambda, cqm=400, rep=50, kratio=6 (rep=25 kNInt=33000 2.5 kratio=4 r1)
//expnum =32 copy of 21
//Expnum =33 NInt=128, dist=8, new lambda, cqm=400, rep=50, kratio=6
//Expnum =34 NInt=128, dist=10, new lambda, cqm=400, rep=50, kratio=6, kNInt=19000, xNInt=3500
//Expnum =35 NInt=128, dist=10, new lambda, cqm=400, rep=50, kratio=6, kNInt=22000, xNInt=3000
//Expnum =36 NInt=128, dist=10, new lambda, cqm=400, rep=50, kratio=6
//Expnum =37 NInt=128, dist=10, lambda =1.9, cqm=400, rep=25, kratio=6
//Expnum =38 NInt=128, dist=10, lambda =1.9, cqm=400, rep=25, kratio=4
//Expnum =39 NInt=128, dist=10, lambda =1.8, cqm=400, rep=25, kratio=6
//Expnum =40 NInt=128, dist=10, lambda =1.8, cqm=400, rep=25, kratio=4
//Expnum =41 NInt=256, dist=12, new lambda, cqm=400, rep=50, kratio=4
//Expnum =42 NInt=128, dist=10, lambda =1.8, cqm=400, rep=20, kratio=4
//Expnum =43 NInt=128, dist=10, lambda =1.9, cqm=400, rep=20, kratio=4
//Expnum =44 NInt=128, dist=10, new lambda, cqm=400, rep=20, kratio=4 (1 kNInt=132000 3.0*1.5)
//Expnum =45 NInt=128, dist=10, new lambda, cqm=450, rep=20, kratio=4
//Expnum =46 NInt=128, dist=10, new lambda, cqm=400, rep=20, kratio=4 (0 kNInt=66000 2.0*1.5) (1 kNInt=55000 3.0) (2 kNInt=33000 7.0 r1)
//Expnum =47 NInt=128, dist=10, new lambda, cqm=400, rep=20, kratio=4 (0 kNInt=33000 2.0) (1 33000 1.5) (2 33000 2.5) (3 33000 2.5 7000) (4 33000 2.5 7000 24) (5 22000 2.5 r1) (6 kNInt=33000 7.0)
//Expnum =48 NInt=128, dist=10, new lambda, cqm=400, rep=20, kratio=4 (0 kNInt=44000 2.0) (1 kNInt=44000 2.5)
//Expnum =49 NInt=128, dist=10, lambda =1.8, cqm=400, rep=20, kratio=4 (0 kNInt=22000 2.5)
//Expnum =50 NInt=128, dist=10, good lambda, cqm=700, rep=20, kratio=4 (0 kNInt=22000 2.5)
//Expnum =51 NInt=128, dist=10, lambda =1.8, cqm=700, rep=20, kratio=4 (0 kNInt=22000 2.5)
//Expnum =52 NInt=128, dist=15, lambda =3.0, cqm=400, rep=30, kratio=4 (0 kNInt=22000 1.5)
//Expnum =53 NInt=128, dist=15, lambda =4.0, cqm=400, rep=30, kratio=4 (0 kNInt=22000 1.5)
//Expnum =54 NInt=128, dist=15, lambda =1.0, cqm=400, rep=30, kratio=4 (0 kNInt=22000 1.5)
//Expnum =55 NInt=128, dist=6, new lambda, cqm=400, rep=25, kratio=4 (0 kNInt=22000 1.5)
//Expnum =56 NInt=128, dist=10, new lambda, cqm=400, rep=20, kratio=4 (0 kNInt=22000 2.5 r1) 
//Expnum =57 NInt=128, dist=8, old lambda, cqm=400, rep=25, kratio=4 (0 kNInt=33000 2.5 r1)
//Expnum =58 NInt=128, dist=6, old lambda, cqm=400, rep=25, kratio=4 (0 kNInt=33000 2.5 r1)
//Expnum =59 NInt=128, dist=8, old lambda, cqm=400, rep=25, kratio=4 (0 kNInt=33000 2.5)
//Expnum =60 NInt=128, dist=6, old lambda, cqm=400, rep=25, kratio=4 (0 kNInt=33000 2.5)
//Expnum =61 NInt=128, dist=8, lambda=1.8, cqm=400, rep=25, kratio=4 (0 kNInt=33000 2.5 r1)
//Expnum =62 NInt=128, dist=6, lambda=1.8, cqm=400, rep=25, kratio=4 (0 kNInt=33000 2.5 r1)
//Expnum =63 NInt=128, dist=8, lambda=1.8, cqm=400, rep=25, kratio=4 (0 kNInt=33000 2.5)
//Expnum =64 NInt=128, dist=6, lambda=1.8, cqm=400, rep=25, kratio=4 (0 kNInt=33000 2.5)
//Expnum =65 NInt=256, dist=15, new lambda, cqm=400, rep=25, kratio=2.5 (0 kNInt=27000 7.0 r1)
//Expnum =66 NInt=256, dist=15, new lambda, cqm=400, rep=25, kratio=2.5 (0 kNInt=27000 7.0)
//Expnum =67 NInt=256, dist=20, new lambda, cqm=400, rep=35, kratio=2.5 (0 kNInt=27000 10.0 r1)
//Expnum =68 NInt=256, dist=20, new lambda, cqm=400, rep=35, kratio=2.5 (0 kNInt=27000 10.0)
//Expnum =69 NInt=256, dist=15, new lambda, cqm=400, rep=25, kratio=4.0 (0 kNInt=33000 2.5 r1)
//Expnum =70 NInt=128, dist=10, new lambda, cqm=400, rep=25, kratio=2.5 (0 kNInt=33000 2.5 r1)
//Expnum =71 NInt=128, dist=10, new lambda, cqm=400, rep=25, kratio=2.5 (0 kNInt=33000 2.5)
//Expnum =72 NInt=128, dist=10, old lambda, cqm=400, rep=25, kratio=2.5 (0 kNInt=33000 2.5 r1)
//Expnum =73 NInt=128, dist=10, old lambda, cqm=400, rep=25, kratio=2.5 (0 kNInt=33000 2.5)
//Expnum =74 NInt=128, dist=10, lambda =1.8, cqm=400, rep=25, kratio=2.5 (0 kNInt=33000 2.5 r1)
//Expnum =75 NInt=128, dist=10, lambda =1.8, cqm=400, rep=25, kratio=2.5 (0 kNInt=33000 2.5)
//Expnum =76 NInt=128, dist=10, new lambda, cqm=400, rep=25, kratio=3.0 (0 kNInt=22000)
//Expnum =77 NInt=128, dist=10, new lambda, cqm=400, rep=25, kratio=2.5 (0 kNInt=22000)
//Expnum =78 NInt=128, dist=12, new lambda, cqm=400, rep=25, kratio=4.0 (0 kNInt=22000)
//Expnum =79 NInt=128, dist=8, new lambda, cqm=400, rep=25, kratio=4.0 (0 kNInt=22000)
//Expnum =80 NInt=128, dist=10, new lambda, cqm=400, rep=25, kratio=4.0 ALT (0 kNInt=22000 alpha=1.0)
//Expnum =81 NInt=128, dist=10, new lambda, cqm=400, rep=25, kratio=4.0 ALT (0 kNInt=22000 alpha=5.0)
//Expnum =82 NInt=128, dist=10, new lambda, cqm=400, rep=25, kratio=4.0 ALT (0 kNInt=22000 alpha=8.0)
//Expnum =83 NInt=128, dist=10, new lambda, cqm=400, rep=25, kratio=4.0 ALT (0 kNInt=22000 alpha=12.0)
//Expnum =84 NInt=128, dist=10, new lambda, cqm=400, rep=25, kratio=4.0 ALT (0 kNInt=22000 alpha=1.0 xdist=4)
//Expnum =85 NInt=128, dist=10, new lambda, cqm=400, rep=25, kratio=4.0 ALT (0 kNInt=22000 alpha=5.0 xdist=4)
//Expnum =86 NInt=128, dist=10, new lambda, cqm=400, rep=25, kratio=4.0 ALT (0 kNInt=22000 alpha=8.0 xdist=4)
//Expnum =87 NInt=128, dist=10, new lambda, cqm=400, rep=25, kratio=4.0 ALT (0 kNInt=22000 alpha=12.0 xdist=4)
//Expnum =88 NInt=128, dist=10, new lambda, cqm=400, rep=25, kratio=4.0 OLD (0 kNInt=22000 xdist=4)
//Expnum =89 NInt=128, dist=10, old lambda, cqm=400, rep=25, kratio=4.0 OLD (0 kNInt=22000 xdist=4)
//Expnum =90 NInt=128, dist=10, lambda=1.9, cqm=400, rep=25, kratio=4.0 OLD (0 kNInt=22000 xdist=4)
//Expnum =91 NInt=128, dist=10, lambda=1.8, cqm=400, rep=25, kratio=4.0 OLD (0 kNInt=22000 xdist=4)
//Expnum =92 NInt=128, dist=12, new lambda, cqm=400, rep=35, kratio=4.0 OLD (0 kNInt=22000 xdist=4)
//Expnum =93 NInt=128, dist=15, new lambda, cqm=400, rep=45, kratio=4.0 OLD (0 kNInt=22000 xdist=4)
//Expnum =94 NInt=128, dist=10, new lambda, cqm=400, rep=25, kratio=4.0 OLD (0 kNInt=22000 xdist=4 xNInt=5000)
//Expnum =95 NInt=128, dist=10, old lambda, cqm=400, rep=25, kratio=4.0 OLD (0 kNInt=22000 xdist=4 xNInt=5000)
//Expnum =96 NInt=128, dist=10, lambda=1.9, cqm=400, rep=25, kratio=4.0 OLD (0 kNInt=22000 xdist=4 xNInt=5000)
//Expnum =97 NInt=128, dist=10, lambda=1.8, cqm=400, rep=25, kratio=4.0 OLD (0 kNInt=22000 xdist=4 xNInt=5000)
//Expnum =98 NInt=128, dist=10, new lambda, cqm=400, rep=25, kratio=4.0 OLD (0 kNInt=22000 xdist=12 r1)
//Expnum =99 NInt=128, dist=10, new lambda, cqm=400, rep=25, kratio=4.0 OLD (0 kNInt=22000 xdist=5 xNInt=5000)
//Expnum =100NInt=128, dist=10, old lambda, cqm=400, rep=25, kratio=4.0 OLD (0 kNInt=22000 xdist=5 xNInt=5000)
//Expnum =101NInt=128, dist=10, lambda=1.8, cqm=400, rep=25, kratio=4.0 OLD (0 kNInt=22000 xdist=5 xNInt=5000)
//Expnum =102NInt=128, dist=10, lambda=1.9, cqm=400, rep=25, kratio=4.0 OLD (0 kNInt=22000 xdist=5 xNInt=5000)
//Expnum =103NInt=128, dist=10, new lambda, cqm=400, rep=25, kratio=4.0 OLD (0 kNInt=22000 xdist=7 xNInt=7000)
//Expnum =104NInt=128, dist=10, old lambda, cqm=400, rep=25, kratio=4.0 OLD (0 kNInt=22000 xdist=7 xNInt=7000)
//Expnum =105NInt=128, dist=10, lambda=1.8, cqm=400, rep=25, kratio=4.0 OLD (0 kNInt=22000 xdist=7 xNInt=7000)
//Expnum =106NInt=128, dist=10, lambda=1.9, cqm=400, rep=25, kratio=4.0 OLD (0 kNInt=22000 xdist=7 xNInt=7000)
//Expnum =107NInt=128, dist=10, new lambda, cqm=400, rep=25, kratio=4.0 OLD (0 kNInt=22000 xdist=5 xNInt=10000)
//Expnum =108NInt=128, dist=10, new lambda, cqm=400, rep=25, kratio=4.0 OLD (0 kNINt=22000 xdist=5 xNInt=5000)
//Expnum =109NInt=128, dist=8 , new lambda, cqm=400, rep=25, kratio=4.0 OLD (0 kNINt=22000 xdist=4 xNInt=5000)
//Expnum =110NInt=128, dist=10, new lambda, cqm=400, rep=25, kratio=4.0 OLD (0 kNInt=22000 xdist=5 xNInt=3500)
//Expnum =111NInt=128, dist=10, new lambda, cqm=400, rep=25, kratio=5.0 OLD (0 kNInt=22000 xdist=4 xNInt=3500)
//Expnum =112NInt=128, dist=10, new lambda, cqm=400, rep=25, kratio=3.0 OLD (0 kNInt=22000 xdist=4 xNInt=3500)
//Expnum =113NInt=128, dist=10, new lambda, cqm=400, rep=25, kratio=7.0 OLD (0 kNInt=15000 xdist=4 xNInt=5000)
//Expnum =114NInt=128, dist=10, new lambda, cqm=400, rep=25, kratio=7.0 OLD (0 kNInt=15000 xdist=3 xNInt=3500)
//Expnum =115NInt=128, dist=10, new lambad, cqm=400, rep=25, kratio=4.0 OLD (0 kNInt=22000 xdist=18xNInt=5000)
//Expnum =116NInt=128, dist=10, new lambda, cqm=400, rep=25, kratio=10.0OLD (0 kNInt=40000 xdist=4 xNInt=3500 r1)
//Expnum =117NInt=128, dist=10, new lambda, cqm=400, rep=25, kratio=10.0OLD (0 kNInt=15000 xdist=4 xNInt=5000 2.2)
//Expnum =118NInt=128, dist=10, new lambda, cqm=400, rep=25, kratio=8.5 OLD (0 kNInt=33000 xdist=4 xNInt=3500)
