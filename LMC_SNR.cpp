/// In the name of GOD: 01/06/99 elahi be omid e to faghat khodet 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <sys/time.h>
#include <sys/timeb.h>
#include "VBBinaryLensingLibrary.h"
using namespace std;


const int Num=10000;
const double MaxD=65.0;///kpc
const double RA=180.0/M_PI;
const double pi= M_PI; 
const double step=MaxD/(double)Num/1.0;///step in kpc
const double Hp= 6.62607004*pow(10.0,-34); 
const double KP=3.08568025*pow(10.,19); // in meter.
const double G= 6.67384*pow(10.,-11.0);// in [m^3/s^2*kg].
const double velocity= 3.0*pow(10.0,8.0);//velosity of light
const double Msun=1.98892*pow(10.,30); //in [kg].
const double Rsun= 6.957*pow(10.0,8.0); ///solar radius [meter]
const double Mjupiter=1.898*pow(10.0,27.0); 
const double Rjupiter= 71.492*pow(10.,6);//radius of jupiter in [m]
const double Mearth=  5.9722*pow(10.0,24.0); 
const double Rearth= 6.3781*pow(10.0,6.0);//meter
const double logT_sun=log10(5778.0);
const double vro_sun=226.0;
const double AU=1.4960*pow(10.0,11.0);
const double year=364.0; 
const double binary_fraction=double(2.0/3.0);
const double df=5.0;
const double Rojupiter= 1326.0 ; //in [kg/m^3]
const double Roearth= 5514.0;// [kg/m^3]
const double epss=0.00000000000038654723954932645; 

const int Nm= 40;//// number of mass range
const int nte=40;//// number of tE range 
const int nrd=10000; 
///============================ Besancon constant ==============///
const double Dsun=8.0;
const double rho0[8]={4.0,7.9,6.2,4.0,5.8,4.9,6.6,3.96};///considering WD
const double d0[8]=  {0.073117,0.0216524,0.0217405,0.0217901,0.0218061,0.0218118,0.0218121,0.0218121};
const double epci[8]={0.014,0.0268,0.0375,0.0551,0.0696,0.0785,0.0791,0.0791};
const double corr[8]={1.0,7.9/4.48419, 6.2/3.52112, 4.0/2.27237, 5.8/3.29525, 4.9/2.78402, 6.6/3.74991, 3.96/2.24994};
const double Rv[4]={3.1,2.5,3.1,3.1};

///============================ WFIRST ===================
const int M=5;//  V  I   K   H   W149
const double FWHM[M]={0.389, 1.0, 1.6, 0.4, 0.33};//V I K H W149
const double mzp[M]={22.0, 20.4, 24.6 , 24.9 , 27.615};
const double texp[M]={200.0, 180.0, 60.0,  54.0, 46.8 };
const double mus[M]= {21.8,  19.7,  13.0,  21.5, 21.48}; 
const double AlAv[M]={1.009, 0.60 , 0.118, 0.184, 0.225};
const double sigma[M]={0.022, 0.022, 0.02, 0.025, 0.025};//MOAاستفاده از مقاله کاردلی
const double satu[M]=  {7.0,  12.0, 15.0, 10.0, 14.8}; //it should be changed
const double detect[M]={20.0, 21.0, 21.0, 21.0, 26.0};
const double Akv=0.118;
////=======================================================

const double bLMC= -32.9; 
const double lLMC= 280.5; 
const double RaLMC= 80.89375;
const double DecLMC=-68.2438888888889;  
const double sLMC= 5.0;///degree (half of size of LMC)
const double DLMC= 49.97;///KPC 
const double ddeg=0.2; 
const int nn=int(455);

///**************  Number constant  ********************///
const int Nz=int(10);
const int Ntem=int(8979);///number of rows in file Fstar_disk.dat
const int Nw=int(123);///sigma_WFIRST.dat
const int YZ=3615;///number of rows in file yzma
const int N1=36224, N2=25000, N3=3818, N4=3500;///CMD_BESANCON, thinD, bulge, thickD, halo
/// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///======================================================
struct source{
    int nums,struc, cl;
    double Ds,TET,FI;
    double lg,bg, right, decli;
    double od_disk,od_ThD,od_bulge,od_halo,opt;///optical  depth
    double od_dlmc, od_blmc, od_hlmc; 
    double rho_disk[Num],rho_ThD[Num],rho_halo[Num],rho_stars[Num],rho_bulge[Num];
    double rho_hlmc[Num], rho_dlmc[Num],rho_blmc[Num]; 
    double Rostar0[Num],Rostari[Num],Nstari[Num];
    double Nstart,Rostart,Romaxs,nstart;
    double Nblend[M],blend[M], Fluxb[M], magb[M], Ai[M], Mab[M], Map[M]; 
    double Mw[Ntem],Mz[Ntem],Tem[Ntem];
    double type, Tstar, logl, col, Rstar, gravity, mass, vs;
    double Romax, u0, ro_star; 
    double SNR[M];   
    double  xv, yv, zv; 
    double Avv; 
};

struct lens{
    double Ml,Dl,vl,vs,Vt,xls,ro_lens;
    double rhomaxl,tE,RE;
    int numl,struc;
    double mul, Rl, Mlj; 
    double pt1, pt2, dt, tstar, t0; 
    int type; 
    int Nlens; 
};
struct CMD{
    double Teff_d[N1],logl_d[N1],Mab_d[M][N1],Rs_d[N1],mass_d[N1],type_d[N1]; int cl_d[N1];  ///thin disk
    double Teff_b[N2],logl_b[N2],Mab_b[M][N2],Rs_b[N2],mass_b[N2],type_b[N2]; int cl_b[N2];  /// bulge
    double Teff_t[N3],logl_t[N3],Mab_t[M][N3],Rs_t[N3],mass_t[N3],type_t[N3]; int cl_t[N3];  ///thick disk
    double Teff_h[N4],logl_h[N4],Mab_h[M][N4],Rs_h[N4],mass_h[N4],type_h[N4]; int cl_h[N4];  /// halo
};
struct galactic{
   double l[nrd],b[nrd], RA[nrd], Dec[nrd]; 
};
struct finiteLens{
double finiteL; 
};
///===================== FUNCTION ==============================================
void read_cmd(CMD & cm);
void func_source(source & s, CMD & cm );
void func_lens(lens & l, source & s, int mstep);
void vrel(source & s,lens & l);
void Disk_model(source & s, int, int);
void optical_depth(source & s);
double RandN(double sigma, double);
double finitelens( finiteLens & fl, double u, double rho, double magn0, double rhol);

/////////////////////////////////////////////
    time_t _timeNow;
      unsigned int _randVal;
    unsigned int _dummyVal;
    FILE * _randStream;
////////////////////////////////////////////
///==============================================================//
///                                                              //
///                  Main program                                //
///                                                              //
///==============================================================//

int main()
{

///****************************************************************************
//gettimeofday(&newTime,NULL);
//ftime(&tp);
	time(&_timeNow);
	_randStream = fopen("/dev/urandom", "r");
	_dummyVal = fread(&_randVal, sizeof(_randVal), 1, _randStream);
	srand(_randVal);
    _dummyVal = fread(&_randVal, sizeof(_randVal), 1, _randStream);
	srand(_randVal);
///****************************************************************************

    source s;
    lens l;
    CMD cm;  
    galactic ga; 
    finiteLens fl; 
    VBBinaryLensing vbb;
    vbb.Tol=1.e-3;
    vbb.LoadESPLTable("./files/ESPL.tbl");

    FILE* convert; 
    convert=fopen("./files/convert_coordinate_2.dat","r");
    int coo=0; 
    for(int i=0; i<100; ++i){
    for(int j=0; j<100; ++j){
    fscanf(convert,"%lf  %lf   %lf  %lf \n",&ga.RA[coo], &ga.Dec[coo], &ga.l[coo], &ga.b[coo]);  coo++; }}   
    fclose(convert); 
    FILE* fil1;  FILE* fil2; 
    fil1=fopen("./files/parl_1.txt","w");
    fil2=fopen("./files/parl_2.txt","w");
    fclose(fil1); 
    fclose(fil2); 
    


    read_cmd(cm);     
    s.FI=bLMC/RA;  
    s.TET=lLMC/RA; 
    s.right=RaLMC;
    s.decli=DecLMC;
    Disk_model(s, 0, 0);
    s.Romax=s.Rostart; 
    cout<<"*********************"<<"\t ro_max:  "<<s.Romax<<endl;

   /*int nnb=100;
   FILE*  mapd; 
   mapd=fopen("./files/Map_density_LMC/mapd.txt","w"); 
   for(int l1=0; l1<nnb; l1++) {
   for(int l2=0; l2<nnb; l2++) {
   s.right = double(RaLMC+  (l1-nnb*0.5)*ddeg); ///degree
   s.decli = double(DecLMC+ (l2-nnb*0.5)*ddeg); ///degree
   //========================================
   hh=-1; mindd=10000.0; 
   for(int i=0; i<10000; ++i){
   dist= (s.right- ga.RA[i] )*(s.right- ga.RA[i]) + (s.decli- ga.Dec[i])*(s.decli- ga.Dec[i]); 
   dist=sqrt(dist); 
   if(dist==0.0 or dist<ddeg) {hh=i;  break;}
   if(dist<= mindd) { mindd=dist; hh= i;   }}
   if(hh<0){
   cout<<"ERROR:right:    "<<s.right<<"\t s.dec:  "<<s.decli<<endl; 
   cout<<"mindd:  "<<mindd<<"\t hh:  "<<hh<<endl;  int uue;  cin>>uue; }
   s.lg=double(ga.l[hh]);   
   s.bg=double(ga.b[hh]); 
   //========================================
   s.TET= s.lg/RA; 
   s.FI=  s.bg/RA;
   Disk_model(s, l1, l2);
   fprintf(mapd,"%.5lf  %.5lf    %.5lf  %.5lf   %.5lf  %.5lf    %.5lf  %.5lf     %.10lf  %.10lf\n",
   s.lg,s.bg,s.right, s.decli,double(s.right-RaLMC),double(s.decli-DecLMC),s.xv,s.yv,log10(s.Rostart),log10(s.Nstart));}
   cout<<"l1:  "<<l1<<"\t lat:  "<<s.bg<<endl; }
   exit(0);  */
  

    

    
//===================== Monte Carlo Simulation ================================
double Magn, magni0, deno; 
double  dist, mindd;
int  Nsi=0,  fla,  hh; 
double dec1= -72.5; 
double dec2= -64.8; 
double rig1=double((4.0+27.0/60.0)*15.0 );
double rig2=double((5.0+50.0/60.0)*15.0 );
cout<<"rig1:  "<<rig1<<"\t rig2:   "<<rig2<<endl;
double ddec= fabs(dec2-dec1)/nn/1.0;
double drig= fabs(rig2-rig1)/nn/1.0;
double exti[455]={0.0}; 
double Nmw=0.0, Nlmc=0.0;
double flagb=0.0, maga[M], flagdd[M], sum; 
int l1a, l2a, Nde; 

    int mstep=0; 
  // for(int mstep=0; mstep<36;  mstep+=5){ 
 //  cout<<"*********************** mstep:   "<<mstep<<endl;

   Nde=0;
   do{
   s.right= double((double)rand()/(double)(RAND_MAX+1.))*4.0 + 79.0;//[79,  83]
   s.decli= double((double)rand()/(double)(RAND_MAX+1.))*7.0 - 72.0;//[-72,-63]
   l1a=int( double(rig2-s.right)/drig +1.0) +1; 
   l2a=int( double(dec2-s.decli)/ddec +1.0) +1; 
   if(l1a>=nn)  l1a=nn-1;
   if(l2a>=nn)  l2a=nn-1; 
   if(l1a<1)    l1a=1;
   if(l2a<1)    l2a=1;  
   //cout<<"*******************************************"<<endl;  
   //cout<<"right:  "<<s.right<<"\t decli:  "<<s.decli<<endl;
   //cout<<"l1a:  "<<l1a<<"\t l2a:    "<<l2a<<endl;


    FILE* Extlmc; 
    Extlmc=fopen("./files/EXTLMC_2.dat","r");
    if(!Extlmc) {cout<<"can not open file EXTLMC:  "<<endl;   int uue ;  cin>>uue; }
    
    for(int i=1; i<=l1a; ++i){
    for(int j=1; j<=nn;  ++j){ 
    fscanf(Extlmc,"%lf  ",&exti[j-1]);}
    fscanf(Extlmc,"\n"); }
    fclose(Extlmc); 
    s.Avv= exti[l2a-1] + 0.2*(double)rand()/(double)(RAND_MAX+1.)-0.1;
    if(s.Avv<0.0) s.Avv=0.0; 
    if(s.Avv>5.0){cout<<"Error Extinction_v:  "<<s.Avv<<"\t exti: "<<exti[l2a-1]<<endl;  int yye; cin>>yye; }
   
   //==================================================
   hh=-1; mindd=10000.0; 
   for(int i=0; i<10000; ++i){
   dist= (s.right- ga.RA[i] )*(s.right- ga.RA[i]) + (s.decli- ga.Dec[i])*(s.decli- ga.Dec[i]); 
   dist=sqrt(dist); 
   if(dist==0.0 or dist<ddeg) {hh=i;  break;}
   if(dist<= mindd) {mindd=dist; hh= i;     }}
   if(hh<0 ){cout<<"ERROR:right:    "<<s.right<<"\t s.dec:  "<<s.decli<<endl;   int uue;  cin>>uue;}
   s.lg=double(ga.l[hh]);   
   s.bg=double(ga.b[hh]); 
   //===============================================
   s.TET= s.lg/RA; 
   s.FI=  s.bg/RA;
   Disk_model(s, l1a, l2a);
   
   


   func_source(s, cm);
   func_lens(l,s, mstep);
   optical_depth(s);
   if(l.struc<4.0) Nmw  +=1.0; 
   else            Nlmc +=1.0; 

   
   
   
   magni0= vbb.ESPLMag2(s.u0,s.ro_star);
   if(magni0<0.999){
   vbb.Tol=1.e-2;
   magni0= vbb.ESPLMag2(s.u0,s.ro_star); vbb.Tol=1.e-3; }
   if(l.ro_lens<0.1)  Magn= magni0;  
   else{Magn=finitelens(fl, s.u0, s.ro_star, magni0, l.ro_lens); }
   if(Magn<1.0 and l.ro_lens<0.1){
   cout<<"Error Magn<1.0:    "<<Magn<<"\t u0:  "<<s.u0<<"\t ro_star:  "<<s.ro_star<<endl;
   cout<<"ro_lens:  "<<l.ro_lens<<"\t magni0:  "<<magni0<<endl;
   int yye;  cin>>yye;}
   
   for(int i=0; i<M; ++i){
   s.SNR[i]=0.0; 
   deno= s.Fluxb[i] + M_PI*pow(FWHM[i]*0.5,2.0)*pow(10.0,-0.4*mus[i]) + (Magn-1.0)*pow(10.0,-0.4*s.Map[i]);
   s.SNR[i] = sqrt(texp[i])*Magn*pow(10.0, -0.4*s.Map[i]+ 0.2*mzp[i]);
   s.SNR[i] = s.SNR[i]/sqrt(deno);
   maga[i]= double(s.magb[i]-2.5*log10(s.blend[i]*Magn + 1.0-s.blend[i]));}
   
   
   sum=0.0;  
   for(int i=0; i<M; ++i){
 //  flagb= double((double)rand()/(double)(RAND_MAX+1.));
   if(s.magb[i]<=detect[i] and  maga[i]>=satu[i] and s.SNR[i]>50.0) flagdd[i]=1.0;
   else flagdd[i]=0.0;
   sum+=flagdd[i];}
   
    
   
   
   if(sum>1.0){
   fil1=fopen("./files/parl_1.txt","a+");
   fil2=fopen("./files/parl_2.txt","a+");
   fprintf(fil1,"%d  %.3lf  %.3lf  %d  %d  %e  %e  %e  %.5lf  %.4lf  %.4lf  %d  %.5lf "
   " %d  %d  %.3lf  %.4lf  %.5lf  %.4lf  %.4lf %.5lf  %.1lf  %.5lf   "
   " %.4lf   %.4lf  %.2lf  %.4lf  %.4lf  %.4lf  %.3lf   %.3lf  %d  " 
   "%.10lf  %.10lf  %.4lf  %.4lf  %.6lf %.3lf  %.2lf  %.1lf  %.1lf  %.5lf  %.5lf  %d\n",
   Nsi, s.bg, s.lg,l.struc, l.numl, log10(l.Ml), l.Mlj, double(l.Rl/Rsun), log10(l.ro_lens), l.Dl, l.vl, l.type,l.tstar,  ///13
   s.struc, s.nums, s.magb[4], s.Rstar,log10(s.ro_star), s.Ds, s.vs, s.blend[4],s.Nblend[4], s.opt*1.0e6, //23
   s.mass, s.gravity, s.type, s.Tstar, s.col, s.logl, s.Mab[4], s.Map[4], s.cl,   ///32
   log10(l.tE), log10(l.RE/AU), l.mul, l.Vt, log10(fabs(s.u0)), l.t0,fl.finiteL,Nmw, Nlmc, s.right, s.decli, mstep);//44
   for(int i=0;i<M; ++i)
   fprintf(fil2,"%d  %.3lf  %.3lf  %.4lf  %.3lf   %.4lf   %.1lf  %.3lf  %.1lf \n",
   Nsi, s.Mab[i], s.Map[i], s.Ai[i], s.magb[i], s.blend[i], s.Nblend[i], s.SNR[i], flagdd[i]); 
   fclose(fil1);  fclose(fil2);
   Nde+=1; }

   if(int(Nsi)%3000==0){
   cout<<"============================================================="<<endl;
   cout<<"right_accention:  "<<s.right<<"\t declination:  "<<s.decli<<endl;
   cout<<"lat:  "<<s.bg<<"\t lon:  "<<s.lg<<endl;
   cout<<"Ds:  "<<s.Ds<<"\t nums:  "<<s.nums<<endl;
   cout<<"Nblend[W149]:  "<<s.Nblend[4]<<"\t belnding[W149]:  "<<s.blend[4]<<endl;
   cout<<"Ml(Sun):  "<<l.Ml<<"\t Ml(Jupiter):  "<<l.Mlj<<endl;
   cout<<"Rl[Rjupiter]:  "<<l.Rl/Rjupiter<<"\t ltype: "<<l.type<<"\t u0:  "<<s.u0<<endl;
   cout<<"t0:  "<<l.t0<<"\t tstar:  "<<l.tstar<<"\t ro_lens:  "<<l.ro_lens<<endl;
   cout<<"vl:  "<<l.vl<<"\t vs:  "<<s.vs<<"\t Vt:  "<<l.Vt<<endl;
   cout<<"struc_lens:  "<<l.struc<<"\t struc_source:  "<<s.struc<<"\t type_lens:  "<<l.type<<endl;
   cout<<"Nmw:  "<<Nmw<<"\t Nlmc:  "<<Nlmc<<"\t Nde:  "<<Nde<<"\t Nsim:  "<<Nsi<<endl;
   cout<<"==============================================================="<<endl;}
   Nsi+=1;
   }while(Nde<200000);
 //  } 
   cout<<"**********************************************************"<<endl;
   return(0);
}
////==========================================================================
///==============================================================//
///                                                              //
///                  finite_lens                                 //
///                                                              //
///==============================================================//
double finitelens( finiteLens & fl,  double u, double rho, double magn0, double rhol)
{
   double u2 = u*u;
   double rho2Tol = rho*rho/0.001;
   double u6 = u2*u2*u2;
   double plus= 0.5*fabs(u+ sqrt(u2+4.0));
   double minu= 0.5*fabs(u- sqrt(u2+4.0)); 
   double A= (u2+2.0)/sqrt(u2*(u2+4.0))+0.00000475349256745; 
   double ratio= (A+1.0)/(A-1.0); 
   double finiteS, msmall, plus1, plus2, minu1, minu2, kkd, magn, Rat1, Rat2;
   
   if (u6*(1+0.003*rho2Tol)>0.027680640625*rho2Tol*rho2Tol) finiteS=0.0;
   else     finiteS=1.0; 
   
   msmall= magn0/(1.0+ratio); 
   magn=magn0;
   fl.finiteL=0.0; 

   if(finiteS<1.0){  
   if(rhol>=plus)     {magn=0.0; fl.finiteL=2.0; } 
   else if(rhol>=minu){magn=0.5*(A+1.0);  fl.finiteL=1.0; }}


   
   if(finiteS>0.0){
   plus1= 0.5*fabs(u+rho + sqrt( (u+rho)*(u+rho)+4.0));
   plus2= 0.5*fabs(u-rho + sqrt( (u-rho)*(u-rho)+4.0));//smaller
   minu1= 0.5*fabs(u+rho - sqrt( (u+rho)*(u+rho)+4.0));
   minu2= 0.5*fabs(u-rho - sqrt( (u-rho)*(u-rho)+4.0));//smaller
   if(minu1<minu2) { kkd=minu2;  minu2= minu1;  minu1= kkd;} 
   if(plus1<plus2) { kkd=plus2;  plus2= plus1;  plus1= kkd;} 
   
   Rat1= (plus1*plus1-rhol*rhol)/(plus1*plus1-plus2*plus2);
   Rat2= (minu1*minu1-rhol*rhol)/(minu1*minu1-minu2*minu2);
   
   
   
if(rhol>=plus1){     magn=0.0;    fl.finiteL=2.0;}
else if(rhol<=minu2){magn=magn0;  fl.finiteL=0.0;}
else {
if(plus2<minu1){
if(rhol>=minu1 and rhol<plus1) {magn=msmall*ratio*Rat1; fl.finiteL=2.0 - Rat1;}
if(rhol>=plus2 and rhol<minu1) {magn=msmall*ratio*Rat1 + msmall*Rat2; fl.finiteL=2.0 - Rat1-Rat2;}
if(rhol>=minu2 and rhol<plus2) {magn=msmall*ratio + msmall*Rat2;  fl.finiteL=1.0 - Rat2;} }

if(plus2>=minu1){
if(rhol>=plus2 and rhol<plus1) {magn=msmall*ratio*Rat1; fl.finiteL=2.0 - Rat1;}
if(rhol>=minu1 and rhol<plus2) {magn=msmall*ratio;      fl.finiteL=1.0;}
if(rhol>=minu2 and rhol<minu1) {magn=msmall*ratio + msmall*Rat2;  fl.finiteL=1.0 - Rat2;}}  }



  
   if(plus1<plus2 or minu1<minu2 or magn>magn0 or fl.finiteL>2.00004 or fl.finiteL<0.0 ){
   cout<<"Rat1:  "<<Rat1<<"\t Rat2:  "<<Rat2<<"\t finiteL:  "<<fl.finiteL<<"\t finite_s:  "<<finiteS<<endl;
   cout<<"ERROR: magn:  "<< magn<<"\t magn0:  "<<magn0<<endl;
   cout<<"ERROR:  plus1:   "<<plus1<<"\t plus2:  "<<plus2<<endl;
   cout<<"rho:  "<<rho<<"\t rhoL:    "<<rhol<<"\t u:  "<<u<<endl;
   cout<<"ERROR:  minu1:  "<<minu1<<"\t minu2:   "<<minu2<<endl; int uue ;  cin>>uue; } }
   
   return(magn); 
}
///==============================================================//6
///                                                              //
///                  optical_depth                               //
///                                                              //
///==============================================================//
void optical_depth(source & s)
{
    double ds =(double)s.nums*step;///kpc
    double CC=4.0*G*M_PI*ds*ds*pow(10.0,9.0)*Msun/(velocity*velocity*KP);
    double dl,x,dx;
    s.od_disk=s.od_ThD=s.od_bulge=s.od_halo=s.opt=0.0;
    s.od_dlmc=s.od_hlmc=s.od_blmc=0.0; 
    for(int k =1;k<s.nums;++k){
    dl =(double)k*step;///kpc
    x=dl/ds;
    dx=(double)step/ds/1.0;
    s.od_disk +=  s.rho_disk[k]*x*(1.0-x)*dx*CC;
    s.od_ThD +=   s.rho_ThD[k]*x*(1.0-x)*dx*CC;
    s.od_bulge += s.rho_bulge[k]*x*(1.0-x)*dx*CC;
    s.od_halo +=  s.rho_halo[k]*x*(1.0-x)*dx*CC;
    s.od_dlmc +=  s.rho_dlmc[k]*x*(1.0-x)*dx*CC;
    s.od_blmc +=  s.rho_blmc[k]*x*(1.0-x)*dx*CC; 
    s.od_hlmc +=  s.rho_hlmc[k]*x*(1.0-x)*dx*CC;  }
    s.opt= fabs(s.od_disk+s.od_ThD+s.od_bulge+s.od_halo+ s.od_dlmc + s.od_blmc + s.od_hlmc );///total
    //cout<<"total_opticalD: "<<s.opd<<"\t od_disk: "<<s.od_disk<<endl;
   // cout<<"od_ThD: "<<s.od_ThD<<"\t od_bulge: "<<s.od_bulge<<"\t od_halo: "<<s.od_halo<<endl;
}
///==============================================================//
///                                                              //
///                  func_source   Initial amounts               //
///                                                              //
///==============================================================//
void func_source(source & s, CMD  &  cm)
{

    int num,struc,nums,yye;
    double rho,rf;
    double Ds,Ai[M],Av;
    double Mab[M],Map[M],temp=0.0;
    double maxn;
    double Radius;//double diff[2]={0.0};


    double gg=RandN(1.0,1.0);
    for(int i=0; i<M; ++i){s.Fluxb[i]=0.0;  s.magb[i]=s.Ai[i]=s.Map[i]=s.Mab[i]=0.0;   }
    
    maxn=0.0; 
    for(int i=0; i<M; ++i){
    s.Nblend[i]  =s.Nstart*pow(FWHM[i]*0.5,2)*M_PI/(3600.0*3600.0);
    s.Nblend[i] += gg*s.Nblend[i];
    if(s.Nblend[i]<1.0 or s.Nblend[i]==1.0) s.Nblend[i]=1.0; 
    if(maxn<s.Nblend[i]) maxn= s.Nblend[i];  }
  //  cout<<"Nblend[i]:  "<<s.Nblend[i]<<"\t filter:  "<<i<<"\t maxn:  "<<maxn<<endl;
    
    

    for(int k=1; k<=int(maxn); ++k){
    do{
    num=int((Num-15.00)*(double)rand()/(double)(RAND_MAX+1.) +10.00);
    rho=fabs((double)rand()/((double)(RAND_MAX+1.))*s.Romaxs);
    Ds=(double)(num*step);
    }while(rho>s.Rostari[num] or Ds<0.1 or Ds>MaxD);///distance larger than 50.0
    //Ds=(double)(num*step);
    if(Ds>MaxD or Ds<0.1){cout<<"ERROR (1): Ds: "<<Ds<<"\t MaxD: "<<MaxD<<"\t step: "<<step<<"\t num: "<<num<<endl; cin>>yye;}
    nums=num;
    if(k==1){s.Ds=Ds;  s.nums=nums;}
    rf=((double)rand()/(double)(RAND_MAX+1.))*s.Rostar0[nums];
     if (rf<= s.rho_disk[nums])                    struc=0;///thin disk
else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums])) struc=1;/// bulge structure
else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums]+s.rho_ThD[nums])) struc=2;///thick disk
else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums]+s.rho_ThD[nums]+ s.rho_halo[nums])) struc=3;///halo
else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums]+s.rho_ThD[nums]+ s.rho_halo[nums]+s.rho_hlmc[nums]))struc=4;///halo_lmc
else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums]+s.rho_ThD[nums]+ s.rho_halo[nums]+s.rho_hlmc[nums]+s.rho_dlmc[nums]))struc=5;//disk_lmc
else struc=6;///bulge_lmc

    if(k==1) s.struc=struc;



///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    if(struc==0 or struc==5 ){///thin disk
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N1-2.0));
    for(int i=0; i<M; ++i){Mab[i]=cm.Mab_d[i][num];}
    temp=cm.Teff_d[num];    Radius= cm.Rs_d[num]; 
    if(k==1){
    s.type=cm.type_d[num];
    s.mass=cm.mass_d[num];
    s.Tstar=temp;
    s.col=Mab[0]-Mab[1]; 
    s.logl=cm.logl_d[num];
    s.cl=  cm.cl_d[num];}
    if(cm.mass_d[num]<0.0 or int(cm.cl_d[num])==6 or temp<0.0 or float(cm.type_d[num])>=8.0 or float(cm.mass_d[num])>100000.0){
    cout<<"Error(thin disk) temp: "<<temp<<"\t mass: "<<cm.mass_d[num]<<"\t counter: "<<num<<endl; cin>>yye; }}


    if(struc==1 or struc==6){///bulge
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N2-2.0));
    for(int i=0; i<M; ++i){Mab[i]=cm.Mab_b[i][num];}
    temp=cm.Teff_b[num]; Radius= cm.Rs_b[num];
    if(k==1){
    s.type=cm.type_b[num];
    s.mass=cm.mass_b[num];
    s.Tstar=temp;
    s.col=Mab[0]-Mab[1];
    s.logl=cm.logl_b[num];
    s.cl= cm.cl_b[num];}
    if(cm.mass_b[num]<0.0 or int(cm.cl_b[num])==6 or temp<0.0 or float(cm.type_b[num])>8.0 or cm.mass_b[num]>100000.0){
    cout<<"Error(bulge) temp: "<<temp<<"\t mass: "<<cm.mass_b[num]<<"\t counter: "<<num<<endl;   cin>>yye; }}


    if(struc==2){///thick disk
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N3-2.0));
    for(int i=0; i<M; ++i){Mab[i]=cm.Mab_t[i][num]; }
    temp=cm.Teff_t[num];  Radius= cm.Rs_t[num];
    if(k==1){
    s.type=cm.type_t[num];
    s.mass=cm.mass_t[num];
    s.Tstar=temp;
    s.col=Mab[0]-Mab[1];
    s.logl=cm.logl_t[num];
    s.cl= cm.cl_t[num];}
    if( cm.mass_t[num]<0.0 or int(cm.cl_t[num])==6 or temp<0.0 or float(cm.type_t[num])>8.0 or cm.mass_t[num]>100000.0){
    cout<<"Error(thick disk) temp: "<<temp<<"\t mass: "<<cm.mass_t[num]<<"\t counter: "<<num<<endl;  cin>>yye;}}



    if(struc==3 or struc==4){///stellar halo
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N4-2.0));
    for(int i=0; i<M; ++i){Mab[i]=cm.Mab_h[i][num];}
    temp=cm.Teff_h[num]; Radius= cm.Rs_h[num];
    if(k==1){
    s.type=cm.type_h[num];
    s.mass=cm.mass_h[num];
    s.Tstar=temp;
    s.col=Mab[0]-Mab[1];
    s.logl=cm.logl_h[num];
    s.cl= cm.cl_h[num];}
    if(cm.mass_h[num]<0.0 or int(cm.cl_h[num])==6 or temp<0.0 or float(cm.type_h[num])>8.0 or cm.mass_h[num]>100000.0){
    cout<<"Error(Galactic halo) temp: "<<temp<<"\t mass: "<<cm.mass_h[num]<<"\t counter: "<<num<<endl;  cin>>yye;}    }
    if(k==1) s.Rstar= Radius;//in solar radius
    if(k==1) s.gravity= log10(G*s.mass*Msun*100.0/(Radius*Radius*Rsun*Rsun)); 

   ///*********************** absolute magnitude in W149, Z087 ****************///
   double W149=(Mab[2]+Mab[3]+Mab[4])/3.0;//absolute magnitude in W149 (K+H+J)/3
   Mab[0]=Mab[0]; ///V-band
   Mab[1]=Mab[1]; ///I-band
   Mab[2]=Mab[2]; ///K-band
   Mab[3]=Mab[3]; ///H-band
   Mab[4]=W149;   ///W149-band
  
  
 

   Av=s.Avv;
   if(Av>6.0 or Av<-0.000000343 or Ds>MaxD  or Ds<0.0){ 
   cout<<"ERROR Ds:  "<<Ds<<" \t Av:  "<<Av<<endl; int yyw;  cin>>yyw;  }
   if(Av<0.0)   Av=0.0;         
   for(int i=0; i<M; ++i){    
   Ai[i]=fabs(Av*AlAv[i])+RandN(sigma[i], 1.5); //extinction in other bands
   if(Ai[i]<0.0) Ai[i]=0.0;
   Map[i]=Mab[i]+5.0*log10(Ds*100.0)+Ai[i];
   if(s.Nblend[i]>=k){ 
   s.Fluxb[i]+=fabs(pow(10.0,-0.4*Map[i])); }}
  
 
   if(k==1){
   for(int i=0; i<M; ++i){s.Ai[i]=Ai[i]; s.Map[i]=Map[i]; s.Mab[i]= Mab[i];}
   s.col=s.col+s.Ai[0]-s.Ai[1]; }}///loop over the stars
   
   
   
    for(int i=0; i<M; ++i){
    s.magb[i]=-2.5*log10(fabs(s.Fluxb[i]));
    s.blend[i]=double(pow(10.0,-0.4*s.Map[i])/s.Fluxb[i]);    
    if(int(s.Nblend[i])<1.0 or (s.Nblend[i]==1.0 and s.blend[i]<1.0) or s.blend[i]>1.0 or s.blend[i]<0.0 or
    s.Fluxb[i]<-0.00000000009543){
    cout<<"BIGG ERRROR nsbl: "<<s.Nblend[i]<<"\t Nblend[i]: "<<s.Nblend[i]<<"\t belnd: "<<s.blend[i]<<endl;
    cout<<"BIG ERROR Flux is negative: "<<s.Fluxb[i]<<"\t No. filter:  "<<i<<endl;   cin>>yye;}  }
    //cout<<"filter:    "<<i<<"\t abs magnitude:  "<<s.Mab[i]<<endl;
    //cout<<"apt magnitude:  "<<s.Map[i]<<"\t extin:  "<<s.Ai[i]<<"\t fluxb:  "<<s.Fluxb[i]<<endl; 
    s.u0= (double)rand()/(double)(RAND_MAX+1.0);   
   
    if(s.type>8.0 or s.type<1.0 or s.mass<=0.0 or s.Tstar<0.0  or s.Rstar<0.0 or s.mass>10000.0 or 
    s.nums>Num or s.nums<=0 or Av<0.0 or s.Fluxb[0]<=0.0 or s.cl<0 ){
    cout<<"ERROR(source):  type: "<<s.type<<"\t struc: "<<struc<<"\t num: "<<num<<endl;   cin>>yye;  }
}
///==============================================================//
///                                                              //
///                  func_lens     Initial amounts              //
///                                                              //
///==================================b============================//
void func_lens(lens & l, source & s, int mstep)
{
    double f,test;
    double rholens[s.nums+2]={0.0};
    l.rhomaxl=0.0;
    double Mmin, Mmax, tt;


    for(int k=1;k<int(s.nums-2) ;++k){
    rholens[k]=0.0;
    l.Dl=double(k*step);
    l.xls=l.Dl/s.Ds;
    if(l.Dl>s.Ds){cout<<"ERROR (Dl>Ds) Dl: "<<l.Dl<<"\t Ds: "<<s.Ds<<endl;  int yye; cin>>yye;}
    rholens[k]= sqrt((s.Ds-l.Dl)*l.Dl/s.Ds)*s.Rostar0[k];
    if(rholens[k]>l.rhomaxl) l.rhomaxl=rholens[k];}



    do{
    l.numl = (int)((double)rand()*1000./((double)(RAND_MAX+1.)*1000.)*(s.nums-2.0)+1.0);
    test = ((double)rand()/(double)(RAND_MAX+1.)*l.rhomaxl);
    if(rholens[l.numl]>l.rhomaxl){cout<<"ERROR: rholens[numl]: "<<rholens[l.numl]<<""<<l.rhomaxl<<endl;
    int ue; cin>>ue;}
    }while(test>rholens[l.numl]);
    l.Dl=l.numl*step;///kpc
    



    double  randflag=((double)rand()/(double)(RAND_MAX+1.))*s.Rostar0[l.numl];
       if (randflag<=fabs(s.rho_disk[l.numl]) ) l.struc=0;//disk
  else if (randflag<=fabs(s.rho_disk[l.numl]+s.rho_bulge[l.numl])) l.struc=1;//bulge
  else if (randflag<=fabs(s.rho_disk[l.numl]+s.rho_bulge[l.numl]+s.rho_ThD[l.numl])) l.struc=2;//thick
  else if (randflag<=fabs(s.rho_disk[l.numl]+s.rho_bulge[l.numl]+s.rho_ThD[l.numl]+s.rho_halo[l.numl]) ) l.struc=3;//halo
  else if (randflag<=fabs( s.Rostar0[l.numl]-s.rho_blmc[l.numl]- s.rho_dlmc[l.numl]) ) l.struc=4;//halo_lmc
  else if (randflag<=fabs( s.Rostar0[l.numl]-s.rho_blmc[l.numl])) l.struc=5;//disk_lmc
  else if (randflag<=fabs( s.Rostar0[l.numl])) l.struc=6;//bar_lmc
  else {  cout<<"randflag: "<<randflag<<"\t rho_star0: "<<s.Rostar0[l.numl]<<endl;  int ye; cin>>ye;}



   
    double mmin=double(0.01);//in Earth mass
    double mmax=double(13.0*Mjupiter/Mearth);//0.08*Msun/Mearth);///0.08 Msun
    double hel, f1, f2; 
    double thre=double(13.0*Mjupiter/Mearth); 
    f1=200.0;
    f2=0.000003;
    if(f1<f2){hel=f1;   f1=f2;  f2= hel;}
    do{ 
    l.Ml=(double)rand()/(double)(RAND_MAX+1.0) *(mmax - mmin) +mmin;//in Earth mass
    test=fabs((double)rand()/(double)(RAND_MAX+1.)*(f1-f2) +f2);
    if(l.Ml<5.2)                  f=2.000*pow(l.Ml, -1.00);
    if(l.Ml>=5.2 and l.Ml<=thre)  f=6.667*pow(l.Ml, -1.73);//pow(l.Ml,-1.73);
    }while(test>f);  
    l.Ml= double(l.Ml*Mearth/Msun); 
    
    
    
    //double mmin=-12.0;
    //double mmax=log10(1.0*Mjupiter/Msun); //-1.8442298477737173;
    //l.Ml= double((double)rand()*1000./((double)(RAND_MAX+1.)*1000.)*5.0 -1.0);
    //l.Ml=pow(10.0,l.Ml)* Mearth/Msun; 
    //double min= pow(10.0,double(-10.0+mstep*0.2));
    //l.Ml= double(min + RandN(min*0.025,4.0) );//in Sun mass
     

    l.xls=l.Dl/s.Ds;
    l.RE=sqrt(4.0*G*l.Ml*Msun*s.Ds*KP)/velocity;
    l.RE=l.RE*sqrt(l.xls*(1.0-l.xls));///meter
    l.Vt=0.0; 
    vrel(s,l);
    l.tE=l.RE/(l.Vt*1000.0*3600.0*24.0);///in day
    s.ro_star=s.Rstar*Rsun*l.xls/l.RE; 
    l.mul=l.Vt*1000.0*180.0*3600.0*1000.0*3600.0*24.0*364.5/(l.Dl*KP*M_PI); 
   
   
   
   
    double MlE; 
    l.Mlj=double(l.Ml*Msun/Mjupiter); 
    MlE=  double(l.Ml*Msun/Mearth);   
    if(l.Ml>0.08 or l.Ml==0.08){//stars
    l.type=0;
    l.Rl=pow(MlE,0.88)*0.00145*Rearth;}
    else if(l.Mlj>0.4 or l.Mlj==0.4){///brown dwarfs
    l.type=1; 
    l.Rl= pow(MlE,-0.04)*16.9*Rearth;}
    else if(MlE>2.0 or MlE==2.0){///giant planets
    l.type=2; 
    l.Rl=pow(MlE,0.59)*0.8*Rearth;}
    else{///rocky planets
    l.type=3; 
    l.Rl=pow(MlE,0.28)*Rearth;}
    
    
    
   
    l.ro_lens= double(l.Rl/l.RE);
    l.tstar= double(l.tE*s.ro_star); 
    l.t0=5.0*year*(double)rand()/(double)(RAND_MAX+1.)+1.0; 
    l.pt1= -3.5*l.tE + l.t0; 
    l.pt2= +3.5*l.tE + l.t0; 
    l.dt= double(l.pt2-l.pt1)/l.Nlens; 

    //cout<<"tstar:  "<<l.tstar<<"\t t0:  "<<l.t0<<"\t dt:  "<<l.dt<<endl;
    //cout<<"pt1:  "<<l.pt1<<"\t pt2:  "<<l.pt2<<endl;

    if(s.ro_star<=0.0 or l.ro_lens>100.0 or l.ro_lens<=0.0 or l.tE<0.0  or l.tE==0.0 or l.Dl>s.Ds or l.Vt<=0.0){
    cout<<"ERROR ro_star:  "<<s.ro_star<<endl;
    cout<<"Vt:  "<<l.Vt<<"\t Rs:  "<<s.Rstar<<endl; 
    cout<<"RE: "<<l.RE/AU<<"\t xls:  "<<l.xls<<"\t tE: "<<l.tE<<endl;
    cout<<"Ml[Sun mass]:  "<<l.Ml<<"\t Dl:  "<<l.Dl<<"\t Ds:  "<<s.Ds<<endl;
    cout<<"ro_lens:  "<<l.ro_lens<<"\t Rl[R_jupiter]:  "<<l.Rl/Rjupiter<<endl;
    cout<<"BIG ERROR te: "<<l.tE<<"\t RE(AU): "<<l.RE/AU<<"\t V_rel: "<<l.Vt<<"\t l.Ml:  "<<l.Ml<<endl;
    int iie; cin>>iie;}
 
}
///==============================================================//
///                                                              //
///                  READ CMD FILE                               //
///                                                              //
///==============================================================//
void read_cmd(CMD & cm)
{



    //mass, Teff, Age, logL,  log(g),  Z,  Rs,  MB, MV, MI, MK, Cl, type (13)
    int yye, uui, k1, k2, h, g;  
    double metal, age,gravity, MB;
    char filename[40];
    FILE *fp2;


    double Age1[YZ]={0.0}; double B1[YZ]={0.0};  double M1[YZ]={0.0};   double mm1[YZ]={0.0}; 
    double Age2[YZ]={0.0}; double B2[YZ]={0.0};  double M2[YZ]={0.0};   double mm2[YZ]={0.0}; 
    int number[70]={0};   int count[70]={0};   double Metal[70]={0.0}; 
    FILE *meta; 
    meta=fopen("./files/CMD_WFIRST/metal.txt","r"); 
    for(int i=0; i<70; ++i){
    fscanf(meta,"%lf   %d  %d\n",&Metal[i],&count[i],&number[i]);    
    if((Metal[i]<Metal[i-1] and i>0) or float(Metal[i])<-0.001 or number[i]==0 or count[i]>YZ or (abs(count[i-1]+number[i-1]-count[i])>2 and i>0)){
    cout<<"ERROR Metal[i]: "<<Metal[i]<<"\t count[i]: "<<count[i]<<"\t number[i]: "<<number[i]<<"\t i: "<<i<<endl; cin>>uui;} }
    fclose(meta); 
    FILE *hks; 
    hks=fopen("./files/CMD_WFIRST/HKS.txt", "r"); 
    for(int i=0; i<YZ; ++i){
    fscanf(hks,"%lf  %lf  %lf  %lf\n",&Age1[i],&mm1[i],&B1[i],&M1[i]); 
    if(Age1[i]<0.0 or mm1[i]<0.0 or fabs(B1[i])>0.3 or M1[i]<0.5 or Age1[i]>18.0){   
    cout<<"ERROR Age(HKS): "<<Age1[i]<<"\t metal: "<<mm1[i]<<"\t B[i]"<<B1[i]<<"\t M[i]: "<<M1[i]<<"\t i: "<<i<<endl; cin>>uui; }}
    fclose(hks);
    FILE *ji; 
    ji=fopen("./files/CMD_WFIRST/JI.txt", "r"); 
    for(int i=0; i<YZ; ++i){
    fscanf(ji,"%lf   %lf   %lf  %lf\n",&Age2[i],&mm2[i],&B2[i],&M2[i]); 
    if(Age2[i]<0.0 or mm2[i]<0.0 or fabs(B2[i])>1.7 or M2[i]<0.5 or Age2[i]>18.0  or Age1[i]!=Age2[i] or mm1[i]!=mm2[i]){   
    cout<<"ERROR Age(JI): "<<Age2[i]<<"\t metal: "<<mm2[i]<<"\t B[i]"<<B2[i]<<"\t M[i]: "<<M2[i]<<"\t i: "<<i<<endl;
    cout<<"Age1[i]:  "<<Age1[i]<<"\t mm1[i]:  "<<mm1[i]<<endl;  cin>>uui;}}
    fclose(ji); 




////=================================== THIN DISK ==============================
    int j=0;
    sprintf(filename,"./files/CMD_WFIRST/%c%c%c%c%c%c.dat",'C','M','D','T','i','W');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDTiW.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %d %lf\n",
    &cm.mass_d[j],&cm.Teff_d[j],&age,&cm.logl_d[j],&gravity,&metal,&cm.Rs_d[j],&MB, 
    &cm.Mab_d[0][j],&cm.Mab_d[1][j],&cm.Mab_d[2][j],&cm.cl_d[j],&cm.type_d[j]);
    ///*******************************************************
    h=-1; 
    if(metal<Metal[0] || metal==Metal[0])         h=0; 
    else if(metal>Metal[69] ||  metal==Metal[69]) h=69;
    else{ 
    for(int i=1; i<70; ++i){
    if((metal-Metal[i-1])*(metal-Metal[i])<0.0 || metal==Metal[i-1]){h=int(i-1); break;}}}
    if(h==-1 or h<0 or h>69){cout<<"ERROR1 h: "<<h<<"\t metal: "<<metal<<endl; cin>>uui;} 
    g=-1; k1=int(count[h]);   k2=int(count[h]+number[h]);
    if(k1==k2 or k1>k2 or k1<0 or k2>YZ or mm1[k1]!=mm1[k2-1] or number[h]==0){
    cout<<"ERROR4: k1: "<<k1<<"\t k2: "<<k2<<"\t h: "<<h<<"\t count[h]: "<<count[h]<<"\t number[h]: "<<number[h]<<endl; 
    cout<<"age: "<<age<<"\t Age1[k1]: "<<Age1[k1]<<"\t Age1[k2-1]: "<<Age1[k2-1]<<endl; 
    cout<<"metal[k1]: "<<mm1[k1]<<"\t metal[k2-1]: "<<mm1[k2-1]<<endl;  cin>>uui;}
    if(age<Age1[k1]  or age==Age1[k1])  g=k1; 
    else if(age>Age1[k2-1] or age==Age1[k2-1]) g=int(k2-1); 
    else{
    for(int k=k1+1;  k<k2;  ++k){
    if(Age1[k-1]>Age1[k] or mm1[k-1]!=mm1[k]){
    cout<<"Bad error: Age1[k-1]: "<<Age1[k-1]<<"\t Age1[k]: "<<Age1[k]<<"\t mm[k-1]:"<<mm1[k-1]<<"\t mm[k]"<<mm1[k]<<endl; cin>>uui;} 
    if((age-Age1[k-1])*(age-Age1[k])<0.0  ||  age==Age1[k-1]) {g=int(k-1);   break;}}}
    if(g==-1 or g<0 or g>(YZ-1) or g<k1 or g>k2){cout<<"ERROR2 g:"<<g<<"\t age:"<<age<<"\t metal: "<<metal<<endl; cin>>uui;} 
    cm.Mab_d[3][j]= double(B1[g]+M1[g]*cm.Mab_d[2][j]*1.05263157894737); ///H-band  Mks/Mk=1.05263157894737
    cm.Mab_d[4][j]= double(B2[g]+M2[g]*cm.Mab_d[1][j]);   ///J-band
    if(fabs(cm.Mab_d[3][j]-cm.Mab_d[2][j])>1.5 or fabs(age-Age1[g])>3.0 or fabs(cm.Mab_d[4][j]-cm.Mab_d[1][j])>1.5){
    cout<<"ERROR: Mab_d(y-band): "<<cm.Mab_d[4][j]<<"\t Mab_d(z-band): "<<cm.Mab_d[3][j]<<"\t metal: "<<metal<<endl; 
    cout<<"age: "<<age<<"\t Age1[g]: "<<Age1[g]<<"\t B[g]: "<<B1[g]<<"\t M[g]: "<<M1[g]<<"\t g: "<<g<<"\t Metal:"<<Metal[h]<<endl;
    cin>>uui;}
    ///********************************************************
    if(cm.mass_d[j]<0.0 or cm.mass_d[j]==0.0 or cm.Teff_d[j]<0.0 or metal>0.1 or age>10 or int(cm.cl_d[j])==6 or 
    float(cm.type_d[j])>8.0 or cm.type_d[j]<1.0 or (cm.cl_d[j]==5 and float(cm.type_d[j])>8) or (cm.cl_d[j]<5 and cm.type_d[j]>=8.0)){
    cout<<"ERROR(reading cmd file) structure thin disk: "<<"\t Nfil: "<<j+1<<endl;
    cout<<"type_d: "<<cm.type_d[j]<<"\t CL(thin disk): "<<cm.cl_d[j]<<endl;
    cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N1){
    cout<<"BIG ERRROR j: "<<j<<"\t N1: "<<N1<<endl;  cin>>yye;}
    //cout<<"End of CMD reading (thin disk):  No. rows file: "<<j<<endl;






////=================================== BULGE ==================================
    j=0;
    sprintf(filename,"./files/CMD_WFIRST/%c%c%c%c%c.dat",'C','M','D','b','W');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDbW.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %d %lf\n",
    &cm.mass_b[j],&cm.Teff_b[j],&age,&cm.logl_b[j],&gravity,&metal,&cm.Rs_b[j],&MB,
    &cm.Mab_b[0][j],&cm.Mab_b[1][j],&cm.Mab_b[2][j],&cm.cl_b[j],&cm.type_b[j]);
    ///*****************************************************
    h=-1; 
    if(metal<Metal[0] || metal==Metal[0]) h=0; 
    else if(metal>Metal[69] ||  metal==Metal[69]) h=69;
    else{ 
    for(int i=1; i<70; ++i){
    if((metal-Metal[i-1])*(metal-Metal[i])<0.0 || metal==Metal[i-1]){h=int(i-1); break;}}}
    if(h==-1 or h<0 or h>69) {cout<<"ERROR1 h: "<<h<<"\t metal: "<<metal<<endl; cin>>uui;} 
    g=-1;    k1=int(count[h]);   k2=int(count[h]+number[h]);
    if(k1==k2 or k1>k2 or k1<0 or k2>YZ or mm1[k1]!=mm1[k2-1] or number[h]==0){
    cout<<"ERROR4: k1: "<<k1<<"\t k2: "<<k2<<"\t h: "<<h<<"\t count[h]: "<<count[h]<<"\t number[h]: "<<number[h]<<endl; 
    cout<<"age: "<<age<<"\t Age1[k1]: "<<Age1[k1]<<"\t Age1[k2-1]: "<<Age1[k2-1]<<endl; 
    cout<<"metal[k1]: "<<mm1[k1]<<"\t metal[k2-1]: "<<mm1[k2-1]<<endl;  cin>>uui;}
    if(age<Age1[k1]  or age==Age1[k1])  g=k1; 
    else if(age>Age1[k2-1] or age==Age1[k2-1]) g=int(k2-1); 
    else{
    for(int k=k1+1;  k<k2;  ++k){
    if(Age1[k-1]>Age1[k] or mm1[k-1]!=mm1[k]){
    cout<<"Bad error: Age1[k-1]: "<<Age1[k-1]<<"\t Age1[k]: "<<Age1[k]<<"\t mm[k-1]:"<<mm1[k-1]<<"\t mm[k]"<<mm1[k]<<endl; cin>>uui;} 
    if((age-Age1[k-1])*(age-Age1[k])<0.0  ||  age==Age1[k-1]) {g=int(k-1);   break;}}}
    if(g==-1 or g<0 or g>(YZ-1) or g<k1 or g>k2){cout<<"ERROR2 g:"<<g<<"\t age:"<<age<<"\t metal: "<<metal<<endl; cin>>uui;} 
    cm.Mab_b[3][j]= double(B1[g]+M1[g]*cm.Mab_b[2][j]*1.05263157894737); ///H-band  ?=Mks/Mk
    cm.Mab_b[4][j]= double(B2[g]+M2[g]*cm.Mab_b[1][j]);   ///J-band
    if(fabs(cm.Mab_b[3][j]-cm.Mab_b[2][j])>1.5 or fabs(age-Age1[g])>3.0 or fabs(cm.Mab_b[4][j]-cm.Mab_b[1][j])>1.5){
    cout<<"ERROR: Mab_b(y-band): "<<cm.Mab_b[4][j]<<"\t Mab_b(z-band): "<<cm.Mab_b[3][j]<<"\t metal: "<<metal<<endl; 
    cout<<"age: "<<age<<"\t Age1[g]: "<<Age1[g]<<"\t B[g]: "<<B1[g]<<"\t M[g]: "<<M1[g]<<"\t g: "<<g<<"\t Metal:"<<Metal[h]<<endl;
    cin>>uui; }
    ///*****************************************************
    if(cm.mass_b[j]<0.0|| cm.mass_b[j]==0.0 ||cm.Teff_b[j]<0.0|| age>10 ||metal>0.9||cm.cl_b[j]==6 ||cm.type_b[j]>=8.0 or
    (cm.cl_b[j]==5 and int(cm.type_b[j])>8.0) or (cm.cl_b[j]==6 and int(cm.type_b[j])!=9) or (cm.cl_b[j]<5 and int(cm.type_b[j])==9)){
    cout<<"ERROR(reading cmd file) structure bulge: "<<"\t Nfil: "<<j+1<<endl;
    cout<<"type_b: "<<cm.type_b[j]<<"\t CL(bulge): "<<cm.cl_b[j]<<endl;
    cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N2){cout<<"BIG ERRROR j: "<<j<<"\t N2: "<<N2<<endl;  cin>>yye;}
   // cout<<"End of CMD reading (bulge):  No. rows file: "<<j<<endl;






////=================================== THICK DISK =============================
    j=0;
    sprintf(filename,"./files/CMD_WFIRST/%c%c%c%c%c%c.dat",'C','M','D','T','k','W');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDTkW.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %d %lf\n",
    &cm.mass_t[j],&cm.Teff_t[j],&age,&cm.logl_t[j],&gravity,&metal,&cm.Rs_t[j],&MB,
    &cm.Mab_t[0][j],&cm.Mab_t[1][j],&cm.Mab_t[2][j],&cm.cl_t[j],&cm.type_t[j]);
    ///*********************************************************
    h=-1; 
    if(metal<Metal[0] || metal==Metal[0]) h=0; 
    else if(metal>Metal[69] ||  metal==Metal[69]) h=69;
    else{ 
    for(int i=1; i<70; ++i){
    if((metal-Metal[i-1])*(metal-Metal[i])<0.0 || metal==Metal[i-1]){h=int(i-1); break;}}}
    if(h==-1 or h<0 or h>69) {cout<<"ERROR1 h: "<<h<<"\t metal: "<<metal<<endl; cin>>uui;} 
    g=-1;    k1=int(count[h]);   k2=int(count[h]+number[h]);
    if(k1==k2 or k1>k2 or k1<0 or k2>YZ or mm1[k1]!=mm1[k2-1] or number[h]==0){
    cout<<"ERROR4: k1: "<<k1<<"\t k2: "<<k2<<"\t h: "<<h<<"\t count[h]: "<<count[h]<<"\t number[h]: "<<number[h]<<endl; 
    cout<<"age: "<<age<<"\t Age1[k1]: "<<Age1[k1]<<"\t Age1[k2-1]: "<<Age1[k2-1]<<endl; 
    cout<<"metal[k1]: "<<mm1[k1]<<"\t metal[k2-1]: "<<mm1[k2-1]<<endl;  cin>>uui;}
    if(age<Age1[k1]  or age==Age1[k1])  g=k1; 
    else if(age>Age1[k2-1] or age==Age1[k2-1]) g=int(k2-1); 
    else{
    for(int k=k1+1;  k<k2;  ++k){
    if(Age1[k-1]>Age1[k] or mm1[k-1]!=mm1[k]){
    cout<<"Bad error: Age1[k-1]: "<<Age1[k-1]<<"\t Age1[k]: "<<Age1[k]<<"\t mm[k-1]:"<<mm1[k-1]<<"\t mm[k]"<<mm1[k]<<endl;  cin>>uui;} 
    if((age-Age1[k-1])*(age-Age1[k])<0.0  ||  age==Age1[k-1]) {g=int(k-1);   break;}}}
    if(g==-1 or g<0 or g>(YZ-1) or g<k1 or g>k2){cout<<"ERROR2 g:"<<g<<"\t age:"<<age<<"\t metal: "<<metal<<endl; cin>>uui;} 
    cm.Mab_t[3][j]= double(B1[g]+M1[g]*cm.Mab_t[2][j]*1.05263157894737); ///H-band  ?=Mks/Mk
    cm.Mab_t[4][j]= double(B2[g]+M2[g]*cm.Mab_t[1][j]);   ///J-band
    if(fabs(cm.Mab_t[3][j]-cm.Mab_t[2][j])>1.5 or fabs(age-Age1[g])>3.0 or fabs(cm.Mab_t[4][j]-cm.Mab_t[1][j])>1.5){
    cout<<"ERROR: Mab_t(y-band): "<<cm.Mab_t[4][j]<<"\t Mab_t(z-band): "<<cm.Mab_t[3][j]<<"\t metal: "<<metal<<endl; 
    cout<<"age: "<<age<<"\t Age1[g]: "<<Age1[g]<<"\t B[g]: "<<B1[g]<<"\t M[g]: "<<M1[g]<<"\t g: "<<g<<"\t Metal:"<<Metal[h]<<endl;
    cin>>uui; }
    ///********************************************************
    if(cm.mass_t[j]<0.0||  cm.mass_t[j]==0.0 or cm.Teff_t[j]<0.0 or metal>0.2||cm.cl_t[j]==6|| cm.type_t[j]>=8.0 or
    (cm.cl_t[j]==5 and int(cm.type_t[j])>8) or (cm.cl_t[j]==6 and int(cm.type_t[j])!=9) or (cm.cl_t[j]<5 and int(cm.type_t[j])==9)){
    cout<<"type_thick: "<<cm.type_t[j]<<"\t CL(thick): "<<cm.cl_t[j]<<endl;
    cout<<"ERROR(reading cmd file) structure thick disk: "<<"\t Nfil: "<<j+1<<endl; cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N3){cout<<"BIG ERRROR j: "<<j<<"\t N3: "<<N3<<endl;  cin>>yye;}
    //cout<<"End of CMD reading (thick disk):  No. rows file: "<<j<<endl;








////=================================== STELLAR HALO ===========================
    j=0;
    sprintf(filename,"./files/CMD_WFIRST/%c%c%c%c%c.dat",'C','M','D','h','W');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDhW.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %d %lf\n",
    &cm.mass_h[j],&cm.Teff_h[j],&age,&cm.logl_h[j],&gravity,&metal,&cm.Rs_h[j],&MB,
    &cm.Mab_h[0][j],&cm.Mab_h[1][j],&cm.Mab_h[2][j],&cm.cl_h[j],&cm.type_h[j]);
    ///************************************************************
    h=-1; 
    if(metal<Metal[0] || metal==Metal[0]) h=0; 
    else if(metal>Metal[69] ||  metal==Metal[69]) h=69;
    else{ 
    for(int i=1; i<70; ++i){
    if((metal-Metal[i-1])*(metal-Metal[i])<0.0 || metal==Metal[i-1]){h=int(i-1); break;}}}
    if(h==-1 or h<0 or h>69) {cout<<"ERROR1 h: "<<h<<"\t metal: "<<metal<<endl; cin>>uui;} 
    g=-1;    k1=int(count[h]);   k2=int(count[h]+number[h]);
    if(k1==k2 or k1>k2 or k1<0 or k2>YZ or mm1[k1]!=mm1[k2-1] or number[h]==0){
    cout<<"ERROR4: k1: "<<k1<<"\t k2: "<<k2<<"\t h: "<<h<<"\t count[h]: "<<count[h]<<"\t number[h]: "<<number[h]<<endl; 
    cout<<"age: "<<age<<"\t Age1[k1]: "<<Age1[k1]<<"\t Age1[k2-1]: "<<Age1[k2-1]<<endl; 
    cout<<"metal[k1]: "<<mm1[k1]<<"\t metal[k2-1]: "<<mm1[k2-1]<<endl;  cin>>uui;}
    if(age<Age1[k1]  or age==Age1[k1])  g=k1; 
    else if(age>Age1[k2-1] or age==Age1[k2-1]) g=int(k2-1); 
    else{
    for(int k=k1+1;  k<k2;  ++k){
    if(Age1[k-1]>Age1[k] or mm1[k-1]!=mm1[k]){
    cout<<"Bad error: Age1[k-1]: "<<Age1[k-1]<<"\t Age1[k]: "<<Age1[k]<<"\t mm[k-1]:"<<mm1[k-1]<<"\t mm[k]"<<mm1[k]<<endl;  cin>>uui;} 
    if((age-Age1[k-1])*(age-Age1[k])<0.0  ||  age==Age1[k-1]) {g=int(k-1);   break;}}}
    if(g==-1 or g<0 or g>(YZ-1) or g<k1 or g>k2){cout<<"ERROR2 g:"<<g<<"\t age:"<<age<<"\t metal: "<<metal<<endl; cin>>uui;} 
    cm.Mab_h[3][j]= double(B1[g]+M1[g]*cm.Mab_h[2][j]*1.05263157894737); ///H-band  ?=Mks/Mk
    cm.Mab_h[4][j]= double(B2[g]+M2[g]*cm.Mab_h[1][j]);   ///J-band
    if(fabs(cm.Mab_h[3][j]-cm.Mab_h[2][j])>1.5 or fabs(age-Age1[g])>3.0 or fabs(cm.Mab_h[4][j]-cm.Mab_h[1][j])>1.5){
    cout<<"ERROR: Mab_h(y-band): "<<cm.Mab_h[4][j]<<"\t Mab_h(z-band): "<<cm.Mab_h[3][j]<<"\t metal: "<<metal<<endl; 
    cout<<"age: "<<age<<"\t Age1[g]: "<<Age1[g]<<"\t B[g]: "<<B1[g]<<"\t M[g]: "<<M1[g]<<"\t g: "<<g<<"\t Metal:"<<Metal[h]<<endl;
    cin>>uui; }
    ///***********************************************************
    if(cm.mass_h[j]<0.0 || cm.mass_h[j]==0.0 || age<0 or cm.cl_h[j]<0  or cm.cl_h[j]==6  or  cm.Teff_h[j]<0.0 or
    metal>0.1 || cm.cl_h[j]>7|| cm.type_h[j]>9 or (cm.cl_h[j]==5 and int(cm.type_h[j])>8) or (cm.cl_h[j]==6 and int(cm.type_h[j])!=9) or
    (cm.cl_h[j]<5 and int(cm.type_h[j])==9)){
    cout<<"type_halo: "<<cm.type_h[j]<<"\t CL(halo): "<<cm.cl_h[j]<<endl;
    cout<<"ERROR(reading cmd file) structure halo: "<<"\t Nfil: "<<j+1<<endl; cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N4){cout<<"BIG ERRROR j: "<<j<<"\t N4: "<<N4<<endl;  cin>>yye;}
   // cout<<"End of CMD reading (halo):  No. rows file: "<<j<<endl;
   // cout<<">>>>>>>>>>>>>>>>> END OF CMD READING <<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;


}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                       Glactic model                            ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void Disk_model(source & s, int l1, int l2)
{
   double x,xb,yb,zb,r4,r2,rdi,Rb;
   double nnf=0.4/0.8;
   double mBarre;///stars/pc^3
   double Rx0, Ry0, Rz0,Rc,cp,cn,rhoE,rhoS;
   double alfa=12.89/RA;
   double xf,yf,zf,rho;
   s.Romaxs=s.Nstart=s.Rostart=0.0;
   double fd=1.0;///see the program mass_averaged.cpp. we do not apply any limitation
   double fb=1.0;///0.657066;////just stars brighter than V=11.5, but we change to consider all stars  
   double fh=1.0;///No limitation 
   double Rdd=2.17;///2.53;///2.17;
   double Rhh=1.33;///1.32;//1.33;


   double alm= 2.0; ////2KPC,  Mancini 2004
   double rlm_0= 1.76*0.01;/// solar mass/ Pc^-3
   double M_dlm=2.6;//[M_sun]mass of disk LMC  from Mancini 2004 paper
   double M_blm=0.15*M_dlm;///[M_sun] mass of bar LMC
   double xb0, yb0, zb0;///KPC
   double frac=0.05; // fraction of halo in the form of compact objects
   double qq=0.688; 
   double Rd0=1.54;
   double xol, yol, zol, x0, y0, z0; 
   double r0,  pos1, inc1, Rlm; 

   
   
   char filename[40];
   FILE *fill;



/*   int flagf=0;
   if(fabs(l1-nn*0.5)<2.0  and fabs(l2-nn*0.5)<2.0 ) {
   flagf=1; 
   sprintf(filename,"./files/density/%c%d%c%d.dat",'d',l1,'_',l2);
   fill=fopen(filename,"w");
   if(!fill){cout<<"cannot open file longtitude : "<<s.lg<<"\t latitude: "<<s.bg<<endl;  exit(0);}
   }*/


for(int i=1;i<Num;++i){
   s.Rostar0[i]=s.Rostari[i]=s.Nstari[i]=0.0;
   s.rho_disk[i]=s.rho_bulge[i]=s.rho_halo[i]=s.rho_ThD[i]=0.0;
   s.rho_hlmc[i]=s.rho_dlmc[i]=s.rho_blmc[i]=0.0; 

   x=i*step;
   zb = sin(s.FI)*x;
   yb = cos(s.FI)*sin(s.TET)*x;
   xb = x*cos(s.FI)*cos(s.TET)-Dsun;
   Rb=sqrt(xb*xb+yb*yb);
   double rsun= sqrt(zb*zb+ yb*yb+ xb*xb); 


///========== Galactic Thin Disk =====================
   for(int ii=0; ii<8; ++ii){
   rdi=Rb*Rb+zb*zb/(epci[ii]*epci[ii]);
   if(ii==0)     rho=exp(-rdi/25.0)-exp(-rdi/9.0);
   else if(ii>0) rho=exp(-sqrt(0.25+rdi/(Rdd*Rdd)))-exp(-sqrt(0.25+rdi/(Rhh*Rhh)));
   s.rho_disk[i]=s.rho_disk[i]+ rho0[ii]*corr[ii]*0.001*rho/d0[ii];}///M_sun/pc^3
///=================================================


///========== Galactic Thick Disk =====================

  double rho00=1.34*0.001+3.04*0.0001;
  if(fabs(zb)<0.4) s.rho_ThD[i]=(rho00/0.999719)*exp(-(Rb-Dsun)/2.5)*(1.0-zb*zb/(0.4*0.8*(2.0+nnf)));
  else s.rho_ThD[i]=(rho00/0.999719)*exp(-(Rb-Dsun)/2.5)*exp(nnf)*exp(-fabs(zb)/0.8)/(1.0+0.5*nnf);///M_sun/pc^3
///=================================================


///========== Galactic Stellar Halo=================
   rdi=sqrt(Rb*Rb+ zb*zb/(0.76*0.76));
   if( rdi <=0.5)  s.rho_halo[i]=frac*(0.932*0.00001/867.067)*pow(0.5/Dsun,-2.44);
   else            s.rho_halo[i]=frac*(0.932*0.00001/867.067)*pow(rdi/Dsun,-2.44);///M_sun/pc^3
///=================================================



///========== Galactic bulge =====================
   xf = xb*cos(alfa) + yb*sin(alfa);
   yf =-xb*sin(alfa) + yb*cos(alfa);
   zf = zb;
   Rx0=1.46, Ry0=0.49, Rz0=0.39; Rc=3.43; cp=3.007;  cn=3.329;  mBarre=35.45/(3.84723);
   r4=pow(pow(fabs(xf/Rx0),cn)+pow(fabs(yf/Ry0),cn),cp/cn)+pow(fabs(zf/Rz0),cp);
   r4=pow(fabs(r4),1.0/cp);
   r2=sqrt(fabs(xf*xf+yf*yf));
   if(r2<=Rc) rhoS= mBarre*1.0/(cosh(-r4)*cosh(-r4));
   else       rhoS= mBarre*1.0/(cosh(-r4)*cosh(-r4))*exp(-4.0*(r2-Rc)*(r2-Rc));

   Rx0=4.44, Ry0=1.31, Rz0=0.80; Rc=6.83; cp=2.786; cn=3.917; mBarre=2.27/87.0;//85.3789;
   r4=pow(fabs(pow(fabs(xf/Rx0),cn)+pow(fabs(yf/Ry0),cn)),cp/cn)+pow(fabs(zf/Rz0),cp);
   r4=pow(r4,1.0/cp);
   r2=sqrt(fabs(xf*xf+yf*yf));
   if(r2<=Rc) rhoE= mBarre*exp(-r4);
   else       rhoE= mBarre*exp(-r4)*exp(-4.0*(r2-Rc)*(r2-Rc));
  s.rho_bulge[i]= fabs(rhoS)+fabs(rhoE);///M_sun/pc^3
///==================================================================





///=====================  LMC density profiles ====================== 
   x0 = -x*cos(s.decli/RA)*sin(s.right/RA-RaLMC/RA);
   y0 =  x*sin(s.decli/RA)*cos(DecLMC/RA)-x*cos(s.decli/RA)*sin(DecLMC/RA)*cos(s.right/RA-RaLMC/RA);
   z0 = DLMC- x*cos(s.decli/RA)*cos(DecLMC/RA)*cos(s.right/RA-RaLMC/RA)-x*sin(s.decli/RA)*sin(DecLMC/RA);
   if(fabs(x-DLMC)<double(step*1.5)){s.xv=x0;  s.yv=y0;   s.zv=z0;}
   r0= sqrt(x0*x0+y0*y0+z0*z0); 
  
///=====================  LMC density HALO ======================
   if(r0<15.0) s.rho_hlmc[i]=fabs(rlm_0*frac/(1.0+r0*r0/alm/alm) );  ///M_sun/PC^3
   else         s.rho_hlmc[i]=0.0;   
  
  
///=====================  LMC density Disk ======================
   pos1=(170.0-90.0)*M_PI/180.0; 
   inc1=34.7*M_PI/180.0;  
   xol= x0*cos(pos1)+ y0* sin(pos1);
   yol=-x0*sin(pos1)*cos(inc1) + y0*cos(pos1)*cos(inc1) -z0*sin(inc1);
   zol=-x0*sin(pos1)*sin(inc1) + y0*cos(pos1)*sin(inc1) +z0*cos(inc1);    
   Rlm= sqrt(xol*xol + yol*yol ); 
   double zd0=0.3;//KP Kim's paper
   double Rd1=1.8;///KPc  Kim's paper
   s.rho_dlmc[i]=fabs(M_dlm/(4.0*M_PI*zd0*Rd1*Rd1)*exp(-Rlm/Rd1)*exp(-fabs(zol/zd0))); ////kim [M_sun/Pc^3]  Kim 2000
   

///=====================  LMC density Bulge ======================
   xb0= 1.2;  yb0=zb0= 0.44;
   pos1=(110.0-90.0)*M_PI/180.0; 
   inc1=0.0*M_PI/180.0;  
   xol= x0* cos(pos1) +  y0* sin(pos1);
   yol=-x0*sin(pos1)*cos(inc1) + y0*cos(pos1)*cos(inc1) -z0*sin(inc1);
   zol=-x0*sin(pos1)*sin(inc1) + y0*cos(pos1)*sin(inc1) +z0*cos(inc1);    
   s.rho_blmc[i]=M_blm*pow(2.0*M_PI,-1.5)/(xb0*yb0*zb0)*exp(-0.5*(xol*xol/xb0/xb0+yol*yol/yb0/yb0+zol*zol/zb0/zb0));///[M_sun/Pc^3]  
///===============================================================================



s.Rostar0[i]=fabs(s.rho_disk[i])+fabs(s.rho_ThD[i])+fabs(s.rho_bulge[i])+fabs(s.rho_halo[i]) +fabs(s.rho_hlmc[i])+fabs(s.rho_dlmc[i])+ fabs(s.rho_blmc[i]);///[M_sun/pc^3]

s.Rostari[i]=s.Rostar0[i]*x*x*step*1.0e9*(M_PI/180.0)*(M_PI/180.0);///[M_sun/deg^2]

s.Nstari[i]= binary_fraction*((s.rho_disk[i]+s.rho_dlmc[i])*fd/0.403445+s.rho_ThD[i]*fh/0.4542+(s.rho_halo[i]+s.rho_hlmc[i])*fh/0.4542+(s.rho_bulge[i]+s.rho_blmc[i])*fb/0.308571);////[Nt/pc^3] 

s.Nstari[i]=s.Nstari[i]*x*x*step*1.0e9*(M_PI/180.0)*(M_PI/180.0);///[Ni/deg^2]

s.Nstart  +=  s.Nstari[i];///[Nt/deg^2]
s.Rostart += s.Rostari[i];///[Mt/deg^2]
if(s.Rostari[i]>s.Romaxs) s.Romaxs=s.Rostari[i];///source selection

//if(flagf>0)
//fprintf(fill,"%e   %e   %e   %e   %e  %e  %e   %e   %e   %e\n",
//  x,s.rho_disk[i],s.rho_bulge[i],s.rho_ThD[i],s.rho_halo[i],s.rho_dlmc[i],s.rho_blmc[i],s.rho_hlmc[i], s.Rostar0[i],s.Nstari[i]); 
 // cout<<"rho_disk(LMC):  "<<s.rho_dlmc[i]<<"\t rho_bar(LMC): "<<s.rho_blmc[i]<<endl;
  //cout<<"rho_halo(LMC):  "<<s.rho_hlmc[i]<<endl;
   }

//if(flagf>0)  fclose(fill);


   //cout<<"xLMC:  "<<xLMC<<"\t yLMC:  "<<yLMC<<"\t zLMC:  "<<zLMC<<endl;
   //cout<<"LMC_distance:  "<<sqrt(xLMC*xLMC + yLMC*yLMC +  zLMC*zLMC)<<endl;
   //cout<<"xol:  "<<xol<<"\t yol:  "<<yol<<"\t  zol:  "<<zol<<endl;
   //cout<<"distance_LMC_Center:  "<<rlm<<"\t projected_distance:  "<<Rlm<<endl;
  // cout<<"Nstart [Nt/deg^2]: "<<s.Nstart<<"\t Ro_star [Mass/deg^2]: "<<s.Rostart<<endl;
   //cout<<">>>>>>>>>>>>>>>>>>>>>>>>> END OF DISK MODLE <<<<<<<<<<<<<<<<<<<<"<<endl;
   //exit(0);
   
}
///>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
///===========================================================================//
double RandN(double sigma, double nnd){
   double rr,f,frand;
   do{
   rr=double(((double)rand()/(double)(RAND_MAX+1.))*2.0-1.0)*sigma*nnd; ///[-N sigma:N sigma]
   f= exp(-0.5*rr*rr/sigma/sigma);
   frand=fabs((double)rand()/((double)(RAND_MAX+1.))*1.0);
   }while(frand>f);
   return rr;
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                       relative velocity                        ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void vrel(source & s, lens & l)
{
if (l.Dl==0.0) l.Dl=0.00034735;
 double Rlc=sqrt(l.Dl*l.Dl*cos(s.FI)*cos(s.FI)+Dsun*Dsun-2.*Dsun*l.Dl*cos(s.TET)*cos(s.FI));
 double Rsc=sqrt(s.Ds*s.Ds*cos(s.FI)*cos(s.FI)+Dsun*Dsun-2.*Dsun*s.Ds*cos(s.TET)*cos(s.FI));
 if(Rlc==0.0) Rlc=0.00034346123;
 if(Rsc==0.0) Rsc=0.0004762654134;  
 ///Source and Lens velocity components in Galactocentric cylindrical coordinates
 double SVT, SVR, SVZ, LVT, LVR, LVZ,SVt,LVt;
 ///Source and Lens velocity components in heliocenteric galactic coordinates
 double SVb, SVl, LVb, LVl;
 double fv, testfv;
 double VSunl,VSunt,VSunb,vls_b,vls_l;
 double betal,betas,deltal,deltas,deltao;

 double NN=3.0;
 double VSunR =-10.3;
 double VSunT =vro_sun*(1.00762+0.00712)+6.3;
 double VSunZ = 5.9;
 double sigma_R_Disk= 43.0,  sigma_T_Disk= 27.8, sigma_Z_Disk=17.5;
 double sigma_R_TDisk=67.0,  sigma_T_TDisk=51.0, sigma_Z_TDisk=42.0;
 double sigma_R_halo= 131.0, sigma_T_halo=106.0, sigma_Z_halo=85.0;
 double sigma_R_Bulge=113.0, sigma_T_Bulge=115.0, sigma_Z_Bulge=100.0;
 double Rho[8]={00.0}; double maxr=0.0;
 for(int i=0; i<8; ++i){  Rho[i]=rho0[i]*corr[i]/d0[i]; maxr=maxr+ Rho[i];}
 
  double v_R_lmc=-57.0;
  double v_T_lmc=-226.0; 
  double v_Z_lmc= 221.0;
  double sigma_LMC=20.2; 
  double err_rlmc= 13.0; ///error of global velocity
  double err_tlmc= 15.0; 
  double err_zlmc= 19.0; 
 

  double test= ((double)rand()/(double)(RAND_MAX+1.))*maxr; ///total ages
     if(test<=Rho[0])     {sigma_R_Disk=16.7; sigma_T_Disk=10.8; sigma_Z_Disk=6.0;}
else if(test<=(Rho[0]+Rho[1])) {sigma_R_Disk=19.8; sigma_T_Disk=12.8; sigma_Z_Disk=8.0;}
else if(test<=(Rho[0]+Rho[1]+Rho[2]))  {sigma_R_Disk=27.2; sigma_T_Disk=17.6; sigma_Z_Disk=10.0;}
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]))  
                           {sigma_R_Disk=30.2; sigma_T_Disk=19.5; sigma_Z_Disk=13.2;}
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]+Rho[4]))  
                           {sigma_R_Disk=36.7; sigma_T_Disk=23.7; sigma_Z_Disk=15.8;}
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]+Rho[4]+Rho[5]))
                           {sigma_R_Disk=43.1; sigma_T_Disk=27.8; sigma_Z_Disk=17.4;}
else if(test<=maxr)        {sigma_R_Disk=43.1; sigma_T_Disk=27.8; sigma_Z_Disk=17.5;}
else  {cout<<"BIG ERROR "<<test<<"\t maxr: "<<maxr<<endl;  int yye; cin>>yye;}  


///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
/// Generate Source velocity components in Glactocenteric cylindrical coordinates(x',y')
    SVR=SVT=SVZ=0.0;
    if(s.struc==0){///Galactic disk
    SVR= RandN(sigma_R_Disk, NN); 
    SVT= RandN(sigma_T_Disk, NN);
    SVZ= RandN(sigma_Z_Disk, NN); 
    SVT =SVT+ vro_sun *(1.00762 * pow(Rsc/Dsun,0.0394) + 0.00712);}

   else if(s.struc==1){///Galactic bulge
   SVZ=RandN(sigma_Z_Bulge, NN); 
   SVR=RandN(sigma_R_Bulge, NN);
   SVT=RandN(sigma_T_Bulge, NN);}

   else if(s.struc==2){///thick disk
   SVR= RandN(sigma_R_TDisk, NN); 
   SVT= RandN(sigma_T_TDisk, NN); 
   SVZ= RandN(sigma_Z_TDisk, NN);
   SVT =SVT+ vro_sun *(1.00762*pow(Rsc/Dsun,0.0394) + 0.00712); }
   
   else if(s.struc==3){///stellar halo
   SVR= RandN(sigma_R_halo, NN); 
   SVT= RandN(sigma_T_halo, NN); 
   SVZ= RandN(sigma_Z_halo, NN);}
   
   else if(s.struc>3){
   SVR= RandN(sigma_LMC, NN); 
   SVT= RandN(sigma_LMC, NN); 
   SVZ= RandN(sigma_LMC, NN); 
   SVZ +=  v_Z_lmc + ((double)rand()/(double)(RAND_MAX+1.)*2.0-1.0)*err_zlmc;
   SVR +=  v_R_lmc + ((double)rand()/(double)(RAND_MAX+1.)*2.0-1.0)*err_rlmc;
   SVT +=  v_T_lmc + ((double)rand()/(double)(RAND_MAX+1.)*2.0-1.0)*err_tlmc;}
   s.vs=sqrt(SVT*SVT+SVZ*SVZ+SVR*SVR);

///======================================================================================
/// Generate Lens velocity components in Glactocenteric cylindrical coordinates(x',y'
    LVR=LVT=LVZ=0.0;
    if(l.struc==0){///Galactic disk
    LVR= RandN(sigma_R_Disk, NN); 
    LVT= RandN(sigma_T_Disk, NN);
    LVZ= RandN(sigma_Z_Disk, NN); 
    LVT =LVT+ vro_sun *(1.00762 * pow(Rlc/Dsun,0.0394) + 0.00712);}

   else if(l.struc==1){///Galactic bulge
   LVZ=RandN(sigma_Z_Bulge, NN); 
   LVR=RandN(sigma_R_Bulge, NN);
   LVT=RandN(sigma_T_Bulge, NN);}

   else if(l.struc==2){///thick disk
   LVR= RandN(sigma_R_TDisk, NN); 
   LVT= RandN(sigma_T_TDisk, NN); 
   LVZ= RandN(sigma_Z_TDisk, NN);
   LVT =LVT+ vro_sun *(1.00762*pow(Rlc/Dsun,0.0394) + 0.00712); }
   
   else if(l.struc==3){///stellar halo
   LVR= RandN(sigma_R_halo, NN); 
   LVT= RandN(sigma_T_halo, NN); 
   LVZ= RandN(sigma_Z_halo, NN);}
   
   else if(l.struc>3){
   LVR= RandN(sigma_LMC, NN); 
   LVT= RandN(sigma_LMC, NN); 
   LVZ= RandN(sigma_LMC, NN); 
   LVZ +=  v_Z_lmc + ((double)rand()/(double)(RAND_MAX+1.)*2.0-1.0)*err_zlmc;
   LVR +=  v_R_lmc + ((double)rand()/(double)(RAND_MAX+1.)*2.0-1.0)*err_rlmc;
   LVT +=  v_T_lmc + ((double)rand()/(double)(RAND_MAX+1.)*2.0-1.0)*err_tlmc;}
   l.vl=sqrt(LVT*LVT+LVZ*LVZ+LVR*LVR);
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH


   if(fabs(l.Dl*cos(s.FI)*sin(s.TET)/Rlc-1.0)<0.05) betal=pi/2.0; 
   else if(fabs(l.Dl*cos(s.FI)*sin(s.TET)/Rlc+1.0)<0.05) betal=-pi/2.0; 
   else  betal=asin(l.Dl*cos(s.FI)*sin(s.TET)/Rlc);///lens[-pi/2,pi/2]
   
   if(fabs(s.Ds*cos(s.FI)*sin(s.TET)/Rsc-1.0)<0.05) betas=pi/2.0; 
   else if(fabs(s.Ds*cos(s.FI)*sin(s.TET)/Rsc+1.0)<0.05) betas=-pi/2.0; 
   else  betas=asin(s.Ds*cos(s.FI)*sin(s.TET)/Rsc);///lens[-pi/2,pi/2]
   
   if(fabs(l.Dl*cos(s.FI)*sin(s.TET)/Rlc)>1.05 || fabs(s.Ds*cos(s.FI)*sin(s.TET)/Rsc)>1.05 || Rlc==0.0 || Rsc==0.0){
   cout<<"ERROR Dl: "<<l.Dl<<"\t Ds: "<<s.Ds<<endl;
   cout<<"FI: "<<s.FI<<"\t TET: "<<s.TET<<endl; 
   cout<<"Rlc: "<<Rlc<<"\t Rsc: "<<Rsc<<endl;
   cout<<"sin(l): "<<l.Dl*cos(s.FI)*sin(s.TET)/Rlc<<"\t sin(s): "<<s.Ds*cos(s.FI)*sin(s.TET)/Rsc<<endl;
   int ew; cin>>ew; }



   
   if(fabs(l.Dl*cos(s.FI))>sqrt(pow(Dsun,2.0)+pow(l.Dl*cos(s.FI)*sin(s.TET),2.0)) ) betal= pi-betal;
   if(fabs(s.Ds*cos(s.FI))>sqrt(pow(Dsun,2.0)+pow(s.Ds*cos(s.FI)*sin(s.TET),2.0)) ) betas= pi-betas;

    if(fabs((Rlc-Dsun*cos(betal))/(l.Dl*cos(s.FI))-1.0)<0.05)   deltal=0.0; 
    else if (fabs((Rlc-Dsun*cos(betal))/(l.Dl*cos(s.FI))+1.0)<0.05) deltal=pi; 
    else    deltal=acos((Rlc-Dsun*cos(betal))/(l.Dl*cos(s.FI)));
    
    if(fabs((Rsc-Dsun*cos(betas))/(s.Ds*cos(s.FI))-1.0)<0.05)   deltas=0.0; 
    else if (fabs((Rsc-Dsun*cos(betas))/(s.Ds*cos(s.FI))+1.0)<0.05) deltas=pi; 
    else    deltas=acos((Rsc-Dsun*cos(betas))/(s.Ds*cos(s.FI)));
    
    if(fabs((Rlc-Dsun*cos(betal))/(l.Dl*cos(s.FI)))>1.05 or fabs((Rsc-Dsun*cos(betas))/(s.Ds*cos(s.FI)))>1.05 or fabs(s.FI)==pi/2.0){
    cout<<"ERROR Dl: "<<l.Dl<<"\t Ds: "<<s.Ds<<endl;
    cout<<"betal: "<<betal<<"\t betas: "<<betas<<endl; 
    cout<<"FI: "<<s.FI<<"\t TET: "<<s.TET<<endl; 
    cout<<"Rlc: "<<Rlc<<"\t Rsc: "<<Rsc<<endl;
    cout<<"cos(dl): "<<(Rlc-Dsun*cos(betal))/(l.Dl*cos(s.FI))<<endl;
    cout<<"cos(ds): "<<(Rsc-Dsun*cos(betas))/(s.Ds*cos(s.FI))<<endl;  int ew; cin>>ew; }


    deltao=pi/2.0;
    SVl=-SVR*    sin(deltas)+ SVT* cos(deltas);
    LVl=-LVR*    sin(deltal)+ LVT* cos(deltal);
    VSunl=-VSunR*sin(deltao)+VSunT*cos(deltao);

    SVt=  1.0*SVR*cos(deltas)+  SVT*sin(deltas);
    LVt=  1.0*LVR*cos(deltal)+  LVT*sin(deltal);
    VSunt=1.0*VSunR*cos(deltao)+VSunT*sin(deltao);

    SVb=-sin(s.FI)*(SVt) + cos(s.FI)* (SVZ);
    LVb=-sin(s.FI)*(LVt) + cos(s.FI)* (LVZ);
    VSunb=-sin(s.FI)*(VSunt)+cos(s.FI)*(VSunZ);

    vls_l= LVl-l.xls*SVl -(1.0-l.xls)*VSunl;
    vls_b= LVb-l.xls*SVb -(1.0-l.xls)*VSunb;
    l.Vt=sqrt(fabs(vls_l*vls_l+ vls_b*vls_b));
    if (l.Vt<0.0 || l.Vt>1.0e6 ){
    cout<<" Vt is very large: "<<l.Vt<<"\t vl: "<<l.vl<<"\t Vs: "<<l.vs<<endl;   int yee; cin>>yee;}
}
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
