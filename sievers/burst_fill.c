#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <omp.h>
#include <cerf.h>

#ifndef M_PI
#define M_PI 3.14159265358979
#endif

//gcc-4.8 -O3 -fPIC --shared --std=c99  burst_fill.c -o libburst_fill.so
//gcc-4.8 -O3 -fopenmp -fPIC --shared --std=c99  burst_fill.c -o libburst_fill.so
// gcc -O3 -fopenmp -fPIC --shared --std=c99  burst_fill.c -o libburst_fill.so 
// gcc -I/home/sievers/local/include -O3 -fopenmp -fPIC --shared --std=c99  burst_fill.c -o libburst_fill.so


#define REF_FREQ 764.2

//#define DM0 4150
#define DM0  4148.808   //from barsdell et al.
void fill_model_general(float dt, float fwhm, float tau, float scat, float alpha, float t0, float DM, float DM_ind, float *freq, int nfreq, int n, void *mat_in)
{

  //printf("dt, fwhm, scat, alpha, t0, and DM are %14.5g %14.5g %14.5g %14.5g %14.5g %14.5g\n",dt,fwhm,scat,alpha,t0,DM);

  //scat=scat*dt;
  int nn=1+n/2;
  
  float complex *mat=(float complex *)mat_in;

  float freq_ref=800.0;
  float sig=fwhm/sqrt(8*log(2))/dt;
  float pi=3.14159265;
  float *dm_lags=(float *)malloc(sizeof(float)*nfreq);
  float *sigvec=(float *)malloc(sizeof(float)*nn);
  float ref_lag=DM0*DM*pow(REF_FREQ,DM_ind);
  for (int i=0;i<nfreq;i++) {
    dm_lags[i]=DM0*DM*pow(freq[i],DM_ind)-ref_lag;
    dm_lags[i]=dm_lags[i]/dt;    
  }

  //float lag_min=dm_lags[0];
  //for (int i=1;i<nfreq;i++)
  //  if (dm_lags[i]<lag_min)
  //    lag_min=dm_lags[i];
  //for (int i=0;i<nfreq;i++)
  //  dm_lags[i]-=lag_min;


  for (int i=0;i<nn;i++)  {
    float fac=2*pi*sig*i/n;
    sigvec[i]=exp(-0.5*fac*fac);
  }
  //sigvec[i]=exp(-0.5dm_lags[i]=DM0*DM/(freq[i]*freq[i]);
  


  float complex *aa2exp=(float complex *)malloc(nn*sizeof(float complex));
  float complex *naa2exp=(float complex *)malloc(nn*sizeof(float complex));
  for (int i=0;i<nn;i++) {
    float complex aa2=-2*pi*i/n*I;
    aa2exp[i]=cexp(aa2);
    naa2exp[i]=cexp(aa2*n);
    
  }

#pragma omp parallel for
  for (int j=0;j<nfreq;j++) {

    
    float mylag=dm_lags[j]+t0;
    float fac=freq[j]/freq_ref;
    float myscat=scat*(fac*fac*fac*fac)+tau;
    
    float mynorm=(1-exp(-myscat*n))/(1-exp(-myscat))/pow( (freq[j]/freq_ref),alpha);
    
    //if (j==1)
    //printf("things are %12.5g %12.5g %12.5g\n",mylag,myscat,mynorm);

    for (int i=0;i<nn;i++) {


      float aa1=-myscat;
      float complex scatvec=(1-exp(n*aa1)*naa2exp[i])/(1-exp(aa1)*aa2exp[i]);
      scatvec=scatvec/mynorm;
      float complex lagvec=cexp(-2*pi*I*i*mylag/n);

      //if ((i==1)&&(j==1)) {
	//printf("aa is %12.5g %12.5g\n",creal(aa),cimag(aa));
	//printf("scatvec is %12.5g %12.5g\n",creal(scatvec),cimag(scatvec));
	//printf("lagvec is %12.5g %12.5g\n",creal(lagvec),cimag(lagvec));
      //}
      mat[j*nn+i]=sigvec[i]*scatvec*lagvec;
    }
  }

  //printf("second element is %12.6g %12.6g\n",creal(mat[1]),cimag(mat[1]));

  //printf("hello from C2!\n");
  //mat[2]=4+5I;
  free(sigvec);
  free(dm_lags);
  free(aa2exp);
  free(naa2exp);

}


/*--------------------------------------------------------------------------------*/
void fill_model(float dt, float fwhm, float scat, float alpha, float t0, float DM, float *freq, int nfreq, int n, void *mat_in)
{

  //printf("dt, fwhm, scat, alpha, t0, and DM are %14.5g %14.5g %14.5g %14.5g %14.5g %14.5g\n",dt,fwhm,scat,alpha,t0,DM);

  //scat=scat*dt;
  int nn=1+n/2;
  
  float complex *mat=(float complex *)mat_in;
  //float DM0=4150;
  float freq_ref=800.0;
  float sig=fwhm/sqrt(8*log(2))/dt;
  float pi=3.14159265;
  float *dm_lags=(float *)malloc(sizeof(float)*nfreq);
  float *sigvec=(float *)malloc(sizeof(float)*nn);
  for (int i=0;i<nfreq;i++) {
    dm_lags[i]=DM0*DM/(freq[i]*freq[i]);
    dm_lags[i]=dm_lags[i]/dt;
  }
  float lag_min=dm_lags[0];
  for (int i=1;i<nfreq;i++)
    if (dm_lags[i]<lag_min)
      lag_min=dm_lags[i];
  for (int i=0;i<nfreq;i++)
    dm_lags[i]-=lag_min;


  for (int i=0;i<nn;i++)  {
    float fac=2*pi*sig*i/n;
    sigvec[i]=exp(-0.5*fac*fac);
  }
  //sigvec[i]=exp(-0.5dm_lags[i]=DM0*DM/(freq[i]*freq[i]);
  


  float complex *aa2exp=(float complex *)malloc(nn*sizeof(float complex));
  float complex *naa2exp=(float complex *)malloc(nn*sizeof(float complex));
  for (int i=0;i<nn;i++) {
    float complex aa2=-2*pi*i/n*I;
    aa2exp[i]=cexp(aa2);
    naa2exp[i]=cexp(aa2*n);
    
  }

#pragma omp parallel for
  for (int j=0;j<nfreq;j++) {

    
    float mylag=dm_lags[j]+t0;
    float fac=freq[j]/freq_ref;
    float myscat=scat*(fac*fac*fac*fac);
    
    float mynorm=(1-exp(-myscat*n))/(1-exp(-myscat))/pow( (freq[j]/freq_ref),alpha);
    
    //if (j==1)
    //printf("things are %12.5g %12.5g %12.5g\n",mylag,myscat,mynorm);

    for (int i=0;i<nn;i++) {


      float aa1=-myscat;
      //float complex aa2=-2*pi*i/n*I;
      //float complex aa=aa1+aa2;
      //float complex aa=-myscat-2*pi*i/n*I;

      //float complex scatvec=(1-cexp(n*aa))/(1-cexp(aa));
      //float complex scatvec=(1-exp(n*aa1)*cexp(n*aa2))/(1-exp(aa1)*cexp(aa2));
      float complex scatvec=(1-exp(n*aa1)*naa2exp[i])/(1-exp(aa1)*aa2exp[i]);
      scatvec=scatvec/mynorm;
      float complex lagvec=cexp(-2*pi*I*i*mylag/n);

      //if ((i==1)&&(j==1)) {
	//printf("aa is %12.5g %12.5g\n",creal(aa),cimag(aa));
	//printf("scatvec is %12.5g %12.5g\n",creal(scatvec),cimag(scatvec));
	//printf("lagvec is %12.5g %12.5g\n",creal(lagvec),cimag(lagvec));
      //}
      mat[j*nn+i]=sigvec[i]*scatvec*lagvec;
    }
  }

  //printf("second element is %12.6g %12.6g\n",creal(mat[1]),cimag(mat[1]));

  //printf("hello from C2!\n");
  //mat[2]=4+5I;
  free(sigvec);
  free(dm_lags);
  free(aa2exp);
  free(naa2exp);

}


/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
void fill_model_double(float dt, float fwhm, float scat, float alpha, float t0, float DM, float *freq, int nfreq, int n, void *mat_in)
{

  //printf("dt, fwhm, scat, alpha, t0, and DM are %14.5g %14.5g %14.5g %14.5g %14.5g %14.5g\n",dt,fwhm,scat,alpha,t0,DM);

  //scat=scat*dt;
  int nn=1+n/2;
  
  float complex *mat=(float complex *)mat_in;
  //double DM0=4150;
  double freq_ref=800.0;
  double sig=fwhm/sqrt(8*log(2))/dt;
  double pi=3.14159265358979;
  double *dm_lags=(double *)malloc(sizeof(double)*nfreq);
  double *sigvec=(double *)malloc(sizeof(double)*nn);
  for (int i=0;i<nfreq;i++) {
    dm_lags[i]=DM0*DM/(freq[i]*freq[i]);
    dm_lags[i]=dm_lags[i]/dt;
  }
  double lag_min=dm_lags[0];
  for (int i=1;i<nfreq;i++)
    if (dm_lags[i]<lag_min)
      lag_min=dm_lags[i];
  for (int i=0;i<nfreq;i++)
    dm_lags[i]-=lag_min;


  for (int i=0;i<nn;i++)  {
    double fac=2*pi*sig*i/n;
    sigvec[i]=exp(-0.5*fac*fac);
  }
  //sigvec[i]=exp(-0.5dm_lags[i]=DM0*DM/(freq[i]*freq[i]);
  


  double complex *aa2exp=(double complex *)malloc(nn*sizeof(double complex));
  double complex *naa2exp=(double complex *)malloc(nn*sizeof(double complex));
  for (int i=0;i<nn;i++) {
    double complex aa2=-2*pi*i/n*I;
    aa2exp[i]=cexp(aa2);
    naa2exp[i]=cexp(aa2*n);
    
  }

#pragma omp parallel for
  for (int j=0;j<nfreq;j++) {

    
    double mylag=dm_lags[j]+t0;
    double fac=freq[j]/freq_ref;
    double myscat=scat*(fac*fac*fac*fac);
    
    double mynorm=(1-exp(-myscat*n))/(1-exp(-myscat))/pow( (freq[j]/freq_ref),alpha);
    
    //if (j==1)
    //printf("things are %12.5g %12.5g %12.5g\n",mylag,myscat,mynorm);

    for (int i=0;i<nn;i++) {


      double aa1=-myscat;
      //float complex aa2=-2*pi*i/n*I;
      //float complex aa=aa1+aa2;
      //float complex aa=-myscat-2*pi*i/n*I;

      //float complex scatvec=(1-cexp(n*aa))/(1-cexp(aa));
      //float complex scatvec=(1-exp(n*aa1)*cexp(n*aa2))/(1-exp(aa1)*cexp(aa2));
      double complex scatvec=(1-exp(n*aa1)*naa2exp[i])/(1-exp(aa1)*aa2exp[i]);
      scatvec=scatvec/mynorm;
      double complex lagvec=cexp(-2*pi*I*i*mylag/n);

      //if ((i==1)&&(j==1)) {
	//printf("aa is %12.5g %12.5g\n",creal(aa),cimag(aa));
	//printf("scatvec is %12.5g %12.5g\n",creal(scatvec),cimag(scatvec));
	//printf("lagvec is %12.5g %12.5g\n",creal(lagvec),cimag(lagvec));
      //}
      mat[j*nn+i]=sigvec[i]*scatvec*lagvec;
    }
  }

  //printf("second element is %12.6g %12.6g\n",creal(mat[1]),cimag(mat[1]));

  //printf("hello from C2!\n");
  //mat[2]=4+5I;
  free(sigvec);
  free(dm_lags);
  free(aa2exp);
  free(naa2exp);

}


/*--------------------------------------------------------------------------------*/
void fill_model_fast(float dt, float fwhm, float scat, float alpha, float t0, float DM, float *freq, int nfreq, int n, void *mat_in)
{

  //printf("dt, fwhm, scat, alpha, t0, and DM are %14.5g %14.5g %14.5g %14.5g %14.5g %14.5g\n",dt,fwhm,scat,alpha,t0,DM);

  //scat=scat*dt;
  int nn=1+n/2;
  
  float complex *mat=(float complex *)mat_in;

  //double DM0=4150;
  double freq_ref=800.0;
  double sig=fwhm/sqrt(8*log(2))/dt;
  double pi=3.14159265;
  double *dm_lags=(double *)malloc(sizeof(double)*nfreq);
  double *sigvec=(double *)malloc(sizeof(double)*nn);
  for (int i=0;i<nfreq;i++) {
    dm_lags[i]=DM0*DM/(freq[i]*freq[i]);
    dm_lags[i]=dm_lags[i]/dt;
  }
  double lag_min=dm_lags[0];
  for (int i=1;i<nfreq;i++)
    if (dm_lags[i]<lag_min)
      lag_min=dm_lags[i];
  for (int i=0;i<nfreq;i++)
    dm_lags[i]-=lag_min;


  for (int i=0;i<nn;i++)  {
    double fac=2*pi*sig*i/n;
    sigvec[i]=exp(-0.5*fac*fac);
  }
  //sigvec[i]=exp(-0.5dm_lags[i]=DM0*DM/(freq[i]*freq[i]);
  


  double complex *aa2exp=(double complex *)malloc(nn*sizeof(double complex));
  double complex *naa2exp=(double complex *)malloc(nn*sizeof(double complex));
  for (int i=0;i<nn;i++) {
    double complex aa2=-2*pi*i/n*I;
    aa2exp[i]=cexp(aa2);
    naa2exp[i]=cexp(aa2*n);
    
  }

#pragma omp parallel for
  for (int j=0;j<nfreq;j++) {

    
    double mylag=dm_lags[j]+t0;
    double fac=freq[j]/freq_ref;
    double myscat=scat*(fac*fac*fac*fac);
    
    double mynorm=(1-exp(-myscat*n))/(1-exp(-myscat))/pow( (freq[j]/freq_ref),alpha);
    
    //if (j==1)
    //printf("things are %12.5g %12.5g %12.5g\n",mylag,myscat,mynorm);
    double aa1=-myscat;
    double aa1exp=exp(aa1);
    double naa1exp=exp(n*aa1);
    double complex lagfac=cexp(-2*pi*I*mylag/n);
    double complex lagvec=1.0;
    for (int i=0;i<nn;i++) {
      naa1exp*=aa1exp;


      double complex scatvec=(1-naa1exp*naa2exp[i])/(1-aa1exp*aa2exp[i]);
      scatvec=scatvec/mynorm;
      //float complex lagvec=cexp(-2*pi*I*i*mylag/n);      
      mat[j*nn+i]=sigvec[i]*scatvec*lagvec;
      lagvec*=lagfac;
    }
  }

  //printf("second element is %12.6g %12.6g\n",creal(mat[1]),cimag(mat[1]));

  //printf("hello from C2!\n");
  //mat[2]=4+5I;
  free(sigvec);
  free(dm_lags);
  free(aa2exp);
  free(naa2exp);

}


/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
void fill_model_fast_scatfit(float dt, float fwhm, float scat, float scatpow, float alpha, float t0, float DM, float *freq, int nfreq, int n, void *mat_in)
{

  //printf("dt, fwhm, scat, alpha, t0, and DM are %14.5g %14.5g %14.5g %14.5g %14.5g %14.5g\n",dt,fwhm,scat,alpha,t0,DM);

  //scat=scat*dt;
  int nn=1+n/2;
  
  float complex *mat=(float complex *)mat_in;

  //double DM0=4150;
  double freq_ref=800.0;
  double sig=fwhm/sqrt(8*log(2))/dt;
  double pi=3.14159265;
  double *dm_lags=(double *)malloc(sizeof(double)*nfreq);
  double *sigvec=(double *)malloc(sizeof(double)*nn);
  for (int i=0;i<nfreq;i++) {
    dm_lags[i]=DM0*DM/(freq[i]*freq[i]);
    dm_lags[i]=dm_lags[i]/dt;
  }
  double lag_min=dm_lags[0];
  for (int i=1;i<nfreq;i++)
    if (dm_lags[i]<lag_min)
      lag_min=dm_lags[i];
  for (int i=0;i<nfreq;i++)
    dm_lags[i]-=lag_min;


  for (int i=0;i<nn;i++)  {
    double fac=2*pi*sig*i/n;
    sigvec[i]=exp(-0.5*fac*fac);
  }
  //sigvec[i]=exp(-0.5dm_lags[i]=DM0*DM/(freq[i]*freq[i]);
  


  double complex *aa2exp=(double complex *)malloc(nn*sizeof(double complex));
  double complex *naa2exp=(double complex *)malloc(nn*sizeof(double complex));
  for (int i=0;i<nn;i++) {
    double complex aa2=-2*pi*i/n*I;
    aa2exp[i]=cexp(aa2);
    naa2exp[i]=cexp(aa2*n);
    
  }

#pragma omp parallel for
  for (int j=0;j<nfreq;j++) {

    
    double mylag=dm_lags[j]+t0;
    double fac=freq[j]/freq_ref;
    double myscat=scat*pow(fac,scatpow);
    
    double mynorm=(1-exp(-myscat*n))/(1-exp(-myscat))/pow( (freq[j]/freq_ref),alpha);
    
    //if (j==1)
    //printf("things are %12.5g %12.5g %12.5g\n",mylag,myscat,mynorm);
    double aa1=-myscat;
    double aa1exp=exp(aa1);
    double naa1exp=exp(n*aa1);
    double complex lagfac=cexp(-2*pi*I*mylag/n);
    double complex lagvec=1.0;
    for (int i=0;i<nn;i++) {
      naa1exp*=aa1exp;


      double complex scatvec=(1-naa1exp*naa2exp[i])/(1-aa1exp*aa2exp[i]);
      scatvec=scatvec/mynorm;
      //float complex lagvec=cexp(-2*pi*I*i*mylag/n);      
      mat[j*nn+i]=sigvec[i]*scatvec*lagvec;
      lagvec*=lagfac;
    }
  }

  //printf("second element is %12.6g %12.6g\n",creal(mat[1]),cimag(mat[1]));

  //printf("hello from C2!\n");
  //mat[2]=4+5I;
  free(sigvec);
  free(dm_lags);
  free(aa2exp);
  free(naa2exp);

}

/*--------------------------------------------------------------------------------*/
static inline float complex cerf_func(double omega, double sig, double tau, double t)
{
  complex double myshift=I*omega*sig*sig;
  complex double erfarg=(sig*sig/tau-(t-myshift))/sqrt(2)/sig;
  complex double exparg=-(t-myshift)/tau;
  double mynorm=exp(-0.5*sig*sig*(omega*omega-1/tau/tau))/2/tau;
  return mynorm*cexp(exparg)*cerfc(erfarg);
}
/*--------------------------------------------------------------------------------*/

void fill_model_real_qu_swing(double dt, double fwhm, double scat, double alpha, double t0, double DM, double omega, double *freq, int nfreq, int n, void *mat_in, int *imin, int *imax)
{

  double tstart=omp_get_wtime();
  float complex *mat=(float complex *)mat_in;
  //printf("input parameters are %10.5f %12.7f %12.7f %8.4f %12.5f %12.4f\n",dt,fwhm,scat,alpha,t0,DM);
  double scat_use=scat/dt;
  double sig=fwhm/sqrt(8*log(2))/dt;
  int sig_width=10;
  int scat_width=20;
  double rt2=sqrt(2.0);


#pragma omp parallel for
  for (int ifreq=0;ifreq<nfreq;ifreq++) {

    double amp=pow(freq[ifreq]/REF_FREQ,alpha);
    double dm_lag=DM0*DM*(1.0/freq[ifreq]/freq[ifreq]-1.0/REF_FREQ/REF_FREQ);
    double myscat=scat_use*pow(freq[ifreq]/REF_FREQ,-4);
    double myarr=(t0+dm_lag)/dt;
    int istart=myarr-sig*sig_width-1;
    int istop=myarr+sig*sig_width+myscat*scat_width+1;
    if (istart<0)
      istart=0;
    if (istop>=n)
      istop=n-1;
    
    //if we want to calculate a fast chi^2, let the calling function know what the non-zero 
    //part of the function are
    if (imin)
      imin[ifreq]=istart;
    if (imax)
      imax[ifreq]=istop;
    //printf("ifreq is %d of %d %d %d\n",ifreq,nfreq,imin[ifreq],imax[ifreq]);

    double omega_use=omega;

    //double norm=0.5/myscat*exp(0.5*sig*sig/myscat/myscat)*amp;
    double erfpart=sig/rt2/myscat;
    double complex myshift=I*omega*sig*sig;
    for (int i=istart;i<=istop;i++) {
      double t=i-myarr;

      float complex f0=cerf_func(omega_use,sig,myscat,t);
      float complex f1=cerf_func(omega_use,sig,myscat,t+0.25);
      float complex f2=cerf_func(omega_use,sig,myscat,t+0.5);
      float complex f3=cerf_func(omega_use,sig,myscat,t+0.75);
      float complex f4=cerf_func(omega_use,sig,myscat,t+1.0);
      //mat[ifreq*n+i]=norm*(7*f0+32*f1+12*f2+32*f3+7*f4)/90.0;
      mat[ifreq*n+i]=amp*(7*f0+32*f1+12*f2+32*f3+7*f4)/90.0;
    }

  }        
  //printf("took %12.5f seconds\n",omp_get_wtime()-tstart);

}



/*--------------------------------------------------------------------------------*/

void fill_model_real(double dt, double fwhm, double scat, double alpha, double t0, double DM, double *freq, int nfreq, int n, float *mat_in, int *imin, int *imax)
{
  //printf("input parameters are %10.5f %12.7f %12.7f %8.4f %12.5f %12.4f\n",dt,fwhm,scat,alpha,t0,DM);
  double scat_use=scat/dt;
  double sig=fwhm/sqrt(8*log(2))/dt;
  int sig_width=10;
  int scat_width=20;
  double rt2=sqrt(2.0);


#pragma omp parallel for
  for (int ifreq=0;ifreq<nfreq;ifreq++) {
    double amp=pow(freq[ifreq]/REF_FREQ,alpha);
    double dm_lag=DM0*DM*(1.0/freq[ifreq]/freq[ifreq]-1.0/REF_FREQ/REF_FREQ);
    double myscat=scat_use*pow(freq[ifreq]/REF_FREQ,-4);
    double myarr=(t0+dm_lag)/dt;
    int istart=myarr-sig*sig_width-1;
    int istop=myarr+sig*sig_width+myscat*scat_width+1;
    if (istart<0)
      istart=0;
    if (istop>=n)
      istop=n-1;
    
    //if we want to calculate a fast chi^2, let the calling function know what the non-zero 
    //part of the function are
    if (imin)
      imin[ifreq]=istart;
    if (imax)
      imax[ifreq]=istop;

    double norm=0.5/myscat*exp(0.5*sig*sig/myscat/myscat)*amp;
    double erfpart=sig/rt2/myscat;
    for (int i=istart;i<=istop;i++) {
      double t=i-myarr;
#if 1
      if (myscat>0.05*sig) {
	double f0=erfc(erfpart-t/(rt2*sig))*exp(-t/myscat);
	double t1=t+0.25;double f1=erfc(erfpart-t1/(rt2*sig))*exp(-t1/myscat);
	double t2=t+0.5;double f2=erfc(erfpart-t2/(rt2*sig))*exp(-t2/myscat);
	double t3=t+0.75;double f3=erfc(erfpart-t3/(rt2*sig))*exp(-t3/myscat);
	double t4=t+1.0;double f4=erfc(erfpart-t4/(rt2*sig))*exp(-t4/myscat);
	mat_in[ifreq*n+i]=norm*(7*f0+32*f1+12*f2+32*f3+7*f4)/90.0;
      }
      else {
	t=t-myscat;
	double f0=exp(-t*t/(2*sig*sig));
	double t1=t+0.25;double f1=exp(-t1*t1/(2*sig*sig));
	double t2=t+0.5;double f2=exp(-t2*t2/(2*sig*sig));
	double t3=t+0.75;double f3=exp(-t3*t3/(2*sig*sig));
	double t4=t+1.0;double f4=exp(-t4*t4/(2*sig*sig));
	norm=amp/sqrt(2*3.141592653589793)/sig;
	mat_in[ifreq*n+i]=norm*(7*f0+32*f1+12*f2+32*f3+7*f4)/90.0;
	
      }
#else
      //if ((ifreq==0)&&(i==istart))
      //	printf("t is %12.5f%12.5f\n",t,myarr);
      double erfarg=erfpart-t/(rt2*sig);
      
      mat_in[ifreq*n+i]=norm*exp(-t/myscat)*erfc(erfarg);
#endif
    }        
  }
}


/*--------------------------------------------------------------------------------*/

void fill_model_oversamp(double dt, double fwhm, double scat, double alpha, double t0, double DM, double *freq_in, int nfreq_in, int n, float *mat_in, int *imin, int *imax)
{

  int nfreq=2*nfreq_in+1;
  double *freq=(double *)malloc(nfreq*sizeof(double));
  double dnu=(freq_in[nfreq_in-1]-freq_in[0])/(nfreq_in-1)/2.0;
  freq[0]=freq_in[0]-dnu;
  //do it this way so if frequency sampling is slightly uneven/roundoff error is an issue
  //the frequencies stay locked to the original ones.
  for (int i=0;i<nfreq_in;i++) {
    freq[2*i]=freq_in[i];
    freq[2*i+1]=freq_in[i]+dnu;
  }

  
  //printf("input parameters are %10.5f %12.7f %12.7f %8.4f %12.5f %12.4f\n",dt,fwhm,scat,alpha,t0,DM);
  double scat_use=scat/dt;
  double sig=fwhm/sqrt(8*log(2))/dt;
  int sig_width=10;
  int scat_width=20;
  double rt2=sqrt(2.0);


#pragma omp parallel for
  for (int ifreq=0;ifreq<nfreq;ifreq++) {  

    double amp=pow(freq[ifreq]/REF_FREQ,alpha);
    double dm_lag=DM0*DM*(1.0/freq[ifreq]/freq[ifreq]-1.0/REF_FREQ/REF_FREQ);
    double myscat=scat_use*pow(freq[ifreq]/REF_FREQ,-4);
    double myarr=(t0+dm_lag)/dt;
    int istart=myarr-sig*sig_width-1;
    int istop=myarr+sig*sig_width+myscat*scat_width+1;
    if (istart<0)
      istart=0;
    if (istop>=n)
      istop=n-1;
    
    //if we want to calculate a fast chi^2, let the calling function know what the non-zero 
    //part of the function are
    if (imin)
      imin[ifreq]=istart;
    if (imax)
      imax[ifreq]=istop;

    double norm=0.5/myscat*exp(0.5*sig*sig/myscat/myscat)*amp;
    double erfpart=sig/rt2/myscat;
    //mat_in[ifreq*n+istart-1]=0;
    //mat_in[ifreq*n+istop+1]=0;
    for (int i=istart;i<=istop;i++) {
      double t=i-myarr;
#if 1
      if (myscat>0.05*sig) {
	double f0=erfc(erfpart-t/(rt2*sig))*exp(-t/myscat);
	double t1=t+0.25;double f1=erfc(erfpart-t1/(rt2*sig))*exp(-t1/myscat);
	double t2=t+0.5;double f2=erfc(erfpart-t2/(rt2*sig))*exp(-t2/myscat);
	double t3=t+0.75;double f3=erfc(erfpart-t3/(rt2*sig))*exp(-t3/myscat);
	double t4=t+1.0;double f4=erfc(erfpart-t4/(rt2*sig))*exp(-t4/myscat);
	mat_in[ifreq*n+i]=norm*(7*f0+32*f1+12*f2+32*f3+7*f4)/90.0;
      }
      else {
	t=t-myscat;
	double f0=exp(-t*t/(2*sig*sig));
	double t1=t+0.25;double f1=exp(-t1*t1/(2*sig*sig));
	double t2=t+0.5;double f2=exp(-t2*t2/(2*sig*sig));
	double t3=t+0.75;double f3=exp(-t3*t3/(2*sig*sig));
	double t4=t+1.0;double f4=exp(-t4*t4/(2*sig*sig));
	norm=amp/sqrt(2*3.141592653589793)/sig;
	mat_in[ifreq*n+i]=norm*(7*f0+32*f1+12*f2+32*f3+7*f4)/90.0;
	
      }
#else
      //if ((ifreq==0)&&(i==istart))
      //	printf("t is %12.5f%12.5f\n",t,myarr);
      double erfarg=erfpart-t/(rt2*sig);
      
      mat_in[ifreq*n+i]=norm*exp(-t/myscat)*erfc(erfarg);
#endif
    }        
  }

  if (1) //do hamming integration over frequency window to second order
    {
      //be careful parallelizing this as it could lead to severe breakage!

      double alpha=0.46*8/M_PI/M_PI + 0.54*2/3;
      double beta=0.23*(1-8/M_PI/M_PI)+0.54/6;

      if (fabs(alpha+2*beta-1)>1e-14)
	fprintf(stderr,"Coeffs don't match up in fill_model_real.\n");
      for (int i=0;i<(nfreq-1)/2;i++) {
	int ii=2*i+1;
	for (int j=imin[ii];j<=imax[ii];j++) {
	  double tmp=alpha*mat_in[ii*n+j]+beta*(mat_in[(ii+1)*n+j]+mat_in[(ii-1)*n+j]);
	  mat_in[i*n+j]=tmp;
	  //mat_in[i*n+j]=mat_in[ii*n+j];
	}
	imin[i]=imin[ii];
	imax[i]=imax[ii];
      }
    }  
  free(freq);
}


/*--------------------------------------------------------------------------------*/

void fill_model_real_general(double dt, double fwhm, double scat, double scatpow, double alpha, double t0, double DM, double DM_pow, double *freq, int nfreq, int n, float *mat_in, int *imin, int *imax)
{
  //printf("input parameters are %10.5f %12.7f %12.7f %8.4f %12.5f %12.4f\n",dt,fwhm,scat,alpha,t0,DM);
  double scat_use=scat/dt;
  double sig=fwhm/sqrt(8*log(2))/dt;
  int sig_width=10;
  int scat_width=20;
  double rt2=sqrt(2.0);
#pragma omp parallel for
  for (int ifreq=0;ifreq<nfreq;ifreq++) {
    double amp=pow(freq[ifreq]/REF_FREQ,alpha);
    double dm_lag=DM0*DM*(pow(freq[ifreq],DM_pow)-pow(REF_FREQ,DM_pow));
    
    double myscat=scat_use*pow(freq[ifreq]/REF_FREQ,scatpow);
    double myarr=(t0+dm_lag)/dt;
    int istart=myarr-sig*sig_width-1;
    int istop=myarr+sig*sig_width+myscat*scat_width+1;
    if (istart<0)
      istart=0;
    if (istop>=n)
      istop=n-1;
    
    //if we want to calculate a fast chi^2, let the calling function know what the non-zero 
    //part of the function are
    if (imin)
      imin[ifreq]=istart;
    if (imax)
      imax[ifreq]=istop;

    double norm=0.5/myscat*exp(0.5*sig*sig/myscat/myscat)*amp;
    double erfpart=sig/rt2/myscat;
    for (int i=istart;i<=istop;i++) {
      double t=i-myarr;
#if 1
      double f0=erfc(erfpart-t/(rt2*sig))*exp(-t/myscat);
      double t1=t+0.25;double f1=erfc(erfpart-t1/(rt2*sig))*exp(-t1/myscat);
      double t2=t+0.5;double f2=erfc(erfpart-t2/(rt2*sig))*exp(-t2/myscat);
      double t3=t+0.75;double f3=erfc(erfpart-t3/(rt2*sig))*exp(-t3/myscat);
      double t4=t+1.0;double f4=erfc(erfpart-t4/(rt2*sig))*exp(-t4/myscat);
      mat_in[ifreq*n+i]=norm*(7*f0+32*f1+12*f2+32*f3+7*f4)/90.0;
#else
      //if ((ifreq==0)&&(i==istart))
      //	printf("t is %12.5f%12.5f\n",t,myarr);
      double erfarg=erfpart-t/(rt2*sig);
      
      mat_in[ifreq*n+i]=norm*exp(-t/myscat)*erfc(erfarg);
#endif
    }        
  }
}

/*--------------------------------------------------------------------------------*/
double cold_plasma_delay(double nu, double nu_p_sqr)
{
  return 1.0/sqrt(1-nu_p_sqr/nu/nu)-1;
}
/*--------------------------------------------------------------------------------*/

void fill_model_real_general_nup(double dt, double fwhm, double scat, double scatpow, double alpha, double t0, double DM, double nu_p_sqr, double *freq, int nfreq, int n, float *mat_in, int *imin, int *imax)
{
  //printf("input parameters are %10.5f %12.7f %12.7f %8.4f %12.5f %12.4f\n",dt,fwhm,scat,alpha,t0,DM);
  double scat_use=scat/dt;
  double sig=fwhm/sqrt(8*log(2))/dt;
  int sig_width=10;
  int scat_width=20;
  double rt2=sqrt(2.0);
  double ref_lag=cold_plasma_delay(REF_FREQ,nu_p_sqr);
  double ref_lag2=cold_plasma_delay(700,nu_p_sqr);
  
  double normal_lag=DM0*DM/REF_FREQ/REF_FREQ;
  double normal_lag2=DM0*DM/700.0/700.0;
  double fac=(normal_lag-normal_lag2)/(ref_lag-ref_lag2);
  
#pragma omp parallel for
  for (int ifreq=0;ifreq<nfreq;ifreq++) {
    double amp=pow(freq[ifreq]/REF_FREQ,alpha);
    //double dm_lag=DM0*DM*(cold_plasma_delay(freq[ifreq],nu_p_sqr)-cold_plasma_delay(REF_FREQ,nu_p_sqr));
    double dm_lag=fac*(cold_plasma_delay(freq[ifreq],nu_p_sqr)-ref_lag);

    
    double myscat=scat_use*pow(freq[ifreq]/REF_FREQ,scatpow);
    double myarr=(t0+dm_lag)/dt;
    int istart=myarr-sig*sig_width-1;
    int istop=myarr+sig*sig_width+myscat*scat_width+1;
    if (istart<0)
      istart=0;
    if (istop>=n)
      istop=n-1;
    
    //if we want to calculate a fast chi^2, let the calling function know what the non-zero 
    //part of the function are
    if (imin)
      imin[ifreq]=istart;
    if (imax)
      imax[ifreq]=istop;

    double norm=0.5/myscat*exp(0.5*sig*sig/myscat/myscat)*amp;
    double erfpart=sig/rt2/myscat;
    for (int i=istart;i<=istop;i++) {
      double t=i-myarr;
#if 1
      double f0=erfc(erfpart-t/(rt2*sig))*exp(-t/myscat);
      double t1=t+0.25;double f1=erfc(erfpart-t1/(rt2*sig))*exp(-t1/myscat);
      double t2=t+0.5;double f2=erfc(erfpart-t2/(rt2*sig))*exp(-t2/myscat);
      double t3=t+0.75;double f3=erfc(erfpart-t3/(rt2*sig))*exp(-t3/myscat);
      double t4=t+1.0;double f4=erfc(erfpart-t4/(rt2*sig))*exp(-t4/myscat);
      mat_in[ifreq*n+i]=norm*(7*f0+32*f1+12*f2+32*f3+7*f4)/90.0;
#else
      //if ((ifreq==0)&&(i==istart))
      //	printf("t is %12.5f%12.5f\n",t,myarr);
      double erfarg=erfpart-t/(rt2*sig);
      
      mat_in[ifreq*n+i]=norm*exp(-t/myscat)*erfc(erfarg);
#endif
    }        
  }
}


/*--------------------------------------------------------------------------------*/

void fill_model_fast_general(float dt, float fwhm, float tau, float scat, float alpha, float t0, float DM, float DM_ind, float *freq, int nfreq, int n, void *mat_in)
{

  //printf("dt, fwhm, scat, alpha, t0, and DM are %14.5g %14.5g %14.5g %14.5g %14.5g %14.5g\n",dt,fwhm,scat,alpha,t0,DM);

  //scat=scat*dt;
  int nn=1+n/2;
  
  float complex *mat=(float complex *)mat_in;

  //double DM0=4150;
  double freq_ref=800.0;
  double sig=fwhm/sqrt(8*log(2))/dt;
  double pi=3.14159265;
  double *dm_lags=(double *)malloc(sizeof(double)*nfreq);
  double *sigvec=(double *)malloc(sizeof(double)*nn);
  for (int i=0;i<nfreq;i++) {
    //dm_lags[i]=DM0*DM/(freq[i]*freq[i]);
    dm_lags[i]=DM0*DM*pow(freq[i],DM_ind);
    dm_lags[i]=dm_lags[i]/dt;
  }
  double lag_min=dm_lags[0];
  for (int i=1;i<nfreq;i++)
    if (dm_lags[i]<lag_min)
      lag_min=dm_lags[i];
  for (int i=0;i<nfreq;i++)
    dm_lags[i]-=lag_min;


  for (int i=0;i<nn;i++)  {
    double fac=2*pi*sig*i/n;
    sigvec[i]=exp(-0.5*fac*fac);
  }
  //sigvec[i]=exp(-0.5dm_lags[i]=DM0*DM/(freq[i]*freq[i]);
  


  double complex *aa2exp=(double complex *)malloc(nn*sizeof(double complex));
  double complex *naa2exp=(double complex *)malloc(nn*sizeof(double complex));
  for (int i=0;i<nn;i++) {
    double complex aa2=-2*pi*i/n*I;
    aa2exp[i]=cexp(aa2);
    naa2exp[i]=cexp(aa2*n);
    
  }

#pragma omp parallel for
  for (int j=0;j<nfreq;j++) {

    
    double mylag=dm_lags[j]+t0;
    double fac=freq[j]/freq_ref;


    double tmp=scat*(fac*fac*fac*fac);
    
    double myscat=1.0/(1.0/(scat*(fac*fac*fac*fac))+1.0/tau);
    if (tau==0)
      myscat=scat*(fac*fac*fac*fac);
    //if (j==1)
    //  printf("scat,tau,tmp, and myscat are %12.7f %12.7f %12.4g %12.4g\n",scat,tau,tmp,myscat);

    double mynorm=(1-exp(-myscat*n))/(1-exp(-myscat))/pow( (freq[j]/freq_ref),alpha);
    
    //if (j==1)
    //printf("things are %12.5g %12.5g %12.5g\n",mylag,myscat,mynorm);
    double aa1=-myscat;
    double aa1exp=exp(aa1);
    double naa1exp=exp(n*aa1);
    double complex lagfac=cexp(-2*pi*I*mylag/n);
    double complex lagvec=1.0;
    for (int i=0;i<nn;i++) {
      naa1exp*=aa1exp;


      double complex scatvec=(1-naa1exp*naa2exp[i])/(1-aa1exp*aa2exp[i]);
      scatvec=scatvec/mynorm;
      //float complex lagvec=cexp(-2*pi*I*i*mylag/n);      
      mat[j*nn+i]=sigvec[i]*scatvec*lagvec;
      lagvec*=lagfac;
    }
  }

  //printf("second element is %12.6g %12.6g\n",creal(mat[1]),cimag(mat[1]));

  //printf("hello from C2!\n");
  //mat[2]=4+5I;
  free(sigvec);
  free(dm_lags);
  free(aa2exp);
  free(naa2exp);

}


/*--------------------------------------------------------------------------------*/

double calculate_chisq(void *dat_in, float *weights, int imin,float dt, float fwhm, float scat, float alpha, float amp, float t0, float DM, float *freq, int nfreq, int n)
{

  int nn=1+n/2;
  

  scat=scat*dt;
  t0=1+t0/dt;  //burst_model is expecting t0 to be in samples
  //float complex *mat=(float complex *)mat_in;
  float complex *mat=(float complex *)malloc(sizeof(float complex)*nn*nfreq);
  //fill_model(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,mat_in);
  double t1=omp_get_wtime();
  fill_model(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,(void *)mat);
  //printf("params are %12.5g %12.5g %12.5g %12.5g %12.5g %12.5g\n",dt,fwhm,scat,alpha,t0,DM);
  //printf("amps are %12.6f %12.6f %12.6f\n",creal(mat[0]),creal(mat[1]),cimag(mat[1]));
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  float complex *dat=(float complex *)dat_in;
  double *all_chisq=(double *)malloc(sizeof(double)*nfreq);
#pragma omp parallel for
  for (int j=0;j<nfreq;j++) {
    all_chisq[j]=0;
    for (int i=imin;i<nn;i++)  {
      double complex delt=amp*mat[j*nn+i]-dat[j*nn+i];
      all_chisq[j]+=(creal(delt)*creal(delt)+cimag(delt)*cimag(delt));
    }
    all_chisq[j]*=weights[j];
  }
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  double chisq=0;
  for (int i=0;i<nfreq;i++)
    chisq+=all_chisq[i];
  
  free(all_chisq);
  free(mat);
  return chisq;
}

/*--------------------------------------------------------------------------------*/

double calculate_chisq_cached(void *dat_in, float *weights, int imin,float dt, float *params, float *freq, int nfreq, int n, void *cached_in)
{

  float fwhm=params[0];
  float scat=1.0/params[1];
  float alpha=params[2];
  float amp=params[3];
  float t0=params[4];
  float DM=params[5];

  int nn=1+n/2;
  

  scat=scat*dt;
  t0=1+t0/dt;  //burst_model is expecting t0 to be in samples
  //float complex *mat=(float complex *)mat_in;
  float complex *mat=(float complex *)cached_in;
  //fill_model(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,mat_in);
  double t1=omp_get_wtime();
  fill_model_fast(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,(void *)mat);
  //printf("params are %12.5g %12.5g %12.5g %12.5g %12.5g %12.5g\n",dt,fwhm,scat,alpha,t0,DM);
  //printf("amps are %12.6f %12.6f %12.6f\n",creal(mat[0]),creal(mat[1]),cimag(mat[1]));
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  float complex *dat=(float complex *)dat_in;
  double *all_chisq=(double *)malloc(sizeof(double)*nfreq);
#pragma omp parallel for
  for (int j=0;j<nfreq;j++) {
    all_chisq[j]=0;
    for (int i=imin;i<nn;i++)  {
      double complex delt=amp*mat[j*nn+i]-dat[j*nn+i];
      all_chisq[j]+=(creal(delt)*creal(delt)+cimag(delt)*cimag(delt));
    }
    all_chisq[j]*=weights[j];
  }
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  double chisq=0;
  for (int i=0;i<nfreq;i++)
    chisq+=all_chisq[i];
  
  free(all_chisq);
  //free(mat);
  return chisq;
}

/*--------------------------------------------------------------------------------*/

double calculate_chisq_real_cached(void *dat_in, double *weights, double dt, double *params, double *freq, int nfreq, int n, void *cached_in)
{

  double fwhm=params[0];
  double scat=params[1];
  double alpha=params[2];
  double amp=params[3];
  double t0=params[4];
  double DM=params[5];

  if (scat<0)
    return 1e20;

  int *imin=(int *)malloc(sizeof(int)*(2*nfreq+1));
  int *imax=(int *)malloc(sizeof(int)*(2*nfreq+1));
  memset(imin,0,nfreq*sizeof(int));
  memset(imax,0,nfreq*sizeof(int));
  float *mat=(float *)cached_in;

  double t1=omp_get_wtime();
  fill_model_real(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,(void *)mat,imin,imax);
  //fill_model_oversamp(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,(void *)mat,imin,imax);

  float *dat=(float *)dat_in;
  double *all_chisq=(double *)malloc(sizeof(double)*nfreq);
  double *all_chisq_ref=(double *)malloc(sizeof(double)*nfreq);
#pragma omp parallel for
  for (int j=0;j<nfreq;j++) {
    all_chisq[j]=0;
    all_chisq_ref[j]=0;
    for (int i=imin[j];i<=imax[j];i++)  {
      double delt=amp*mat[j*n+i]-dat[j*n+i];
      all_chisq[j]+=delt*delt;
      all_chisq_ref[j]+=dat[j*n+i]*dat[j*n+i];
    }
    //all_chisq[j]*=weights[j];
    all_chisq[j]=weights[j]*(all_chisq[j]-all_chisq_ref[j]);

  }
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  double chisq=0;
  for (int i=0;i<nfreq;i++)
    chisq+=(all_chisq[i]);
  free(imin);
  free(imax);
  free(all_chisq);
  free(all_chisq_ref);
  //free(mat);
  return chisq;
}

/*--------------------------------------------------------------------------------*/

double calculate_chisq_real_cached_oversamp(void *dat_in, double *weights, double dt, double *params, double *freq, int nfreq, int n, void *cached_in)
{

  double fwhm=params[0];
  double scat=params[1];
  double alpha=params[2];
  double amp=params[3];
  double t0=params[4];
  double DM=params[5];

  if (scat<0)
    return 1e20;

  int *imin=(int *)malloc(sizeof(int)*(nfreq));
  int *imax=(int *)malloc(sizeof(int)*(nfreq));

  float *mat=(float *)cached_in;

  double t1=omp_get_wtime();
  fill_model_real(dt,fwhm,scat,alpha,t0,DM,freq,2*nfreq+1,n,(void *)mat,imin,imax);


  float *dat=(float *)dat_in;
  double *all_chisq=(double *)malloc(sizeof(double)*nfreq);
  double *all_chisq_ref=(double *)malloc(sizeof(double)*nfreq);
#pragma omp parallel for
  for (int j=0;j<nfreq;j++) {
    all_chisq[j]=0;
    all_chisq_ref[j]=0;
    for (int i=imin[j];i<=imax[j];i++)  {
      double delt=amp*mat[j*n+i]-dat[j*n+i];
      all_chisq[j]+=delt*delt;
      all_chisq_ref[j]+=dat[j*n+i]*dat[j*n+i];
    }
    //all_chisq[j]*=weights[j];
    all_chisq[j]=weights[j]*(all_chisq[j]-all_chisq_ref[j]);

  }
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  double chisq=0;
  for (int i=0;i<nfreq;i++)
    chisq+=(all_chisq[i]);
  free(imin);
  free(imax);
  free(all_chisq);
  free(all_chisq_ref);
  //free(mat);
  return chisq;
}

/*--------------------------------------------------------------------------------*/

double calculate_chisq_linfit_real_cached(void *dat_in, double *weights, double dt, double *params, double *freq, int nfreq, int n, void *cached_in)
{

  double fwhm=params[0];
  double scat=params[1];
  //double alpha=params[2];
  //double amp=params[3];
  double alpha=0.0;
  double slope=params[2];
  double offset=params[3];
  double t0=params[4];
  double DM=params[5];

  int *imin=(int *)malloc(sizeof(int)*(nfreq));
  int *imax=(int *)malloc(sizeof(int)*(nfreq));

  float *mat=(float *)cached_in;

  double t1=omp_get_wtime();
  fill_model_real(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,(void *)mat,imin,imax);

  float *dat=(float *)dat_in;
  double *all_chisq=(double *)malloc(sizeof(double)*nfreq);
  double *all_chisq_ref=(double *)malloc(sizeof(double)*nfreq);
#pragma omp parallel for
  for (int j=0;j<nfreq;j++) {
    all_chisq[j]=0;
    all_chisq_ref[j]=0;
    double nu_width=100.0;
    double amp=offset+slope*(freq[j]-REF_FREQ)/nu_width;
    for (int i=imin[j];i<=imax[j];i++)  {
      double delt=amp*mat[j*n+i]-dat[j*n+i];
      all_chisq[j]+=delt*delt;
      all_chisq_ref[j]+=dat[j*n+i]*dat[j*n+i];
    }
    //all_chisq[j]*=weights[j];
    all_chisq[j]=weights[j]*(all_chisq[j]-all_chisq_ref[j]);

  }
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  double chisq=0;
  for (int i=0;i<nfreq;i++)
    chisq+=(all_chisq[i]);
  free(imin);
  free(imax);
  free(all_chisq);
  free(all_chisq_ref);
  //free(mat);
  return chisq;
}

/*--------------------------------------------------------------------------------*/

double calculate_chisq_qu_real_cached(void *q_in, void *u_in,double *weights, double dt, double *params, double *freq, int nfreq, int n, void *cached_in)
{

  double fwhm=params[0];
  double scat=params[1];
  double alpha=params[2];
  double amp=params[3];
  double t0=params[4];
  double DM=params[5];
  double RM=params[6];
  double phi=params[7];


  if (scat<0)
    return 1e20;


  int *imin=(int *)malloc(sizeof(int)*(nfreq));
  int *imax=(int *)malloc(sizeof(int)*(nfreq));

  float *mat=(float *)cached_in;

  double t1=omp_get_wtime();
  fill_model_real(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,(void *)mat,imin,imax);

  float *q=(float *)q_in;
  float *u=(float *)u_in;
  double *all_chisq=(double *)malloc(sizeof(double)*nfreq);
  double *all_chisq_ref=(double *)malloc(sizeof(double)*nfreq);
#pragma omp parallel for
  for (int j=0;j<nfreq;j++) {
    all_chisq[j]=0;
    all_chisq_ref[j]=0;
    double lambda=299.79/freq[j];
    double lambda_ref=299.79/REF_FREQ;
    double fac=lambda_ref*lambda_ref;
    double mycos=cos(RM*(lambda*lambda-fac)+phi);
    double mysin=sin(RM*(lambda*lambda-fac)+phi);
    for (int i=imin[j];i<=imax[j];i++)  {
      double myq=mycos*q[j*n+i]+mysin*u[j*n+i];
      double myu=-mysin*q[j*n+i]+mycos*u[j*n+i];
      double delt=amp*mat[j*n+i]-myq;
      
      all_chisq[j]+=delt*delt+myu*myu;
      all_chisq_ref[j]+=myq*myq+myu*myu;
    }
    //all_chisq[j]*=weights[j];
    all_chisq[j]=weights[j]*(all_chisq[j]-all_chisq_ref[j]);

  }
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  double chisq=0;
  for (int i=0;i<nfreq;i++)
    chisq+=(all_chisq[i]);
  free(imin);
  free(imax);
  free(all_chisq);
  free(all_chisq_ref);
  //free(mat);
  return chisq;
}

/*--------------------------------------------------------------------------------*/

double calculate_chisq_qu_ramp_real_cached(void *q_in, void *u_in,double *weights, double dt, double *params, double *freq, int nfreq, int n, void *cached_in)
{

  double fwhm=params[0];
  double scat=params[1];
  double alpha=params[2];
  double amp=params[3];
  double t0=params[4];
  double DM=params[5];
  double RM=params[6];
  double phi=params[7];
  double phigrad=params[8];


  if (scat<0)
    return 1e20;


  int *imin=(int *)malloc(sizeof(int)*(nfreq));
  int *imax=(int *)malloc(sizeof(int)*(nfreq));

  float *mat=(float *)cached_in;

  double t1=omp_get_wtime();
  fill_model_real(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,(void *)mat,imin,imax);

  float *q=(float *)q_in;
  float *u=(float *)u_in;
  double *all_chisq=(double *)malloc(sizeof(double)*nfreq);
  double *all_chisq_ref=(double *)malloc(sizeof(double)*nfreq);
#pragma omp parallel for
  for (int j=0;j<nfreq;j++) {
    all_chisq[j]=0;
    all_chisq_ref[j]=0;
    double lambda=299.79/freq[j];
    double lambda_ref=299.79/REF_FREQ;
    double fac=lambda_ref*lambda_ref;
    double t_arr=t0+DM*DM0*(1/freq[j]/freq[j]-1/REF_FREQ/REF_FREQ);
    t_arr=t_arr/dt;  //do this in samples.  
    //if (j==0)
    // printf("t0 and t_arr are %12.5g %12.5g\n",t0,t_arr);
    for (int i=imin[j];i<=imax[j];i++)  {




      double mycos=cos(RM*(lambda*lambda-fac)+phigrad*(i-t_arr)+phi);
      double mysin=sin(RM*(lambda*lambda-fac)+phigrad*(i-t_arr)+phi);

      //if ((i==imin[j])&&(j==0))
      //if (j==0)
      //	printf("tarr and imin/imax are %12.6g %d %d, with terms %14.5g %14.5g\n",t_arr,imin[j],imax[j],mycos,mysin);


      double myq=mycos*q[j*n+i]+mysin*u[j*n+i];
      double myu=-mysin*q[j*n+i]+mycos*u[j*n+i];
      double delt=amp*mat[j*n+i]-myq;
      
      all_chisq[j]+=delt*delt+myu*myu;
      all_chisq_ref[j]+=myq*myq+myu*myu;
    }
    //all_chisq[j]*=weights[j];
    all_chisq[j]=weights[j]*(all_chisq[j]-all_chisq_ref[j]);

  }
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  double chisq=0;
  for (int i=0;i<nfreq;i++)
    chisq+=(all_chisq[i]);
  free(imin);
  free(imax);
  free(all_chisq);
  free(all_chisq_ref);
  //free(mat);
  return chisq;
}

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/

double calculate_chisq_qu_ramp_conv_real_cached(void *q_in, void *u_in,double *weights, double dt, double *params, double *freq, int nfreq, int n, void *cached_in)
{

  double fwhm=params[0];
  double scat=params[1];
  double alpha=params[2];
  double amp=params[3];
  double t0=params[4];
  double DM=params[5];
  double RM=params[6];
  double phi=params[7];
  double phigrad=params[8];
  

  if (scat<0)
    return 1e20;


  int *imin=(int *)malloc(sizeof(int)*(nfreq));
  int *imax=(int *)malloc(sizeof(int)*(nfreq));

  float complex *mat=(float complex *)cached_in;
  //printf("filling model.\n");
  double t1=omp_get_wtime();
  fill_model_real_qu_swing(dt,fwhm,scat,alpha,t0,DM,phigrad,freq,nfreq,n,(void *)mat,imin,imax);
  //printf("filled.\n");

  float *q=(float *)q_in;
  float *u=(float *)u_in;
  double *all_chisq=(double *)malloc(sizeof(double)*nfreq);
  double *all_chisq_ref=(double *)malloc(sizeof(double)*nfreq);
#pragma omp parallel for
  for (int j=0;j<nfreq;j++) {
    all_chisq[j]=0;
    all_chisq_ref[j]=0;
    double lambda=299.79/freq[j];
    double lambda_ref=299.79/REF_FREQ;
    double fac=lambda_ref*lambda_ref;
    double t_arr=t0+DM*DM0*(1/freq[j]/freq[j]-1/REF_FREQ/REF_FREQ);
    t_arr=t_arr/dt;  //do this in samples.  
    //if (j==0)
    // printf("t0 and t_arr are %12.5g %12.5g\n",t0,t_arr);

    double mycos=cos(RM*(lambda*lambda-fac)+phi);
    double mysin=sin(RM*(lambda*lambda-fac)+phi);

    for (int i=imin[j];i<=imax[j];i++)  {

      //double mycos=cos(RM*(lambda*lambda-fac)+phigrad*(i-t_arr)+phi);
      //double mysin=sin(RM*(lambda*lambda-fac)+phigrad*(i-t_arr)+phi);

      double myq=mycos*q[j*n+i]+mysin*u[j*n+i];
      double myu=-mysin*q[j*n+i]+mycos*u[j*n+i];
      double qdelt=amp*creal(mat[j*n+i])-myq;
      double udelt=amp*cimag(mat[j*n+i])-myu;
      
      all_chisq[j]+=qdelt*qdelt+udelt*udelt;
      all_chisq_ref[j]+=myq*myq+myu*myu;
    }
    //all_chisq[j]*=weights[j];
    all_chisq[j]=weights[j]*(all_chisq[j]-all_chisq_ref[j]);

  }
  //printf("got chi^2.\n");
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  double chisq=0;
  for (int i=0;i<nfreq;i++)
    chisq+=(all_chisq[i]);
  free(imin);
  free(imax);
  free(all_chisq);
  free(all_chisq_ref);
  //free(mat);
  return chisq;
}

/*--------------------------------------------------------------------------------*/

double calculate_chisq_qu_rmpow_real_cached(void *q_in, void *u_in,double *weights, double dt, double *params, double *freq, int nfreq, int n, void *cached_in)
{

  double fwhm=params[0];
  double scat=params[1];
  double alpha=params[2];
  double amp=params[3];
  double t0=params[4];
  double DM=params[5];
  double RM=params[6];
  double phi=params[7];
  double rmpow=params[8];

  if (scat<0)
    return 1e20;


  int *imin=(int *)malloc(sizeof(int)*(nfreq));
  int *imax=(int *)malloc(sizeof(int)*(nfreq));

  float *mat=(float *)cached_in;

  double t1=omp_get_wtime();
  fill_model_real(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,(void *)mat,imin,imax);

  float *q=(float *)q_in;
  float *u=(float *)u_in;
  double *all_chisq=(double *)malloc(sizeof(double)*nfreq);
  double *all_chisq_ref=(double *)malloc(sizeof(double)*nfreq);
  double lambda_ref=299.79/REF_FREQ;
  double lambda_0=299.79/700;
  double RM_use=RM*( (lambda_ref*lambda_ref-lambda_0*lambda_0)/(pow(lambda_ref,rmpow)-pow(lambda_0,rmpow)));
#pragma omp parallel for
  for (int j=0;j<nfreq;j++) {
    all_chisq[j]=0;
    all_chisq_ref[j]=0;
    double lambda=299.79/freq[j];
    double lambda_ref=299.79/REF_FREQ;
    double fac=pow(lambda,rmpow)-pow(lambda_ref,rmpow);
    double mycos=cos(RM_use*fac+phi);
    double mysin=sin(RM_use*fac+phi);
    for (int i=imin[j];i<=imax[j];i++)  {
      double myq=mycos*q[j*n+i]+mysin*u[j*n+i];
      double myu=-mysin*q[j*n+i]+mycos*u[j*n+i];
      double delt=amp*mat[j*n+i]-myq;
      
      all_chisq[j]+=delt*delt+myu*myu;
      all_chisq_ref[j]+=myq*myq+myu*myu;
    }
    //all_chisq[j]*=weights[j];
    all_chisq[j]=weights[j]*(all_chisq[j]-all_chisq_ref[j]);

  }
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  double chisq=0;
  for (int i=0;i<nfreq;i++)
    chisq+=(all_chisq[i]);
  free(imin);
  free(imax);
  free(all_chisq);
  free(all_chisq_ref);
  //free(mat);
  return chisq;
}

/*--------------------------------------------------------------------------------*/

double calculate_chisq_general_real_cached(void *dat_in, double *weights, double dt, double *params, double *freq, int nfreq, int n, void *cached_in)
{

  double fwhm=params[0];
  double scat=params[1];
  double scatpow=params[2];
  double alpha=params[3];
  double amp=params[4];
  double t0=params[5];
  double DM=params[6];
  double DM_pow=params[7];
  

  int *imin=(int *)malloc(sizeof(int)*(nfreq));
  int *imax=(int *)malloc(sizeof(int)*(nfreq));

  float *mat=(float *)cached_in;

  double t1=omp_get_wtime();
  fill_model_real_general(dt,fwhm,scat,scatpow,alpha,t0,DM,DM_pow,freq,nfreq,n,(void *)mat,imin,imax);

  float *dat=(float *)dat_in;
  double *all_chisq=(double *)malloc(sizeof(double)*nfreq);
  double *all_chisq_ref=(double *)malloc(sizeof(double)*nfreq);
#pragma omp parallel for
  for (int j=0;j<nfreq;j++) {
    all_chisq[j]=0;
    all_chisq_ref[j]=0;
    for (int i=imin[j];i<=imax[j];i++)  {
      double delt=amp*mat[j*n+i]-dat[j*n+i];
      all_chisq[j]+=delt*delt;
      all_chisq_ref[j]+=dat[j*n+i]*dat[j*n+i];
    }
    //all_chisq[j]*=weights[j];
    all_chisq[j]=weights[j]*(all_chisq[j]-all_chisq_ref[j]);

  }
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  double chisq=0;
  for (int i=0;i<nfreq;i++)
    chisq+=(all_chisq[i]);
  free(imin);
  free(imax);
  free(all_chisq);
  free(all_chisq_ref);
  //free(mat);
  return chisq;
}

/*--------------------------------------------------------------------------------*/


double calculate_chisq_general_nup_real_cached(void *dat_in, double *weights, double dt, double *params, double *freq, int nfreq, int n, void *cached_in)
{

  double fwhm=params[0];
  double scat=params[1];
  double scatpow=params[2];
  double alpha=params[3];
  double amp=params[4];
  double t0=params[5];
  double DM=params[6];
  double nu_p_sqr=params[7];
  

  int *imin=(int *)malloc(sizeof(int)*(nfreq));
  int *imax=(int *)malloc(sizeof(int)*(nfreq));

  float *mat=(float *)cached_in;

  double t1=omp_get_wtime();
  fill_model_real_general_nup(dt,fwhm,scat,scatpow,alpha,t0,DM,nu_p_sqr,freq,nfreq,n,(void *)mat,imin,imax);

  float *dat=(float *)dat_in;
  double *all_chisq=(double *)malloc(sizeof(double)*nfreq);
  double *all_chisq_ref=(double *)malloc(sizeof(double)*nfreq);
#pragma omp parallel for
  for (int j=0;j<nfreq;j++) {
    all_chisq[j]=0;
    all_chisq_ref[j]=0;
    for (int i=imin[j];i<=imax[j];i++)  {
      double delt=amp*mat[j*n+i]-dat[j*n+i];
      all_chisq[j]+=delt*delt;
      all_chisq_ref[j]+=dat[j*n+i]*dat[j*n+i];
    }
    //all_chisq[j]*=weights[j];
    all_chisq[j]=weights[j]*(all_chisq[j]-all_chisq_ref[j]);

  }
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  double chisq=0;
  for (int i=0;i<nfreq;i++)
    chisq+=(all_chisq[i]);
  free(imin);
  free(imax);
  free(all_chisq);
  free(all_chisq_ref);
  //free(mat);
  return chisq;
}

/*--------------------------------------------------------------------------------*/
double calculate_chisq_scatfit_cached(void *dat_in, float *weights, int imin,float dt, float *params, float *freq, int nfreq, int n, void *cached_in)
{

  float fwhm=params[0];
  float scat=1.0/params[1];
  float alpha=params[2];
  float amp=params[3];
  float t0=params[4];
  float DM=params[5];
  float scatpow=params[6];
  
  int nn=1+n/2;
  

  scat=scat*dt;
  t0=1+t0/dt;  //burst_model is expecting t0 to be in samples
  //float complex *mat=(float complex *)mat_in;
  float complex *mat=(float complex *)cached_in;
  //fill_model(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,mat_in);
  double t1=omp_get_wtime();
  fill_model_fast_scatfit(dt,fwhm,scat,scatpow,alpha,t0,DM,freq,nfreq,n,(void *)mat);
  //printf("params are %12.5g %12.5g %12.5g %12.5g %12.5g %12.5g\n",dt,fwhm,scat,alpha,t0,DM);
  //printf("amps are %12.6f %12.6f %12.6f\n",creal(mat[0]),creal(mat[1]),cimag(mat[1]));
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  float complex *dat=(float complex *)dat_in;
  double *all_chisq=(double *)malloc(sizeof(double)*nfreq);
#pragma omp parallel for
  for (int j=0;j<nfreq;j++) {
    all_chisq[j]=0;
    for (int i=imin;i<nn;i++)  {
      double complex delt=amp*mat[j*nn+i]-dat[j*nn+i];
      all_chisq[j]+=(creal(delt)*creal(delt)+cimag(delt)*cimag(delt));
    }
    all_chisq[j]*=weights[j];
  }
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  double chisq=0;
  for (int i=0;i<nfreq;i++)
    chisq+=all_chisq[i];
  
  free(all_chisq);
  //free(mat);
  return chisq;
}

/*--------------------------------------------------------------------------------*/

double calculate_chisq_linfit_cached(void *dat_in, float *weights, int imin,float dt, float *params, float *freq, int nfreq, int n, void *cached_in)
{

  float fwhm=params[0];
  float scat=1.0/params[1];
  //float alpha=params[2];
  //float amp=params[3];
  float slope=params[2];
  float offset=params[3];
  float t0=params[4];
  float DM=params[5];

  int nn=1+n/2;
  

  scat=scat*dt;
  t0=1+t0/dt;  //burst_model is expecting t0 to be in samples
  //float complex *mat=(float complex *)mat_in;
  float complex *mat=(float complex *)cached_in;
  //fill_model(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,mat_in);
  double t1=omp_get_wtime();
  float alpha=0.0;
  fill_model_fast(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,(void *)mat);
  //printf("params are %12.5g %12.5g %12.5g %12.5g %12.5g %12.5g\n",dt,fwhm,scat,alpha,t0,DM);
  //printf("amps are %12.6f %12.6f %12.6f\n",creal(mat[0]),creal(mat[1]),cimag(mat[1]));
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  float complex *dat=(float complex *)dat_in;
  double *all_chisq=(double *)malloc(sizeof(double)*nfreq);
#pragma omp parallel for
  for (int j=0;j<nfreq;j++) {
    all_chisq[j]=0;
    for (int i=imin;i<nn;i++)  {
      float nu_ref=800.0;
      float nu_width=100.0;
      double amp=offset+slope*(freq[j]-nu_ref)/nu_width;
      double complex delt=amp*mat[j*nn+i]-dat[j*nn+i];
      all_chisq[j]+=(creal(delt)*creal(delt)+cimag(delt)*cimag(delt));
    }
    all_chisq[j]*=weights[j];
  }
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  double chisq=0;
  for (int i=0;i<nfreq;i++)
    chisq+=all_chisq[i];
  
  free(all_chisq);
  //free(mat);
  return chisq;
}

/*--------------------------------------------------------------------------------*/


double calculate_chisq_cached_dmpow(void *dat_in, float *weights, int imin,float dt, float *params, float *freq, int nfreq, int n, void *cached_in)
{

  float fwhm=params[0];
  float scat=1.0/params[1];
  float alpha=params[2];
  float amp=params[3];
  float t0=params[4];
  float DM=params[5];
  float DM_pow=params[6];
  float tau=0;
  int nn=1+n/2;
  
  double nu2=700;  //need to adjust t0 and DM to decorrelate the DM index.  Requires a second frequency
  double a=pow(REF_FREQ,DM_pow+2);
  double b=REF_FREQ*REF_FREQ/DM0;
  double c=pow(nu2,DM_pow+2);
  double d=nu2*nu2/DM0;
  double e=DM;
  double f=DM;
  double DM_new=(f*b-e*d)/(b*c-a*d);
  double delta_t=(e*c-a*f)/(b*c-a*d);

  //printf("DM, eps, are %12.4f %12.6g, DM_new and dt are %12.4f %12.6g\n",DM,DM_pow+2,DM_new,delta_t);


  scat=scat*dt;
  tau=tau*dt;
  t0=1+(t0-0*delta_t)/dt;  //burst_model is expecting t0 to be in samples
  //float complex *mat=(float complex *)mat_in;
  float complex *mat=(float complex *)cached_in;
  //fill_model(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,mat_in);
  double t1=omp_get_wtime();
  //fill_model_fast(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,(void *)mat);
  printf("filling model.\n");
  //fill_model_fast_general(dt,fwhm,tau,scat,alpha,t0,DM,DM_pow,freq,nfreq,n,(void *)mat);
  fill_model_fast_general(dt,fwhm,tau,scat,alpha,t0,DM_new,DM_pow,freq,nfreq,n,(void *)mat);
  //printf("params are %12.5g %12.5g %12.5g %12.5g %12.5g %12.5g\n",dt,fwhm,scat,alpha,t0,DM);
  //printf("amps are %12.6f %12.6f %12.6f\n",creal(mat[0]),creal(mat[1]),cimag(mat[1]));
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  float complex *dat=(float complex *)dat_in;
  double *all_chisq=(double *)malloc(sizeof(double)*nfreq);
#pragma omp parallel for
  for (int j=0;j<nfreq;j++) {
    all_chisq[j]=0;
    for (int i=imin;i<nn;i++)  {
      double complex delt=amp*mat[j*nn+i]-dat[j*nn+i];
      all_chisq[j]+=(creal(delt)*creal(delt)+cimag(delt)*cimag(delt));
    }
    all_chisq[j]*=weights[j];
  }
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  double chisq=0;
  for (int i=0;i<nfreq;i++)
    chisq+=all_chisq[i];
  
  free(all_chisq);
  //free(mat);
  return chisq;
}

/*--------------------------------------------------------------------------------*/


double calculate_chisq_cached_tau(void *dat_in, float *weights, int imin,float dt, float *params, float *freq, int nfreq, int n, void *cached_in)
{

  float fwhm=params[0];
  float scat=1.0/params[1];
  float alpha=params[2];
  float amp=params[3];
  float t0=params[4];
  float DM=params[5];
  float DM_pow=-2.0;
  float tau=1.0/params[6];
  
  if (tau<0) 
    return 1e20;
  if (scat<0)
    return 1e20;
  if (fwhm<0)
    return 1e20;
    
  int nn=1+n/2;
  

  scat=scat*dt;
  tau=tau*dt;
  t0=1+t0/dt;  //burst_model is expecting t0 to be in samples
  //float complex *mat=(float complex *)mat_in;
  float complex *mat=(float complex *)cached_in;
  //fill_model(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,mat_in);
  double t1=omp_get_wtime();
  //fill_model_fast(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,(void *)mat);
  fill_model_fast_general(dt,fwhm,tau,scat,alpha,t0,DM,DM_pow,freq,nfreq,n,(void *)mat);
  //printf("params are %12.5g %12.5g %12.5g %12.5g %12.5g %12.5g\n",dt,fwhm,scat,alpha,t0,DM);
  //printf("amps are %12.6f %12.6f %12.6f\n",creal(mat[0]),creal(mat[1]),cimag(mat[1]));
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  float complex *dat=(float complex *)dat_in;
  double *all_chisq=(double *)malloc(sizeof(double)*nfreq);
#pragma omp parallel for
  for (int j=0;j<nfreq;j++) {
    all_chisq[j]=0;
    for (int i=imin;i<nn;i++)  {
      double complex delt=amp*mat[j*nn+i]-dat[j*nn+i];
      all_chisq[j]+=(creal(delt)*creal(delt)+cimag(delt)*cimag(delt));
    }
    all_chisq[j]*=weights[j];
  }
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  double chisq=0;
  for (int i=0;i<nfreq;i++)
    chisq+=all_chisq[i];
  
  free(all_chisq);
  //free(mat);
  return chisq;
}

/*--------------------------------------------------------------------------------*/
//double calculate_chisq(void *dat_in, float *weights, int imin,float dt, float fwhm, float scat, float alpha, float amp, float t0, float DM, float *freq, int nfreq, int n, void *mat_in)
double calculate_chisq_v2(void *dat_in, float *weights, int imin,float dt, float *params, float *freq, int nfreq, int n)
{

  float fwhm=params[0];
  float scat=params[1];
  float alpha=params[2];
  float amp=params[3];
  float t0=params[4];
  float DM=params[5];
  float DM_pow=params[6];

  int nn=1+n/2;
  

  scat=scat*dt;
  t0=1+t0/dt;  //burst_model is expecting t0 to be in samples
  //float complex *mat=(float complex *)mat_in;
  float complex *mat=(float complex *)malloc(sizeof(float complex)*nn*nfreq);
  //fill_model(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,mat_in);
  double t1=omp_get_wtime();
  fill_model(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,(void *)mat);
  //printf("params are %12.5g %12.5g %12.5g %12.5g %12.5g %12.5g\n",dt,fwhm,scat,alpha,t0,DM);
  //printf("amps are %12.6f %12.6f %12.6f\n",creal(mat[0]),creal(mat[1]),cimag(mat[1]));
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  float complex *dat=(float complex *)dat_in;
  double *all_chisq=(double *)malloc(sizeof(double)*nfreq);
#pragma omp parallel for
  for (int j=0;j<nfreq;j++) {
    all_chisq[j]=0;
    for (int i=imin;i<nn;i++)  {
      double complex delt=amp*mat[j*nn+i]-dat[j*nn+i];
      all_chisq[j]+=(creal(delt)*creal(delt)+cimag(delt)*cimag(delt));
    }
    all_chisq[j]*=weights[j];
  }
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  double chisq=0;
  for (int i=0;i<nfreq;i++)
    chisq+=all_chisq[i];
  
  free(all_chisq);
  free(mat);
  return chisq;
}

/*--------------------------------------------------------------------------------*/
//double calculate_chisq(void *dat_in, float *weights, int imin,float dt, float fwhm, float scat, float alpha, float amp, float t0, float DM, float *freq, int nfreq, int n, void *mat_in)
double calculate_chisq_ampfit(void *dat_in, float *weights, int imin,float dt, float fwhm, float scat, float alpha, float *amp_fit, float t0, float DM, float *freq, int nfreq, int n)
{

  int nn=1+n/2;
  

  scat=scat*dt;
  t0=1+t0/dt;  //burst_model is expecting t0 to be in samples
  //float complex *mat=(float complex *)mat_in;
  float complex *mat=(float complex *)malloc(sizeof(float complex)*nn*nfreq);
  //fill_model(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,mat_in);
  double t1=omp_get_wtime();
  fill_model(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,(void *)mat);
  //printf("params are %12.5g %12.5g %12.5g %12.5g %12.5g %12.5g\n",dt,fwhm,scat,alpha,t0,DM);
  //printf("amps are %12.6f %12.6f %12.6f\n",creal(mat[0]),creal(mat[1]),cimag(mat[1]));

  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  float complex *dat=(float complex *)dat_in;
  double *vec1=(double *)malloc(sizeof(double)*nfreq);
  double *vec2=(double *)malloc(sizeof(double)*nfreq);
  double *all_chisq=(double *)malloc(sizeof(double)*nfreq);

#pragma omp parallel for
  for (int j=0;j<nfreq;j++) {
    vec1[j]=0;
    vec2[j]=0;
    for (int i=imin;i<nn;i++) {
      float r1=creal(mat[j*nn+i]);
      float i1=cimag(mat[j*nn+i]);
      float r2=creal(dat[j*nn+i]);
      float i2=cimag(dat[j*nn+i]);
      vec1[j]+=r1*r1+i1*i1;
      vec2[j]+=r1*r2+i1*i2;

    }
    
  }
  double denom=0;
  double num=0;
  for (int j=0;j<nfreq;j++) {
    denom+=vec1[j]*weights[j];
    num+=vec2[j]*weights[j];
  }
  float amp=num/denom;
  *amp_fit=amp;

#pragma omp parallel for
  for (int j=0;j<nfreq;j++) {
    all_chisq[j]=0;
    for (int i=imin;i<nn;i++)  {
      double complex delt=amp*mat[j*nn+i]-dat[j*nn+i];
      all_chisq[j]+=(creal(delt)*creal(delt)+cimag(delt)*cimag(delt));
    }
    all_chisq[j]*=weights[j];
  }
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  double chisq=0;
  for (int i=0;i<nfreq;i++)
    chisq+=all_chisq[i];
  
  
  free(vec1);
  free(vec2);
  free(all_chisq);
  free(mat);
  return chisq;
}

/*--------------------------------------------------------------------------------*/
double calculate_chisq_qu(void *dat_q_in, void *dat_u_in, float *weights, int imin,float dt,float *params, float *freq, int nfreq, int n)

{

  float fwhm=params[0];
  float scat=1.0/params[1];
  float alpha=params[2];
  float amp=params[3];
  float t0=params[4];
  float DM=params[5];
  float RM=params[6];
  float phi=params[7]; 
  int nn=1+n/2;

  //printf("fwhm, scat, amp, RM, and phi are %12.4g %12.4g %12.4g %12.4g %12.4g\n",fwhm,scat,amp,RM,phi);
  

  scat=scat*dt;
  t0=1+t0/dt;  //burst_model is expecting t0 to be in samples
  //float complex *mat=(float complex *)mat_in;
  float complex *mat=(float complex *)malloc(sizeof(float complex)*nn*nfreq);
  //fill_model(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,mat_in);
  double t1=omp_get_wtime();
  fill_model(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,(void *)mat);
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  float complex *dat_q=(float complex *)dat_q_in;
  float complex *dat_u=(float complex *)dat_u_in;
  double *all_chisq=(double *)malloc(sizeof(double)*nfreq);

  
  //find that if RM goes up by 1, phi should go down by .143
  //so, to keep phi constant, should be RM*(lambda^2-.143)
  //implies lambda_ref is .3777,and nu_ref is 793 



#pragma omp parallel for
  for (int j=0;j<nfreq;j++) {
    float lambda=299.79/freq[j];
    float nu_ref=793;
    float lambda_ref =299.79/nu_ref;
    float fac=lambda_ref*lambda_ref+0.010973; //.010973 comes from looking at some likelihoods
    if (j==0)
      printf("fac is %12.6f\n",fac);
    float amp_q=amp*cos(RM*(lambda*lambda-fac)+phi);
    float amp_u=amp*sin(RM*(lambda*lambda-fac)+phi);
    float mycos=cos(RM*(lambda*lambda-fac)+phi);
    float mysin=sin(RM*(lambda*lambda-fac)+phi);
   
    if (j==0) 
      printf("lambda is %12.4f and amps are %12.4f %12.4f %12.6f %12.6f %12.6f\n",lambda,amp_q,amp_u,creal(mat[0]),creal(mat[1]),cimag(mat[1]));

    all_chisq[j]=0;
    for (int i=imin;i<nn;i++)  {
#if 0
      double complex deltq=amp_q*mat[j*nn+i]-dat_q[j*nn+i];
      double complex deltu=amp_u*mat[j*nn+i]-dat_u[j*nn+i];
#else
      float complex myq=mycos*dat_q[j*nn+i]+mysin*dat_u[j*nn+i];
      float  complex myu=0-mysin*dat_q[j*nn+i]+mycos*dat_u[j*nn+i];
      double complex deltq=myq-amp*mat[j*nn+i];
      double complex deltu=myu;
#endif
      all_chisq[j]+=(creal(deltq)*creal(deltq)+cimag(deltq)*cimag(deltq)+creal(deltu)*creal(deltu)+cimag(deltu)*cimag(deltu));
    }
    all_chisq[j]*=weights[j];
  }
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  double chisq=0;
  for (int i=0;i<nfreq;i++)
    chisq+=all_chisq[i];
  
  free(all_chisq);
  free(mat);
  return chisq;
}



/*--------------------------------------------------------------------------------*/

double calculate_chisq_qu_cached(void *dat_q_in, void *dat_u_in, float *weights, int imin,float dt,float *params, float *freq, int nfreq, int n, void *myscratch)

{

  float fwhm=params[0];
  float scat=1.0/params[1];
  float alpha=params[2];
  float amp=params[3];
  float t0=params[4];
  float DM=params[5];
  float RM=params[6];
  float phi=params[7]; 
  int nn=1+n/2;

  //printf("fwhm, scat, amp, RM, and phi are %12.4g %12.4g %12.4g %12.4g %12.4g\n",fwhm,scat,amp,RM,phi);
  

  scat=scat*dt;
  t0=1+t0/dt;  //burst_model is expecting t0 to be in samples
  //float complex *mat=(float complex *)mat_in;
  float complex *mat=(float complex *)myscratch;
  //fill_model(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,mat_in);
  double t1=omp_get_wtime();
  fill_model_fast_general(dt,fwhm,0,scat,alpha,t0,DM,-2.0,freq,nfreq,n,(void *)mat);
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  float complex *dat_q=(float complex *)dat_q_in;
  float complex *dat_u=(float complex *)dat_u_in;
  double *all_chisq=(double *)malloc(sizeof(double)*nfreq);

  
  //find that if RM goes up by 1, phi should go down by .143
  //so, to keep phi constant, should be RM*(lambda^2-.143)
  //implies lambda_ref is .3777,and nu_ref is 793 



#pragma omp parallel for
  for (int j=0;j<nfreq;j++) {
    float lambda=299.79/freq[j];
    float nu_ref=793;
    float lambda_ref =299.79/nu_ref;
    float fac=lambda_ref*lambda_ref+0.010973; //.010973 comes from looking at some likelihoods
    //if (j==0)
    //  printf("fac is %12.6f\n",fac);
    float amp_q=amp*cos(RM*(lambda*lambda-fac)+phi);
    float amp_u=amp*sin(RM*(lambda*lambda-fac)+phi);
    float mycos=cos(RM*(lambda*lambda-fac)+phi);
    float mysin=sin(RM*(lambda*lambda-fac)+phi);
   
    //if (j==0)  
    //  printf("lambda is %12.4f and amps are %12.4f %12.4f %12.6f %12.6f %12.6f\n",lambda,amp_q,amp_u,creal(mat[0]),creal(mat[1]),cimag(mat[1]));

    all_chisq[j]=0;
    for (int i=imin;i<nn;i++)  {
#if 0
      double complex deltq=amp_q*mat[j*nn+i]-dat_q[j*nn+i];
      double complex deltu=amp_u*mat[j*nn+i]-dat_u[j*nn+i];
#else
      float complex myq=mycos*dat_q[j*nn+i]+mysin*dat_u[j*nn+i];
      float  complex myu=0-mysin*dat_q[j*nn+i]+mycos*dat_u[j*nn+i];
      double complex deltq=myq-amp*mat[j*nn+i];
      double complex deltu=myu;
#endif
      all_chisq[j]+=(creal(deltq)*creal(deltq)+cimag(deltq)*cimag(deltq)+creal(deltu)*creal(deltu)+cimag(deltu)*cimag(deltu));
    }
    //if (j==0)
    //  printf("all_chisq[%d] is %12.5f, data are %12.4f %12.4f\n",j,all_chisq[j],creal(dat_u[0]),cimag(dat_u[0]));
      
    all_chisq[j]*=weights[j];
  }
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  double chisq=0;
  for (int i=0;i<nfreq;i++)
    chisq+=all_chisq[i];
  
  free(all_chisq);
  //free(mat);
  return chisq;
}



/*--------------------------------------------------------------------------------*/

double calculate_chisq_qu_rmpow_cached(void *dat_q_in, void *dat_u_in, float *weights, int imin,float dt,float *params, float *freq, int nfreq, int n, void *myscratch)

{

  float fwhm=params[0];
  float scat=1.0/params[1];
  float alpha=params[2];
  float amp=params[3];
  float t0=params[4];
  float DM=params[5];
  float RM=params[6];
  float phi=params[7]; 
  float rmpow=params[8];
  int nn=1+n/2;

  //printf("fwhm, scat, amp, RM, and phi are %12.4g %12.4g %12.4g %12.4g %12.4g\n",fwhm,scat,amp,RM,phi);
  

  scat=scat*dt;
  t0=1+t0/dt;  //burst_model is expecting t0 to be in samples
  //float complex *mat=(float complex *)mat_in;
  float complex *mat=(float complex *)myscratch;
  //fill_model(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,mat_in);
  double t1=omp_get_wtime();
  fill_model_fast_general(dt,fwhm,0,scat,alpha,t0,DM,-2.0,freq,nfreq,n,(void *)mat);
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  float complex *dat_q=(float complex *)dat_q_in;
  float complex *dat_u=(float complex *)dat_u_in;
  double *all_chisq=(double *)malloc(sizeof(double)*nfreq);

  
  //find that if RM goes up by 1, phi should go down by .143
  //so, to keep phi constant, should be RM*(lambda^2-.143)
  //implies lambda_ref is .3777,and nu_ref is 793 

  float lambda_ref=299.79/REF_FREQ;
  float lambda_0=299.79/700;
  float RM_use=RM*( (lambda_ref*lambda_ref-lambda_0*lambda_0)/(pow(lambda_ref,rmpow)-pow(lambda_0,rmpow)));

#pragma omp parallel for
  for (int j=0;j<nfreq;j++) {
    float lambda=299.79/freq[j];
    float lambda_ref =299.79/REF_FREQ;
    float fac=pow(lambda,rmpow)-pow(lambda_ref,rmpow);

    float mycos=cos(RM_use*fac+phi);
    float mysin=sin(RM_use*fac+phi);
   
    //if (j==0)  
    //  printf("lambda is %12.4f and amps are %12.4f %12.4f %12.6f %12.6f %12.6f\n",lambda,amp_q,amp_u,creal(mat[0]),creal(mat[1]),cimag(mat[1]));

    all_chisq[j]=0;
    for (int i=imin;i<nn;i++)  {
#if 0
      double complex deltq=amp_q*mat[j*nn+i]-dat_q[j*nn+i];
      double complex deltu=amp_u*mat[j*nn+i]-dat_u[j*nn+i];
#else
      float complex myq=mycos*dat_q[j*nn+i]+mysin*dat_u[j*nn+i];
      float  complex myu=0-mysin*dat_q[j*nn+i]+mycos*dat_u[j*nn+i];
      double complex deltq=myq-amp*mat[j*nn+i];
      double complex deltu=myu;
#endif
      all_chisq[j]+=(creal(deltq)*creal(deltq)+cimag(deltq)*cimag(deltq)+creal(deltu)*creal(deltu)+cimag(deltu)*cimag(deltu));
    }
    //if (j==0)
    //  printf("all_chisq[%d] is %12.5f, data are %12.4f %12.4f\n",j,all_chisq[j],creal(dat_u[0]),cimag(dat_u[0]));
      
    all_chisq[j]*=weights[j];
  }
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  double chisq=0;
  for (int i=0;i<nfreq;i++)
    chisq+=all_chisq[i];
  
  free(all_chisq);
  //free(mat);
  return chisq;
}



/*--------------------------------------------------------------------------------*/
double calculate_chisq_qu_rmpow(void *dat_q_in, void *dat_u_in, float *weights, int imin,float dt,float *params, float *freq, int nfreq, int n)

{

  float fwhm=params[0];
  float scat=1.0/params[1];
  float alpha=params[2];
  float amp=params[3];
  float t0=params[4];
  float DM=params[5];
  float RM=params[6];
  float phi=params[7]; 
  float rmpow=params[8];
  int nn=1+n/2;

  //printf("fwhm, scat, amp, RM, and phi are %12.4g %12.4g %12.4g %12.4g %12.4g\n",fwhm,scat,amp,RM,phi);
  

  scat=scat*dt;
  t0=1+t0/dt;  //burst_model is expecting t0 to be in samples
  //float complex *mat=(float complex *)mat_in;
  float complex *mat=(float complex *)malloc(sizeof(float complex)*nn*nfreq);
  //fill_model(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,mat_in);
  double t1=omp_get_wtime();
  fill_model(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,(void *)mat);
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  float complex *dat_q=(float complex *)dat_q_in;
  float complex *dat_u=(float complex *)dat_u_in;
  double *all_chisq=(double *)malloc(sizeof(double)*nfreq);

  
  //find that if RM goes up by 1, phi should go down by .143
  //so, to keep phi constant, should be RM*(lambda^2-.143)
  //implies lambda_ref is .3777,and nu_ref is 793 



#pragma omp parallel for
  for (int j=0;j<nfreq;j++) {
    float lambda=299.79/freq[j];
    float lambda_ref =299.79/REF_FREQ;
    float fac=lambda_ref*lambda_ref;
    float myarg=RM*(pow(lambda,rmpow)-pow(lambda_ref,rmpow))+phi; //get the d

    float amp_q=amp*cos(myarg);
    float amp_u=amp*sin(myarg);
    float mycos=cos(myarg);
    float mysin=sin(myarg);
   
    //if (j==0) 
    // printf("lambda is %12.4f and amps are %12.4f %12.4f %12.6f %12.6f %12.6f\n",lambda,amp_q,amp_u,creal(mat[0]),creal(mat[1]),cimag(mat[1]));

    all_chisq[j]=0;
    for (int i=imin;i<nn;i++)  {
      float complex myq=mycos*dat_q[j*nn+i]+mysin*dat_u[j*nn+i];
      float  complex myu=0-mysin*dat_q[j*nn+i]+mycos*dat_u[j*nn+i];
      double complex deltq=myq-amp*mat[j*nn+i];
      double complex deltu=myu;
      
      all_chisq[j]+=(creal(deltq)*creal(deltq)+cimag(deltq)*cimag(deltq)+creal(deltu)*creal(deltu)+cimag(deltu)*cimag(deltu));
    }
    all_chisq[j]*=weights[j];
  }
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  double chisq=0;
  for (int i=0;i<nfreq;i++)
    chisq+=all_chisq[i];
  
  free(all_chisq);
  free(mat);
  return chisq;
}


/*--------------------------------------------------------------------------------*/
//double calculate_chisq(void *dat_in, float *weights, int imin,float dt, float fwhm, float scat, float alpha, float amp, float t0, float DM, float *freq, int nfreq, int n, void *mat_in)
void calculate_many_chisq_qu(void *dat_q_in, void *dat_u_in, float *weights, int imin,float dt,float *params, float *freq, int nfreq, int n, float *amps, float *RMs, float *phis, int nlike, double *chisq_vec)
// float fwhm, float scat, float alpha, float *amp_fit, float t0, float DM, float *freq, int nfreq, int n)
{

  float fwhm=params[0];
  float scat=1.0/params[1];
  float alpha=params[2];
  float amp=params[3];
  float t0=params[4];
  float DM=params[5];
  float RM=params[6];
  float phi=params[7]; 
  int nn=1+n/2;

  //printf("fwhm, scat, amp, RM, and phi are %12.4g %12.4g %12.4g %12.4g %12.4g\n",fwhm,scat,amp,RM,phi);
  
  
  scat=scat*dt;
  t0=1+t0/dt;  //burst_model is expecting t0 to be in samples
  //float complex *mat=(float complex *)mat_in;
  float complex *mat=(float complex *)malloc(sizeof(float complex)*nn*nfreq);
  //fill_model(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,mat_in);
  double t1=omp_get_wtime();
  printf("filling model.\n");
  fill_model(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,(void *)mat);
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  float complex *dat_q=(float complex *)dat_q_in;
  float complex *dat_u=(float complex *)dat_u_in;
  double *all_chisq=(double *)malloc(sizeof(double)*nfreq);
  double *ref_chisq=(double *)malloc(sizeof(double)*nfreq);

  #pragma omp parallel for
  for (int i=0;i<nfreq;i++) {
    ref_chisq[i]=0;
    for (int j=imin;j<nn;j++) {
      ref_chisq[i]+=creal(dat_q[i*nn+j])*creal(dat_q[i*nn+j])+cimag(dat_q[i*nn+j])*cimag(dat_q[i*nn+j])+creal(dat_u[i*nn+j])*creal(dat_u[i*nn+j])+cimag(dat_u[i*nn+j])*cimag(dat_u[i*nn+j]);
    }
    ref_chisq[i]*=weights[i];
    if (i==1500)
      printf("ref_chisq is %14.4g\n",ref_chisq[i]);
  
  }
  

  //float lambda_ref=299.79/750.0;
  float lambda_ref=299.79/793.0;
  float fac=lambda_ref*lambda_ref+0.010973;
  printf("lambda_ref is %12.4f %12.5f\n",lambda_ref,fac);
  for (int ii=0;ii<nlike;ii++) {
    
    amp=amps[ii];
    RM=RMs[ii];
    phi=phis[ii];
#pragma omp parallel for
    for (int j=0;j<nfreq;j++) {
      float lambda=299.79/freq[j];
      float amp_q=amp*cos(RM*(lambda*lambda-fac)+phi);
      float amp_u=amp*sin(RM*(lambda*lambda-fac)+phi);
      float mycos=cos(RM*(lambda*lambda-fac)+phi);
      float mysin=sin(RM*(lambda*lambda-fac)+phi);
      
      //if (j==0) 
      //printf("lambda is %12.4f and amps are %12.4f %12.4f %12.6f %12.6f %12.6f\n",lambda,amp_q,amp_u,creal(mat[0]),creal(mat[1]),cimag(mat[1]));
      
      all_chisq[j]=0;
      for (int i=imin;i<nn;i++)  {
#if 0
	double complex deltq=amp_q*mat[j*nn+i]-dat_q[j*nn+i];
	double complex deltu=amp_u*mat[j*nn+i]-dat_u[j*nn+i];
#else
	float complex myq=mycos*dat_q[j*nn+i]+mysin*dat_u[j*nn+i];
	float  complex myu=0-mysin*dat_q[j*nn+i]+mycos*dat_u[j*nn+i];
	double complex deltq=myq-amp*mat[j*nn+i];
	double complex deltu=myu;
#endif
	all_chisq[j]+=(creal(deltq)*creal(deltq)+cimag(deltq)*cimag(deltq)+creal(deltu)*creal(deltu)+cimag(deltu)*cimag(deltu));
      }
      all_chisq[j]*=weights[j];
    }
    //printf("time here is %12.5f\n",omp_get_wtime()-t1);
    double chisq=0;
    for (int i=0;i<nfreq;i++)
      chisq+=all_chisq[i]-ref_chisq[i];
    chisq_vec[ii]=chisq;
  }
  free(all_chisq);
  free(ref_chisq);
  free(mat);
  return;
}

/*--------------------------------------------------------------------------------*/

void calculate_many_chisq_qu_ampphase_fit(void *dat_q_in, void *dat_u_in, float *weights, int imin,float dt,float *params, float *freq, int nfreq, int n, float *amps, float *RMs, float *phis, int nlike, double *chisq_vec)
// float fwhm, float scat, float alpha, float *amp_fit, float t0, float DM, float *freq, int nfreq, int n)
{

  float fwhm=params[0];
  float scat=1.0/params[1];
  float alpha=params[2];
  float amp=params[3];
  float t0=params[4];
  float DM=params[5];
  float RM=params[6];
  float phi=params[7]; 
  int nn=1+n/2;

  //printf("fwhm, scat, amp, RM, and phi are %12.4g %12.4g %12.4g %12.4g %12.4g\n",fwhm,scat,amp,RM,phi);
  
  
  scat=scat*dt;
  t0=1+t0/dt;  //burst_model is expecting t0 to be in samples
  //float complex *mat=(float complex *)mat_in;
  float complex *mat=(float complex *)malloc(sizeof(float complex)*nn*nfreq);
  //fill_model(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,mat_in);
  double t1=omp_get_wtime();
  printf("filling model.\n");
  fill_model(dt,fwhm,scat,alpha,t0,DM,freq,nfreq,n,(void *)mat);
  //printf("time here is %12.5f\n",omp_get_wtime()-t1);
  float complex *dat_q=(float complex *)dat_q_in;
  float complex *dat_u=(float complex *)dat_u_in;
  double *all_chisq=(double *)malloc(sizeof(double)*nfreq);

  

  float lambda_ref=299.79/750.0;
  printf("lambda_ref is %12.4f\n",lambda_ref);
  for (int ii=0;ii<nlike;ii++) {
    
    amp=amps[ii];
    RM=RMs[ii];
    phi=phis[ii];
#pragma omp parallel for
    for (int j=0;j<nfreq;j++) {
      float lambda=299.79/freq[j];
      float amp_q=amp*cos(RM*(lambda*lambda-lambda_ref*lambda_ref)+phi);
      float amp_u=amp*sin(RM*(lambda*lambda-lambda_ref*lambda_ref)+phi);
      float mycos=cos(RM*(lambda*lambda-lambda_ref*lambda_ref)+phi);
      float mysin=sin(RM*(lambda*lambda-lambda_ref*lambda_ref)+phi);
      
      //if (j==0) 
      //printf("lambda is %12.4f and amps are %12.4f %12.4f %12.6f %12.6f %12.6f\n",lambda,amp_q,amp_u,creal(mat[0]),creal(mat[1]),cimag(mat[1]));
      
      all_chisq[j]=0;
      for (int i=imin;i<nn;i++)  {
#if 0
	double complex deltq=amp_q*mat[j*nn+i]-dat_q[j*nn+i];
	double complex deltu=amp_u*mat[j*nn+i]-dat_u[j*nn+i];
#else
	float complex myq=mycos*dat_q[j*nn+i]+mysin*dat_u[j*nn+i];
	float  complex myu=0-mysin*dat_q[j*nn+i]+mycos*dat_u[j*nn+i];
	double complex deltq=myq-amp*mat[j*nn+i];
	double complex deltu=myu;
#endif
	all_chisq[j]+=(creal(deltq)*creal(deltq)+cimag(deltq)*cimag(deltq)+creal(deltu)*creal(deltu)+cimag(deltu)*cimag(deltu));
      }
      all_chisq[j]*=weights[j];
    }
    //printf("time here is %12.5f\n",omp_get_wtime()-t1);
    double chisq=0;
    for (int i=0;i<nfreq;i++)
      chisq+=all_chisq[i];
    chisq_vec[ii]=chisq;
  }
  free(all_chisq);
  free(mat);
  return;
}
