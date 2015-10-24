#include <octave/oct.h>
#include <iostream>
#include <stdio.h>
#include <omp.h>
#include <complex>

//export DL_LD=g++-4.8
//export CXX=g++-4.8
//export CXXFLAGS="-O3 -fopenmp"
//mkoctfile-3.8.0 burst_routines.cpp -lgomp -v
//mkoctfile-3.8.0 burst_routines.cpp -I. -L. -lburst_fill  -lgomp -v
//  mkoctfile burst_routines.cpp -I. -L. -lburst_fill  -lgomp -v
//  mkoctfile burst_routines.cpp -I. -L. -L/home/sievers/local/lib -lburst_fill -lcerf -lgomp -v

#ifdef __cplusplus
extern "C"
{
#endif
#include <burst_fill.h>
  //#include <dirfile.h>
#ifdef __cplusplus
}  /* end extern "C" */
#endif

/*--------------------------------------------------------------------------------*/
void *get_pointer(octave_value val)
{
  int64NDArray myptr=val.array_value();
  long myptr2=myptr(0,0);
  return (void *)myptr2;

}



/*--------------------------------------------------------------------------------*/

double get_value(octave_value val)
{
  NDArray myptr=val.array_value();
  double myval=(double)myptr(0,0);
  return myval;
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (freeptr,args,nargout,"Free a pointer.\n")
{
  void *ptr=get_pointer(args(0));
  free(ptr);
  return octave_value_list();
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (floatcomplex2ptr, args, nargout, "Turn a float complex matrix into a pointer.\n")
{

  FloatComplexMatrix dat=args(0).float_complex_matrix_value();
  dim_vector dm=dat.dims();
  long nelem=dm(0)*dm(1);
  printf("have %ld elements in matrix.\n",nelem);
  void *datptr=malloc(sizeof(float complex)*nelem);
  memcpy(datptr,dat.fortran_vec(),nelem*sizeof(float complex));
  int64NDArray myptr(1);
  
  myptr(0)=(long)datptr;
  //long *ptr=myptr.fortran_vec();
  //ptr[0]=(long)mytod;
  return octave_value(myptr);
}


/*--------------------------------------------------------------------------------*/
DEFUN_DLD (float2ptr, args, nargout, "Turn a float matrix into a pointer.\n")
{

  FloatMatrix dat=args(0).float_matrix_value();
  dim_vector dm=dat.dims();
  long nelem=dm(0)*dm(1);
  printf("have %ld elements in matrix.\n",nelem);
  void *datptr=malloc(sizeof(float)*nelem);
  memcpy(datptr,dat.fortran_vec(),nelem*sizeof(float));
  int64NDArray myptr(1);
  
  myptr(0)=(long)datptr;
  //long *ptr=myptr.fortran_vec();
  //ptr[0]=(long)mytod;
  return octave_value(myptr);
}


/*--------------------------------------------------------------------------------*/


DEFUN_DLD (octave_pulsars_test, args, nargout, "Say hi!\n")
{
  if (args.length()<1) {
    fprintf(stderr,"Need arguments.\n");
    return octave_value_list();
  }
  Matrix mat=args(0).matrix_value();
  
  double *ptr=mat.fortran_vec();
  printf("hello!  First element is %12.4f\n",ptr[0]);
#pragma omp parallel
  printf("I am %d of %d\n",omp_get_thread_num(),omp_get_num_threads());
  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (make_burst_model_real_c, args, nargout, "Calculate a model for an FRB, brute-force\n")
{
  
  if (args.length()<8) {
    printf("need at least 8 arguments to make_burst_model_c.\n");
    return octave_value_list();
  }
  double dt=get_value(args(0));
  int n=(int)get_value(args(1));
  ColumnVector freq=args(2).column_vector_value();
  double fwhm=get_value(args(3));
  double scat=get_value(args(4));
  double alpha=get_value(args(5));
  double t0=get_value(args(6));
  double DM=get_value(args(7));
  int nfreq=freq.length();
  dim_vector dm;
  dm(0)=n;
  dm(1)=nfreq;
  FloatMatrix mat(dm);
  double t1=omp_get_wtime();
  memset(mat.fortran_vec(),0,n*nfreq*sizeof(float));
  fill_model_real(dt,fwhm,scat,alpha,t0,DM,freq.fortran_vec(),nfreq,n,mat.fortran_vec(),NULL,NULL);
  double t2=omp_get_wtime();
  printf("model filling took %12.4f seconds.\n",t2-t1);
  
  return octave_value(mat);

}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (make_burst_model_real_qu_swing_c, args, nargout, "Calculate a model for an FRB, brute-force\n")
{
  
  if (args.length()<9) {
    printf("need at least 8 arguments to make_burst_model_c.\n");
    return octave_value_list();
  }
  double dt=get_value(args(0));
  int n=(int)get_value(args(1));
  ColumnVector freq=args(2).column_vector_value();
  double fwhm=get_value(args(3));
  double scat=get_value(args(4));
  double alpha=get_value(args(5));
  double t0=get_value(args(6));
  double DM=get_value(args(7));
  double omega=get_value(args(8));
  int nfreq=freq.length();
  dim_vector dm;
  dm(0)=n;
  dm(1)=nfreq;
  FloatComplexMatrix mat(dm);
  memset(mat.fortran_vec(),0,2*n*nfreq*sizeof(float));
  double t1=omp_get_wtime();
  fill_model_real_qu_swing(dt,fwhm,scat,alpha,t0,DM,omega,freq.fortran_vec(),nfreq,n,mat.fortran_vec(),NULL,NULL);
  double t2=omp_get_wtime();
  printf("model filling took %12.4f seconds.\n",t2-t1);
  
  return octave_value(mat);

}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (make_burst_model_c, args, nargout, "Calculate a model for an FRB, brute-force\n")
{
  if (args.length()<8) {
    printf("need at least 8 arguments to make_burst_model_c.\n");
    return octave_value_list();
  }
  float pi=3.141592653;
  float freq_ref=800;
  float dt=get_value(args(0));
  int n=(int)get_value(args(1));
  FloatColumnVector freq=args(2).float_column_vector_value();
  //printf("have %d frequencies\n",freq.length());
  float fwhm=get_value(args(3));
  float scat=get_value(args(4));
  float alpha=get_value(args(5));
  float t0=get_value(args(6));
  float DM=get_value(args(7));
  float DM0=4150;
  scat=scat*dt;
  float sig=fwhm/sqrt(8*log(2))/dt;
  t0=1+t0/dt;

  int nfreq=freq.length();
  int nn=1+n/2;
#if 0
  FloatColumnVector dm_lags(freq.length());
  for (int i=0;i<nfreq;i++)
    dm_lags(i)=DM0*DM/(freq(i)*freq(i));
  printf("min lag is %12.5f\n",dm_lags.min());


  FloatColumnVector sigvec(nn);
  for (int i=0;i<nn;i++) {
    float fac=2*pi*sig*i/n;
    sigvec(i)=exp(-0.5*fac*fac);
  }
#endif


  dim_vector dm;
  dm(0)=nn;
  dm(1)=freq.length();

  FloatComplexMatrix mat(dm);
  //float _Complex *a=mat.fortran_vec();
  //std::complex<float> *ptr=mat.fortran_vec();
  void *ptr=(void *)mat.fortran_vec();


  //a[0]=1+2I;


  //fill_model(float dt, float fwhm, float scat, float alpha, float t0, float DM, float *freq, int nfreq, int n, void *mat_in)

  double t1=omp_get_wtime();
  fill_model_fast(dt,fwhm,scat,alpha,t0,DM,freq.fortran_vec(),nfreq,n,mat.fortran_vec());
  double t2=omp_get_wtime();
  printf("model filling took %12.4f seconds.\n",t2-t1);

  
#if 0
  for (int j=0;j<nfreq;j++) {
    float mylag=dm_lags(j)+t0;
    float fac=freq(j)/freq_ref;
    float myscat=scat/(fac*fac*fac*fac);
    std::complex<float> *myptr=ptr+nn*j;
    for (int i=0;i<nn;i++) {
      
      std::complex<float> aa=(-myscat,-2*pi*i/n);
      //std::complex<float> scatvec=(1-exp(n*aa))/(1-exp(aa));
      //std::complex<float> scatvec=std::exp(n*aa);
      //scatvec=1.0-scatvec;
	
    }
  }
#endif

  //double complex ptr=mat.fortran_vec()[0];
  //std::complex<double> c(1,2);
  //printf("c is %g %g\n",c.real(),c.imag());
  //cout << c << "\n";

  printf("nn is %d\n",nn);



  return octave_value(mat);
  
  
  //ComplexMatrix mat(dm);
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_burst_chisq_c, args, nargout, "Calculate chisq for an FRB, brute-force\n")
{
  if (args.length()<7) {
    printf("Need 7 args to get_burst_chisq_c.\n");
    return octave_value_list();
  }
  FloatColumnVector guess=args(0).float_column_vector_value();

  float fwhm=guess(0);
  float scat=1.0/guess(1);
  float alpha=guess(2);
  float amp=guess(3);
  float t0=guess(4);
  float DM=guess(5);

  
  FloatComplexMatrix dat=args(1).float_complex_matrix_value();;

  //dim_vector dm=dat.dims();
  //FloatComplexMatrix mat(dm);
  FloatColumnVector freqs=args(2).float_column_vector_value();
  int n=(int)get_value(args(3));
  FloatColumnVector nvec=args(4).float_column_vector_value();
  float dt=(float)get_value(args(5));
  int imin=(int)get_value(args(6));
  
  //double chisq=calculate_chisq(dat.fortran_vec(),nvec.fortran_vec(),imin,dt,fwhm,scat,alpha,amp,t0,DM,freqs.fortran_vec(),freqs.length(),n,mat.fortran_vec());
  double chisq=calculate_chisq(dat.fortran_vec(),nvec.fortran_vec(),imin,dt,fwhm,scat,alpha,amp,t0,DM,freqs.fortran_vec(),freqs.length(),n);



  return octave_value(chisq);
    
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_burst_chisq_cached_c, args, nargout, "Calculate chisq for an FRB, brute-force, data should already be cached.\n")
{
  if (args.length()<8) {
    printf("Need 8 args to get_burst_chisq_cached_c.\n");
    return octave_value_list();
  }
  FloatColumnVector guess=args(0).float_column_vector_value();


  void *datft=get_pointer(args(1));
  FloatColumnVector freqs=args(2).float_column_vector_value();
  int n=(int)get_value(args(3));
  FloatColumnVector nvec=args(4).float_column_vector_value();
  float dt=(float)get_value(args(5));
  int imin=(int)get_value(args(6));
  void *myscratch=get_pointer(args(7));

  //double chisq=calculate_chisq(dat.fortran_vec(),nvec.fortran_vec(),imin,dt,fwhm,scat,alpha,amp,t0,DM,freqs.fortran_vec(),freqs.length(),n,mat.fortran_vec());
  double chisq=calculate_chisq_cached(datft,nvec.fortran_vec(),imin,dt,guess.fortran_vec(),freqs.fortran_vec(),freqs.length(),n,myscratch);



  return octave_value(chisq);
    
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_burst_real_chisq_cached_c, args, nargout, "Calculate chisq for an FRB, brute-force, data should already be cached.\n")
{
  if (args.length()<7) {
    printf("Need 7 args to get_burst_real_chisq_cached_c.\n");
    return octave_value_list();
  }
  ColumnVector guess=args(0).column_vector_value();
  

  void *datft=get_pointer(args(1));
  ColumnVector freqs=args(2).column_vector_value();
  int n=(int)get_value(args(3));
  ColumnVector nvec=args(4).column_vector_value();
  double dt=get_value(args(5));
  void *myscratch=get_pointer(args(6));


  //double chisq=calculate_chisq(dat.fortran_vec(),nvec.fortran_vec(),imin,dt,fwhm,scat,alpha,amp,t0,DM,freqs.fortran_vec(),freqs.length(),n,mat.fortran_vec());
  double chisq=calculate_chisq_real_cached(datft,nvec.fortran_vec(),dt,guess.fortran_vec(),freqs.fortran_vec(),freqs.length(),n,myscratch);
  


  return octave_value(chisq);
    
}
/*--------------------------------------------------------------------------------*/


DEFUN_DLD (get_burst_real_chisq_linfit_cached_c, args, nargout, "Calculate chisq for an FRB, brute-force, data should already be cached.\n")
{
  if (args.length()<7) {
    printf("Need 7 args to get_burst_real_chisq_linfit_cached_c.\n");
    return octave_value_list();
  }
  ColumnVector guess=args(0).column_vector_value();
  

  void *datft=get_pointer(args(1));
  ColumnVector freqs=args(2).column_vector_value();
  int n=(int)get_value(args(3));
  ColumnVector nvec=args(4).column_vector_value();
  double dt=get_value(args(5));
  void *myscratch=get_pointer(args(6));


  //double chisq=calculate_chisq(dat.fortran_vec(),nvec.fortran_vec(),imin,dt,fwhm,scat,alpha,amp,t0,DM,freqs.fortran_vec(),freqs.length(),n,mat.fortran_vec());
  double chisq=calculate_chisq_linfit_real_cached(datft,nvec.fortran_vec(),dt,guess.fortran_vec(),freqs.fortran_vec(),freqs.length(),n,myscratch);
  


  return octave_value(chisq);
    
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_burst_real_chisq_qu_cached_c, args, nargout, "Calculate polarized chisq for an FRB, brute-force, data should already be cached.\n")
{
  if (args.length()<8) {
    printf("Need 8 args to get_burst_real_chisq_qu_cached_c.\n");
    return octave_value_list();
  }
  ColumnVector guess=args(0).column_vector_value();
  

  void *q=get_pointer(args(1));
  void *u=get_pointer(args(2));
  ColumnVector freqs=args(3).column_vector_value();
  int n=(int)get_value(args(4));
  ColumnVector nvec=args(5).column_vector_value();
  double dt=get_value(args(6));
  void *myscratch=get_pointer(args(7));


  //double chisq=calculate_chisq(dat.fortran_vec(),nvec.fortran_vec(),imin,dt,fwhm,scat,alpha,amp,t0,DM,freqs.fortran_vec(),freqs.length(),n,mat.fortran_vec());
  double chisq=calculate_chisq_qu_real_cached(q,u,nvec.fortran_vec(),dt,guess.fortran_vec(),freqs.fortran_vec(),freqs.length(),n,myscratch);
  


  return octave_value(chisq);
    
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_burst_real_chisq_qu_ramp_cached_c, args, nargout, "Calculate polarized chisq for an FRB, brute-force, with a phase gradient w.r.t. time.  data should already be cached.\n")
{
  if (args.length()<8) {
    printf("Need 8 args to get_burst_real_chisq_qu_ramp_cached_c.\n");
    return octave_value_list();
  }
  ColumnVector guess=args(0).column_vector_value();
  

  void *q=get_pointer(args(1));
  void *u=get_pointer(args(2));
  ColumnVector freqs=args(3).column_vector_value();
  int n=(int)get_value(args(4));
  ColumnVector nvec=args(5).column_vector_value();
  double dt=get_value(args(6));
  void *myscratch=get_pointer(args(7));


  //double chisq=calculate_chisq(dat.fortran_vec(),nvec.fortran_vec(),imin,dt,fwhm,scat,alpha,amp,t0,DM,freqs.fortran_vec(),freqs.length(),n,mat.fortran_vec());
  double chisq=calculate_chisq_qu_ramp_real_cached(q,u,nvec.fortran_vec(),dt,guess.fortran_vec(),freqs.fortran_vec(),freqs.length(),n,myscratch);
  


  return octave_value(chisq);
    
}
/*--------------------------------------------------------------------------------*/


DEFUN_DLD (get_burst_real_chisq_qu_ramp_conv_cached_c, args, nargout, "Calculate polarized chisq for an FRB, brute-force, with a phase gradient w.r.t. time.  data should already be cached.\n")
{
  if (args.length()<8) {
    printf("Need 8 args to get_burst_real_chisq_qu_ramp_conv_cached_c.\n");
    return octave_value_list();
  }
  ColumnVector guess=args(0).column_vector_value();
  
  //printf("getting ready to get inputs.\n");
  void *q=get_pointer(args(1));
  void *u=get_pointer(args(2));
  ColumnVector freqs=args(3).column_vector_value();
  int n=(int)get_value(args(4));
  ColumnVector nvec=args(5).column_vector_value();
  double dt=get_value(args(6));
  void *myscratch=get_pointer(args(7));
  //printf("got em.\n");


  //double chisq=calculate_chisq(dat.fortran_vec(),nvec.fortran_vec(),imin,dt,fwhm,scat,alpha,amp,t0,DM,freqs.fortran_vec(),freqs.length(),n,mat.fortran_vec());
  double chisq=calculate_chisq_qu_ramp_conv_real_cached(q,u,nvec.fortran_vec(),dt,guess.fortran_vec(),freqs.fortran_vec(),freqs.length(),n,myscratch);
  //printf("back with chisq %12.5f\n",chisq);


  return octave_value(chisq);
    
}
/*--------------------------------------------------------------------------------*/


DEFUN_DLD (get_burst_real_chisq_qu_rmpow_cached_c, args, nargout, "Calculate polarized chisq for an FRB, brute-force, data should already be cached.\n")
{
  if (args.length()<8) {
    printf("Need 8 args to get_burst_real_chisq_qu_rmpow_real_cached_c.\n");
    return octave_value_list();
  }
  ColumnVector guess=args(0).column_vector_value();
  

  void *q=get_pointer(args(1));
  void *u=get_pointer(args(2));
  ColumnVector freqs=args(3).column_vector_value();
  int n=(int)get_value(args(4));
  ColumnVector nvec=args(5).column_vector_value();
  double dt=get_value(args(6));
  void *myscratch=get_pointer(args(7));


  //double chisq=calculate_chisq(dat.fortran_vec(),nvec.fortran_vec(),imin,dt,fwhm,scat,alpha,amp,t0,DM,freqs.fortran_vec(),freqs.length(),n,mat.fortran_vec());
  double chisq=calculate_chisq_qu_rmpow_real_cached(q,u,nvec.fortran_vec(),dt,guess.fortran_vec(),freqs.fortran_vec(),freqs.length(),n,myscratch);
  


  return octave_value(chisq);
    
}
/*--------------------------------------------------------------------------------*/


DEFUN_DLD (get_burst_real_general_chisq_cached_c, args, nargout, "Calculate chisq for an FRB, brute-force, data should already be cached.\n")
{
  if (args.length()<7) {
    printf("Need 7 args to get_burst_real_chisq_cached_c.\n");
    return octave_value_list();
  }
  ColumnVector guess=args(0).column_vector_value();
  

  void *datft=get_pointer(args(1));
  ColumnVector freqs=args(2).column_vector_value();
  int n=(int)get_value(args(3));
  ColumnVector nvec=args(4).column_vector_value();
  double dt=get_value(args(5));
  void *myscratch=get_pointer(args(6));


  //double chisq=calculate_chisq(dat.fortran_vec(),nvec.fortran_vec(),imin,dt,fwhm,scat,alpha,amp,t0,DM,freqs.fortran_vec(),freqs.length(),n,mat.fortran_vec());
  double chisq=calculate_chisq_general_real_cached(datft,nvec.fortran_vec(),dt,guess.fortran_vec(),freqs.fortran_vec(),freqs.length(),n,myscratch);
  


  return octave_value(chisq);
    
}
/*--------------------------------------------------------------------------------*/


DEFUN_DLD (get_burst_real_general_nup_chisq_cached_c, args, nargout, "Calculate chisq for an FRB, brute-force, data should already be cached.\n")
{
  if (args.length()<7) {
    printf("Need 7 args to get_burst_real_chisq_cached_c.\n");
    return octave_value_list();
  }
  ColumnVector guess=args(0).column_vector_value();
  

  void *datft=get_pointer(args(1));
  ColumnVector freqs=args(2).column_vector_value();
  int n=(int)get_value(args(3));
  ColumnVector nvec=args(4).column_vector_value();
  double dt=get_value(args(5));
  void *myscratch=get_pointer(args(6));


  //double chisq=calculate_chisq(dat.fortran_vec(),nvec.fortran_vec(),imin,dt,fwhm,scat,alpha,amp,t0,DM,freqs.fortran_vec(),freqs.length(),n,mat.fortran_vec());
  double chisq=calculate_chisq_general_nup_real_cached(datft,nvec.fortran_vec(),dt,guess.fortran_vec(),freqs.fortran_vec(),freqs.length(),n,myscratch);
  


  return octave_value(chisq);
    
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_burst_chisq_scatfit_cached_c, args, nargout, "Calculate chisq for an FRB with arbitrary scattering power law, brute-force, data should already be cached.\n")
{
  if (args.length()<8) {
    printf("Need 8 args to get_burst_chisq_scatfit_cached_c.\n");
    return octave_value_list();
  }
  FloatColumnVector guess=args(0).float_column_vector_value();


  void *datft=get_pointer(args(1));
  FloatColumnVector freqs=args(2).float_column_vector_value();
  int n=(int)get_value(args(3));
  FloatColumnVector nvec=args(4).float_column_vector_value();
  float dt=(float)get_value(args(5));
  int imin=(int)get_value(args(6));
  void *myscratch=get_pointer(args(7));

  //double chisq=calculate_chisq(dat.fortran_vec(),nvec.fortran_vec(),imin,dt,fwhm,scat,alpha,amp,t0,DM,freqs.fortran_vec(),freqs.length(),n,mat.fortran_vec());
  double chisq=calculate_chisq_scatfit_cached(datft,nvec.fortran_vec(),imin,dt,guess.fortran_vec(),freqs.fortran_vec(),freqs.length(),n,myscratch);



  return octave_value(chisq);
    
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_burst_chisq_linfit_cached_c, args, nargout, "Calculate chisq for an FRB, brute-force, data should already be cached.\n")
{
  if (args.length()<8) {
    printf("Need 8 args to get_burst_chisq_linfit_cached_c.\n");
    return octave_value_list();
  }
  FloatColumnVector guess=args(0).float_column_vector_value();


  void *datft=get_pointer(args(1));
  FloatColumnVector freqs=args(2).float_column_vector_value();
  int n=(int)get_value(args(3));
  FloatColumnVector nvec=args(4).float_column_vector_value();
  float dt=(float)get_value(args(5));
  int imin=(int)get_value(args(6));
  void *myscratch=get_pointer(args(7));

  //double chisq=calculate_chisq(dat.fortran_vec(),nvec.fortran_vec(),imin,dt,fwhm,scat,alpha,amp,t0,DM,freqs.fortran_vec(),freqs.length(),n,mat.fortran_vec());
  double chisq=calculate_chisq_linfit_cached(datft,nvec.fortran_vec(),imin,dt,guess.fortran_vec(),freqs.fortran_vec(),freqs.length(),n,myscratch);



  return octave_value(chisq);
    
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_burst_chisq_dmpow_cached_c, args, nargout, "Calculate chisq for an FRB, brute-force, data should already be cached.\n")
{
  if (args.length()<8) {
    printf("Need 8 args to get_burst_chisq_cached_c.\n");
    return octave_value_list();
  }
  FloatColumnVector guess=args(0).float_column_vector_value();


  void *datft=get_pointer(args(1));
  FloatColumnVector freqs=args(2).float_column_vector_value();
  int n=(int)get_value(args(3));
  FloatColumnVector nvec=args(4).float_column_vector_value();
  float dt=(float)get_value(args(5));
  int imin=(int)get_value(args(6));
  void *myscratch=get_pointer(args(7));

  //double chisq=calculate_chisq(dat.fortran_vec(),nvec.fortran_vec(),imin,dt,fwhm,scat,alpha,amp,t0,DM,freqs.fortran_vec(),freqs.length(),n,mat.fortran_vec());
  double chisq=calculate_chisq_cached_dmpow(datft,nvec.fortran_vec(),imin,dt,guess.fortran_vec(),freqs.fortran_vec(),freqs.length(),n,myscratch);



  return octave_value(chisq);
    
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_burst_chisq_tau_cached_c, args, nargout, "Calculate chisq for an FRB, brute-force, data should already be cached.\n")
{
  if (args.length()<8) {
    printf("Need 8 args to get_burst_chisq_cached_c.\n");
    return octave_value_list();
  }
  FloatColumnVector guess=args(0).float_column_vector_value();


  void *datft=get_pointer(args(1));
  FloatColumnVector freqs=args(2).float_column_vector_value();
  int n=(int)get_value(args(3));
  FloatColumnVector nvec=args(4).float_column_vector_value();
  float dt=(float)get_value(args(5));
  int imin=(int)get_value(args(6));
  void *myscratch=get_pointer(args(7));

  //double chisq=calculate_chisq(dat.fortran_vec(),nvec.fortran_vec(),imin,dt,fwhm,scat,alpha,amp,t0,DM,freqs.fortran_vec(),freqs.length(),n,mat.fortran_vec());
  double chisq=calculate_chisq_cached_tau(datft,nvec.fortran_vec(),imin,dt,guess.fortran_vec(),freqs.fortran_vec(),freqs.length(),n,myscratch);



  return octave_value(chisq);
    
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_burst_chisq_qu_cached_c, args, nargout, "Calculate polarized chisq for an FRB, brute-force, data should already be cached.\n")
{
  if (args.length()<9) {
    printf("Need 9 args to get_burst_chisq_qu_cached_c.\n");
    return octave_value_list();
  }
  FloatColumnVector guess=args(0).float_column_vector_value();


  void *qft=get_pointer(args(1));
  void *uft=get_pointer(args(2));
  FloatColumnVector freqs=args(3).float_column_vector_value();
  int n=(int)get_value(args(4));
  FloatColumnVector nvec=args(5).float_column_vector_value();
  float dt=(float)get_value(args(6));
  int imin=(int)get_value(args(7));
  void *myscratch=get_pointer(args(8));

  //double chisq=calculate_chisq(dat.fortran_vec(),nvec.fortran_vec(),imin,dt,fwhm,scat,alpha,amp,t0,DM,freqs.fortran_vec(),freqs.length(),n,mat.fortran_vec());
  double chisq=calculate_chisq_qu_cached(qft,uft,nvec.fortran_vec(),imin,dt,guess.fortran_vec(),freqs.fortran_vec(),freqs.length(),n,myscratch);



  return octave_value(chisq);
    
}
/*--------------------------------------------------------------------------------*/


DEFUN_DLD (get_burst_chisq_qu_rmpow_cached_c, args, nargout, "Calculate polarized chisq for an FRB with non-standard rotation measure response, brute-force.  data should already be cached.\n")
{
  if (args.length()<9) {
    printf("Need 9 args to get_burst_chisq_qu_cached_c.\n");
    return octave_value_list();
  }
  FloatColumnVector guess=args(0).float_column_vector_value();


  void *qft=get_pointer(args(1));
  void *uft=get_pointer(args(2));
  FloatColumnVector freqs=args(3).float_column_vector_value();
  int n=(int)get_value(args(4));
  FloatColumnVector nvec=args(5).float_column_vector_value();
  float dt=(float)get_value(args(6));
  int imin=(int)get_value(args(7));
  void *myscratch=get_pointer(args(8));

  double chisq=calculate_chisq_qu_rmpow_cached(qft,uft,nvec.fortran_vec(),imin,dt,guess.fortran_vec(),freqs.fortran_vec(),freqs.length(),n,myscratch);

  return octave_value(chisq);
    
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_burst_chisq_ampfit_c, args, nargout, "Calculate chisq for an FRB, brute-force\n")
{
  if (args.length()<7) {
    printf("Need 7 args to get_burst_chisq_c.\n");
    return octave_value_list();
  }
  FloatColumnVector guess=args(0).float_column_vector_value();

  float fwhm=guess(0);
  float scat=1.0/guess(1);
  float alpha=guess(2);
  float amp=guess(3);
  float t0=guess(4);
  float DM=guess(5);

  
  FloatComplexMatrix dat=args(1).float_complex_matrix_value();;

  //dim_vector dm=dat.dims();
  //FloatComplexMatrix mat(dm);
  FloatColumnVector freqs=args(2).float_column_vector_value();
  int n=(int)get_value(args(3));
  FloatColumnVector nvec=args(4).float_column_vector_value();
  float dt=(float)get_value(args(5));
  int imin=(int)get_value(args(6));
  
  //double chisq=calculate_chisq(dat.fortran_vec(),nvec.fortran_vec(),imin,dt,fwhm,scat,alpha,amp,t0,DM,freqs.fortran_vec(),freqs.length(),n,mat.fortran_vec());
  double chisq=calculate_chisq_ampfit(dat.fortran_vec(),nvec.fortran_vec(),imin,dt,fwhm,scat,alpha,&amp,t0,DM,freqs.fortran_vec(),freqs.length(),n);


  octave_value_list retval;
  retval(0)=chisq;
  retval(1)=amp;
  return retval;
  //return octave_value(chisq);
    
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_burst_chisq_qu_c, args, nargout, "Calculate chisq for an FRB, brute-force\n")
{
  if (args.length()<8) {
    printf("Need 8 args to get_burst_chisq_qu_c.\n");
    return octave_value_list();
  }
  FloatColumnVector guess=args(0).float_column_vector_value();

  
  FloatComplexMatrix dat_q=args(1).float_complex_matrix_value();;
  FloatComplexMatrix dat_u=args(2).float_complex_matrix_value();;

  //dim_vector dm=dat.dims();
  //FloatComplexMatrix mat(dm);
  FloatColumnVector freqs=args(3).float_column_vector_value();
  int n=(int)get_value(args(4));
  FloatColumnVector nvec=args(5).float_column_vector_value();
  float dt=(float)get_value(args(6));
  int imin=(int)get_value(args(7));
  
  //double chisq=calculate_chisq(dat.fortran_vec(),nvec.fortran_vec(),imin,dt,fwhm,scat,alpha,amp,t0,DM,freqs.fortran_vec(),freqs.length(),n,mat.fortran_vec());
  double chisq=calculate_chisq_qu(dat_q.fortran_vec(),dat_u.fortran_vec(),nvec.fortran_vec(),imin,dt,guess.fortran_vec(),freqs.fortran_vec(),freqs.length(),n);


  octave_value_list retval;
  retval(0)=chisq;
  //retval(1)=amp;
  return retval;
  //return octave_value(chisq);
    
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_burst_chisq_qu_rmpow_c, args, nargout, "Calculate chisq for an FRB, brute-force, power law of rotation measure can be !=2\n")
{
  if (args.length()<8) {
    printf("Need 8 args to get_burst_chisq_qu_rmpow_c.\n");
    return octave_value_list();
  }
  FloatColumnVector guess=args(0).float_column_vector_value();

  
  FloatComplexMatrix dat_q=args(1).float_complex_matrix_value();;
  FloatComplexMatrix dat_u=args(2).float_complex_matrix_value();;

  //dim_vector dm=dat.dims();
  //FloatComplexMatrix mat(dm);
  FloatColumnVector freqs=args(3).float_column_vector_value();
  int n=(int)get_value(args(4));
  FloatColumnVector nvec=args(5).float_column_vector_value();
  float dt=(float)get_value(args(6));
  int imin=(int)get_value(args(7));
  
  //double chisq=calculate_chisq(dat.fortran_vec(),nvec.fortran_vec(),imin,dt,fwhm,scat,alpha,amp,t0,DM,freqs.fortran_vec(),freqs.length(),n,mat.fortran_vec());
  double chisq=calculate_chisq_qu_rmpow(dat_q.fortran_vec(),dat_u.fortran_vec(),nvec.fortran_vec(),imin,dt,guess.fortran_vec(),freqs.fortran_vec(),freqs.length(),n);

  return octave_value(chisq);
    
}


/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_burst_many_chisq_qu_c, args, nargout, "Calculate pol chisq for an FRB over fast parameters.  \n")
{
  if (args.length()<10) {
    printf("Need 10 args to get_burst_many_chisq_qu_c.\n");
    return octave_value_list();
  }
  FloatColumnVector guess=args(0).float_column_vector_value();

  
  FloatComplexMatrix dat_q=args(1).float_complex_matrix_value();;
  FloatComplexMatrix dat_u=args(2).float_complex_matrix_value();;

  //dim_vector dm=dat.dims();
  //FloatComplexMatrix mat(dm);
  FloatColumnVector freqs=args(3).float_column_vector_value();
  int n=(int)get_value(args(4));
  FloatColumnVector nvec=args(5).float_column_vector_value();
  float dt=(float)get_value(args(6));
  int imin=(int)get_value(args(7));
  FloatColumnVector ampvec=args(8).float_column_vector_value();
  FloatColumnVector RMvec=args(9).float_column_vector_value();
  FloatColumnVector phivec=args(10).float_column_vector_value();
  int nchi=ampvec.length();
  ColumnVector chivec(nchi);

  //double chisq=calculate_chisq(dat.fortran_vec(),nvec.fortran_vec(),imin,dt,fwhm,scat,alpha,amp,t0,DM,freqs.fortran_vec(),freqs.length(),n,mat.fortran_vec());
  //double chisq=calculate_chisq_qu(dat_q.fortran_vec(),dat_u.fortran_vec(),nvec.fortran_vec(),imin,dt,guess.fortran_vec(),freqs.fortran_vec(),freqs.length(),n);
  calculate_many_chisq_qu(dat_q.fortran_vec(),dat_u.fortran_vec(),nvec.fortran_vec(),imin,dt,guess.fortran_vec(),freqs.fortran_vec(),freqs.length(),n,ampvec.fortran_vec(),RMvec.fortran_vec(),phivec.fortran_vec(),nchi,chivec.fortran_vec());
  
  
  
  return octave_value(chivec);

  //octave_value_list retval;
  //retval(0)=chisq;
  //retval(1)=amp;
  //return retval;
  //return octave_value(chisq);
    
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (circshift_mat,args,nargout,"Circularly shift each column of a matrix.\n")
{
  FloatMatrix dat=args(0).float_matrix_value();
  dim_vector dm=dat.dims();
  ColumnVector shifts=args(1).column_vector_value();
  int ncol=dm(1);
  int nrow=dm(0);
  //printf("Have %d %d rows and columns.\n",nrow,ncol);
  int *shft=(int *)malloc(sizeof(int)*ncol);
  double *shiftptr=shifts.fortran_vec();
  for (int i=0;i<ncol;i++) {
    shft[i]=shiftptr[i];
    if (fabs(shft[i]-shiftptr[i])>1e-10) {
      printf("non-integer shifting going on in circshift_mat.\n");
      return octave_value_list();
    }
  }
  
  float *tmp=(float *)malloc(sizeof(float)*nrow);
  float *dd=dat.fortran_vec();
  
  
  for (int i=0;i<ncol;i++) {    
    int myshift=shft[i];
    while (myshift<0)
      myshift+=nrow;
    while (myshift>=nrow)
      myshift -=nrow;
    if (myshift>0) {  //don't worry about doing anything if shift=0;
      for (int j=0;j<nrow-myshift;j++) {
	//printf("assigning %d to %d from shift %d\n",j,j+myshift,myshift);
	tmp[j+myshift]=dd[nrow*i+j];
      }
      for (int j=nrow-myshift;j<nrow;j++) {
	//printf("assigning %d to %d from shift %d\n",j,j+myshift-nrow,myshift);
	tmp[j+myshift-nrow]=dd[nrow*i+j];
      }
      memcpy(dd+nrow*i,tmp,nrow*sizeof(float));
      
    }

  }
  
  free(tmp);
  free(shft);
  return octave_value(dat);

}
/*--------------------------------------------------------------------------------*/
