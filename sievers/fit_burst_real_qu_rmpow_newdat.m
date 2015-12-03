try
  mpi_init;
end
more off
format short g

if ~exist('dat')
  dat_raw=read_npy('filtered.npy');
  %dat_raw=read_npy('calibrated.npy');
  freq=read_npy('freq.npy');
  dat=squeeze(dat_raw(:,1,:));
  dat_q=squeeze(dat_raw(:,2,:));
  dat_u=squeeze(dat_raw(:,3,:));
  clear dat_raw
  tvec=read_npy('time.npy');
  dt=median(diff(tvec));
  good_chan=sum(isnan(dat))==0;
  dat=dat(:,good_chan);
  dat_q=-1*dat_q(:,good_chan);
  dat_u=dat_u(:,good_chan);
  freq_use=freq(good_chan);


  best_guess=load('initial_conditions/rm_pow_mean.txt');
  mycov=load('initial_conditions/rm_pow_cov.txt');

   
end

myopts.nstep=500000;
myopts.outroot='chains/chain_qu_rmpow_real_newdat.txt';
myopts.noise_scale=1.025;
myopts.func=@get_burst_real_chisq_qu_rmpow_cached_c;
[pp,ll,nvec]=mcmc_burst_real_qu(dat_q,dat_u,myopts,freq_use,best_guess,0.7*mycov,dt);
