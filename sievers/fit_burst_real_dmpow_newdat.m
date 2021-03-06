try
  mpi_init;
end
more off
format short g

if ~exist('dat')
  dat_raw=read_npy('filtered.npy');
  freq=read_npy('freq.npy');
  dat=squeeze(dat_raw(:,1,:));
  dat_q=squeeze(dat_raw(:,2,:));
  dat_u=squeeze(dat_raw(:,3,:));
  clear dat_raw
  tvec=read_npy('time.npy');
  dt=median(diff(tvec));
  good_chan=sum(isnan(dat))==0;
  dat=dat(:,good_chan);
  dat_q=dat_q(:,good_chan);
  dat_u=dat_u(:,good_chan);
  freq_use=freq(good_chan);

  mycov=load('initial_conditions/dm_pow_cov.txt');
  best_guess=load('initial_conditions/dm_pow_mean.txt');
   
end



myopts.nstep=500000;
myopts.outroot='chains/chain_tt_real_dmpow_newdat.txt';
myopts.func=@get_chisq_real_dmpow;
myopts.noise_scale=1.025;
[pp,ll,nvec]=mcmc_burst_real(dat,myopts,freq_use,best_guess,0.5*mycov,dt,2.0);

