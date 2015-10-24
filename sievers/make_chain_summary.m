if ~exist('vanilla')
  dr=['../chains/'];
  nfac=1.025;
  nu_ref=764.2;
  tsamp=1.024;
  polchains=read_chains([dr 'chain_qu_real_newdat.txt'],0.2,'');polchains=polchains(:,3:end);
  %polchains(:,7)=polchains(:,7)/2;
  dmpow=read_chains([dr 'chain_tt_real_dmpow_newdat.txt'],0.2,'');dmpow=dmpow(:,3:end);
  vanilla=read_chains([dr 'chain_tt_real_newdat.txt'],0.2,'');vanilla=vanilla(:,3:end);
  circ=read_chains([dr '/chain_v_real_newdat_final.txt'],0.2,'');circ=circ(:,3:end);
  ramp=read_chains([dr '/chain_qu_ramp_real_newdat.txt'],0.2,'');ramp=ramp(:,3:end);
  ramp2=read_chains([dr '/chain_qu_ramp_conv_real_newdat_good.txt'],0.2,'');ramp2=ramp2(:,3:end);
  rmpow=read_chains([dr '/chain_qu_rmpow_real_newdat.txt'],0.2,'');rmpow=rmpow(:,3:end);
  scatpow=read_chains([dr '/chain_tt_real_scatpow_newdat.txt'],0.2,'');scatpow=scatpow(:,3:end);
end


polmean=mean(polchains);polerr=std(polchains)*nfac;
mycov=cov(polchains);save('-ascii','initial_conditions/vanilla_pol_cov.txt','mycov');
mm=mean(polchains);save('-ascii','initial_conditions/vanilla_pol_mean.txt','mm');
disp('pol fit:')
disp([polmean' polerr'])


disp('dm power: ');
disp([mean(dmpow)' std(dmpow)'*nfac]);
mycov=cov(dmpow);save('-ascii','initial_conditions/dm_pow_cov.txt','mycov');
mm=mean(dmpow);save('-ascii','initial_conditions/dm_pow_mean.txt','mm');


disp('vanilla :');
disp([mean(vanilla)' std(vanilla)'*nfac]);
scat=(vanilla(:,2));
mycov=cov(vanilla);save('-ascii','initial_conditions/vanilla_cov.txt','mycov');
mm=mean(vanilla);save('-ascii','initial_conditions/vanilla_mean.txt','mm');

disp(['scattering at 800 is ' num2str([mean(scat) std(scat)*nfac]*((nu_ref/800)^4)*tsamp*1000)]);

f800=vanilla(:,4).*(800/nu_ref).^vanilla(:,3);f800=f800/tsamp;
disp(['f800 is ' num2str([mean(f800) std(f800)*nfac])])
disp(['peak at 800 is ' num2str([mean(f800) std(f800)*nfac]/sqrt(8*log(2)))])

nn=min(length(polchains),length(vanilla));
rat=polchains(1:nn,4)./vanilla(1:nn,4);
disp(['pol fraction:'])
disp(num2str([mean(rat) std(rat)*nfac]*100))

disp('circular params:')
disp([mean(circ)' std(circ)'*nfac])
nn=min(length(circ),length(vanilla));
rat=circ(1:nn,4)./vanilla(1:nn,4);
disp(['circ fraction:'])
disp(num2str([mean(rat) std(rat)*nfac]*100))
mycov=cov(circ);save('-ascii','initial_conditions/circ_cov.txt','mycov');
mm=mean(circ);save('-ascii','initial_conditions/circ_mean.txt','mm');

disp(['phase ramp params:'])
disp([mean(ramp)' std(ramp)'*nfac])
disp(['phase ramp is ' num2str([mean(ramp(:,end)) std(ramp(:,end))*nfac]/tsamp)])
disp(['convolved phase ramp params:'])
disp([mean(ramp2)' std(ramp2)'*nfac])
mycov=cov(ramp);save('-ascii','initial_conditions/phase_ramp_cov.txt','mycov');
mycov=cov(ramp2);save('-ascii','initial_conditions/phase_ramp_convolved_cov.txt','mycov');
mm=mean(ramp);save('-ascii','initial_conditions/phase_ramp_mean.txt','mm');
mm=mean(ramp2);save('-ascii','initial_conditions/phase_ramp_convolved_mean.txt','mm');


disp(['rmpow params:'])
disp([mean(rmpow)' std(rmpow)'*nfac])
mycov=cov(rmpow);save('-ascii','initial_conditions/rm_pow_cov.txt','mycov');
mm=mean(rmpow);save('-ascii','initial_conditions/rm_pow_mean.txt','mm');
disp('scatpow params:')
disp([mean(scatpow)' std(scatpow)'*nfac])
mycov=cov(scatpow);save('-ascii','initial_conditions/scat_pow_cov.txt','mycov');
mm=mean(scatpow);save('-ascii','initial_conditions/scat_pow_mean.txt','mm');