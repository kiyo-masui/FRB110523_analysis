function[chisq]=get_chisq_real_dmpow(guess,dataptr,freqs,n,nvec,dt,model);
big_guess=zeros(8,1);
big_guess(1:2)=guess(1:2);
big_guess(3)=-4; %scattering index
big_guess(4:end)=guess(3:end);

freq1=764.2;
freq2=700;
%dm0=4148.808;

%work out what DM should be to preserve lag between central and edge frequencies.
%significantly helps chains converge.

dm=guess(end-1);
dmpow=guess(end);
%lag=dm0*dm*(1/freq2^2-1/freq1^2);
%lag2=dm0*(freq2^dmpow-freq1^dmpow);
%dm_use=lag/lag2;
dm_use=dm*(freq1^-2-freq2^-2)/(freq1^dmpow-freq2^dmpow);
%dm_use=dm_use^2/dm;
%disp([dm_use dm_use2])

big_guess(end-1)=dm_use;

chisq=get_burst_real_general_chisq_cached_c(big_guess,dataptr,freqs,n,nvec,dt,model);