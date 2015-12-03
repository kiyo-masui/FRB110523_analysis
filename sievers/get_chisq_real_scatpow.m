function[chisq]=get_chisq_real_scatpow(guess,dataptr,freqs,n,nvec,dt,model);
big_guess=zeros(8,1);
%big_guess(1:length(guess))=guess;
big_guess(1:2)=guess(1:2);
big_guess(3)=guess(end);
big_guess(4:7)=guess(3:end-1);
big_guess(end)=-2.0;  %DM index
chisq=get_burst_real_general_chisq_cached_c(big_guess,dataptr,freqs,n,nvec,dt,model);