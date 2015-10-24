function[pp,ll,nvec]=mcmc_burst_real_qu_complex(q,u,myopts,freqs,guess,sigs,dt,fix_params,true_params)

outroot=get_struct_mem(myopts,'outroot','chain_tt.txt');
nstep=get_struct_mem(myopts,'nstep',3);
chifun=get_struct_mem(myopts,'func',@get_burst_real_chisq_qu_ramp_conv_cached_c);


if ~exist('dt')
  dt=1e-3;
end
n=size(q,1);

pp=zeros(nstep,numel(guess));
ll=zeros(nstep,1);
cur=guess;

noises=double(std(q));
nvec=1.0./noises.^2;

unwind_protect

qptr=float2ptr(q);
uptr=float2ptr(u);
model=floatcomplex2ptr(q+i*u);




if exist('true_params')
  guess(fix_params)=true_params(fix_params);
end


try
  myid=mpi_comm_rank+1;
  fid=fopen([outroot '_'  num2str(myid)],'w');
catch
  fid=fopen(outroot,'w');
end


tic;chisq=feval(chifun,guess,qptr,uptr,freqs,n,nvec,dt,model);toc
disp(['initial chisq is ' num2str(chisq)]);



pp(1,:)=guess;
ll(1)=chisq;
nrep=1;
aaa=now;
for j=2:nstep

  if rem(j,1000)==0
    bbb=now;
    disp([j 86400*(bbb-aaa)]);
    fflush(fid);
  end
  guess=cur+mc_step(sigs);
  if exist('true_params')
    guess(fix_params)=true_params(fix_params);
  end

  chisq=feval(chifun,guess,qptr,uptr,freqs,n,nvec,dt,model);

  if (exp(-0.5*(chisq-ll(j-1)))>rand(1))
    %disp('accepting');
    fprintf(fid,'%d ',nrep);
    fprintf(fid,'%14.7f ',ll(j-1));
    fprintf(fid,'%15.9g ',pp(j-1,:));
    fprintf(fid,'\n');
    fflush(fid);


    ll(j)=chisq;
    pp(j,:)=guess;
    cur=guess;
    nrep=1;
    %disp(['accepting with chisq ' num2str(chisq)]);
  else
    %disp('not accepting')
    %disp(['not accempting with chisq ' num2str(chisq)])
    ll(j)=ll(j-1);
    pp(j,:)=cur;
    nrep=nrep+1;
  end

end
unwind_protect_cleanup
disp('cleaning up')
freeptr(model);
freeptr(qptr);
freeptr(uptr);
fclose(fid);
end_unwind_protect

