function[dat]=read_npy(fname)

fid=fopen(fname,'r');
if (1)

  finfo=stat(fname);
  nn=fread(fid,2,'int');
  hh=fread(fid,1,'int16');
  hh=fread(fid,[1 hh],'char=>char');
  nb=8+2+numel(hh);  %size of header in bytes
  ii=find(hh=='<');
  tt=hh(ii+1:end);
  ii=find(tt=='''');
  tt=tt(1:ii(1)-1);
  i1=find(hh=='(');
  i2=find(hh==')');
  sz=str2num([ '[' hh(i1+1:i2-1) ']' ]);
  if numel(sz)==1
    sz=[1 sz ];
  end



  fmt='';
  if strcmp(tt,'f4')
    fmt='float=>float';
    nb=nb+4*prod(sz);
  end
  if strcmp(tt,'f8')
    fmt='double';
    nb=nb+8*prod(sz);
  end
  assert(~isempty(fmt));  %if this fails, we didn't find a suitable data type
  assert(nb==finfo.size)
  
  %dat=fread(fid,fliplr(sz),fmt);
  dat=fread(fid,inf,fmt);
  assert(numel(dat)==prod(sz));
  dat=reshape(dat,fliplr(sz));
  
  

else
  
  hh=repmat(' ',200);
  for j=1:length(hh)
    hh(j)=fread(fid,1,'char=>char');
    if hh(j)=='}'
      hh=hh(1:j);
      break;
    end
  end
  whos hh
end



fclose(fid);