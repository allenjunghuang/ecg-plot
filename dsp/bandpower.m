%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pwr rpwr dbpwr dbrpwr] =dspstat_bandpwr(X,fs,N0, freqband)
if size(X,2)>1
    rho=mean(X,2);
else
    rho=X;
end
E_total=0;
for kk=1:size(freqband,1)
    pwr(kk)=sum((rho(freqband(kk,1)*N0/fs+1:freqband(kk,2)*N0/fs).^2)*2);
    E_total=E_total+pwr(kk);
end

for ss=1:length(pwr)
    rpwr(ss)=pwr(ss)/E_total;
end

dbs=20*log10(rho);
E_total=0;
for kk=1:size(freqband,1)
    dbpwr(kk)=sum((dbs(freqband(kk,1)*N0/fs+1:freqband(kk,2)*N0/fs).^2)*2);
    E_total=E_total+dbpwr(kk);
end

for ss=1:length(pwr)
    dbrpwr(ss)=dbpwr(ss)/E_total;
end
