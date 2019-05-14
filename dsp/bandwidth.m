% Function:
%   Compute bandwidth with different energy criteria
% Inputs:
%   Xn   -   input signal
%   fs   -   sample rate
%   spc  -   energy specification
% Outputs:
%   BW   -   bandwidth
% Author: J.W. Huang. NSYSU, Medical Mechatronics Lab

function [BW]=dspstat_bandwidth(Xn, fs, spc)

Xn_t=length(Xn)/fs; 
Xn=Xn-mean(Xn);
Xn_fft=fft(Xn,Xn_t*fs);
[theta, rho]=cart2pol(real(Xn_fft),imag(Xn_fft));
total_pwr=sum(rho(2:end/2+1).^2);

for pp=1:length(spc)
    pwr=0;
    for kk=1:length(rho)
        if pwr<=spc(pp)*total_pwr
            pwr=pwr+rho(1+kk)^2;
            BW(pp)=kk/Xn_t;
        end
    end
end

