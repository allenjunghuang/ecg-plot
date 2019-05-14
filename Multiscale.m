%Multiscale entropy analysys
%This program content three decompostion method: EMD, DWT and
%coarse-graining.
%Author: J.W Huang. NSYSU.

clear all; close all; clc;
path(path,'H:\Lab\m_files\toolbox\complexity_toolbox');
path(path,'H:\Lab\m_files\toolbox\ncu_emd'); % EMD
path(path,'H:\Lab\m_files\m_file\functions'); %wavelet

Fs=1000;
%color noise
% path(path,'H:\Lab\m_files\toolbox\noise_toolbox\'); % color noise
% emgdata=powernoise(0, 12000, 'normalize');  %0:white noise, 1:pink noise
% cosine function
% T0=12;
% T=1/Fs;
% N0=T0/T;
% t=0:T:T*(N0-1); t=t';
% w0=[8*pi 20*pi 30*pi 50*pi];
% for hh=1:length(w0)
%     emgdata(:,hh)=cos(w0(1,hh).*t);
% end

musclelabel{1}='triceps brachii'; musclelabel{2}='biceps brachii';
musclelabel{3}='anterior deltoid'; musclelabel{4}='posterior deltoid';
musclelabel{5}='flexor radialis'; musclelabel{6}='extensor radialis';
musclelabel{7}='flexor digitorum'; musclelabel{8}='extensor digitorum';

group{1}='young';group{2}='old';
mvc=[100];
for gg=1
    for mm=1:length(mvc)
        mvcgrad=num2str(mvc(mm));
        for pp=7
            qq=num2str(pp);
            fprintf(['\nNo.' qq ' '])
            load(['H:\Lab\data\2015_data\grip_control\healthy_' group{gg} '(new)\' qq '\emg' mvcgrad 'mvc1.mat']);
            emgdata=emgdata(1:Fs*13,:);
            for kk=4;%1:size(emgdata,2)
                xn=emgdata(:,kk);
                
                pass_hz=[48 52]; stop_hz=[49 51];
                [xn,fc, num, den] = fn_buterworth(xn,pass_hz,stop_hz,Fs);
                pass_hz=[5 400]; stop_hz=[4 412];
                [xn,fc, num, den] = fn_buterworth(xn,pass_hz,stop_hz,Fs);
                xn=xn(1+0.5*Fs:12.5*Fs);
                %         dn=abs(fft(xn));
                %         figure(101)
                %         plot((1:length(dn)/2)/length(dn)*Fs,dn(2:length(dn)/2+1)); xlabel('Hz'); ylabel('Power Spectral Density');
                
                xn=(xn-mean(xn))/std(xn);
                %Empirical mode decomposition
                fprintf('EMD')
                allmode=emd(xn,0,1,9); %column 1: original data, columns 2,3... m: IMFs from high to low frequency, m+1 (last) column: residual
                allmode(:,1)=[];
                %IMF fine2coarse
                for ii=1:size(allmode,2)
                    imf_f2c=sum(allmode(end-500:end,ii:end),2);
%                     IMEn_f2c(kk,ii)=sampEn(imf_f2c,2,0.15,0,0);
                    IMEn_f2c(kk,ii)=fuzzyEn(imf_f2c,0.2*std(imf_f2c));
                    fprintf('.')
                    clear imf_f2c
                end
                %IMF coarse2fine
                for ii=1:size(allmode,2)
                    imf_c2f=sum(allmode(end-500:end,1:end+1-ii),2);
%                     IMEn_c2f(kk,ii)=sampEn(imf_c2f,2,0.15,0,0);
                    IMEn_c2f(kk,ii)=fuzzyEn(imf_c2f,0.2*std(imf_c2f));
                    fprintf('.')
                    clear imf_c2f
                end
                
                %Wavelet decomposition
                fprintf('DWT')
                wname = 'haar';
                dyadic = 4;
                [wlet passamp]= fn_1Ddwt(xn,dyadic,wname,'padding');
                
                dump=(passamp-length(xn))/(2^dyadic);
                wtmp=wlet{1,1}(1+ceil(dump/2):end-ceil(dump/2));
%                 wletEn(kk,1)=sampEn(wtmp,2,0.15,0,0);
%                 wletEn(kk,1)=fuzzyEn(wtmp(end-500:end),0.2*std(wtmp(end-500:end)));
                clear wtmp
                jj=1;
                for ii=length(wlet):-1:2
                    dump=(passamp-length(xn))/(2*jj);
                    wtmp=wlet{1,ii}(1+ceil(dump/2):end-ceil(dump/2));
                    
%                     wletEn(kk,ii)=sampEn(wtmp,2,0.15,0,0);
                    wletEn(kk,ii)=fuzzyEn(wtmp(end-500:end),0.2*std(wtmp(end-500:end)));
                    jj=jj*2;
                    fprintf('.')
                    clear wtmp
                end
                
                %Coarse-graining decomposition
                fprintf('CGD')
                for ii=1:8  %scale range
                    ss=0;
                    yn=zeros(1,fix(length(xn)/ii)); %initialize scaled vector
                    for jj=1:length(yn)
                        yn(jj)=sum(xn(ss+1:jj*ii))/ii;
                        ss=ss+ii;
                    end
%                     cgEn(kk,ii)=sampEn(yn,2,0.15,0,0); %m=2,r=0.15
%                     cgEn(kk,ii)=fuzzyEn(yn(end-500:end),0.2*std(yn(end-500:end))); %m=2,r=0.2
                    fprintf('.')
                    clear yn
                end
                clearvars -except IMEn_f2c IMEn_c2f wletEn cgEn emgdata Fs qq gg mm pp mvc mvcgrad group
            end
%             save(['H:\Lab\data\2015_data\grip_control\entropy\healthy_' group{gg} '_' mvcgrad '_' qq '.mat'], 'wletEn', 'cgEn', 'IMEn_f2c', 'IMEn_c2f');
            clearvars -except Fs mvc mvcgrad group gg mm pp
        end
        clearvars -except Fs mvc group gg mm pp
    end
end