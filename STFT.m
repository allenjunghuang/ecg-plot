%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: RRI time-domain and frequency-domain featrues.
% The algorithm is refer to Piia Kaikkonen, (2007) Heart rate variabilty
% dynamics during early recovery after different endurance exercises
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all;

sample_rate=512;
rri_fs=5;
testee={'fat_1','james_1','how','allen','wang_1','hero','pu'};
datnum=[8 7 9 9 9 9 9];
% testee={'fat','james','how','allen','wang','hero','pu'};
% datnum=[11 10 10 10 10 9 9];
yearold=[24 24 24 27 23 24 24];
raw=struct; MA=struct;

for ss=1:length(testee)
    bb=num2str(ss);
    name=testee{1,ss};
    age=yearold(ss);
    
    for pp=1:datnum(ss)
        qq=num2str(pp);
        rri_data=[]; rri_interp=[]; ac_rri=[]; rri_seq=[]; stft=[]; dn_mag=[];
        
        load(['H:\Lab\temp\stepinfo\' name '\' qq '_stepinfo.mat']);
        load(['H:\Lab\lab_data\2013_data\ECG\' name '\constant_speed\processed_data\5hz(linear)\' qq '_inter_RRI.mat']);
        
        %%%%%%%%%% HRV %%%%%%%%%%
        % frequency-domain
        rri_data=rri_interp*1000; %unit: msec
        ac_rri=rri_data-mean(rri_data);
        [rri_seq,dumy_hz, numX,denX]=fn_Buterworth(ac_rri,[0.04 1],[0.03 1.06],rri_fs);      
        % % rri and filter spectrum
        % rri_fft1=fft(ac_rri)./length(ac_rri);
        % [theta1, rho1]=cart2pol(real(rri_fft1),imag(rri_fft1));
        % rri_fft2=fft(rri_seq)./length(rri_seq);
        % [theta2, rho2]=cart2pol(real(rri_fft2),imag(rri_fft2));
        % [hX, wX]=freqz(numX,denX,length(ac_rri),rri_fs);
        % figure(900)
        % plotyy((1:length(rho1)/2)/length(rho1)*rri_fs,rho1(2:length(rho1)/2+1),wX, (abs(hX)));hold on;
        % plot((1:length(rho2)/2)/length(rho2)*rri_fs,rho2(2:length(rho2)/2+1),'r');
        % legend('before','after','filter');
        % figure(901)
        % plot(ac_rri,'b'); hold on; plot(rri_seq,'r');              
        wlen=512;
        N0=wlen;
        [stft, fhz, tsec] = fn_stft(rri_seq, wlen, 1, wlen, rri_fs);
        dn_mag=log(abs(stft));
        
        for rr=1:size(stft,2);
            HF(rr)=0; LF(rr)=0; TP(rr)=0; VLF(rr)=0;
            for kk=1:N0/2
                if(rri_fs*(kk-1)/N0 >= 0.15 && rri_fs*(kk-1)/N0 < 1) % high frequrncy power
                    HF(rr)=HF(rr)+dn_mag(kk,rr)^2*2;
                end
                if(rri_fs*(kk-1)/N0 >= 0.04 && rri_fs*(kk-1)/N0 < 0.15) % low frequrncy power
                    LF(rr)=LF(rr)+dn_mag(kk,rr)^2*2;
                end
                if(rri_fs*(kk-1)/N0 < 0.04) % very low frequrncy power
                    VLF(rr)=VLF(rr)+dn_mag(kk,rr)^2*2;
                end
                if(rri_fs*(kk-1)/N0 >= 0.04 && rri_fs*(kk-1)/N0 < 1) % total power
                    TP(rr)=TP(rr)+dn_mag(kk,rr)^2*2;
                end
            end
            LFn(rr)=LF(rr)/TP(rr);
            HFn(rr)=HF(rr)/TP(rr);
            LFHFratio(rr)=LFn(rr)/HFn(rr);
        end
       
        % time-domain
        for rr=0:length(rri_data)-wlen
            rri_window=[]; rri_diff=[]; 
            rri_window=rri_data(rr+1:rr+wlen);
            rri_diff=diff(rri_window);
            AVNN(rr+1)=mean(rri_window); %RRI的平均
            rAVNN(rr+1)=1000*60./AVNN(rr+1); %RRI的平均
            SDNN(rr+1)=std(rri_window); %HRV的整體變異
            rMSSD(rr+1)=sqrt(sum(rri_diff.^2)); %HRV中高頻變異
            SDANN=[]; %HRV中低頻變異
            SDANNI=[];
            
            NN50=[]; NN20=[]; NN10=[]; NN05=[]; NN03=[]; NN01=[];
%             NN50=find(rri_diff>50); %相鄰RRI差異>50ms
%             NN20=find(rri_diff>20); %相鄰RRI差異>20ms
            NN10=find(rri_diff>10); %相鄰RRI差異>10ms
            NN05=find(rri_diff>5); %相鄰RRI差異>5ms
            NN03=find(rri_diff>3); %相鄰RRI差異>3ms
            NN01=find(rri_diff>1); %相鄰RRI差異>1ms
%             pNN50(rr+1)=length(NN50)/length(rri_diff);
%             pNN20(rr+1)=length(NN20)/length(rri_diff);
            pNN10(rr+1)=length(NN10)/length(rri_diff);
            pNN05(rr+1)=length(NN05)/length(rri_diff);
            pNN03(rr+1)=length(NN03)/length(rri_diff);
            pNN01(rr+1)=length(NN01)/length(rri_diff);         
        end
        
        Entire.TP{ss,pp}=TP;
        Entire.HF{ss,pp}=HF;
        Entire.LF{ss,pp}=LF;
        Entire.HFn{ss,pp}=HFn;
        Entire.VLF{ss,pp}=VLF;
        Entire.LFn{ss,pp}=LFn;
        Entire.LFHFratio{ss,pp}=LFHFratio;
        Entire.AVNN{ss,pp}=AVNN;
        Entire.rMSSD{ss,pp}=rMSSD;
        Entire.SDNN{ss,pp}=SDNN;
        Entire.pNN10{ss,pp}=pNN10;
        Entire.pNN05{ss,pp}=pNN05;
        Entire.pNN03{ss,pp}=pNN03;
        Entire.pNN01{ss,pp}=pNN01;
        
        dumy=[]; dumy=find(time_interp<3.5*60);
        flg1=dumy(end);
        dumy=[]; dumy=find(time_interp<33.5*60);
        flg2=dumy(end);
        Run.TP{ss,pp}=fn_movingavg(TP(flg1+1-fix(wlen/2):flg2-fix(wlen/2)), time_interp(flg1+1:flg2)-(3.5*60), 30*60, 60, 30, 1);
        Run.HF{ss,pp}=fn_movingavg(HF(flg1+1-fix(wlen/2):flg2-fix(wlen/2)), time_interp(flg1+1:flg2)-(3.5*60), 30*60, 60, 30, 1);
        Run.LF{ss,pp}=fn_movingavg(LF(flg1+1-fix(wlen/2):flg2-fix(wlen/2)), time_interp(flg1+1:flg2)-(3.5*60), 30*60, 60, 30, 1);
        Run.HFn{ss,pp}=fn_movingavg(HFn(flg1+1-fix(wlen/2):flg2-fix(wlen/2)), time_interp(flg1+1:flg2)-(3.5*60), 30*60, 60, 30, 1);
        Run.VLF{ss,pp}=fn_movingavg(VLF(flg1+1-fix(wlen/2):flg2-fix(wlen/2)), time_interp(flg1+1:flg2)-(3.5*60), 30*60, 60, 30, 1);
        Run.LFn{ss,pp}=fn_movingavg(LFn(flg1+1-fix(wlen/2):flg2-fix(wlen/2)), time_interp(flg1+1:flg2)-(3.5*60), 30*60, 60, 30, 1);
        Run.LFHFratio{ss,pp}=fn_movingavg(LFHFratio(flg1+1-fix(wlen/2):flg2-fix(wlen/2)), time_interp(flg1+1:flg2)-(3.5*60), 30*60, 60, 30, 1);
        Run.AVNN{ss,pp}=fn_movingavg(AVNN(flg1+1-fix(wlen/2):flg2-fix(wlen/2)), time_interp(flg1+1:flg2)-(3.5*60), 30*60, 60, 30, 1);
        Run.rMSSD{ss,pp}=fn_movingavg(rMSSD(flg1+1-fix(wlen/2):flg2-fix(wlen/2)), time_interp(flg1+1:flg2)-(3.5*60), 30*60, 60, 30, 1);
        Run.SDNN{ss,pp}=fn_movingavg(SDNN(flg1+1-fix(wlen/2):flg2-fix(wlen/2)), time_interp(flg1+1:flg2)-(3.5*60), 30*60, 60, 30, 1);
        Run.pNN10{ss,pp}=fn_movingavg(pNN10(flg1+1-fix(wlen/2):flg2-fix(wlen/2)), time_interp(flg1+1:flg2)-(3.5*60), 30*60, 60, 30, 1);
        Run.pNN05{ss,pp}=fn_movingavg(pNN05(flg1+1-fix(wlen/2):flg2-fix(wlen/2)), time_interp(flg1+1:flg2)-(3.5*60), 30*60, 60, 30, 1);
        Run.pNN03{ss,pp}=fn_movingavg(pNN03(flg1+1-fix(wlen/2):flg2-fix(wlen/2)), time_interp(flg1+1:flg2)-(3.5*60), 30*60, 60, 30, 1);
        Run.pNN01{ss,pp}=fn_movingavg(pNN01(flg1+1-fix(wlen/2):flg2-fix(wlen/2)), time_interp(flg1+1:flg2)-(3.5*60), 30*60, 60, 30, 1);
        
        dumy=[]; dumy=find(time_interp<48*60);
        flg3=dumy(end);
        Recovery.TP{ss,pp}=fn_movingavg(TP(flg2+1-fix(wlen/2):flg3-fix(wlen/2)), time_interp(flg2+1:flg3)-(33.5*60), 15*60, 60, 30, 1);
        Recovery.HF{ss,pp}=fn_movingavg(HF(flg2+1-fix(wlen/2):flg3-fix(wlen/2)), time_interp(flg2+1:flg3)-(33.5*60), 15*60, 60, 30, 1);
        Recovery.LF{ss,pp}=fn_movingavg(LF(flg2+1-fix(wlen/2):flg3-fix(wlen/2)), time_interp(flg2+1:flg3)-(33.5*60), 15*60, 60, 30, 1);
        Recovery.VLF{ss,pp}=fn_movingavg(VLF(flg2+1-fix(wlen/2):flg3-fix(wlen/2)), time_interp(flg2+1:flg3)-(33.5*60), 15*60, 60, 30, 1);
        Recovery.HFn{ss,pp}=fn_movingavg(HFn(flg2+1-fix(wlen/2):flg3-fix(wlen/2)), time_interp(flg2+1:flg3)-(33.5*60), 15*60, 60, 30, 1);
        Recovery.LFn{ss,pp}=fn_movingavg(LFn(flg2+1-fix(wlen/2):flg3-fix(wlen/2)), time_interp(flg2+1:flg3)-(33.5*60), 15*60, 60, 30, 1);
        Recovery.LFHFratio{ss,pp}=fn_movingavg(LFHFratio(flg2+1-fix(wlen/2):flg3-fix(wlen/2)), time_interp(flg2+1:flg3)-(33.5*60), 15*60, 60, 30, 1);
        Recovery.AVNN{ss,pp}=fn_movingavg(AVNN(flg2+1-fix(wlen/2):flg3-fix(wlen/2)), time_interp(flg2+1:flg3)-(33.5*60), 15*60, 60, 30, 1);
        Recovery.rMSSD{ss,pp}=fn_movingavg(rMSSD(flg2+1-fix(wlen/2):flg3-fix(wlen/2)), time_interp(flg2+1:flg3)-(33.5*60), 15*60, 60, 30, 1);
        Recovery.SDNN{ss,pp}=fn_movingavg(SDNN(flg2+1-fix(wlen/2):flg3-fix(wlen/2)), time_interp(flg2+1:flg3)-(33.5*60), 15*60, 60, 30, 1);
        Recovery.pNN10{ss,pp}=fn_movingavg(pNN10(flg2+1-fix(wlen/2):flg3-fix(wlen/2)), time_interp(flg2+1:flg3)-(33.5*60), 15*60, 60, 30, 1);
        Recovery.pNN05{ss,pp}=fn_movingavg(pNN05(flg2+1-fix(wlen/2):flg3-fix(wlen/2)), time_interp(flg2+1:flg3)-(33.5*60), 15*60, 60, 30, 1);
        Recovery.pNN03{ss,pp}=fn_movingavg(pNN03(flg2+1-fix(wlen/2):flg3-fix(wlen/2)), time_interp(flg2+1:flg3)-(33.5*60), 15*60, 60, 30, 1);
        Recovery.pNN01{ss,pp}=fn_movingavg(pNN01(flg2+1-fix(wlen/2):flg3-fix(wlen/2)), time_interp(flg2+1:flg3)-(33.5*60), 15*60, 60, 30, 1);
        
        clear HF LF VLF TP LFn HFn LFHFratio AVNN rAVNN SDNN rMSSD pNN20 pNN10 pNN05 pNN03 pNN01
    end
    save(['H:\Lab\temp\HRV\HRV.mat'],'Entire','Run','Recovery');
    
end