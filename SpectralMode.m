clear all; clc; close all;
path(path,'H:\Lab\m_files\eeglab13_3_2b\functions\sigprocfunc\'); 
path(path,'H:\Lab\m_files\m_file\functions\'); 

Fs=250;
motion='f';
for pp=14
    qq=num2str(pp);   
    
%     Fs_grf=512;
%     net_grf=load(['H:\Lab\data\2014_data\eeg_ambulating\Cz\raw\naf\' qq '.txt']);
%     raw_grf=load(['H:\Lab\data\2014_data\eeg_ambulating\Cz\raw\grf\' qq '' motion '_grf.txt']);
%     for kk=1:4
%         cell_netgrf(kk)=mean(net_grf(1+Fs_grf:6*Fs_grf,kk));
%         cell_loadgrf(:,kk)=raw_grf(:,kk)-cell_netgrf(kk);
%     end
%     sum_grf=cell_loadgrf(:,1)+cell_loadgrf(:,2)+cell_loadgrf(:,3)+cell_loadgrf(:,4);
%     %FitPlus-7355 FFT modeling
%     signal_start=0;  
%     signal_time=60*10;
%     bufer=0;
%     sample_start=signal_start-bufer;
%     sample_time=signal_time+2*bufer;   
%     sample_rate=512;
%     daq_fs=512;
%     fft_fs=512;
%     N0=signal_time*sample_rate;
%     
%     load(['H:\Lab\data\calibration_data\impulse_130906\130906_black_alen10615-0054.mat']);    
%     if sample_rate~=fft_fs
%         Imp_data=resample(data_B, sample_rate, fft_fs);
%     else
%         Imp_data=data_B;
%     end
%     clear data_B;
%     input_F=Imp_data(1+sample_rate*2:sample_rate*28,1)+Imp_data(1+sample_rate*2:sample_rate*28,2)+Imp_data(1+sample_rate*2:sample_rate*28,3)+Imp_data(1+sample_rate*2:sample_rate*28,4);
%     output_F=Imp_data(1+sample_rate*2:sample_rate*28,5)+Imp_data(1+sample_rate*2:sample_rate*28,6)+Imp_data(1+sample_rate*2:sample_rate*28,7)+Imp_data(1+sample_rate*2:sample_rate*28,8);        
%     [Max_value Max_index]=max(input_F); % find out maximum value and which value. Note the fft_input data range.
%     initial_input_F=mean(input_F(1:Max_index-2*sample_rate));
%     initial_output_F=mean(output_F(1:Max_index-2*sample_rate));
%     input_F=input_F(Max_index:Max_index+sample_rate*10-1); %collect 10 seconds data after maximun value happened 
%     output_F=output_F(Max_index:Max_index+sample_rate*10-1); %find out loadcell value when maximum impulse value happened and 10 seconds later
% 
%     input_F=(input_F-initial_input_F);
%     output_F=(output_F-initial_output_F);            
%     zero(1:sample_rate*(sample_time),1)=0; %evaluate zero matrix            
%     input_F=[input_F ; zero]; % zero padding
%     output_F=[output_F ; zero];     
%     
%     input_F_fft=fft(input_F,sample_time*sample_rate);   
%     output_F_fft=fft(output_F,sample_time*sample_rate);    
%     
%     H_A_fft=output_F_fft./input_F_fft; %FFT model  
%     [H_A_phi H_A_rho]=cart2pol(real(H_A_fft),imag(H_A_fft));  
%     H_A_hz=(1:sample_time*sample_rate/2)/sample_time;
%     H_A_Dnmag=20*log10(H_A_rho);    
%     
%     output_B_fft=fft(sum_grf,sample_time*sample_rate);
%     input_B_fft=output_B_fft./H_A_fft;
%     input_B=ifft(input_B_fft);  % obtain time doamin input force
%     input_B(find(input_B<0))=0;
%     pass_hz=25; stop_hz=30;
%     [inv_grf,cutoff_hz, dumy_num,dumy_den]=fn_buterworth(input_B,pass_hz,stop_hz,sample_rate);
%     
%     w_grf=inv_grf/mean(inv_grf(1+Fs_grf*60:Fs_grf*120)); % unit: subject's body weight
%     walk_grf=w_grf(1+Fs_grf*173:Fs_grf*415);
%     
%     % spectral analysis
%     walk_grf_fft=fft(walk_grf)/length(walk_grf);
%     [walk_grf_phi walk_grf_rho]=cart2pol(real(walk_grf_fft),imag(walk_grf_fft));
%     w0_pwr(pp)=sum((walk_grf_rho(0.5*N0/Fs_grf:2.5*N0/Fs_grf).^2)*2);
%     w0_rpwr(pp)=w0_pwr(pp)/sum(walk_grf_rho(2:end).^2);
%     w1_pwr(pp)=sum((walk_grf_rho(0.5*N0/Fs_grf:4.5*N0/Fs_grf).^2)*2);
%     w1_rpwr(pp)=w1_pwr(pp)/sum(walk_grf_rho(2:end).^2);
%     
%     walk_grf_dB=20*log10(walk_grf_rho);    
%     walk_tvisum(pp)=sum(walk_grf);
%     walk_std(pp)=std(walk_grf);
%     [grf_E99BW(pp)]=fn_bandwidth(walk_grf, Fs_grf, 0.99);
%     [grf_E95BW(pp)]=fn_bandwidth(walk_grf, Fs_grf, 0.95);
%     [grf_E90BW(pp)]=fn_bandwidth(walk_grf, Fs_grf, 0.90);
%       
%     hz_axis=(1:length(walk_grf_rho)/2)/length(walk_grf_rho)*Fs_grf;
%     [pmv,pmi]=max(walk_grf_rho(2:length(walk_grf_rho)/2+1));
%     w_0(pp)=hz_axis(pmi);
%     w_1(pp)=2*w_0(pp);
%     w_2(pp)=3*w_0(pp);
%     w_3(pp)=4*w_0(pp);
%     
%     % time-domain analysis
%     [cycle_index,cycle_interval,peak_index,peak_interval,trough_index,trough_interval]=grf_walkcycle(walk_grf);
%    
% %     figure(3001)
% %     plot(walk_grf); hold on;
% %     plot(cycle_index,walk_grf(cycle_index),'g*'); hold on;
% %     plot(peak_index,walk_grf(peak_index),'r*'); hold on;
% %     plot(trough_index,walk_grf(trough_index),'r*');
% %     figure(3002)
% %     plot(cycle_interval/Fs_grf,'r'); hold on;
% %     plot(peak_interval/Fs_grf,'g'); hold on;
% %     plot(trough_interval/Fs_grf,'b'); hold on;
% %     ylabel('interval (sec)'); xlabel('step')
%     
%     mode_seq=walk_grf(cycle_index(1):cycle_index(end));
%     T_0=length(mode_seq)/Fs_grf/(length(cycle_index)-1);
%     W_0=1/T_0;
%     P_c=1000;
%     N_c=fix(length(mode_seq)*P_c/(length(cycle_index)-1));
%     if rem(N_c,2) ==1; N_c=N_c+1; end
%     mode_fft=fft(mode_seq-1,N_c)/N_c;
%     [mode_phi mode_rho]=cart2pol(real(mode_fft),imag(mode_fft));  
%     w0_axis=(1:length(mode_rho)/2)/P_c;
%     mode1_pwr=sum((mode_rho(0.5*P_c+1:1.5*P_c).^2)*2);
%     mode2_pwr=sum((mode_rho(1.5*P_c+1:2.5*P_c).^2)*2);
%     mode3_pwr=sum((mode_rho(2.5*P_c+1:3.5*P_c).^2)*2);
%     mode4_pwr=sum((mode_rho(3.5*P_c+1:4.5*P_c).^2)*2);
%     
%     figure(3003)
%     plot(0.5, [0:0.0001:0.03], 'b'); hold on; plot(1.5, [0:0.0001:0.03], 'b'); hold on; plot(2.5, [0:0.0001:0.03], 'b'); hold on; plot(3.5, [0:0.0001:0.03], 'b'); hold on; plot(4.5, [0:0.0001:0.03], 'b'); hold on
%     plot(w0_axis,mode_rho(2:length(mode_rho)/2+1),'k'); xlabel('\omega_0'); ylabel('Body Weight^2/Hz'); title('Power Spectral Density'); axis([0 5 0 0.03]);
%     
%     clf; fig1=figure(1);
%     plot((1:length(walk_grf))/Fs_grf, walk_grf); xlabel('time (sec)'); ylabel('body weight'); axis([120 123 0.75 1.5])
%     saveas(fig1 , ['H:\Lab\figures\eeg\grf\grf_sig' qq '_' motion '.png'])
%     
%     clf; fig11=figure(11);    
%     plot(grf_E90BW(pp), [0:0.0001:0.03], 'r'); hold on; plot(2.5, [0:0.0001:0.03], 'g'); hold on; plot(4.5, [0:0.0001:0.03], 'b'); hold on
%     plot(hz_axis,walk_grf_pwr(2:length(walk_grf_pwr)/2+1),'k'); 
%     xlabel('Hz'); ylabel('Body Weight^2/Hz'); title('PSD'); axis([0 12 0 0.03]);
%     saveas(fig11 , ['H:\Lab\figures\eeg\grf\grf_psd' qq '_' motion '.png'])  
    
    [X, header] = edf2mat(['H:\Lab\data\2014_data\eeg_ambulating\Cz\edf\edf\' qq '' motion '.edf']);
    epch(1,:)=[27 147]; epch(2,:)=[178 410]; epch(3,:)=[425 593];
    grf_stand=X(10,epch(1,1)*Fs+1:epch(1,2)*Fs);
    grf_walk=X(10,epch(2,1)*Fs+1:epch(2,2)*Fs);
    
    w_grf=grf_walk/mean(grf_stand(1+Fs*60:Fs*120));
    w_grf_fft=fft(w_grf)/length(w_grf);
    [w_grf_phi w_grf_rho]=cart2pol(real(w_grf_fft),imag(w_grf_fft));
    
    w_ecg=X(9,epch(2,1)*Fs+1:epch(2,2)*Fs);
    w_ecg_fft=fft(w_ecg)/length(w_ecg);
    [w_ecg_phi w_ecg_rho]=cart2pol(real(w_ecg_fft),imag(w_ecg_fft));
    
    w_eeg=X(1:8,epch(2,1)*Fs+1:epch(2,2)*Fs);
     
    cut_off=8.5;    
    for kk=1:8
        [wl_eeg(kk,:),fc, num, den]=fn_buterworth(w_eeg(kk,:),cut_off,cut_off+1,Fs);
        for ss=1:Fs
            x1(ss,kk)=corr(wl_eeg(kk,1+ss:end)',w_grf(1:end-ss)');
        end
%         w_eeg_fft(kk,:)=fft(w_eeg(kk,:))/length(w_eeg(kk,:));
%         [w_eeg_phi(kk,:) w_eeg_rho(kk,:)]=cart2pol(real(w_eeg_fft(kk,:)),imag(w_eeg_fft(kk,:)));        
    end
    
%     hz_axis=(1:length(w_grf_rho)/2)/length(w_grf_rho)*Fs;
%     figure(3003)
%     plot(hz_axis,w_grf_rho(2:length(w_grf_rho)/2+1),'k'); hold on;
%     plot(hz_axis,w_ecg_rho(2:length(w_ecg_rho)/2+1),'r'); hold on;
%     plot(hz_axis,w_eeg_rho(2,2:length(w_ecg_rho)/2+1),'b'); hold on;
%     plot(hz_axis,w_eeg_rho(3,2:length(w_ecg_rho)/2+1),'g'); 
    sec_axis=(1:length(x1))/Fs;
    clf;fig3004=figure(3004);
    plot(sec_axis,x1); ylabel('correlation coefficient'); xlabel('time (sec)')
    saveas(fig3004 , ['H:\Lab\figures\eeg\' qq '_' motion '.png'])  
    
    clearvars -except Fs motion grf_E99BW grf_E95BW grf_E90BW w_0 w_1 w_2 w_3 w0_pwr w1_pwr w0_rpwr w1_rpwr


end

