clear all; clc; close all;
path(path,'H:\Lab\m_files\eeglab13_3_2b\functions\sigprocfunc\'); 
path(path,'H:\Lab\m_files\m_file\functions\'); 

Fs=250;
motion='b';
options = { 'TYPE' 'EDF' };

for pp=1:19
    qq=num2str(pp);   
    [X, header] = edf2mat(['H:\Lab\data\2014_data\eeg_ambulating\Cz\edf\edf\' qq '' motion '.edf']);             
    epch=[174-1 414+1];
    x2=X(1:10,epch(1)*Fs+1:epch(2)*Fs);
    fact=std(x2,0,2);
    grf=(x2(10,:)-mean(x2(10,:)))*mean(fact(1:8))/fact(10);
    [E99BW_grf]=fn_bandwidth(grf, Fs, 0.99)
    [E95BW_grf]=fn_bandwidth(grf, Fs, 0.95)
    
    %LPF+LMS
    cut_off=8.5;
    fres=100;
    Rn=resample(grf,fres,Fs);
    Sn=zeros(8,length(Rn));
    for hh=1:8
        xn=x2(hh,:);
        [xn,fc, num, den] = fn_buterworth(xn,cut_off,cut_off+1,Fs);
        Sn(hh,:)=resample(xn,fres,Fs);
    end
    
    % Adaption filtering
    % --------------------------------------------------------
    % delay_0=3;
    % wv_0=(cn_0-1)/delay_0;
    cn_0=3;
    mu_0=0.0003;
    
    recnum=1;
    for kk=1:length(cn_0)
        for ss=1:length(mu_0)
            %Decimation LMS
            %         tsec=(1:length(grf))/Fs;
            %         mu_n(kk,ss)=mu_0(ss)/(cn_0(kk)*var(grf));
            %         [E,Y,Wn] = adaptfilt_decilms(x2(1:8,:), mu_n(kk,ss), cn_0(kk), grf, delay_0);
            %         Wn1=mean(Wn(:,:,fix(length(grf)/1.5:end)),3);
            %         [E2,Y2,Wn2] = adaptfilt_decilms2(x2(1:8,:), mu_n(kk,ss), cn_0(kk), grf, Wn1, 1, delay_0);
            
            %         Z=[E2(:,5*Fs+1:end-5*Fs);x2(9,5*Fs+1:end-5*Fs)*mean(fact(1:8))/fact(9);];
            %         EEG = struct('data',Z,'srate',Fs,'Label',{{'O2';'P4';'C4';'F4';'O1';'P3';'C3';'F3';'ECG';'GRF'}});
            %         filename=['decilms_' num2str(pp) motion '_' num2str(cn_0(kk)) 'th_' num2str(delay_0) 'td.edf'];
            %         writeeeg(filename, EEG.data(:,:), EEG.srate,'Label',EEG.Label, options{:});
            
            %LPF+LMS
            tsec=(1:size(Rn,2))/fres;
            mu_n(kk,ss)=mu_0(ss)/(cn_0(kk)*var(Rn))/3;
            [E,Y,Wn] = adaptfilt_lms(Sn, mu_n(kk,ss), cn_0(kk), Rn);
            Wh1=mean(Wn(:,:,fix(size(Rn,2)/1.5:end)),3);
            [E2,Y2,Wn2] = adaptfilt_lms2(Sn, mu_n(kk,ss), cn_0(kk), Rn, Wh1, recnum);
            
            if length(cn_0)==1 && length(mu_0)==1
                d0=resample(Y2',Fs,fres);
                d1=x2(1:8,:)-d0';
                Z=[d1(:,5*Fs+1:end-5*Fs);x2(9,5*Fs+1:end-5*Fs)*mean(fact(1:8))/fact(9);];
                EEG = struct('data',Z,'srate',Fs,'Label',{{'O2';'P4';'C4';'F4';'O1';'P3';'C3';'F3';'ECG';'GRF'}});
                filename=['lms_' num2str(pp) motion '_' num2str(cn_0) 'th_fs' num2str(fres) '_fc' num2str(mu_0) '.edf'];
                writeeeg(filename, EEG.data(:,:), EEG.srate,'Label',EEG.Label, options{:});
            else            
                W{kk,ss}=Wn;
                W2{kk,ss}=Wn2;
                clear E Y Wn E2 Y2 Wn2 Wh1
            end
        end
    end
    
end

% Check filter parameters
% --------------------------------------------------------
% isp_ch=3;
% if length(mu_0)==1 && length(cn_0)~=1%Find filter order
%     for kk=1:cn_0(1) %weight vector
%         for pp=1:length(cn_0) %filter order
%             wh_1(kk,:,pp)=W{pp,1}(kk,isp_ch,:);          
%         end
%         figure(kk+100)
%         plot(tsec, wh_1(kk,:,1),'r'); hold on;
%         plot(tsec, wh_1(kk,:,2),'g'); hold on;
%         plot(tsec, wh_1(kk,:,3),'b'); hold on;
%         plot(tsec, wh_1(kk,:,4),'k'); hold on;
%         plot(tsec, wh_1(kk,:,5),'y'); 
%         axis([-inf inf -inf inf]); xlabel('Time (sec)'); ylabel(['w_' num2str(kk-1)]);
%         h1=legend(['M=' num2str(cn_0(1))],['M=' num2str(cn_0(2))],['M=' num2str(cn_0(3))],['M=' num2str(cn_0(4))],['M=' num2str(cn_0(5))]); set(h1,'box','off');
%         title(['\mu_0=' num2str(mu_0(1))])
%     end
%     
% elseif length(cn_0)==1 && length(mu_0)~=1%Find learning rate
%     for pp=1:cn_0
%         for kk=1:6
%             mu_1(pp,:,kk)=W{1,kk}(pp,isp_ch,:);   
%             mu_2(pp,:,kk)=W2{1,kk}(pp,isp_ch,:);
%         end        
%         fig5xx=figure(pp+500);
%         subplot(3,2,1); 
%         plot(tsec, mu_1(pp,:,1),'r'); hold on; plot(tsec, mu_2(pp,:,1),'g'); 
%         h1=legend('1st loop',[num2str(recnum+1) 'nd loop']); set(h1,'box','off','fontsize',10); title(['\mu_0=' num2str(mu_0(1)) ' M=' num2str(cn_0)],'fontsize',11)
%         axis([1 inf -inf inf]); xlabel('time (sec)','fontsize',10); ylabel(['w_' num2str(pp-1)],'fontsize',10);
%         subplot(3,2,2); 
%         plot(tsec, mu_1(pp,:,2),'r'); hold on; plot(tsec, mu_2(pp,:,2),'g');
%         h1=legend('1st loop',[num2str(recnum+1) 'nd loop']); set(h1,'box','off','fontsize',10); title(['\mu_0=' num2str(mu_0(2)) ' M=' num2str(cn_0)],'fontsize',11)
%         axis([1 inf -inf inf]); xlabel('time (sec)','fontsize',10); ylabel(['w_' num2str(pp-1)],'fontsize',10);
%         subplot(3,2,3); 
%         plot(tsec, mu_1(pp,:,3),'r'); hold on; plot(tsec, mu_2(pp,:,3),'g'); 
%         h1=legend('1st loop',[num2str(recnum+1) 'nd loop']); set(h1,'box','off','fontsize',10); title(['\mu_0=' num2str(mu_0(3)) ' M=' num2str(cn_0)],'fontsize',11)
%         axis([1 inf -inf inf]); xlabel('time (sec)','fontsize',10); ylabel(['w_' num2str(pp-1)],'fontsize',10);
%         subplot(3,2,4); 
%         plot(tsec, mu_1(pp,:,4),'r'); hold on; plot(tsec, mu_2(pp,:,4),'g'); 
%         h1=legend('1st loop',[num2str(recnum+1) 'nd loop']); set(h1,'box','off','fontsize',10); title(['\mu_0=' num2str(mu_0(4)) ' M=' num2str(cn_0)],'fontsize',11)
%         axis([1 inf -inf inf]); xlabel('time (sec)','fontsize',10); ylabel(['w_' num2str(pp-1)],'fontsize',10);
%         subplot(3,2,5); 
%         plot(tsec, mu_1(pp,:,5),'r'); hold on; plot(tsec, mu_2(pp,:,5),'g'); 
%         h1=legend('1st loop',[num2str(recnum+1) 'nd loop']); set(h1,'box','off','fontsize',10); title(['\mu_0=' num2str(mu_0(5)) ' M=' num2str(cn_0)],'fontsize',11)
%         axis([1 inf -inf inf]); xlabel('time (sec)','fontsize',10); ylabel(['w_' num2str(pp-1)],'fontsize',10);
%         subplot(3,2,6); 
%         plot(tsec, mu_1(pp,:,6),'r'); hold on; plot(tsec, mu_2(pp,:,6),'g'); 
%         h1=legend('1st loop',[num2str(recnum+1) 'nd loop']); set(h1,'box','off','fontsize',10); title(['\mu_0=' num2str(mu_0(6)) ' M=' num2str(cn_0)],'fontsize',11)
%         axis([1 inf -inf inf]); xlabel('time (sec)','fontsize',10); ylabel(['w_' num2str(pp-1)],'fontsize',10);
%         set(fig5xx, 'Units', 'centimeters','Position', [2 2 20 16]); 
%         clear fig5xx
%     end          
% 
% elseif length(cn_0)==1 && length(mu_0)==1
%     for pp=1:size(Wn,1)
%         wt1(pp,:)=Wn(pp,isp_ch,cn_0:end);
%         wt2(pp,:)=Wn2(pp,isp_ch,cn_0:end);
%     end
%     fig700=figure(700);
%     subplot(3,1,1)
%     plot((1:size(wt1,2))/Fs, wt1(1,:),'r'); hold on; plot((1:size(wt2,2))/Fs, wt2(1,:),'g');
%     h1=legend('1st loop',[num2str(recnum+1) 'nd loop']); set(h1,'box','off','fontsize',10);
%     title(['\mu_0=' num2str(mu_0) ' M=' num2str(cn_0)],'fontsize',11)
%     xlabel('time (sec)','fontsize',10); ylabel('w_0','fontsize',10);
%     subplot(3,1,2)
%     plot((1:size(wt1,2))/Fs, wt1(2,:),'r'); hold on; plot((1:size(wt2,2))/Fs, wt2(2,:),'g');
%     h1=legend('1st loop',[num2str(recnum+1) 'nd loop']); set(h1,'box','off','fontsize',10);
%     xlabel('time (sec)','fontsize',10); ylabel('w_1','fontsize',10);
%     title(['\mu_0=' num2str(mu_0) ' M=' num2str(cn_0)],'fontsize',11)
%     subplot(3,1,3)
%     plot((1:size(wt1,2))/Fs, wt1(3,:),'r'); hold on; plot((1:size(wt2,2))/Fs, wt2(3,:),'g'); 
%     h1=legend('1st loop',[num2str(recnum+1) 'nd loop']); set(h1,'box','off','fontsize',10);
%     xlabel('time (sec)','fontsize',10); ylabel('w_2','fontsize',10);
%     title(['\mu_0=' num2str(mu_0) ' M=' num2str(cn_0)],'fontsize',11)
%     set(fig700, 'Units', 'centimeters','Position', [2 2 10 16]);
% end