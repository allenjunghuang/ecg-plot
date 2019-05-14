%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  MAIN PROGRAM: COP wavelet analysis 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;clc;close all;

for pp=2
    p=num2str(pp);
    load(['H:\Lab\lab_data\figures\' p '.mat']); % input parameters: weight_cycle,copx_cycle,Tvi_copxN,Tvit_copxN,Tvi_weight,Tvit_weight,binary_step,InvFzs
    %%%%%%%%%%%%%%% 將兩腳step訊號 內插成相同點數 %%%%%%%%%%%%%
    step_norm=512;
    for kk=1:2:length(copx_cycle)-1
        step1=InvFzs(copx_cycle(kk):copx_cycle(kk+1));
        step2=InvFzs(copx_cycle(kk+1):copx_cycle(kk+2));
        xi=linspace(1,length(step1),step_norm);
        x=1:length(step1); 
        s1(:,1)=interp1(x,step1,xi);%stride 內插後
        clear xi x step1;
        xi=linspace(1,length(step2),step_norm);
        x=1:length(step2);
        s2(:,1)=interp1(x,step2,xi);%stride 內插後
        clear xi x step2;
        rr=(kk+1)/2;
        corre_step(rr)=corr(s1,s2);
        corre_step_all{pp,rr}=corr(s1,s2);
        
        wname = 'haar';
        [dwt_data1 dwt_length1]= fn_1Ddwt(s1,4,wname,'padding');%低頻至高頻 A3 D3 D2 D1...
        [dwt_data2 dwt_length2]= fn_1Ddwt(s2,4,wname,'padding');%低頻至高頻 A3 D3 D2 D1...
        % [C,L] = wavedec(F1_open(7169:23552,1),4,'db12');
        
%         figure(1)
%         subplot(5,1,1);plot(dwt_data1{1,1});
%         xlabel('index');ylabel('magnitude');legend('A4');
%         title(['Wavelet Analysis (' wname ')']);
%         subplot(5,1,2);plot(dwt_data1{1,2});
%         xlabel('index');ylabel('magnitude');legend('D4');
%         subplot(5,1,3);plot(dwt_data1{1,3});
%         xlabel('index');ylabel('magnitude');legend('D3');
%         subplot(5,1,4);plot(dwt_data1{1,4});
%         xlabel('index');ylabel('magnitude');legend('D2');
%         subplot(5,1,5);plot(dwt_data1{1,5});
%         xlabel('index');ylabel('magnitude');legend('D1');        
%         figure(2)
%         subplot(5,1,1);plot(dwt_data2{1,1});
%         xlabel('index');ylabel('magnitude');legend('A4');
%         title(['Wavelet Analysis (' wname ')']);
%         subplot(5,1,2);plot(dwt_data2{1,2});
%         xlabel('index');ylabel('magnitude');legend('D4');
%         subplot(5,1,3);plot(dwt_data2{1,3});
%         xlabel('index');ylabel('magnitude');legend('D3');
%         subplot(5,1,4);plot(dwt_data2{1,4});
%         xlabel('index');ylabel('magnitude');legend('D2');
%         subplot(5,1,5);plot(dwt_data2{1,5});
%         xlabel('index');ylabel('magnitude');legend('D1');

        A4_1=dwt_data1{1,1};%小波拆解完訊號 單腳
        D1_1=dwt_data1{1,5};
        D2_1=dwt_data1{1,4};
        D3_1=dwt_data1{1,3};
        D4_1=dwt_data1{1,2};
        
        Total_E1= sum(A4_1.^2)+sum(D1_1.^2)+sum(D2_1.^2)+sum(D3_1.^2)+sum(D4_1.^2);
        A4_ratioE1(rr)=sum(A4_1.^2)/Total_E1;%energy_ratio
        D1_ratioE1(rr)=sum(D1_1.^2)/Total_E1;
        D2_ratioE1(rr)=sum(D2_1.^2)/Total_E1;
        D3_ratioE1(rr)=sum(D3_1.^2)/Total_E1;
        D4_ratioE1(rr)=sum(D4_1.^2)/Total_E1;
        
        A4_2=dwt_data2{1,1};%小波拆解完訊號 單腳
        D1_2=dwt_data2{1,5};
        D2_2=dwt_data2{1,4};
        D3_2=dwt_data2{1,3};
        D4_2=dwt_data2{1,2};
        
        Total_E2=sum(A4_2.^2)+sum(D1_2.^2)+sum(D2_2.^2)+sum(D3_2.^2)+sum(D4_2.^2);
        A4_ratioE2(rr)=sum(A4_2.^2)/Total_E2;
        D1_ratioE2(rr)=sum(D1_2.^2)/Total_E2;
        D2_ratioE2(rr)=sum(D2_2.^2)/Total_E2;
        D3_ratioE2(rr)=sum(D3_2.^2)/Total_E2;
        D4_ratioE2(rr)=sum(D4_2.^2)/Total_E2;
        
        corr_A4(rr)=corr(A4_2',A4_1');
        corr_D1(rr)=corr(D1_2',D1_1');
        corr_D2(rr)=corr(D2_2',D2_1');
        corr_D3(rr)=corr(D3_2',D3_1');
        corr_D4(rr)=corr(D4_2',D4_1');
        %%%%% asymmetry index %%%%%
        SI_A4(rr)=sum(abs(A4_2-A4_1)./(0.5*abs(A4_2+A4_1))); 
        SI_D1(rr)=sum(abs(D1_2-D1_1)./(0.5*abs(D1_2+D1_1))); 
        SI_D2(rr)=sum(abs(D2_2-D2_1)./(0.5*abs(D2_2+D2_1)));
        SI_D3(rr)=sum(abs(D3_2-D3_1)./(0.5*abs(D3_2+D3_1)));
        SI_D4(rr)=sum(abs(D4_2-D4_1)./(0.5*abs(D4_2+D4_1)));

        clear A4_1 A4_2 D1_1 D1_2 D2_1 D2_2 D3_1 D3_2 D4_1 D4_2 dwt_data1 dwt_data2 s1 s2
    end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
%     t_d_bar= corr_D3; % 可依據需求更換
%     t_median=median(t_d_bar); %  尋找Correlation 中值
%     t_temp=find(t_d_bar>t_median); % 尋找大於中值的點(index)
%     t_temp_diff=[diff(t_temp) 0]; % index 間距
%     add=1;bb=0;
%     for cc=1:length(t_temp_diff)
%         if t_temp_diff(cc) ==1 %index 間距=1
%             add=add+1;
%         else %index 間距 > 1
%             bb=bb+1;
%             count_tt(bb)=add; %若間距 >1, 紀錄index間距=1連續的次數; 若無連續次數, 紀錄1
%             add=1;
%         end
%     end
%     
%     ALL_X(pp)=max(count_tt); %間距=1連續的最大時間
% %     the_ratio(pp)=sum(count_tt)/length(t_d_bar); %t_d_bar中超過中值的比例???
% %     count_tt=count_tt(find(count_tt(1:end-1)>0));
%     x=1:35;%max(count_tt);%%%
%     [n1,xout1]=hist(count_tt,x);
%     n1=(n1./sum(n1))*100; %正規化
% %     figure(102)
% %     bar(xout1,n1)
%     clear count_tt t_temp_diff t_temp t_d_bar t_median
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     t_d_bar= corr_D3;  % 可依據需求更換
%     t_median=median(t_d_bar);
%     t_temp=find(t_d_bar < t_median);
%     t_temp_diff=[diff(t_temp) 0];
%     add=1;bb=0;
%     for cc=1:length(t_temp_diff)
%         if t_temp_diff(cc) ==1
%             add=add+1;
%         else
%             bb=bb+1;
%             count_ss(bb)=add;
%             add=1;
%         end
%     end
%     
%     ALL_X_less(pp)=max(count_ss);
% %    the_ratio_less(pp)=sum(count_ss)/length(t_d_bar);
% %    count_tt=count_tt(find(count_tt(1:end-1)>0));
%     xs=1:35;% max(count_tt);
%     [n2,xout2]=hist(count_ss,xs);
%     n2=(n2./sum(n2))*100;
% %     figure(102)
% %     bar(xout2,n2)
%     
%     %%%%%%%%%% entropy %%%%%%%%%%
%     for rr=1:xs(end)
%         templog2_1(rr)=n1(rr)*(log2(n1(rr)))*-1;
%         templog2_2(rr)=n2(rr)*(log2(n2(rr)))*-1;
%     end
%     nanlog2_1=isnan(templog2_1);
%     nanlog2_2=isnan(templog2_2);
%     temp=find(nanlog2_1==0);
%     greater_log2=sum(templog2_1(temp));clear temp
%     temp=find(nanlog2_2==0);
%     less_log2=sum(templog2_2(temp));clear temp
%     
%     clear templog2_1 templog2_2 templog10_1 templog10_2
%     clear count_tt t_temp_diff t_temp t_d_bar xs count_ss t_median
%     %%%%%%%%%% entropy end %%%%%%%%%%
%  
    mean_A4_ratioE1(pp)=mean(A4_ratioE1);
    mean_D1_ratioE1(pp)=mean(D1_ratioE1);
    mean_D2_ratioE1(pp)=mean(D2_ratioE1);
    mean_D3_ratioE1(pp)=mean(D3_ratioE1);
    mean_D4_ratioE1(pp)=mean(D4_ratioE1);
    
    mean_A4_ratioE2(pp)=mean(A4_ratioE2);
    mean_D1_ratioE2(pp)=mean(D1_ratioE2);
    mean_D2_ratioE2(pp)=mean(D2_ratioE2);
    mean_D3_ratioE2(pp)=mean(D3_ratioE2);
    mean_D4_ratioE2(pp)=mean(D4_ratioE2);
%     
%     if pp<72.5
%         
%         mean_A4_corr_Co(pp,:)=mean(corr_A4);
%         mean_D1_corr_Co(pp,:)=mean(corr_D1);
%         mean_D2_corr_Co(pp,:)=mean(corr_D2);
%         mean_D3_corr_Co(pp,:)=mean(corr_D3);
%         mean_D4_corr_Co(pp,:)=mean(corr_D4);
%         
%         hist_P_co(pp,:)=n1;
%         hist_P_less_co(pp,:)=n2;
%         hist_P_diff_co(pp,:)=(((n1-n2)));
%         hist_P_diff_abssum_co(pp,:)=sum(abs((n1-n2)));
%         hist_P_corr_co(pp,:)=(corr(n1',n2'));
%         
%         greater_log2_co(pp,:)=greater_log2;
%         less_log2_co(pp,:)=less_log2;
%         log2_diff_co(pp,:)=(greater_log2-less_log2);
%         log2_mean_co(pp,:)=(greater_log2+less_log2)/2;
%         log2_add_co(pp,:)=(greater_log2+less_log2);
%         %         log2_min_co(pp,:)=max([greater_log2 less_log2]);
%         corre_step_co(pp,:)=mean(corre_step);
%         
%         SI_A4_co(pp,:)=mean(SI_A4);
%         SI_D1_co(pp,:)=mean(SI_D1);
%         SI_D2_co(pp,:)=mean(SI_D2);
%         SI_D3_co(pp,:)=mean(SI_D3);
%         SI_D4_co(pp,:)=mean(SI_D4);
%         
%         clear n1 n2 greater_log2  less_log2 greater_log10 less_log10
%         clear SI_A4 SI_D1 SI_D2 SI_D3 SI_D4
%     else
%         kk=pp-72;
%         mean_A4_corr_Pt(kk,:)=mean(corr_A4);
%         mean_D1_corr_Pt(kk,:)=mean(corr_D1);
%         mean_D2_corr_Pt(kk,:)=mean(corr_D2);
%         mean_D3_corr_Pt(kk,:)=mean(corr_D3);
%         mean_D4_corr_Pt(kk,:)=mean(corr_D4);
%         
%         hist_P_pt(kk,:)=n1;
%         hist_P_less_pt(kk,:)=n2;
%         hist_P_diff_pt(kk,:)=((n1-n2));
%         hist_P_diff_abssum_pt(kk,:)=sum(abs(n1-n2));
%         hist_P_corr_pt(kk,:)=(corr(n1',n2'));
%         greater_log2_pt(kk,:)=greater_log2;
%         less_log2_pt(kk,:)=less_log2;
%         log2_diff_pt(kk,:)=(greater_log2-less_log2);
%         log2_mean_pt(kk,:)=(greater_log2+less_log2)/2;
%         log2_add_pt(kk,:)=(greater_log2+less_log2);
%         %         log2_min_pt(kk,:)=max([greater_log2 less_log2]);
%         corre_step_pt(kk,:)=mean(corre_step);
%         
%         SI_A4_pt(kk,:)=mean(SI_A4);
%         SI_D1_pt(kk,:)=mean(SI_D1);
%         SI_D2_pt(kk,:)=mean(SI_D2);
%         SI_D3_pt(kk,:)=mean(SI_D3);
%         SI_D4_pt(kk,:)=mean(SI_D4);
%         
%         clear n1 n2 greater_log2  less_log2 greater_log10 less_log10
%         clear SI_A4 SI_D1 SI_D2 SI_D3 SI_D4
%     end
%     clear A4_ratioE2 D1_ratioE2 D2_ratioE2 D3_ratioE2   D4_ratioE2 GRF_norm
%     clear A4_ratioE1 D1_ratioE1 D2_ratioE1 D3_ratioE1   D4_ratioE1
%     clear  corr_A4 corr_D1  corr_D2  corr_D3    corr_D4 data GRF_norm Stand_pahse_index2 Stand_pahse_index
%     clear nanlog10_1 nanlog10_2 nanlog2_1 nanlog2_2 stride stand_phase ZeroCrossing_index2 ZeroCrossing_index1
%     clear F1_open F2_open corre_step
end
% % final_Mean_Co=[mean_A4_corr_Co mean_D4_corr_Co mean_D3_corr_Co mean_D2_corr_Co mean_D1_corr_Co hist_P_co(:,1:4) ((sum(hist_P_co(:,5:end)'))') hist_P_less_co(:,1:4) ((sum(hist_P_less_co(:,5:end)'))') hist_P_diff_co(:,1:4) ((sum(hist_P_diff_co(:,5:end)'))') greater_log2_co less_log2_co log2_add_co log2_mean_co corre_step_co hist_P_diff_abssum_co SI_A4_co SI_D4_co SI_D3_co SI_D2_co SI_D1_co];%低頻至高頻 A3 D3 D2 D1...
% % final_Mean_Pt=[mean_A4_corr_Pt mean_D4_corr_Pt mean_D3_corr_Pt mean_D2_corr_Pt mean_D1_corr_Pt hist_P_pt(:,1:4) ((sum(hist_P_pt(:,5:end)'))') hist_P_less_pt(:,1:4) ((sum(hist_P_less_pt(:,5:end)'))') hist_P_diff_pt(:,1:4) ((sum(hist_P_diff_pt(:,5:end)'))') greater_log2_pt less_log2_pt log2_add_pt log2_mean_pt corre_step_pt hist_P_diff_abssum_pt SI_A4_pt SI_D4_pt SI_D3_pt SI_D2_pt SI_D1_pt];
% % final_Mean_Pt(23,:)=[];final_Mean_Pt(4,:)=[];final_Mean_Pt(3,:)=[];final_Mean_Pt(44,:)=[];final_Mean_Co(60,:)=[];
% % PTT_1=final_Mean_Co;
% % PTT_2=final_Mean_Pt;
% % 
% % for kk=1:size(PTT_2')
% %     std_PT1=std(PTT_1(:,kk));
% %     std_PT2=std(PTT_2(:,kk));
% %     mean_PT1=mean(PTT_1(:,kk));
% %     mean_PT2=mean(PTT_2(:,kk));
% %     temp3=std(PTT_1(:,kk))^2/length(PTT_1(:,kk));
% %     temp4=std(PTT_2(:,kk))^2/length(PTT_2(:,kk));
% %     t_star(kk)=(mean_PT1-mean_PT2)/sqrt(temp3+temp4);
% % end

% Signal=dwt_data1{1,1};
% SF=new_SF;
% disPeak=diff(Signal)./SF;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gg=find(disPeak>1.4);
% disPeak(gg)=[];
%%%
% dispeak_mean=mean(disPeak)
% n_fs=1/mean(disPeak);          %Peak間距頻率
% disPeak_1=disPeak-mean(disPeak);

%     n=65536;%FFT點數
%     p=1;
%     temp_ppg=Signal-mean(Signal);
%     temp_power(:,p)=sum(temp_ppg.^2)*(1/SF);      %% 計算每週期時域能量
%     ppg_fft1(:,p)=fft(temp_ppg,n);
%     [x2,y2]=cart2pol(real(ppg_fft1(:,p)),imag(ppg_fft1(:,p)));
%     ppg_fft_x2(:,p)=x2;
%     ppg_fft_y2(:,p)=y2;
%     figure(1)
%     plot((1:n/2)/n*SF,y2(2:n/2+1,p),'color',[0 0.5 0])
%     xlabel('Hz')
% ylabel('Amplitude')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%頻域能量
%     ppg_fft_y2_pow(:,p)=(ppg_fft_y2(:,p).^2);
%     ppg_fft_y2_pow(1,p)=ppg_fft_y2_pow(1,p)/2;
%     P_t(:,p)=sum(ppg_fft_y2_pow(:,p))/n/SF; 
%     %     P_t(:,p)=(ppg_fft_y2_pow(129,p)+(sum(ppg_fft_y2_pow(1:128,p))*2))/SF;
%     Pyy(:,p)=(ppg_fft_y2_pow(:,p))*2/n/SF;
%     Pyy(n/2+1,p)=ppg_fft_y2(n/2+1,p)/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%計算95%實質頻寬項次%%%%%%%%%%%%%
%      for j=1:length(P_t)
%           p_95(:,j)=P_t(1,j)*95/100;
%           temp_p95=0;
%           c=0;
%          for q=1:n/2+1
%              if Pyy(q,j)+temp_p95 <= p_95(:,j)
%                  temp(q,j)=Pyy(q,j);
%                  a1=sum(temp(:,j));
%                  temp_p95=a1;
%                  c=c+1;
% %                  break;
%              else
%                  temp(q,j)=0;
%                  temp_p95=a1+100;
%                  total(:,j)=c-1;
%              end
%          end
%      end
%      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 計算頻寬
%      sample=SF;    % 取樣頻率
% t0=n/sample;        %資料筆數/取樣頻率
% f_fre=2*pi/t0;      %基本頻率
% for w=1:length(total)
%     bw(:,w)=total(:,w)*f_fre;
% end
% fre=bw./(2*pi)                  %%%%% 頻寬