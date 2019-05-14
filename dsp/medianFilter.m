%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: 
% Remove baseline-drifting by median filter. The WindowLength depends
% on the COP features, hop-size is one point.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output_seq,driftline]=Fn_COPmedian_filter(input_seq,Fs,signal_time,WindowLength)
if (size(input_seq,1) > size(input_seq,2)) %correct data dimension to row vector
    input_seq = input_seq';
end    
%%%%%%%%%%%%%%%%%%%%%%% median filter 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WindowNum=(signal_time-WindowLength)*Fs+1;
% Signal_temp=zeros(length(input_seq),1);
driftline=zeros(1,fix(WindowNum));

for kk=1:WindowNum
    Signal_Window=input_seq(kk:fix(kk+WindowLength*Fs-1));
    driftline(kk)=median(Signal_Window);
%     Signal_MeanWindow=Signal_Window-median(Signal_Window); % zero mean
%     Signal_temp(temp:temp+length(Signal_Window)-1)=Signal_MeanWindow;
end
t=fix(WindowLength*Fs/2);
output_seq=input_seq(t+1:t+length(driftline))-driftline;
%%%%%%%%%%%%%%%%%%%%%%%% median filter 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% WindowNum2=(fix(length(driftline)/Fs)-WindowLength2)*Fs;
% driftline2=zeros(1,fix(WindowNum2));
% 
% for kk=1:WindowNum2
%     Signal_Window=driftline(kk:fix(kk+WindowLength2*Fs-1));  
%     driftline2(kk)=median(Signal_Window);
% %     Signal_MeanWindow2=Signal_Window-median(Signal_Window); % zero mean
% %     signal_output(temp:temp+length(Signal_Window)-1)=Signal_MeanWindow2;
% end
% t=(WindowLength*Fs+WindowLength2*Fs)/2;
% output_seq=input_seq(t+1:t+length(driftline2))-driftline2;