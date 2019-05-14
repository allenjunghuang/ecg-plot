% Function: 1-D wavelet decomposition 
% Inputs:
%   xn: ��J�T��
%   level_num: �p�i���Ѧ���
%   mother_wave: �p�i�����
% Reference: 
function [wlet, wlen]= dspwlet_1Ddwt(xn,level_num, mother_wave, style)

if (size(xn,1) > size(xn,2)) %correct data dimension to row vector
    xn = xn';
end       

ss=1;
if strcmp(style, 'truncate')==1 % �N��J�T���h�Y�h�����̤j��"2������"����
    while length(xn) > 2^ss
        wlen = 2^(ss);
        ss=ss+1;
    end
    data_range = fix((length(xn) - wlen)/2);
    xn = xn(:,data_range+1:data_range+wlen);
elseif strcmp(style, 'padding')==1 % �N��J�T���[�Y�[�����̤j��"2������"����
    while length(xn) > 2^ss
        wlen = 2^(ss+1);
        ss=ss+1;
    end  
    data_range = fix((wlen-length(xn))/2);
    xn = [zeros(1,data_range) xn zeros(1,(wlen-length(xn)-data_range))];
else
    error('Please select truncte or padding mode');
end

% �p�i����
for kk = 1 : level_num
    [cA,cD] = dwt(xn,mother_wave); 
    clear xn
    matrix{kk} = [cD];
    if (kk == level_num)
        matrix{kk+1} = [cA];
    end
    xn = cA;
    clear cA cD
end
clear xn

% �N�Ҧ��W�a���Y�Ƹɨ쵥��
% �U�W�a�T�����׭���(�Y�Ƹɪ��ץ�)
% aa = [];
% for kk = 1 : level_num-1
%     aa = [aa 2^kk];
% end
% freq_matrix{1} = matrix{1};
% for jj = 2 : level_num
%     for ii = 1 : length(matrix{jj})
%         freq_matrix{jj}((ii*aa(jj-1))-aa(jj-1)+1:(ii*aa(jj-1))) = matrix{jj}(ii);
%     end
%     clear ii
%     freq_matrix{jj}(length(matrix{1})+1:end) = [];
%     if (jj == level_num)
%         for ii = 1 : length(matrix{jj})
%             freq_matrix{jj+1}((ii*aa(jj-1))-aa(jj-1)+1:(ii*aa(jj-1))) = matrix{jj+1}(ii);
%         end
%         clear ii
%         freq_matrix{jj+1}(length(matrix{1})+1:end) = [];
%     end
% end
% clear jj matrix aa

% �N�U�W�a�ѧC�W�ƨ찪�W
wlet{1} = matrix{level_num+1};
for ii = 1 : level_num
    wlet{ii+1} = matrix{level_num+1-ii};
end




