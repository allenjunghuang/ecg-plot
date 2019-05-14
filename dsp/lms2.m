% Function:
%   Artifact removal using LMS (least mean squares) adaptive regression.
% Inputs:
%   D           - Input data matrix (ddim x N)
%   U           - Reference noise signal (udim x N) 
%   M           - Order of the adaptive filter 
%   mu          - Learning rate 
% Outputs:
%   Y   - Output data matrix (artifact corrected)
%   W   - Final filter weights (M*udim x ddim)
%   Wh  - filter weights evolution (M*udim x ddim x N)
% References:
% [1] S. Haykin. Adaptive Filter Theory, (2002), Prentice Hall
% [2] G. Gomez-Herrero. Automatic Artifact Removal (AAR) toolbox (2007)

function [E2,Y2,Wh2] = dspfilt_lms2(D, mu, M, U, W, num)

OVERFLOW = 1E12; 
[ddim,dlen] = size(D);
[udim,ulen] = size(U);
if dlen~=ulen 
    error('adaptive_lms2: input signal and reference noise must have the same length'); 
end

% initialization of the adaptation loop
% ---------------------------------------------
E2 = zeros(ddim,dlen);
Wh2 = zeros(udim*M,ddim,dlen);
Wh2(:,:,1:M-1) = repmat(W,[1,1,M-1]);
Y2 = zeros(ddim,dlen);


% adaptation loop
% ---------------------------------------------
fprintf('\ndspfilt_lms2: repetitive loop');
for kk=1:num
     fprintf('>');
    for n = M:dlen
        r = U(:,n:-1:(n-M+1)); % adaptive filter input vector
        r = reshape(r', M*udim,1);
        yn = r'*W; % adaptive filter output vector
        en = D(:,n)'-yn; % recovered signals
        if ~isempty(find(abs(yn(:))>OVERFLOW, 1)),
            error('dspfilt_lms2: LMS algorithm became unstable');
        end
        W = W+mu*r*en;
        if kk==num
            E2(:,n) = en;
            Wh2(:,:,n) = W;
            Y2(:,n) = yn;
        end
        if ~mod(n,floor(dlen/10)), fprintf('.'); end
    end
end
fprintf('[done]\n');

return


