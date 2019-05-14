% Function:
%   Artifact removal using LMS (least mean squares) adaptive regression.
% 
% Inputs:
%   D   - Input data matrix (ddim x N)
%   U   - Reference noise signal (udim x N) 
%   M   - Order of the adaptive filter 
%   mu  - Learning rate 
% 
% Outputs:
%   E   - Output data matrix (artifact corrected)
%   Y   - Adaptive filter output
%   Wh  - filter weights evolution (M*udim x ddim x N)
% 
% References:
% [1] S. Haykin. Adaptive Filter Theory, (2002), Prentice Hall
% [2] G. Gomez-Herrero. Automatic Artifact Removal (AAR) toolbox (2007)

function [E,Y,Wh] = dspfilt_lms(D, mu, M, U)

OVERFLOW = 1E12; 
[ddim,dlen] = size(D);
[udim,ulen] = size(U);
if dlen~=ulen 
    error('dspfilt_lms: input signal and reference noise must have the same length'); 
end

% Initialization of the adaptation loop
% ---------------------------------------------
W = zeros(udim*M,ddim); 
E = zeros(ddim,dlen);
Wh = zeros(udim*M,ddim,dlen);
Wh(:,:,1:M-1) = repmat(W,[1,1,M-1]);
Y = zeros(ddim,dlen);


% Adaptation loop1
% ---------------------------------------------
fprintf('\nadaptive_lms: 1st loop');
for n = M:dlen
    r = U(:,n:-1:(n-M+1)); % adaptive filter input vector
    r = reshape(r', M*udim,1); 
    yn = r'*W; % adaptive filter output vector
    en = D(:,n)'-yn; % recovered signals   
    if ~isempty(find(abs(yn(:))>OVERFLOW, 1)),
        error('dspfilt_lms: LMS algorithm became unstable');
    end
    W = W+mu*r*en; 
    E(:,n) = en;   
    Wh(:,:,n) = W;
    Y(:,n) = yn;
   
    if ~mod(n,floor(dlen/10)), fprintf('.'); end
end
fprintf('[done]\n');



