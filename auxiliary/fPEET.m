function W=fPEET(type,Delta_t,Delta_ts)
% -------------------------------------------------------------------------
%
% IQClab:      Version 3.03
% Copyright:   This is copyrighted material owned by Novantec B.V.
% Terms:       IQClab is available for non-commercial usage under a
%              Creative Commons (Attribution-NonCommercial-NoDerivatives
%              4.0 International (CC BY-NC-ND 4.0))license: 
%              https://creativecommons.org/licenses/by-nc-nd/4.0/
%              Commercial usage is only permitted with a commercial
%              license. For further information please visit iqclab.eu
% Author:      J.Veenman
% Date:        15-01-2021
% 
% -------------------------------------------------------------------------
%
% Description: This function generates the pointing error performance
%              analysis weighting functions W_APE, W_MPE, W_RPE, W_PDE.
%
% Syntax:      W = fPEET(type,Delta_t,Delta_ts)
%
% Usage:       Provide the inputs:
%              - "type": There are four options:
%                 1.) 'ape' (Absolute performance error)
%                 2.) 'mpe' (Mean performance error)
%                 3.) 'rpe' (Relative performance error)
%                 4.) 'pde' (Performance drift error)
%              - "Delta_t":  Time window for type = 'mpe', 'rpe', 'pde'
%              - "Delta_ts": Time window for type = 'pde'
%
% -------------------------------------------------------------------------
if nargin == 1
    if ~strcmp(type,'ape')
       error('Not enough input arguments.'); 
    end
elseif nargin == 2
    if strcmp(type,'pde')
       error('Not enough input arguments.'); 
    end
end

s = tf('s');
switch type
    case 'ape'
        W  = 1;
    case 'mpe'
        W  = 2*(s*Delta_t+6)/(s^2*Delta_t^2+6*s*Delta_t+12);
    case 'rpe'
        W  = s*Delta_t*(s*Delta_t+sqrt(12))/(s^2*Delta_t^2+6*s*Delta_t+12);
    case 'pde'
        W1 = 2*(s*Delta_t+6 )/(s^2*Delta_t^2 +6*s*Delta_t +12);
        W2 = 2*(s*Delta_ts+6)/(s^2*Delta_ts^2+6*s*Delta_ts+12);
        W  = W1*W2;
end
W = balreal(ss(W));
end


