function m = fss2m(sys)
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
% Date:        08-12-2019
% 
% -------------------------------------------------------------------------
%
% Description: This function creates a KYP outer-factor from a state space.
%
% Syntax:      m = fss2m(sys)
%
% Usage:       m = fss2m(sys) creates the matrix
%
%                  [I    ,0    ]
%              m = [sys.a,sys.b]
%                  [sys.c,sys.d]
% -------------------------------------------------------------------------

sA1 = size(sys.a,1);
sB2 = size(sys.b,2);
m   = [eye(sA1),zeros(sA1,sB2);sys.a,sys.b;sys.c,sys.d];
end

