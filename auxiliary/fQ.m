function B = fQ(A,P)
% -------------------------------------------------------------------------
%
% IQClab:      Version 3.4.0
% Copyright:   This is copyrighted material owned by Novantec B.V.
% Terms:       IQClab is available for non-commercial usage under a
%              Creative Commons (Attribution-NoDerivatives 4.0
%              International (CC BY-ND 4.0)) license:  
%              https://creativecommons.org/licenses/by-nd/4.0/
%              For further information please visit iqclab.eu
% Author:      J.Veenman
% Date:        26-04-2012
% 
% -------------------------------------------------------------------------
%
% Description: This function generates, depending on the number of inputs,
%              the matrix:
%
%                1.) B = A'*A,   or   2.) B = A'*P*A
%
% Syntax:      B = fQ(A)
%              B = fQ(A,P)
%
% -------------------------------------------------------------------------
if nargin == 1
    B = A'*A;
elseif nargin == 2
    B = A'*P*A;
end