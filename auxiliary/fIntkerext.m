function [a,b,c,d] = fIntkerext(A,B)
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
% Date:        20-04-2011
% 
% -------------------------------------------------------------------------
%
% Description: This function computes orthonormal bases of 
% 
%                - c         = ker(A) \cap ker(B)
%                - [a c]     = ker(A)
%                - [b c]     = ker(B)
%                - [a b c d] = full space
%
% Syntax:      [a,b,c,d] = fIntkerext(A,B)
%
% -------------------------------------------------------------------------
u = orth(A');
v = orth(B');
c = null([u';v']);
a = null([u';c']);
b = null([v';c']);
d = null([a';b';c']);
end