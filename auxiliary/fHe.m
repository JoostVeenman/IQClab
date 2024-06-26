function B = fHe(A)
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
% Date:        07-04-2011
% 
% -------------------------------------------------------------------------
%
% Description: This function generates the matrix B = A' + A
%
% Syntax:      B = fHe(A)
%
% Usage:       For square matrices A, this function generates the symmetric
%              or Hermitian matrix B = A' + A depending on whether A is
%              real or complex respectively.
%
% -------------------------------------------------------------------------
[noa,nia] = size(A);
if noa ~= nia
    error('A should be a square matrix');
end   
B = A + A';
end
