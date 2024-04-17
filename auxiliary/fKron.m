function C = fKron(A,B)
% -------------------------------------------------------------------------
%
% IQClab:      Version 3.04
% Copyright:   This is copyrighted material owned by Novantec B.V.
% Terms:       IQClab is available for non-commercial usage under a
%              Creative Commons (Attribution-NoDerivatives 4.0
%              International (CC BY-ND 4.0)) license:  
%              https://creativecommons.org/licenses/by-nd/4.0/
%              For further information please visit iqclab.eu
% Author:      J.Veenman
% Date:        19-03-2020
% 
% -------------------------------------------------------------------------
%
% Description: This function computes the Kronecker product C = kron(A,B).
%              The dfference with the Matlab command "kron" is that this
%              function is compatible with uncertainty objects.
%
% Syntax:      C = fKron(A,B)
%
% Usage:       The function C = fKron(A,B) returns the Kronecker product
%
%                               [a11*B, ... ,a1n*B]
%               C = kron(A,B) = [ ... , ... , ... ]
%                               [am1*B, ... ,amn*B]
%
%              where aij, i\in[1,n], j\in[1,m] are the elements of the
%              matrix A.
%
% -------------------------------------------------------------------------
[k,l] = size(A);
[m,n] = size(B);

for i = 1:k
   for j = 1:l
       C(m*(i-1)+1:i*m,n*(j-1)+1:j*n)=A(i,j)*B;
   end
end
end