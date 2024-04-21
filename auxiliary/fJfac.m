function [M,J] = fJfac(P,ind)
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
% Description: This function factorizes the symmetric matrix P = [Q,S;S',R]
%              where Q > 0 and R - S'*Q^-1*S < 0 as P = M'*J*M, where
%
%                P =M'*J*M = [*]^T[I_ind,0 ][M11,M12]
%                            [*]  [0    ,-I][0  ,M22]
%
% Syntax:      [M,J] = fJfac(P,ind)
%
% Usage:       As input one should provide:
%                - P = P' is a symmetric matrix with the structure
%
%                  P = [Q,S;S',R],
%
%                  where: Q > 0 and R -S'*Q^-1*S < 
%
%                - ind = size of Q
%
%              As output one obtains:
%                - M   = upper triangular matrix 
%                - J   = [I_ind,0;0,-I]
%
% -------------------------------------------------------------------------
if size(P,1) > size(P,2) || size(P,1) < size(P,2)
    error('The matrix A should be square');
elseif norm(P-P',inf)>1e-16
    error('The matrix A should be symmetric');   
end
P          = 0.5*(P+P');
n          = size(P,1);
com        = setdiff(1:n,ind);
Q          = P(ind,ind);
S          = P(ind,com);
R          = P(com,com);
M11        = chol(-Q);
M12        = -eye(size(M11,1))/(M11')*S;
M22        = chol(R-S'/Q*S);

M(ind,ind) = M11;
M(ind,com) = M12;
M(com,com) = M22;
J          = zeros(1,n);
J(com)     = 1;
J(ind)     = -1;
J          = diag(J);
end