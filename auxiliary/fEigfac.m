function [M,J,v,in] = fEigfac(P)
% -------------------------------------------------------------------------
%
% Copyright:   This is copyrighted material owned by Novantec B.V.
% Terms:       IQClab is available under a Creative Commons
%              (Attribution-NoDerivatives 4.0 International (CC BY-ND 4.0))
%              license: https://creativecommons.org/licenses/by-nd/4.0/
%              For further information please visit iqclab.eu
%
% Author:      J.Veenman
% Date:        20-04-2011
% 
% -------------------------------------------------------------------------
%
% Description: Given the symmetric matrix P, this function factorizes the
%              symmetric matrix P as P = M'*J*M = M'*diag(I,0,-I)*M
%
% Syntax:      [M,J,v,in] = fEigfac(P)
%
% Usage:       The function [M,J,v,in] = fEigfac(P) factorizes the
%              symmetric matrix P as:
%
%                              [I_npos, 0      , 0      ]
%              P = M'*J*M = M'*[0     , 0_nzero, 0      ]*M
%                              [0     , 0      , I_npos ]
%
%              by using the command "eig".
%
%              As input one sould provide:
%
%                - the symmetric matrix P.
% 
%              As output one obtains:
%
%                - the outer factor M
%                - the middle matrix J=diag(I,0,-I)
%                - the eigenvectors of v
%                - the number of positive (npos), negative (nneg), and zero
%                  (nzero) eigenvalues in = [npos;nneg;nzero]
%
% -------------------------------------------------------------------------

if size(P,1) > size(P,2) || size(P,1) < size(P,2)
    error('The matrix A should be square');
elseif norm(P-P',inf)>1e-16
    error('The matrix A should be symmetric');
end
Ps       = 0.5*(P + P');
[v,d]    = eig(Ps);
J        = (d>0)-(d<0);
in       = [sum(diag(d>0));sum(diag(d<0));size(J,1)-sum(diag(d>0))-sum(diag(d<0))];
T        = blkdiag(eye(in(1)),[zeros(in(2),in(3)),eye(in(2));eye(in(3)),zeros(in(3),in(2))]);
[~,~,v1] = svd(double(d>0));
d        = T'*v1'*d*v1*T;
v        = v*v1*T;
J        = T'*v1'*J*v1*T;
D        = diag(J*sqrt(J*diag(d)));
M        = (v*D)';
end