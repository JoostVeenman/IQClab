function B = fSchur(A,T)
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
% Date:        26-04-2012
% 
% -------------------------------------------------------------------------
%
% Description: This fuction computes the Schur complement of the block
%              structured symmetric matrix A with respect to the pioking
%              matrix T.
%
% Syntax:      B = fSchur(A,T)
%
% Usage:       Let T be a given picking matrix
%
%                    [0  ,0  ,...]
%                    [...,0  ,...]
%                    [0  ,...,...]
%                    [I  ,...,...]
%                T = [0  ,0  ,...]
%                    [...,I  ,...]
%                    [...,0  ,...]
%                    [...,...,...]
%                    [0  ,0  ,...]
%
%              and let Tc be the complementary picking matrix (i.e. the
%              orthogonal complement) such that [T,Tc]'*[T,Tc] = I.
%
%              For a given symmetric matrix A, the permutation operation
%              [T,Tc]'*A*[T,Tc] would then yield the structured matrix
%
%                [T,Tc]'*A*[T,Tc] = [T'*A*T ,T'*A*Tc ] = [Q ,S]
%                                   [Tc'*A*T,Tc'*A*Tc]   [S',R]
%
%              As output, the function subsequently computes the Schur
%              complement
%
%                B=R-S'*Q^-1*S
%
% -------------------------------------------------------------------------

if size(A,1) > size(A,2) || size(A,1) < size(A,2)
    error('The matrix A should be square');
elseif norm(A-A',inf) > 0
    error('The matrix A should be symmetric');   
elseif (1 > norm(T)) && (norm(T) > 1)
    error('The matrix T is not properly defined');
elseif norm(T - abs(T)) > 0
    error('The matrix T is not properly defined');
elseif size(T,2) - rank(T) > 0
    error('The matrix T is not properly defined');
elseif size(A,2) < size(T,1) || size(A,2) > size(T,1)
    error('The dimensions of A and T are not compatible');
end

To          = abs(null(T'));
A22         = To'*A*To;
[U,S,V]     = svd(A22);
s           = diag(S);
r           = sum(s>1e-15);

if length(s)-r == 0
    if min(real(eig(A22))) > 0
       T2   = chol(A22);
       A12t = A*To/T2;
       B    = T'*(A-A12t*A12t')*T;
    else
        L   = A*To*V/sqrt(S);
        R   = sqrt(S)\U'*To'*A;
        B   = T'*(A-0.5*fHe(L*R))*T;
    end
else
    error('The matrix To^T*A*To is ill-conditioned.');
end
end