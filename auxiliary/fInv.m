function B = fInv(A)
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
% Description: This function computes the inverse of A while exploiting
%              additional properties
%
% Syntax:      B = fInv(A)
%
% -------------------------------------------------------------------------

% Check if A is square
[noa,nia]     = size(A);
if abs(noa-nia) > 0
    error('A should be a square matrix');
end
if norm(triu(A)-A) < 1e-15 && abs(prod(diag(A))) > 1e-15
    A         = triu(A);
    B         = eye(noa)/A;
elseif norm(tril(A)-A) < 1e-15 && abs(prod(diag(A))) > 1e-15
    A         = tril(A);
    B         = eye(noa)/A;
else
    [U,S,V]   = svd(A);
    s         = diag(S);
    r         = sum(s > 1e-15);
    % Check if A is invertible
    if length(s)-r > 0
       error('The matrix A is ill-conditioned');
    end
    if norm(U-V) < 1e-15 % check if A is positive definite
        T     = U/sqrt(S);
        B     = T*T';
    elseif norm(U+V) < 1e-15 % check if A is negative definite
        T     = U/sqrt(S);
        B     = -T*T';
    else
        L     = V/sqrt(S);
        R     = sqrt(S)\U';
        if norm(abs(U)-abs(V)) < 1e-15 % Check if A is symmetric
            B = 0.5*fHe(L*R);
        else
            B = L*R;
        end
    end
end
end
