function Y = fOblkdiag(X,O1,O2)
% -------------------------------------------------------------------------
%
% Copyright:   This is copyrighted material owned by Novantec B.V.
% Terms:       IQClab is available under a Creative Commons
%              (Attribution-NoDerivatives 4.0 International (CC BY-ND 4.0))
%              license: https://creativecommons.org/licenses/by-nd/4.0/
%              For further information please visit iqclab.eu
%
% Author:      J.Veenman
% Date:        07-04-2011
% 
% -------------------------------------------------------------------------
%
% Description: This function generates anti-block diagonal matrices
%
% Syntax:      Y = fOblkdiag(X)
%              Y = fOblkdiag(X,O1)
%              Y = fOblkdiag(X,O1,O2)
%
% Usage:       If nargin = 1: Generates the matrix
%  
%              Y = [ 0    X ]
%                  [ X^T  0 ]
%
%              If nargin = 2: Generates the matrix
%  
%                  [ 0    0   X ] 
%              Y = [ 0    O1  0 ],
%                  [ X^T  0   0 ]
%
%              where O1 is a square and symmetric matrix
%
%              If nargin=3: Generates the matrix
%  
%                  [ 0   0    0   0   0 ]
%                  [ 0   0    0   X   0 ] 
%              Y = [ 0   0    O1  0   0 ]
%                  [ 0   X^T  0   0   0 ]
%                  [ 0   0    0   0   0 ]
%
%              where: 
%               - O2 = [a,b]
%               - the left-upper zero block is "a times a"
%               - the right-lower zero block is "b times b"
%
% -------------------------------------------------------------------------
if nargin==1
    [nox,nix]=size(X);
    Y=[zeros(nox,nox),X;X',zeros(nix,nix)];
elseif nargin==2
    [nox,nix]=size(X);
    [noo,nio]=size(O1);
    if abs(noo-nio)>0
        error('O1 must be a square matrix');
    end
    Y=[zeros(nox,nox+nio),X;zeros(noo,nox),O1,zeros(noo,nix);X',zeros(nix,nio+nix)];
elseif nargin==3
    [nox,nix]=size(X);
    [noo,nio]=size(O1);
    if abs(noo-nio)>0
        error('O1 must be a square matrix');
    elseif size(O2,2)-size(O2,1)>1 || size(O2,2)-size(O2,1)<1
        error('O2 must be a 1 x 2 matrix');
    end
    Y=blkdiag(zeros(O2(1),O2(1)),[zeros(nox,nox+nio),X;zeros(noo,nox),O1,zeros(noo,nix);X',zeros(nix,nio+nix)],zeros(O2(2),O2(2)));
end
end