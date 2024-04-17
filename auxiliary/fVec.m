function x = fVec(A,type)
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
% Date:        11-04-2012
% 
% -------------------------------------------------------------------------
%
% Description: This function vectorizes the matrices in different fashions
%
% Syntax:      x = fVec(A,type)
%
% Usage:       The functiuon x = fVec(A,type) vectorizes the matrix A as
%              follows:
%
%              1.) For type = 'col' (default) the function performs a
%                  column-wise vectorization:
%
%                        [a1,a4]        [a1]
%                    A = [a2,a5] -> x = [..]
%                        [a3,a6]        [a6]
%
%              2.) For type = 'row' the function performs a column-wise
%                  vectorization:
%
%                        [a1,a2]        [a1]
%                    A = [a3,a4] -> x = [..]
%                        [a5,a6]        [a6]
%
%              3.) For type = 'sym' the function performs a symmetry-wise
%                  vectorization:
%
%                        [a1,a2,a3]        [a1]
%                    A = [a2,a4,a5] -> x = [..]
%                        [a3,a5,a6]        [a6]
%
%              4.) For type = 'skw' the function performs a skew-wise
%                  vectorization:
%
%                        [0 ,-a1,-a2]        [a1]
%                    A = [a1,  0,-a3] -> x = [..]
%                        [a2, a3,  0]        [a3]
%
% -------------------------------------------------------------------------
[a1,a2] = size(A);
if nargin == 1
    x = reshape(A,a1*a2,1);
elseif nargin == 2
    switch type
        case 'col'
            x = reshape(A,a1*a2,1);
        case 'row'
            x = reshape(A',a1*a2,1);
        case 'sym'
            if a1-a2 == 0
                if norm(A-A')<1e-13
                    x = [];
                    for ik = 1:a1
                        for ki = 1:ik
                            x = [x;A(ik,ki)];
                        end
                    end
                else
                    error('A must be square');
                end
            else
                error('A must be a square matrix');
            end    
        case 'skw'
            if a1-a2 == 0
                if norm(A+A') < 1e-13
                    x = [];
                    for ik = 2:a1
                        for ki = 1:ik-1
                            x = [x;A(ik,ki)];
                        end
                    end
                else
                    error('A must be square');
                end
            else
                error('A must be a skew symmetric matrix');
            end    
    end
end