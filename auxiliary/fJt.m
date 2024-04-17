function d = fJt(a,b,c)
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
% Date:        07-04-2011
% 
% -------------------------------------------------------------------------
%
% Description: This function generates the matrix
%                  [0] dim b x a
%              d = [I] dim a x a
%                  [0] dim c x a
%
% Syntax:      d = fCutFig(n,m)
%
% Usage:       1.) if nargin = 1:  d = [ 0_{a x a} ]
%                                      [ I_a       ]
%
%              2.) if nargin = 2:  d = [ 0_{b x a} ]
%                                      [ I_a       ]
%
%              3.) if nargin = 3:  d = [ 0_{b x a} ]
%                                      [ I_a       ]
%                                      [ 0_{c x a} ]
%
% -------------------------------------------------------------------------
if nargin == 1
    d = [zeros(a,a);eye(a)];
elseif nargin == 2
    d = [zeros(b,a);eye(a)];
elseif nargin == 3
    d = [zeros(b,a);eye(a);zeros(c,a)];
end
end