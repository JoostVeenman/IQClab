function d = fJ(a,b,c)
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
% Description: This function generates the matrix
%                  [0] dim c x a
%              d = [I] dim a x a
%                  [0] dim b x a
%
% Syntax:      d = fJ(a,b,c)
%
% Usage:       1.) if nargin = 1:  d = [ I_a       ]
%                                      [ 0_{a x a} ]
%
%              2.) if nargin = 2:  d = [ I_a       ]
%                                      [ 0_{b x a} ]
%
%              3.) if nargin = 3:  d = [ 0_{c x a} ]
%                                      [ I_a       ]
%                                      [ 0_{b x a} ]
%
% -------------------------------------------------------------------------
if nargin == 1
    d = [eye(a);zeros(a,a)];
elseif nargin==2
    d = [eye(a);zeros(b,a)];
elseif nargin==3
    d = [zeros(c,a);eye(a);zeros(b,a)];
end
end