function A=ftriang(A1,A2,A3,options)
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
% Date:        29-03-2012
% 
% -------------------------------------------------------------------------
%
% Description: This function generates triangular matrices
%
% Syntax:      A = ftriang(A1,A2,A3,options)
%
% Usage:       This function generates one of the four the matrices
%
%              1.) With options = 1:
%
%                  A = [A1 A2]
%                      [0  A3]
%
%              2.) With options = 2:
%
%                  A = [0  A1]
%                      [A2 A3]
%
%              3.) With options = 3:
%
%                  A = [A1 0 ]
%                      [A2 A3]
%
%              4.) With options = 4:
%
%                  A = [A1 A2]
%                      [A3 0 ]
%
% -------------------------------------------------------------------------
switch options
    case 1
        A = [A1,A2;zeros(size(A3,1),size(A1,2)),A3];
    case 2
        A = [zeros(size(A1,1),size(A2,2)),A1;A2,A3];
    case 3
        A = [A1,zeros(size(A1,1),size(A3,2));A2,A3];
    case 4
        A = [A1,A2;A3,zeros(size(A3,1),size(A2,2))];
end 
end