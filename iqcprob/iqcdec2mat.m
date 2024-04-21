function var = iqcdec2mat(iqcprob,var)
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
% Date:        26-11-2019
%
% -------------------------------------------------------------------------
%
% Description: iqcdec2mat returns the numerical value of an LMI variable
%
% Syntax:      var = iqcdec2mat(iqcprob,var)
%
% Usage:       This function returns the numerical value of the specified
%              LMI variable "var". 
%
%              This function can be used once iqcsolve has resulted in a
%              feasible solution of the corresponding LMI problem.
%
% -------------------------------------------------------------------------
if strcmp(iqcprob.Parser,'Yalmip')
    var = value(var.var);
elseif strcmp(iqcprob.Parser,'LMIlab')
    if isempty(iqcprob.xsol)
        disp('Error: The LMIs are not feasible.');
    else
        var = dec2mat(iqcprob.lmi,iqcprob.xsol,var.var);
    end
end
end

