function [Bound,muinfo,RSmargin] = fRS(N,om)
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
% Date:        13-12-2019
%
% -------------------------------------------------------------------------
%
% Description: fRS computes the robust stability margins using the mussv
%              command from the Robust Control Toolbox.
%
% Syntax:      [Bound,muinfo,RSmargin] = fRS(N,om)
%
% Usage:       fRS(N,om) computes the robust stability margins for the
%              uncertain plant N and for the frequency grid om =
%              [om1,...,omN]. It returns the mu ssv upper- and lower-bounds
%              by means of the mussv command from the Robust Control
%              Toolbox. The outputs are: 
%              1.) Bound:    A structure with the lower- and upper-bounds
%                            for each frequency om in{om1,...,omN}.
%              2.) muinfo:   The mu-information obtained during the
%                            computations.
%              3.) RSmargin: The robust stability margin
% -------------------------------------------------------------------------

[so,si]        = size(N.NominalValue.d);
[M,~,Blk]      = lftdata(N);
M              = M(1:end-so,1:end-si);
Mfrd           = frd(M,om);
BS             = [];

for i = 1:length(Blk)
    if strcmp(Blk(i).Type,'ureal')
        BS     = [BS;-Blk(i).Occurrences,0];
    elseif strcmp(Blk(i).Type,'ultidyn')
        BS     = [BS;Blk(i).Occurrences,0];
    end
end

[Bound,muinfo] = mussv(Mfrd,BS);
RSmargin       = 1/max(squeeze(muinfo.bnds(1).Response));
end