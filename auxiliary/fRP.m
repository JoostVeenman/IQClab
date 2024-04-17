function [Bound,muinfo,RPmargin] = fRP(N,om,Usc)
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
% Description: fRP computes the robust performnace margins using the mussv
%              command from the Robust Control Toolbox.
%
% Syntax:      [Bound,muinfo,RSmargin] = fRP(N,om,Usc)
%
% Usage:       fRP(N,om) computes the robust performance margins for the
%              uncertain plant N and for the frequency grid om =
%              [om1,...,omN]. It returns the mu ssv upper- and lower-bounds
%              by means of the mussv command from the Robust Control
%              Toolbox. The outputs are: 
%              1.) Bound:    A structure with the lower- and upper-bounds
%                            for each frequency om in{om1,...,omN}.
%              2.) muinfo:   The mu-information obtained during the
%                            comptations.
%              3.) RPmargin: The robust performance margin
%
%              If specifying the 3rd input (i.e. fRP(N,om,Usc)), one can
%              rescale the analysis as follows:
%
%              Given the uncertain plant N = [M,N12;N21,N22] the
%              uncertainty scaling factor Usc scales the size of the
%              uncertainty such that
%
%              1.) mu_Delta(M(i\omega) <= gamma/Usc
%              2.) mu_Delta(N(i\omega) <= gamma
%
%              for all ||Delta|| < Usc/gamma
%
%              If not specifies, Usc is set to the default value of 1.
%
% -------------------------------------------------------------------------
[M,~,Blk]      = lftdata(N);
BS             = [];
DS             = 0;

if nargin < 3
    Usc = 1;
end
for i = 1:length(Blk)
    if strcmp(Blk(i).Type,'ureal')
        BS     = [BS;-Blk(i).Occurrences,0];
    elseif strcmp(Blk(i).Type,'ultidyn')
        BS     = [BS;Blk(i).Occurrences,0];
    end
    DS         = DS + Blk(i).Occurrences;
end

nBS            = sum(abs(BS(:,1)));
[so,si]        = size(N.NominalValue.d);
BS             = [BS;si,so];
sc             = blkdiag(Usc*eye(nBS),eye(si));
Mfrd           = frd(M*sc,om);
[Bound,muinfo] = mussv(Mfrd,BS);
RPmargin       = 1/max(squeeze(muinfo.bnds(1).Response));
end