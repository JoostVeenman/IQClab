function [ga,info] = mu2iqcanalysis(uobject,order,varargin)
% -------------------------------------------------------------------------
%
% Copyright:   This is copyrighted material owned by Novantec B.V.
% Terms:       IQClab is available under a Creative Commons
%              (Attribution-NoDerivatives 4.0 International (CC BY-ND 4.0))
%              license: https://creativecommons.org/licenses/by-nd/4.0/
%              For further information please visit iqclab.eu
%
% Author:      J.Veenman
% Date:        21-05-2022
% 
% -------------------------------------------------------------------------
%
% Description: Given the uncertain system LFT(Delta,M), being an
%              uncertainty object, with M and Delta being LFT this function
%              performs a worst-case IQC-based perfromance analysis
%              corresponding to the command wcgain.
%
%              The code automatically sets up the IQC-analysis. The only
%              necessary additional input is "order", being an integer:
%
%                              order\in[0,1,...,N].
%
%              Note: This value corresponds to the McMillan degree of the
%              basis function Psi, which equals the Length + 1.
%
%              The advantage of this function is that it is not necessary
%              to thing about setting up the IQC analysis. This is all done
%              in the code.
%
%              In addition, similar to the function iqcanalysis, it is
%              possible to specify several properties, which should come in
%              pairs
%
%                1.) For the option 'Parser', one can specify 'LMIlab' or
%                    'Yalmip'.
%
%                2.) For the option 'Solver', one can specify 'mincx', when
%                    using LMIlab as parser, or any solver compatible with
%                    the parser 'Yalmip'.
%
%                3.) For the option 'SolChk', one can specify 'on' or 'off'
%
% -------------------------------------------------------------------------
I.Parser                 = 'LMIlab';
I.Solver                 = 'mincx';
I.SolChk                 = 'on';

if nargin > 2
    m                    = length(varargin);
    if mod(m,2) ~= 0
        error('Error: The input arguments should come in pairs')
    else
        j                = linspace(2,m,m/2);
    end
    for i = 1:length(j)
        prop             = varargin{j(i)-1};
        switch prop
            case 'Parser'
                I.Parser = varargin{j(i)};
            case 'Solver'
                I.Solver = varargin{j(i)};
            case 'SolChk'
                I.SolChk = varargin{j(i)};
        end
    end
end
[M,d]                    = fLFTdata(uobject);
n                        = length(d.type);
for i = 1:n
    if strcmp(d.type{i},'ureal')
        d.unc{i}         = iqcassign(d.unc{i},'ultis','Length',order+1);
    elseif strcmp(d.type{i},'ultidyn')
        d.unc{i}         = iqcassign(d.unc{i},'ultid','Length',order+1);
    end
end

if isfield(d,'perf') && isobject(d.perf)
    d.unc{n+1}           = d.perf;
end
prob                     = iqcanalysis(M,d.unc,'SolChk',I.SolChk,'Parser',I.Parser,'Solver',I.Solver,'eps',1e-7);
ga                       = prob.gamma;

if nargout == 2
    info                 = prob;
end









