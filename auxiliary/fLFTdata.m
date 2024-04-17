function [M,d] = fLFTdata(uobject)
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
% Date:        12-08-2021
% 
% -------------------------------------------------------------------------
%
% Description: This function extracts the LFT data from an uncertain object
%              and generates the following information as input for the
%              IQC-analysis: in the structure d
%                 - The open-loop LTI plant M
%                 - The structure d with the following fields:
%                       - d.name        (Name of the uncertainty)
%                       - d.size        (Size of the uncertainty)
%                       - d.type        (Uncertainty type: ltis or ltid)
%                       - d.occurrences (Number of repetitions)
%                       - d.bounds      (Bounds on the uncertainty)
%                       - d.in          (Channels ids. that connect the
%                                        output of the uncertainty to the
%                                        input of M)
%                       - d.out         (Channels ids. that connect the
%                                        input of the uncertainty to the
%                                        output of M)
%                       - d.pin         (Performance output channel ids.)
%                       - d.pout        (Performance input channel ids.)
%                       - d.unc         (IQC uncertainty objects)
%                       - d.perf        (IQC performance objects)
% -------------------------------------------------------------------------

[M,D,B]               = lftdata(uobject);

for i = 1:length(B)
    d.name{i}        = B(i).Name;
    d.size(i,1:2)    = B(i).Size;
    d.type{i}        = B(i).Type;
    d.occurrences(i) = B(i).Occurrences;
    if strcmp(B(i).Type,'ureal')
        d.bounds{i}  = eval(['D.Uncertainty.',B(i).Name,'.Range']);
    elseif strcmp(B(i).Type,'ultidyn')
        d.bounds{i}  = eval(['D.Uncertainty.',B(i).Name,'.Bound']);
    end
end

jin                  = 1;
jout                 = 1;
for i = 1:length(d.occurrences)
    d.in{i}          = [jin:jin+d.occurrences(i)*d.size(i,1)-1];
    d.out{i}         = [jout:jout+d.occurrences(i)*d.size(i,2)-1];
    jin              = jin+d.occurrences(i)*d.size(i,1);
    jout             = jout+d.occurrences(i)*d.size(i,1);
end

d.pin                = d.in{end}(end)+1:size(M,2);
d.pout               = d.out{end}(end)+1:size(M,1);

for i = 1:length(B)
    if strcmp(d.type{i},'ureal')
        d.unc{i}     = iqcdelta(d.name{i},'InputChannel',d.in{i} ,'OutputChannel',d.out{i},'Bounds',d.bounds{i});
    elseif strcmp(d.type{i},'ultidyn')
        d.unc{i}     = iqcdelta(d.name{i},'InputChannel',d.in{i} ,'OutputChannel',d.out{i},'StaticDynamic','D','NormBounds',d.bounds{i});
    end
end

if ~isempty(d.pin) && ~isempty(d.pout)
    d.perf           = iqcdelta('perf','InputChannel',d.pin,'OutputChannel',d.pout,'ChannelClass','P','PerfMetric','L2');
end
end