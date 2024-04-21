function [To,Ti] = fT(mo,mi,co,ci)
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
% Date:        23-12-2019
% 
% -------------------------------------------------------------------------
%
% Description: This function creates permutation matrices which allow to
%              permute the in- and output channels of the uncertain plant M
%              in accordance with the given uncertainty blocks.
%
% Syntax:      [To,Ti] = fT(mo,mi,co,ci)
%
% Usage:       ci and co are vectors that contain the in- and output
%              channel numbers of the uncertain plant that are associated
%              with the uncertainty blocks.
%
%              mi and mo are the permutation matrices.
%
% -------------------------------------------------------------------------
lo = length(co);
if lo == 0
    loi = 0;
else
    for i = 1:lo
        loi(i) = length(co{i});
    end    
end

li = length(ci);
if li == 0
    lij = 0;
else 
    for j = 1:li
        lij(j) = length(ci{j});
    end
end

To = zeros(sum(loi),mo);
Ti = zeros(mi,sum(lij));

k  = 1;
for i = 1:lo
    for j = 1:loi(i)
        To(k,co{i}(j)) = 1;
        k = k + 1;
    end
end

k  = 1;
for i = 1:li
    for j = 1:lij(i)
        Ti(ci{i}(j),k) = 1;
        k = k + 1;
    end
end

end