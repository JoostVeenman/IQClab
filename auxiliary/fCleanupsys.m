function Gcleaned = fCleanupsys(G,threshold)
% -------------------------------------------------------------------------
%
% Copyright:   This is copyrighted material owned by Novantec B.V.
% Terms:       IQClab is available under a Creative Commons
%              (Attribution-NoDerivatives 4.0 International (CC BY-ND 4.0))
%              license: https://creativecommons.org/licenses/by-nd/4.0/
%              For further information please visit iqclab.eu
%
% Author:      J.Veenman
% Date:        22-08-2011
% 
% -------------------------------------------------------------------------
%
% Description: This function replaces all elements from the state-space
%              realization matrices of the system "G" that have an absolute
%              value smaller or equal than "threshold" with zero.
%
% Syntax:      Gcleaned = fCleanupsys(G,threshold)
%
% -------------------------------------------------------------------------
A = G.a;
B = G.b;
C = G.c;
D = G.d;

for i = 1:size(A,1)
    for j = 1:size(A,2)
        if abs(A(i,j)) <= threshold
            A(i,j) = 0;
        end
    end
end
for i = 1:size(B,1)
    for j = 1:size(B,2)
        if abs(B(i,j)) <= threshold
            B(i,j) = 0;
        end
    end
end
for i = 1:size(C,1)
    for j = 1:size(C,2)
        if abs(C(i,j)) <= threshold
            C(i,j) = 0;
        end
    end
end
for i = 1:size(D,1)
    for j = 1:size(D,2)
        if abs(D(i,j)) <= threshold
            D(i,j) = 0;
        end
    end
end
Gcleaned = ss(A,B,C,D);
end