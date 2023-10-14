function G=fmultss(G1,G2)
% -------------------------------------------------------------------------
%
% IQClab:      Version 3.03
% Copyright:   This is copyrighted material owned by Novantec B.V.
% Terms:       IQClab is available for non-commercial usage under a
%              Creative Commons (Attribution-NonCommercial-NoDerivatives
%              4.0 International (CC BY-NC-ND 4.0))license: 
%              https://creativecommons.org/licenses/by-nc-nd/4.0/
%              Commercial usage is only permitted with a commercial
%              license. For further information please visit iqclab.eu
% Author:      J.Veenman
% Date:        18-04-2012
% 
% -------------------------------------------------------------------------
%
% Description: Multiply two state-space objects G1 and G2 as G=G1*G2.
%
% Syntax:      G = fmultss(G1,G2)
%
% Usage:       Suppose that G1 and G2 admit the realization:
%
%                   [ A1 | B1 ]       [ A2 | B2 ]
%              G1 = [----|----], G2 = [----|----]
%                   [ C1 | D1 ]       [ C2 | D2 ]
%
%              Then we have:
%
%                          [ A1 B1*C2 | B1*D2 ]
%              G = G1*G2 = [ 0  A2    | B2    ]    
%                          [------------------]
%                          [ C1 D1*C2 | D1*D2 ]
%
% -------------------------------------------------------------------------

if abs(size(G1.d,2)-size(G2.d,1))>0
    error('G1 should have as many inputs as G2 has outputs');
end
if G1.Ts ~= G2.Ts
    error('Both systems should have the same sampling time');
end

Ts = G1.Ts;
A  = ftriang(G1.a,G1.b*G2.c,G2.a,1);
B  = [G1.b*G2.d;G2.b];
C  = [G1.c,G1.d*G2.c];
D  = G1.d*G2.d;
G  = ss(A,B,C,D,Ts);
end