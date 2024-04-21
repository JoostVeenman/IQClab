function Gr = fRFZP(G,w)
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
% Date:        13-03-2020
% 
% -------------------------------------------------------------------------
%
% Description: Given the (possibly multi-variable) system G, this function
%              removes all poles and zeros after a given frequency w.
%
% Syntax:      Gr = fRFZP(G)
%              Gr = fRFZP(G,w)
%
% Usage:       As input one should specify:
%
%                - The plant G
%                - The cut-off frequency w (rad/s)
%
%              If nargin is 1, the cut-off frequency w must be entered in
%              the command window
%
%              As output one obtains the reduced plant Gr
%
% -------------------------------------------------------------------------

[n1,n2] = size(G);
[z,p,k] = zpkdata(G,'v');

if n1 == 1 && n2 == 1
    zstr{1,1} = z;
    pstr{1,1} = p;
    z         = zstr;
    p         = pstr;
end

if nargin == 1
    x = input('Enter the desired cut-off frequency w (rad/s): ');
%     [x,~] = ginput(1);
elseif nargin > 1
    x = w;
end

for i = 1:n1
    for j = 1:n2
        zr{i,j} = z{i,j}(abs(z{i,j})<x);
        pr{i,j} = p{i,j}(abs(p{i,j})<x);     
    end
end

Gr     = zpk(zr,pr,k);
Kr     = abs(freqresp(Gr,1e-4*x));
K      = abs(freqresp(G,1e-4*x));
Ki     = Kr == 0;
Kr(Ki) = 1e-3;
Ka     = K./Kr;
Gr     = ss(Ka.*Gr);
end