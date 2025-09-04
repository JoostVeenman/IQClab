function Gr = fSmooth(G,w)
% -------------------------------------------------------------------------
%
% Copyright:   This is copyrighted material owned by Novantec B.V.
% Terms:       IQClab is available under a Creative Commons
%              (Attribution-NoDerivatives 4.0 International (CC BY-ND 4.0))
%              license: https://creativecommons.org/licenses/by-nd/4.0/
%              For further information please visit iqclab.eu
%
% Author:      J.Veenman
% Date:        13-03-2020
% 
% -------------------------------------------------------------------------
%
% Description:  Given the (possibly multi-variable) system G, this function
%               removes all poles and zeros in a given frequency band
%               [w1,w2].
%
% Syntax:       Gr = fSmooth(G)
%               Gr = fSmooth(G,w)
%
% Usage:        As input one should specify:
%
%                 - The plant G
%                 - The interval in which the poles or zeros are located
%                   (i.e. [w1,w2]). 
%                 - The gain-matching frequency w3 (i.e. the frequency
%                   where the gain of the old system G and the new system
%                   Gr is the same).
%
%               The second last two inputs must be:
%
%                 - entered manually through a dialog if nargin is 1. 
%                 - given in a vector w = [w1,w2,w3] if  nargin is 2. 
%
%               As output one obtaines the reduced order plant Gr.
%
%-------------------------------------------------------------------------

[n1,n2] = size(G);
[z,p,k] = zpkdata(G,'v');

if n1 == 1 && n2 == 1
    zstr{1,1} = z;
    pstr{1,1} = p;
    z         = zstr;
    p         = pstr;
end

if nargin == 1
    disp('Enter smoothing window [w1,w2] (rad/s): ');
    w1 = input('w1 = ');
    w2 = input('w2 = ');
%     [x,~] = ginput(2);
%     w1 = x(1);
%     w2 = x(2);
elseif nargin > 1
    w1 = w(1);
    w2 = w(2);
end

for i = 1:n1
    for j = 1:n2
        z1      = z{i,j}(abs(z{i,j})<w1);
        z2      = z{i,j}(abs(z{i,j})>w2);
        zr{i,j} = [z1;z2];
        p1      = p{i,j}(abs(p{i,j})<w1);
        p2      = p{i,j}(abs(p{i,j})>w2);
        pr{i,j} = [p1;p2];
    end
end

Gr     = zpk(zr,pr,k);

if nargin == 1
    w3 = input('Enter the gain matching frequency w3 (rad/s): ');
%     disp('Enter the gain matching frequency w3 (rad/s): ');
%     [w3,~] = ginput(1);
elseif nargin > 1
    w3 = w(3);
end

Kr     = abs(freqresp(Gr,w3));
K      = abs(freqresp(G,w3));
Ki     = Kr == 0;
Kr(Ki) = 1e-3;
Ka     = Kr.\K;
Gr     = ss(Ka.*Gr);
end