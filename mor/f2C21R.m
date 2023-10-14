function Gr = f2C21R(G,w)
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
% Date:        13-03-2020
% 
% -------------------------------------------------------------------------
%
% Description:  Given the (possibly multi-variable) system G, this function
%               either allows you to replace: 
%
%                 1.) a complex pole or zero pair with one real pole or
%                     zero respectively. 
%
%                     x_old1 = a1 + i * b1
%                     x_old2 = a1 - i * b1
%
%                     x_new  = abs(x_old1)
%
%                 2.) two real poles or zeros in close proximity with one
%                     real pole or zero respectively.
%
%                     x_old1 = a1
%                     x_old2 = a2
%
%                     x_new  = sign(x_old1)*mean(abs(x_old1+x_old2))
%
%               Restrictions: x_old1 and x_old2 must have the same sign.
%
% Syntax:       Gr = fSmooth(G)
%               Gr = fSmooth(G,w)
%
% Usage:        As input one should specify:
%
%                 - G      : The plant G
%                 - [w1,w2]: The interval in which the poles or zeros are
%                            located (i.e. [w1,w2]). 
%                 - w3     : The gain-matching frequency w3 (i.e. the
%                            frequency where the gain of the old system G
%                            and the new system Gr is the same). 
%
%               The last two inputs must be:
%
%                 - entered manually through a dialog if nargin is 1. 
%                 - given in a vector w = [w1,w2,w3] if  nargin is 2. 
%
%               As output one obtains the reduced plant Gr.
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
    disp('Enter the (tight) interval in which the pole or zero pairs are located [w1,w2]  (rad/s): ');
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
        zt  = z{i,j};
        pt  = p{i,j};
        zc1 = logical((abs(zt)<w2)-(abs(zt)<w1));
        pc1 = logical((abs(pt)<w2)-(abs(pt)<w1));
        zc2 = logical((abs(zt)>w2)-(abs(zt)<w1));
        pc2 = logical((abs(pt)>w2)-(abs(pt)<w1));
        if sum(zc1) == 2
            x1 = find(zc1==1);
%             x2 = find(zc2==1);
            x3 = sign(real(z{i,j}(x1)));
            if x3(1)+x3(2) == 0
                zr{i,j} = z{i,j};
            else
                x4 = mean(abs(z{i,j}(x1)));
                x5 = x3(1)*x4;
                zr{i,j} = [z{i,j}(zc2);x5];
            end
        else
            zr{i,j} = z{i,j};
        end    
        if sum(pc1) == 2
            x1 = find(pc1==1);
%             x2 = find(pc2==1);
            x3 = sign(real(p{i,j}(x1)));
            if x3(1)+x3(2) == 0
                pr{i,j} = p{i,j};
            else
                x4 = mean(abs(p{i,j}(x1)));
                x5 = x3(1)*x4;
                pr{i,j} = [p{i,j}(pc2);x5];
            end
        else
            pr{i,j} = p{i,j};
        end
    end
end

Gr     = zpk(zr,pr,k);

if nargin == 1
    w3     = input('Enter the gain matching frequency w3 (rad/s): ');
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