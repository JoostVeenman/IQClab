function [Cl,Cr] = fCoprime(G)
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
% Date:        05-04-2020
% 
% -------------------------------------------------------------------------
%
% Description: This function computes the coprime factors of the state
%              continuous or discrete time space realization
%              G = ss(A,B,C,D)
%
% Syntax:      [Tl,Tr] = fCoprime(G)
%
% Usage:       For each proper real-rational transfer matrix G there exist
%              eight RH_inf transfer matrices satisfying the equations
%
%              G = N*M^-1 = Mt^-1*Nt  and  [Xt , -Yt]*[M,Y] = [I, 0]    (1)
%                                          [-Nt,  Mt] [N,X]   [0, I]
%
%              The pairs (N,M) and (Nt,Mt) constitute the right and left
%              coprime factors of the transfer matrix G, respectively.
%
%              They are said to be normalized when each of the transfer
%              matrices [M;N], [Nt,Mt] is norm preserving.
%
%              Equation (1) is known as the Bezout identity. The state
%              space formulae for the normalized coprime factors, and their
%              certificates of coprimeness are given by:
%
%              Let G = ss(A,B,C,D) be a stabilizable and detectable
%              realization and choose R, S, At as:
%
%                # R'*R = (I+D'*D)^-1
%                # S*S' = (I+D*D')^-1
%                # At   = A-B*(I+D'*D)^-1*D'*C
%
%              Let P, Z be the unique stabilizing solutions for each of the
%              following algebraic Riccati equations (AREs), respectively:
%
%                # At'*P+P*At - P*B*(I+D'*D)^-1*B'*P + C'*(I+D*D')^-1*C = 0
%                # At*Z+Z*At' - Z*C'*(I+D*D')^-1*C*Z + B*(I+D'*D)^-1*B' = 0
%
%              Define the state feedback and observer gains F, H as:
%
%                # F = -(I+D'*D)^-1*(B'*P+D'*C)
%                # H = -(B*D'+Z*C')(I+D*D')^-1
%
%              Then
%
%                                [A+B*F  | BR , -HS^-1 ]
%                Tr = [M,Y]    = [-------|-------------]                (2)
%                     [N,X]      [F      | R  , 0      ]
%                                [C+D*F  | D*R, S^-1   ]
%
%                                [A+H*C  | -(B+H*D), H ]
%                Tl = [Xt,-Yt] = [-------|-------------]                (3)
%                     [-Nt,Mt]   [R^-1*F | R^-1    , 0 ]
%                                [S*C    | -S*D    , S ]
%
%              Satisfy (1) and (N,M), (Nt,Mt) are normalized
%
%              As input one should provide:
%
%                # The stabilizable and detectable realization G =
%                  ss(A,B,C,D)
%
%              As output one obains:
%
%                # The structures Cl and Cr with the realizations
%
%                  Cr: Tr, M,  N,  Y,  X
%                  Cl: Tl, Mt, Nt, Yt, Xt
%
% -------------------------------------------------------------------------

G           = ss(G);
A           = G.a;
B           = G.b;
C           = G.c;
D           = G.d;
Ts          = G.Ts;
[no,ni]     = size(D);

R           = (chol(eye(ni)+D'*D)')^-1;
S           = (chol(eye(no)+D*D'))^-1;
At          = A-B/(eye(ni)+D'*D)*D'*C;

if G.Ts == 0
    % Solve ARE: At'*P + P*At - P*B*(I+D'*D)^-1*B'*P + C'*(I+D*D')^-1*C = 0
    [P,~,~] = icare(At,B,C'*(eye(no)+D*D')^-1*C,eye(ni)+D'*D);
    
    % Solve ARE: At*Z + Z*At' - Z*C'*(I+D*D')^-1*C*Z + B*(I+D'*D)^-1*B' = 0
    [Z,~,~] = icare(At',C',B*(eye(ni)+D'*D)^-1*B',eye(no)+D*D');
else
    [P,~,~] = idare(At,B,C'*(eye(no)+D*D')^-1*C,eye(ni)+D'*D);
    [Z,~,~] = idare(At',C',B*(eye(ni)+D'*D)^-1*B',eye(no)+D*D');
end

% Define state feedback and observer gains F, H
F           = -(eye(ni)+D'*D)\(B'*P+D'*C);
H           = -(B*D'+Z*C')/(eye(no)+D*D');

% Define the right coprime factors
Cr.Tr       = ss(A+B*F,[B*R,-H/S],[F;C+D*F],[R,zeros(size(R,1),size(S,2));D*R,S^-1],Ts);
Cr.M        = ss(A+B*F,B*R,F,R,Ts);
Cr.N        = ss(A+B*F,B*R,C+D*F,D*R,Ts);
Cr.Y        = ss(A+B*F,-H/S,F,0,Ts);
Cr.X        = ss(A+B*F,-H/S,C+D*F,S^-1,Ts);

% Define the left coprime factors
Cl.Tl       = ss(A+H*C,[-(B+H*D),H],[R\F;S*C],[R^-1,zeros(size(R,1),size(S,2));-S*D,S],Ts);
Cl.Xt       = ss(A+H*C,-(B+H*D),R\F,R^-1,Ts);
Cl.Yt       = -ss(A+H*C,H,R\F,0,Ts);
Cl.Nt       = -ss(A+H*C,-(B+H*D),S*C,-S*D,Ts);
Cl.Mt       = ss(A+H*C,H,S*C,S,Ts);
end