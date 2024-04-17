function [Gb,Tb,Tbi] = fHinfbalreal(G)
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
% Date:        18-04-2012
% 
% -------------------------------------------------------------------------
%
% Description: This function computes the H-infinity balanced realization
%              matrices for the stable plant G
%
% Syntax:      [Gb,Tb,Tbi] = fHinfbalreal(G)
%
% Usage:       The function [Gb,Tb,Tbi] = fHinfbalreal(G) requires the
%              stable state space realization G = ss(A,B,C,D).
%
%              As output one obtains the H-infinity balanced realization
%
%              Gb = ss(Tb*A*Tbi,Tb*B,C*Tbi,D)
%
%              where Tbi = Tb^-1
%
%              This function follows the following proceedure:
%
%              Compute the H-infinity norm of G and consider the LMI
%              (bounded Real Lemma) 
%
%              [X*A + A'*X,X*B       ] + [C,D]'*gamma^-1*[C,D] < 0      (1)
%              [B'*X      , -gamma*I ]
%
%              and its dual
%
%              [-Y*A' - Y*A, -X*C'   ] - [B',D']'*gamma^-1*[B',D'] > 0  (2)
%              [-C*Y       , gamma I ]  
%
%              Define: Dh'*Dh = -(gamma^-1*D'*D-gamma I) < 0
%              Define: Dc*Dc' = -(gamma^-1*D*D'-gamma I) < 0
%
%              We conclude that (1) and (2) imply the AREs:
%
%              Zx*A + A'*Zx + gamma^-1*C'*C +
%                     + (*)'(Dh'*Dh)^-1*(Zx*B + gamma^-1*C'*D)'  = 0
%
%              Zy*A' + A*Zy + gamma^-1*B*B' +
%                     + (Zy*C' + gamma^-1*B*D')*(Dc*Dc')^-1*(*)' = 0
%
%              With the solutions Zx and Zy it is now possible to perform
%              standard balancing.
%
% -------------------------------------------------------------------------
if max(real(eig(G))) > -1e-10
    error('G must be a stable state-space object');
end

gamma   = norm(G,inf);
gamma   = 1.01*gamma;

Qh      = -G.c'*G.c/gamma;
Qc      = -G.b*G.b'/gamma;
Rh      = -G.d'*G.d/gamma+gamma*eye(size(G.d'*G.d,1));
Rc      = -G.d*G.d'/gamma+gamma*eye(size(G.d*G.d',1));
Sh      = G.c'*G.d/gamma;
Sc      = G.b*G.d'/gamma;

P       = care(-G.a,G.b,Qh,Rh,Sh,eye(size(G.a,1)));
Q       = care(-G.a',G.c',Qc,Rc,Sc,eye(size(G.a,1)));

if min(eig(P)) < 1e-10 || min(eig(Q)) < 1e-10
    error('P and Q should be positive definite');
end

Qsqrt   = sqrtm(Q);
[U,S,~] = svd(Qsqrt*P*Qsqrt);
Sigma   = sqrt(S);
Tb      = fInv(sqrt(Sigma))*U'*Qsqrt;
Tbi     = fInv(Tb);

if norm(Tb*P*Tb'-Tbi'*Q*Tbi,'fro') > 1e-12 || norm(Tb*Tbi-Tbi*Tb,'fro') > 1e-12
    error('Accuracy is less than 1e-12')
end
Gb      = ss(Tb*G.a*Tbi,Tb*G.b,G.c*Tbi,G.d);
end