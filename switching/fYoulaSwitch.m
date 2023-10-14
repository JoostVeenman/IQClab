function [L,Q] = fYoulaSwitch(G,K,varargin)
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
% Date:        06-04-2020 
%              04-07-2023 (modified code to be able to process more that
%                          two controllers)
% 
% -------------------------------------------------------------------------
%
% Description: This function computes a Youla based switching scheme
%              between two continuous or dicrete time controllers K1, K2,
%              ..., KN.
%
% Syntax:      [L,Q] = fYoulaSwitch(G,K1,K2,KN)
%
% Usage:       Given the plant, G, seen by the controllers K1, K2, ..., KN,
%              where K1, K2, ..., KN are stabilizing G. This function
%              generates a Youla based switching scheme that quarantees
%              stability during the switching.
%
%              Indeed, switching between controllers can be implemented as
%              K = alpha1*K1 + alpha2*K2 + ... + alphaN*KN for 
%              alphai\in[0,1], \Sigma_1^N alphai = 1. However, for this
%              configuration it cannot be quaranteed that the plant remains
%              stable for alphai\in(0,1). 
%
%              Instead, one can consider the Youla parametrization based
%              approach, where the switching between the controllers 
%              K1, K2, ..., KN is implemented as
%
%                                  ---------
%                           y  --->|       |--->  u
%                                  |   L   |
%                        |-------->|       |--------|
%                        |         ---------        |
%                     yq |                          | uq
%                        |   -----------   ------   |
%                        O<--| alpha1  |<--| Q1 |<--|
%                        ^   -----------   ------   |
%                        |                          |
%                        |   -----------   ------   |
%                        O<--| alpha2  |<--| Q2 |<--|
%                        ^   -----------   ------
%                        |                          |
%                        |   -----------   ------   |
%                        |---| alphaN  |<--| QN |<--|
%                            -----------   ------
%
%              For this switching scheme, stability of the closed loop
%              interconnection is guaranteed for all alphai\in[0,1] and
%              \Sigma_1^N alphai = 1.
%
%              This proceeds as follows:
%
%              Given the coprime factorization G = N*M^-1 = Mt^-1*Nt and
%              the Bezout identity
%
%                 [Xt ,-Yt][M,Y] = I
%                 [-Nt, Mt][N,X]
%
%              one can parameterize all stabilizing controllers as
%
%                K = (Y-M*Q)*(X-N*Q)^-1=(Xt-Q*Nt)^-1*(Yt-Q*Mt)
%           
%              Here Q is any stable transfer matrix in RHinf such that X-NQ
%              is invertible. The transfer matrix Q is the celebrated Youla
%              parameter.
%
%              The stabilizing controller can also be written as a lower
%              linear fractional representation of a fixed transfer matrix
%              L with the free Youla parameter Q:
%
%                K = (Y-M*Q)*(X-N*Q)^-1 = lft(L,Q)
%
%              Given a controller and its coprime factors K = U*V^-1 =
%              Vt^-1*Ut that stabilizes G, the corresponding Youla
%              parameter Q can be taken to be:
%
%                Q = -(Xt*U - Yt*V)*(-Nt*U + Mt*V)^-1 = 
%                              = -(Xt*K - Y)*(-Nt*K + Mt)^-1\in RHinf
%
%
%              Performing this factorization for K1 and K2 we obtain the
%              switching scheme
%
%                lft(L,Q) with Q = alpha1*Q1 + alpha2*Q2  + ... + alphaN*QN
%
%              This parametrization of the controller guarantees that the
%              closed-loop plant remains stable for all alphai(t)\in[0,1].
%
%              As input one should provide:
%
%                # The continuous or discrete time plant G with
%                  stabilizable and detectable realization G = ss(A,B,C,D).
%                # The continuous or discrete time controllers Ki =
%                  ss(Aki,Bki,Cki,Dki) which (internally) stabilize G.
%
%              As output one obtaines the realizations of L, Q1, ..., QN.
%
% -------------------------------------------------------------------------

Nk = nargin - 1;
if nargin <= 1
    error('This function requires the plant G and at least one controller as input');
end

G               = ss(G);
Ts              = G.Ts;

% put all controllers in a cell
Kc{1}           = ss(K);
for i = 1:length(varargin)
    Kc{i+1}     = ss(varargin{i});
end

% Check if all sampling times are the same
for i = 1:Nk
    if G.Ts ~= Kc{i}.Ts
        error('All systems should have the same sample time.');
    end
end

% Check if Ki, i = 1, ..., N are stabilizing the plant G
[nGo,nGi]       = size(G.d);
nx              = size(G.a,1);

for i = 1:Nk
    G           = ss(G);
    eCL         = eig((eye(nGo)-G*Kc{i})^-1);
    if Ts == 0
        if sum(real(eCL) > -1e-16) > 0
            error('One of the controllers is not stabilizing the plant G.');
        end
    else
        if sum(abs(eCL) > 1-1e-16) > 0
            error('One of the controllers is not stabilizing the plant G.');
        end
    end
end

% check the eigenvalues of G
if Ts == 0
    if  max(real(eig(G))) < -1e-16
        Tr      = ss(G.a,[G.b,zeros(nx,nGo)],[zeros(nGi,nx);G.c],[eye(nGi),zeros(nGi,nGo);G.d,eye(nGo)]);
        Tl      = ss(G.a,[G.b,zeros(nx,nGo)],[zeros(nGi,nx);-G.c],[eye(nGi),zeros(nGi,nGo);-G.d,eye(nGo)]);
        nmi     = nGi;
        nmo     = nGi;
        nxto    = nGi;
    else
        % compute the coprime factors and the Bezout identity of G
        [Cl,Cr] = fCoprime(G);
        Tr      = Cr.Tr;
        Tl      = Cl.Tl;
        nmi     = size(Cr.M.d,2);
        nmo     = size(Cr.M.d,1);
        nxto    = size(Cl.Xt.d,1);
    end
else
    if  max(abs(eig(G))) < 1-1e-16
        Tr      = ss(G.a,[G.b,zeros(nx,nGo)],[zeros(nGi,nx);G.c],[eye(nGi),zeros(nGi,nGo);G.d,eye(nGo)],Ts);
        Tl      = ss(G.a,[G.b,zeros(nx,nGo)],[zeros(nGi,nx);-G.c],[eye(nGi),zeros(nGi,nGo);-G.d,eye(nGo)],Ts);
        nmi     = nGi;
        nmo     = nGi;
        nxto    = nGi;
    else
        % compute the coprime factors and the Bezout identity of G
        [Cl,Cr] = fCoprime(G);
        Tr      = Cr.Tr;
        Tl      = Cl.Tl;
        nmi     = size(Cr.M.d,2);
        nmo     = size(Cr.M.d,1);
        nxto    = size(Cl.Xt.d,1);
    end 
end

% Compute L
A               = Tr.a;
B1              = -Tr.b(:,1:nmi);
B2              = Tr.b(:,nmi+1:end);
C1              = Tr.c(1:nmo,:);
C2              = Tr.c(nmo+1:end,:);
D11             = -Tr.d(1:nmo,1:nmi);
D12             = Tr.d(1:nmo,nmi+1:end);
D21             = -Tr.d(nmo+1:end,1:nmi);
D22             = Tr.d(nmo+1:end,nmi+1:end);
D22i            = D22^-1;
AL              = A-B2*D22i*C2;
BL              = [B2*D22i,B1-B2*D22i*D21];
CL              = [C1-D12*D22i*C2;-D22i*C2];
DL              = [D12*D22i,D11-D12*D22i*D21;D22i,-D22i*D21];
L               = ss(AL,BL,CL,DL,Ts);

A               = Tl.a;
B               = Tl.b;
C1              = Tl.c(1:nxto,:);
C2              = Tl.c(nxto+1:end,:);
D1              = Tl.d(1:nxto,:);
D2              = Tl.d(nxto+1:end,:);

for i = 1:Nk
    % Compute Youla parametrization of first controller
    [~,Kr]      = fCoprime(Kc{i});
    
    Ak          = Kr.N.a;
    Bk          = Kr.N.b;
    Ck          = [Kr.N.c;Kr.M.c];
    Dk          = [Kr.N.d;Kr.M.d];

    Di          = (D2*Dk)^-1;
    
    Aq          = [A-B*Dk*Di*C2,B*Ck-B*Dk*Di*D2*Ck;-Bk*Di*C2,Ak-Bk*Di*D2*Ck];
    Bq          = [B*Dk*Di;Bk*Di];
    Cq          = -[C1-D1*Dk*Di*C2,D1*Ck-D1*Dk*Di*D2*Ck];
    Dq          = -D1*Dk*Di;
    
    Q{i}        = ss(Aq,Bq,Cq,Dq,Ts);

    % Test if the parametrizations are stable
    if Ts == 0
        if sum(real(eig(Q{i})) > -1e-16) > 0
            error('One of the Youla parametrizations is not stable.');
        end
    else
        if sum(abs(eig(Q{i})) > 1-1e-16) > 0
            error('One of the Youla parametrizations is not stable.');
        end
    end
end
end