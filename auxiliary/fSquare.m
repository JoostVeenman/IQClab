function [Psih,Mh,Z,Rm] = fSquare(Psi1,Psi2,M,options)
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
% Description: This function squares the outer-factors of the
%              iqc-mulitpliers 
%
%              Pi = [Psi1,Psi2]^*M[Psi1,Psi2]
%
%              where:
%                1.) M = M'
%
%                2.) [Psi1^* M Psi1] > 0
%
%                3.) [Psi2^* M Psi2]-[Psi2^* M Psi1][Psi1^* M Psi1]^-1*...
%                                                       [Psi1^* M Psi2] < 0
%
% Syntax:      [Psih,Mh,Z,Rm]=fSquare(Psi1,Psi2,M,options)
%
% Usage:       As input one should provide: 
%
%                                                   [ Ai | Bi ]
%                - the minimal realizations: Psii = [----|----],  i = 1,2
%                                                   [ Ci | Di ]
%                - the symmetric matrix:     M
%                - the option:               options = 'eq' or 'ineq'
%
%              As output (in case options = 'eq') one obtains:
%
%                - the factorization
%
%                  Pi=[Psi1,Psi2]^*M[Psi1,Psi2] = [*]^*[I, 0][Psi1h,Psi3h]
%                                                 [*]  [0,-I][  0  ,Psi2h]
%
%                  where Psi1h, Psi2h, Psi3h, Psi1h^-*, Psi2h^-1 are
%                  stable.
%
%                - the matrix Mh=[I,0;0,-I]
%
%                - the controllable realization
%
%                                         [A1 ,0   | B1 ,0  ]
%                  Psih = [Psi1h,Psi3h] = [0  ,A2c | 0  ,B2c]
%                         [  0  ,Psi2h]   [--------|--------]
%                                         [C1h,C3c | D1h,D3c]
%                                         [0  ,C2c | 0  ,D2c]
%
%                - the matrix Z=[Z11,Z12;Z12',Z22] that certifies the
%                  linear matrix equation:
%
%                      [0   ,0  ,Z11 ,Z12,0, 0, 0][I  ,0     ,0,  0  ]
%                      [0   ,0  ,Z12',Z22,0, 0, 0][0  ,I     ,0,  0  ]
%                      [Z11 ,Z12,0   ,0  ,0, 0, 0][A1 ,0     ,B1 ,0  ]
%                  [x]'[Z12',Z22,0   ,0  ,0, 0, 0][0  ,A2c   ,0  ,B2c] = 0 (1)
%                      [0   ,0  ,0   ,0  ,I, 0, 0][C1h,C3c   ,D1h,D3c]
%                      [0   ,0  ,0   ,0  ,0,-I, 0][0  ,C2c   ,0  ,D2c]
%                      [0   ,0  ,0   ,0  ,0, 0,-M][C1 ,[0,C2],D1 ,D2 ]
%
%
%
%                - The realization matrices (as a structure RM):
%                  C1h, D1h, D1hi=D1h^-1, A2c, B2c, C2c, D2c, D2ci=D2c^-1,
%                  C3c, D3c 
%
%              As output (in case options = 'ineq') one obtains:
%
%                - the new factorization
%
%                  Pi = [*]^*M[Psi1,Psi2] < [*]^*[I, 0][Psi1h,Psi3h]
%                                           [*]  [0,-I][  0  ,Psi2h]
%
%                  where Psi1h, Psi2h, Psi3h, Psi1h^-*, Psi2h^-1 are
%                  stable. 
%
%                - the matrix Z=[Z11,Z12;Z12',Z22]+R that renders (1)
%                  strictly negative (where R is a perfurbation matrix)
%
% -------------------------------------------------------------------------

% solve first Algebraic Riccati Equation
Z11            = care(-Psi1.a,-Psi1.b,Psi1.c'*M*Psi1.c,Psi1.d'*M*Psi1.d,Psi1.c'*M*Psi1.d);

% Construct new realization matrices
D1h            = chol(Psi1.d'*M*Psi1.d);
D1hi           = fInv(D1h);
C1h            = (D1h')\(Psi1.d'*M*Psi1.c-Psi1.b'*Z11);
Pt             = blkdiag(M,-eye(size(Psi1.b,2)));
P1h            = (Psi1.c'-C1h'/(D1h')*Psi1.d')*M;
P2h            = (D1h')\Psi1.d'*M;

% remove uncontrolable modes
K              = ctrb(-Psi1.a'+C1h'/(D1h')*Psi1.b',[P1h*Psi2.c,P1h*Psi2.d]);
[L,U]          = lu(K);
s              = svd(U);
r              = sum(s>1e-10);
L1             = L(:,1:r);
Li             = fInv(L);
L1i            = Li(1:r,:);
L1e            = blkdiag(L1,eye(size(Psi2.a,1)));
L1ie           = blkdiag(L1i,eye(size(Psi2.a,1)));

% Construct new realization matrices
A2c            = L1ie*[-Psi1.a'+C1h'/(D1h')*Psi1.b',P1h*Psi2.c;zeros(size(Psi1.a,2),size(Psi2.a,1)),Psi2.a]*L1e;
B2c            = L1ie*[P1h*Psi2.d;Psi2.b];
C3c            = [-(D1h')\Psi1.b',P2h*Psi2.c]*L1e;
D3c            = P2h*Psi2.d;
C2t            = [zeros(size(Psi2.c,1),size(Psi1.a,1)),Psi2.c;-(D1h')\Psi1.b',P2h*Psi2.c]*L1e;
D2t            = [Psi2.d;D3c];

% Solve second Algebraic Riccati Equation
[sbal,T]       = ssbal(ss(A2c,B2c,C2t,0));
Z22            = -care(sbal.a,sbal.b,sbal.c'*Pt*sbal.c,D2t'*Pt*D2t,sbal.c'*Pt*D2t);
Z22            = T'*Z22*T;

% Generate new factorization 
D2c            = chol(-D2t'*Pt*D2t);
D2ci           = fInv(D2c);
C2c            = (D2c')\(B2c'*Z22-D2t'*Pt*C2t);
Psih           = ss(blkdiag(Psi1.a,A2c),blkdiag(Psi1.b,B2c),[C1h,C3c;zeros(size(C2c,1),size(C1h,2)),C2c],...
                 [D1h,D3c;zeros(size(D2c,1),size(D1h,2)),D2c]);
Mh             = blkdiag(eye(size(C1h,1)),-eye(size(C2c,1)));

% Construct the solution of the general indefinite ARE
Z12            = [L1,zeros(size(L1,1),size(Psi2.a,1))];
Z              = [Z11,Z12;Z12',Z22];

% Improve the solution by Newton iteration
Apsi           = Psih.a;
Bpsi           = Psih.b;
Cpsi           = [Psi1.c,zeros(size(Psi1.c,1),size(A2c,1)-size(Psi2.a,1)),Psi2.c];
Dpsi           = [Psi1.d,Psi2.d];
R              = fInv(Dpsi'*M*Dpsi);

for i = 1
    Fx         = fHe(-Z*Apsi)+Cpsi'*M*Cpsi-(Bpsi'*-Z+Dpsi'*M*Cpsi)'*R*(Bpsi'*-Z+Dpsi'*M*Cpsi);
    Acl        = Apsi-Bpsi*R*(Dpsi'*M*Cpsi-Bpsi'*Z);
    H          = lyap(Acl',Fx); % solve (Apsi+R*-Z)'*H+H*(A+R*-Z)+Fx=0 => Z=Z+H
    Z          = Z-0.5*fHe(H);
end

% Redefine system matrices
Z11            = Z(1:size(Psi1.a),1:size(Psi1.a));
Z22            = Z(size(Psi1.a)+1:end,size(Psi1.a)+1:end);

C1h            = (D1h')\(Psi1.d'*M*Psi1.c-Psi1.b'*Z11);
P1h            = (Psi1.c'-C1h'/(D1h')*Psi1.d')*M;

% remove uncontrolable modes
K              = ctrb(-Psi1.a'+C1h'/(D1h')*Psi1.b',[P1h*Psi2.c,P1h*Psi2.d]);
[L,~]          = lu(K);
L1             = L(:,1:r);
Li             = fInv(L);
L1i            = Li(1:r,:);
L1e            = blkdiag(L1,eye(size(Psi2.a,1)));
L1ie           = blkdiag(L1i,eye(size(Psi2.a,1)));

% Construct new realization matrices
A2c            = L1ie*[-Psi1.a'+C1h'/(D1h')*Psi1.b',P1h*Psi2.c;zeros(size(Psi1.a,2),size(Psi2.a,1)),Psi2.a]*L1e;
B2c            = L1ie*[P1h*Psi2.d;Psi2.b];
C3c            = [-(D1h')\Psi1.b',P2h*Psi2.c]*L1e;
D3c            = P2h*Psi2.d;
C2t            = [zeros(size(Psi2.c,1),size(Psi1.a,1)),Psi2.c;-(D1h')\Psi1.b',P2h*Psi2.c]*L1e;
D2t            = [Psi2.d;D3c];

C2c            = (D2c')\(B2c'*Z22-D2t'*Pt*C2t);
Psih           = ss(blkdiag(Psi1.a,A2c),blkdiag(Psi1.b,B2c),[C1h,C3c;zeros(size(C2c,1),size(C1h,2)),C2c],...
                 [D1h,D3c;zeros(size(D2c,1),size(D1h,2)),D2c]);

Rm.C1h         = C1h;
Rm.D1h         = D1h;
Rm.D1hi        = D1hi;
Rm.A2c         = A2c;
Rm.B2c         = B2c;
Rm.C2c         = C2c;
Rm.D2c         = D2c;
Rm.D2ci        = D2ci;
Rm.C3c         = C3c;
Rm.D3c         = D3c;

% Check the results
        Middle = blkdiag(fOblkdiag(Z),Mh,-M);
        Outer  = [fTriang(eye(size(Psih.a,1)),Psih.a,Psih.b,3);Psih.c,Psih.d;Cpsi,Dpsi];
        sqerr  = norm(Outer'*Middle*Outer);
switch options
    case 'eq'
        disp(['Squaring error  <<1: ', num2str(sqerr)]);
    case 'ineq'
        scl    = 10;
        Psi    = ss(blkdiag(Psi1.a,Psi2.a),blkdiag(Psi1.b,Psi2.b),[Psi1.c,Psi2.c],[Psi1.d,Psi2.d]);
        X      = care(Psi.a,Psi.b,scl*Psi.c'*Psi.c-scl*eye(size(Psi.c'*Psi.c,1)),scl*Psi.d'*Psi.d,scl*Psi.c'*Psi.d);
        Middle = blkdiag(fOblkdiag(X),scl*eye(size(Psi.d,1)));
        Outer  = [fTriang(eye(size(Psi.a,1)),Psi.a,Psi.b,3);Psi.c,Psi.d];
        eig(Outer'*Middle*Outer);

        % Compute the perturbation matrix W
        nx4    = size(A2c,1)-size(Psi2.a,1);
        Au     = A2c(1:nx4,1:nx4);
        Auo    = [zeros(nx4,size(Psi1.a,1)),A2c(1:nx4,nx4+1:end)];
        Bu     = [zeros(nx4,size(Psi1.b,2)),B2c(1:nx4,:)];

        T      = chol(Outer'*Middle*Outer);
        R      = lyapchol(Au,eye(nx4));
        Xu     = -R'*R;
        Outer  = Xu*[Auo,Bu]/T;
        U      = Outer'*fInv(fHe(Xu*Au))*Outer;
        eps    = 0.99999/max(eig(U));

        T1     = blkdiag(fJ(size(Psi1.a,1),nx4)',eye(size(Psi2.a,1)));
        T2     = blkdiag(zeros(size(Psi1.a,1)),eps*Xu,zeros(size(Psi2.a,1)));
        Xe     = T1'*X*T1+T2;

        del    = sqerr;
        Middle = blkdiag(fOblkdiag(del*Xe),del*eye(size(M,1)));
        Outer  = [fTriang(eye(size(Psih.a,1)),Psih.a,Psih.b,3);Cpsi,Dpsi];
        eig(Outer'*Middle*Outer);

        % Check the results
        Middle = blkdiag(fOblkdiag(Z-del*Xe),Mh,-M-del*scl*eye(size(M,1)));
        Outer  = [fTriang(eye(size(Psih.a,1)),Psih.a,Psih.b,3);Psih.c,Psih.d;Cpsi,Dpsi];
        X      = Outer'*Middle*Outer;
        disp(['Squaring error      <<1: ', num2str(norm(X))]);
        disp(['Largest Eigenvalues  <0: ', num2str(max(eig(X)))]);
end
end