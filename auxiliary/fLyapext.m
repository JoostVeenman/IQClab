function Xe = fLyapext(X,Y,type)
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
% Date:        26-04-2012
% 
% -------------------------------------------------------------------------
%
% Description: Given the symmetric and nonsingular matrices X and Y, this
%              function constructs the extended Lyapunov matrix Xe in
%              different fashions. 
%
% Syntax:      Xe = fLyapext(X,Y,type)
%
% Usage:       Given the symmetric and nonsingular matrices X and Y, the
%              function Xe = fLyapext(X,Y,type) constructs the extended
%              Lyapunov matrix Xe in different fashions such that
%
%                1.) Given X, Y, Xe^-1 yields the extended Lyapunov
%
%                    Xe^-1 = [ Y , * ]
%                            [ * , * ]
%
%                    This is obtained for type = 1 - 6
%
%                2.) Given X, Y, Xe^-1 yields the extended Lyapunov
%
%                    Xe^-1 = [ Y^-1 , * ]
%                            [ *    , * ]
%
%                    This is obtained for type = 7 - 12
%
%                3.) Given X = [X11 , X12], Y = [Y11 , Y12]
%                              [X12', X22]      [Y12', Y22]
%                             
%                    and 
%
%                    [Y11 + X11, Y12, -X12]
%                    [Y12'     , Y22, I   ] > 0
%                    [-X12'    , I  , X22 ]
%
%                    Xe^-1 yields the extended Lyapunov matrix
%
%                         [X11,  X12 ,X13] 
%                    Xe = [X12', X22 ,X23], with
%                         [X13', X23',X33]
%
%                         [X11, X12 , X13 ]^-1 [ Y11 ,Y12 , * ]
%                    Xe = [X12',X22 , X23 ]   =[ Y12',Y22 , * ]
%                         [X13',X23', X33 ]    [  *  , *  , * ]
%
%                    This is obtained for type = 13 - 15
%
%
%               Types 1-14 are constructed as follows:
%
%                - type = 1:  Xe = [X, I          ]
%                                  [I, (X-Y^-1)^-1]
%
%                - type = 2:  Xe = [X     , X-Y^-1]
%                                  [X-Y^-1, X-Y^-1]
% 
%                - type = 3:  If Y is positive definite, let Y=Yc'*Yc > 0
%                           ( Choleski factorization) and Yci = Yc^-1, then 
%
%                             Xe = [X        , X*Yc'-Yci]
%                                  [Yc*X-Yci', Yc*X*Yc'-I]
%
%                             If Y is indefinite, let Y = M'*J*M, with
%                             J = blkdiag(I,-I) and Mi = M^-1, then
%
%                             Xe = [X        , X*M'-Mi*J]
%                                  [M*X-J*Mi', M*X*M'-J ]
%
%                - type = 4:  Xe = [X    , X*Y-I  ]
%                                  [Y*X-I, Y*X*Y-Y]
%
%                - type = 5:  If X and Y are positive definite, let X =
%                             M'*M, Y = N'*N, (Choleski factorization) Mi =
%                             M^-1, Ni = N^-1 and u*s*v' = Ni'*Mi-N*M'.
%                             Then
%
%                             Xe = [Mi',0        ]^-1[M  , V*sqrt(s)]
%                                  [N  ,U*sqrt(s)]   [Ni', 0        ]
%
%                             If X or Y is indefinite, let X = M'*Jx*M, Y =
%                             N'*Jy*N, Jx = blkdiag(I,-I), Jy =
%                             blkdiag(I,-I), Mi = M^-1, Ni = N^-1 and 
%                             usv' = (Ni'*Mi - Jy*N*M'*Jx). Then
%
%                             Xe = [Mi' , 0        ]^-1[Jx*M, V*sqrt(s)]
%                                  [Jy*N, U*sqrt(s)]   [Ni' , 0        ]
%
%                - type = 6:  If X is positive definite, let X = Xc'Xc>0
%                             (Choleski factorization) and Xci=Xc^-1, then 
%
%                             Xe = [X            ,   Tx'-Y^-1*Txi ]
%                                  [Tx'-Txi'*Y^-1, I-Txi'*Y^-1*Txi ]
%
%                             If X is indefinite, let X = M'*J*M, with
%                             J = blkdiag(I,-I) and Mi = M^-1, then
% 
%                             Xe = [X           ,  M'*J-Y^-1*Mi ]
%                                  [J*M-Mi'*Y^-1, J-Mi'*Y^-1*Mi ]
%
%                - type = 7:  Xe = [X,I       ]
%                                  [I,(X-Y)^-1]
%
%                - type = 8:  Xe = [X  ,X-Y]
%                                  [X-Y,X-Y]
%
%                - type = 9:  If Y is positive definite, let Y=Yc'Yc>0
%                             (Choleski factorization) and Yci = Yc^-1, then
%
%                             Xe = [X        , X*Yci-Yc'   ]
%                                  [Yci'*X-Yc, Yci'*X*Yci-I]
%
%                             If Y is indefinite, let Y = M'*J*M, with
%                             J=blkdiag(I,-I) and Mi = M^-1, then
%
%                             Xe = [X        , X*Mi-M'*J]
%                                  [Mi'*X-J*M, Mi'*X*M-J]
%
%                - type = 10: Xe = [X       , X*Y^-1-I         ]
%                                  [Y^-1*X-I , Y^-1*X*Y^-1-Y^-1]
%
%                - type = 11: If X and Y are positive definite, let X =
%                             M'*M, Y = N'*N, (Choleski factorization) Mi =
%                             M^-1, Ni = N^-1 and u*s*v' = N*Mi-Ni'*M',
%                             then
%
%                             Xe = [Mi', 0        ]^-1[M, V*sqrt(s)]
%                                  [Ni', U*sqrt(s)]   [N, 0        ]
%
%                             If X or Y is indefinite, let X = M'*Jx*M, Y =
%                             N'*Jy*N, Jx = blkdiag(I,-I), Jy =
%                             blkdiag(I,-I), Mi = M^-1, Ni = N^-1 and 
%                             u*s*v' = (N*Mi - Jy*Ni'*M'*Jx). Then
%
%                             Xe = [Mi'   , 0        ]^-1[Jx*M, V*sqrt(s)]
%                                  [Jy*Ni', U*sqrt(s)]   [N   , 0        ]
%
%                - type = 12: If X is positive definite, let X = Xc'*Xc > 0
%                             (Choleski factorization) and Xci = Xc^-1,
%                             then 
%
%                             Xe = [X        , Xc'-Y*Xci   ]
%                                  [Xc-Xci'*Y, I-Xci'*Y*Xci]
%
%                             If X is indefinite, Let X = M'*J*M, with
%                             J = blkdiag(I,-I) and Mi = M^-1, then
%
%                             Xe = [X        , M'*J-Y*Mi]
%                                  [J*M-Mi'*Y, J-Mi'*Y*M]
%
%                - type = 13: Xe = [X11 ,X12 ,X13]
%                                  [X12',X22 ,X23]
%                                  [X13',X23',X13]
%
%                             where:
%
%                               - X13 = [M12;I]
%                               - X23 = [M22;0]
%                               - X33 = blkdiag(M22,M33)
%                               - M12 = X12 + Y12/Y22
%                               - M22 = X22-Y22^-1
%                               - M33 = (X11+Y11-Y12/Y22*Y12'- ... 
%                                                        M12/M22*M12')^-1
%
%                - type = 14: Xe = [X11 ,X12 ,X13]
%                                  [X12',X22 ,X23]
%                                  [X13',X23',X13]
%
%                             where:
%
%                               - X13 = [A,C]
%                               - X23 = [B,A']
%                               - X33 = [B,A';A,C]
%                               - A   = X12+Y12/Y22
%                               - B   = X22-Y22^-1
%                               - C   = X11+Y11-Y12/Y22*Y12'
%
%                - type = 13: Xe = [X11 ,X12 ,X13]
%                                  [X12',X22 ,X23]
%                                  [X13',X23',X13]
%
%                             where:
%
%                               - X13 = [A;B]
%                               - X23 = [B;C]
%                               - X33 = [E,F;F',G]
%                               - A   = X12*Y22+Y12
%                               - B   = Y11+X12*Y12'
%                               - C   = X22*Y22-I
%                               - C   = X22*Y12'
%                               - E   = Y22*X22*Y22-Y22
%                               - F   = (Y22*X22-I)*Y12'
%                               - G   = Y12*X22*Y12'-Y11
% -------------------------------------------------------------------------
nx                    = size(X,1);
switch type
    case 1
        Xe            = [X,eye(nx);eye(nx),fInv(X-fInv(Y))];
    case 2
        Yinv          = fInv(Y);
        Xe            = [X,X-Yinv;X-Yinv,X-Yinv];
    case 3
        [M,J,~,in]    = fEigfac(Y);
        if in(2) == 0 && in(3) == 0
            Yc        = chol(Y);
            Yci       = fInv(Yc);
            Yt        = X*Yc'-Yci;
            Yh        = Yc*X*Yc';
            Xe        = [X,Yt;Yt',0.5*fHe(Yh)-eye(nx)];
        elseif in(3) == 0
            Mi        = fInv(M);
            Yt        = X*M'-Mi*J;
            Yh        = M*X*M';
            Xe        = [X,Yt;Yt',0.5*fHe(Yh)-J];
        else
            error('The matrix Y should be nonsingular');
        end
    case 4
        Yt            = X*Y-eye(nx);
        Yh            = Y*X*Y;
        Xe            = [X,Yt;Yt',0.5*fHe(Yh)-Y];
    case 5
        [M,Jx,~,inx]  = fEigfac(X);
        [N,Jy,~,iny]  = fEigfac(Y);
        if inx(2) == 0 && inx(3) == 0 && iny(2) == 0 && iny(3) == 0
            M         = chol(X);
            Mi        = fInv(M);
            N         = chol(Y);
            Ni        = fInv(N);
            [U,S,V]   = svd(Ni'*Mi-N*M');
            Xe1       = [Mi',zeros(nx);N,U*sqrt(S)];
            Xe2       = [M,V*sqrt(S);Ni',zeros(nx)];
            Xe        = 0.5*fHe(Xe1\Xe2);
        elseif inx(3) == 0 && iny(3) == 0
            Mi        = fInv(M);
            Ni        = fInv(N);
            [U,S,V]   = svd(Ni'*Mi-Jy*N*M'*Jx);
            Xe1       = [Mi',zeros(nx);Jy*N,U*sqrt(S)];
            Xe2       = [Jx*M,V*sqrt(S);Ni',zeros(nx)];
            Xe        = 0.5*fHe(Xe1\Xe2);
        else
            error('The matrix Y should be nonsingular');
        end
    case 6
        [Mx,Jx,~,inx] = fEigfac(X);
        [My,Jy,~,iny] = fEigfac(Y);
        if inx(2) == 0 && inx(3) == 0 && iny(2) == 0 && iny(3) == 0
            Xc        = chol(X);
            Yc        = chol(Y);
            Xci       = fInv(Xc);
            Xt        = Xc'-Y\Xci;
            Xh        = Xci'/Yc;
            Xe        = [X,Xt;Xt',eye(nx)-Xh*Xh'];
        elseif inx(2) == 0 && iny(2) == 0
            Mxi       = fInv(Mx);
            Myi       = fInv(My);
            Xt        = Mx'*Jx-Y\Mxi;
            Xh        = Mxi'*Myi*Jy*Myi'*Mxi;
            Xe        = [X,Xt;Xt',Jx-Xh];
        else
            error('Both the matrices X and Y should be nonsingular');
        end
    case 7
        Xe            = [X,eye(nx);eye(nx),0.5*fHe(fInv(X-Y))];
    case 8
        Xe            = [X,X-Y;X-Y,X-Y];
    case 9
        [M,J,~,in]    = fEigfac(Y);
        if in(2) == 0 && in(3) == 0
            Yc        = chol(Y);
            Yci       = fInv(Yc);
            Yt        = X*Yci-Yc';
            Yh        = Yci'*(X)*Yci;
            Xe        = [X,Yt;Yt',0.5*fHe(Yh)-eye(nx)];
        elseif in(3) == 0
            Mi        = fInv(M);
            Yt        = X*Mi-M'*J;
            Yh        = Mi'*X*Mi;
            Xe        = [X,Yt;Yt',0.5*fHe(Yh)-J];
        else
            error('The matrix Y should be nonsingular');
        end
    case 10
        Yi            = fInv(Y);
        Yt            = X*Yi-eye(nx);
        Yh            = Yi*X*Yi;
        Xe            = [X,Yt;Yt',0.5*fHe(Yh)-Yi];
    case 11
        [M,Jx,~,inx]  = fEigfac(X);
        [N,Jy,~,iny]  = fEigfac(Y);
        if inx(2) == 0 && inx(3) == 0 && iny(2) == 0 && iny(3) == 0
            M         = chol(X);
            Mi        = fInv(M);
            N         = chol(Y);
            Ni        = fInv(N);
            [U,S,V]   = svd(N*Mi-Ni'*M');
            Xe1       = [Mi',zeros(nx);Ni',U*sqrt(S)];
            Xe2       = [M,V*sqrt(S);N,zeros(nx)];
            Xe        = 0.5*fHe(Xe1\Xe2);
        elseif inx(3) == 0 && iny(3) == 0
            Mi        = fInv(M);
            Ni        = fInv(N);
            [U,S,V]   = svd(N*Mi-Jy*Ni'*M'*Jx);
            Xe1       = [Mi',zeros(nx);Jy*Ni',U*sqrt(S)];
            Xe2       = [Jx*M,V*sqrt(S);N,zeros(nx)];
            Xe        = 0.5*fHe(Xe1\Xe2);
        else
            error('The matrix Y should be nonsingular');
        end
    case 12
        [M,J,~,in]    = fEigfac(X);
        if in(2) == 0 && in(3) == 0
            Xc        = chol(X);
            Xt        = Xc'-Y/Xc;
            Yh        = Xc'\Y/Xc;
            Xe        = [X,Xt;Xt',eye(nx)-0.5*fHe(Yh)];
        elseif in(3) == 0
            Xt        = M'*J-Y/M;
            Xh        = M'\Y/M;
            Xe        = [X,Xt;Xt',J-0.5*fHe(Xh)];
        else
            error('The matrix Y should be nonsingular');
        end
    case 13
        n1            = X.dimX(1);
        n2            = X.dimX(2);
        X11           = X.X(1:n1,1:n1);
        X12           = X.X(1:n1,n1+1:end);
        X22           = X.X(n1+1:end,n1+1:end);
        Y11           = Y.Y(1:n1,1:n1);
        Y12           = Y.Y(1:n1,n1+1:end);
        Y22           = Y.Y(n1+1:end,n1+1:end);
        Y22i          = fInv(Y22);
        XY22          = X22-Y22i;
        XY12          = X12+Y12*Y22i;
        X23i          = fInv(XY22);
        X13           = [XY12,eye(n1)];
        X23           = [XY22,zeros(n2,n1)];
        X33           = blkdiag(XY22,fInv(X11+Y11-Y12*(Y22i)*Y12'-XY12*(X23i)*XY12'));
        Xe            = [X11,X12,X13;X12',X22,X23;X13',X23',X33];
    case 14
        n1            = X.dimX(1);
        X11           = X.X(1:n1,1:n1);
        X12           = X.X(1:n1,n1+1:end);
        X22           = X.X(n1+1:end,n1+1:end);
        Y11           = Y.Y(1:n1,1:n1);
        Y12           = Y.Y(1:n1,n1+1:end);
        Y22           = Y.Y(n1+1:end,n1+1:end);
        Y22i          = fInv(Y22);
        A             = X12+Y12/Y22;
        B             = X22-Y22i;
        C             = X11+Y11-Y12*Y22i*Y12';
        X13           = [A,C];
        X23           = [B,A'];
        X33           = [B,A';A,C];
        Xe            = [X11,X12,X13;X12',X22,X23;X13',X23',X33];
    case 15
        n1            = X.dimX(1);
        n2            = X.dimX(2);
        X11           = X.X(1:n1,1:n1);
        X12           = X.X(1:n1,n1+1:end);
        X22           = X.X(n1+1:end,n1+1:end);
        Y11           = Y.Y(1:n1,1:n1);
        Y12           = Y.Y(1:n1,n1+1:end);
        Y22           = Y.Y(n1+1:end,n1+1:end);
        A             = X12*Y22+Y12;
        B             = Y11+X12*Y12';
        C             = X22*Y22-eye(n2);
        D             = X22*Y12';
        E             = 0.5*fHe(Y22*X22*Y22-Y22);
        F             = (Y22*X22-eye(n2))*Y12';
        G             = 0.5*fHe(Y12*X22*Y12')-Y11;
        X13           = [A,B];
        X23           = [C,D];
        X33           = [E,F;F',G];
        Xe            = [X11,X12,X13;X12',X22,X23;X13',X23',X33];    
end
end