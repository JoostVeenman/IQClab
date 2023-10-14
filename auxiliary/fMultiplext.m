function [Pe,Qe,Se,Re] = fMultiplext(Pp,Pd,type)
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
% Description: This function constructs the extended Full-block Multiplier
%              Pe in different fashions.
%
% Syntax:      [Pe,Qe,Se,Re] = fMultiplext(Pp,Pd,npos,nneg,type)  
%
% Usage:       Let:
%                1.) Pp = [Qp ,Sp] with Qp > 0
%                         [Sp',Rp]
%
%                2.) Pp = [Qd ,Sd] with Rd < 0
%                         [Sd',Rd]
%
%                3.) npos denote the number of positive eigenvalues of Pp
%                4.) nneg denote the number of negative eigenvalues of Pp
%                5.) J = (I_npos,0_nneg)' and Jt = (0_npos,I_nneg)'
%
%              Then one can construct the extended multiplier 
%
%                Pe = [ Qe  , Se ]
%                     [ Qe' , Re ]
%
%              as follows:
%
%              Use type = 1-3 to obtain Pe^-1 = [ Qd  , * , Sd , * ]
%                                               [ *   , * , *  , * ]
%                                               [ Sd' , * , Rd , * ]
%                                               [ *   , * , *  , * ]
%
%              Use type = 4-6 to obtain Pe^-1 = [ Qdi  , * , Sdi , * ]
%                                               [ *    , * , *   , * ]
%                                               [ Sdi' , * , Rdi , * ]
%                                               [ *    , * , *   , * ]
%
%              where Pdi = [Qdi,Sdi;Sdi',Rdi] = Pd^-1.
%
%              These are defined as follows:
%
%                - type=1: - Define Pi = (Pp-Pd^-1)^-1
%                          - Let T1 denote the collection of eigenvectors
%                            corresponding to the positive eigenvalues of
%                            eig(Pi-J/(J'*Pp*J)*J') 
%                          - Let T2 denote the collection of eigenvectors
%                            corresponding to the negative eigenvalues of
%                            eig(Pi-Jt/(Jt'*Pp*Jt)*Jt') 
%                          - Let T = [T1,T2]
%                          - Then Pe = [ Qe , Se ] = [ Pp , T       ]
%                                      [ Se', Re ]   [ T' , T'*Pi*T ]
%
%                - type=2: - Define Pi = (Pp-Pd^-1)
%                          - Let T1 denote the collection of eigenvectors
%                            corresponding to the positive eigenvalues of
%                            eig(Pi-PiJ/(J'*Pp*J)*J'Pi) 
%                          - Let T2 denote the collection of eigenvectors
%                            corresponding to the negative eigenvalues of
%                            eig(Pi-PiJt/(Jt'*Pp*Jt)*Jt'Pi) 
%                          - Let T = [T1,T2]
%                          - Then Pe = [ Qe , Se ] = [ Pp   , PiT     ]
%                                      [ Se', Re ]   [ T'Pi , T'*Pi*T ]
%
%                - type=3: - Factorize Pp = Mp'*Jp*Mp, Jp = blkdiag(I,-I)
%                            and define Mpi = Mp^-1 
%                          - Factorize Pd = Md'*Jd*Md, Jd = blkdiag(I,-I)
%                            and define Mdi=Md^-1 
%                          - Define u*s*v' = Mdi'*Mpi - Jd*Md*Mp'*Jp, Sd =
%                            u*sqrt(s), Sp = v*sqrt(s) 
%                          - Let T1 and T2 denote the collection of
%                            eigenvectors corresponding to the positive and
%                            negative eigenvalues of 
%
%                            eig(-(Sd\Jd*Md)*Mp'*Sp - ...
%                                            Sp'*Mp*J'/(J'*Pp*J)*J'*Mp'*Sp) 
%                            and
%
%                            eig(-(Sd\Jd*Md)*Mp'*Sp - ...
%                                         Sp'*Mp*Jt/(Jt'*Pp*Jt)*Jt'*Mp'*Sp)
%
%                            respectively.
%                          - Let T = [T1,T2]
%                          - Then Pe = [ Qe , Se ] = ...
%                                      [ Se', Re ]
%
%                                  =  [ Mpi  , 0       ]^-1[ Jp*Mp , Sp*T ]
%                                     [ Jd*Md, Sd/(T') ]   [ Mdi   , 0    ]
%
%                - type=4: - Define Pi = (Pp-Pd)^-1
%                          - Let T1 be the collection of eigenvectors
%                            corresponding to the positive eigenvalues of
%                            eig(Pi-J/(J'*Pp*J)*J') 
%                          - Let T2 be the collection of eigenvectors
%                            corresponding to the negative eigenvalues of
%                            eig(Pi-Jt/(Jt'*Pp*Jt)*Jt') 
%                          - Let T = [T1,T2]
%                          - Then Pe = [ Qe , Se ] = [ Pp , T       ]
%                                      [ Se', Re ]   [ T' , T'*Pi*T ]
%
%                - type=5: - Define Pi =(Pp-Pd)^-1
%                          - Let T1 be the collection of eigenvectors
%                            corresponding to the positive eigenvalues of
%                            eig(Pi-PiJ/(J'*Pp*J)*J'Pi) 
%                          - Let T2 be the collection of eigenvectors
%                            corresponding to the negative eigenvalues of
%                            eig(Pi-PiJt/(Jt'*Pp*Jt)*Jt'Pi) 
%                          - Let T = [T1,T2]
%                          - Then Pe = [ Qe , Se ] = [ Pp   , PiT     ]
%                                      [ Se', Re ]   [ T'Pi , T'*Pi*T ]
%
%                - type=6: - Factorize Pp = Mp'*Jp*Mp, Jp = blkdiag(I,-I)
%                            and define Mpi = Mp^-1 
%                          - Factorize Pd = Md'*Jd*Md, Jd = blkdiag(I,-I)
%                            and define Mdi = Md^-1 
%                          - Define u*s*v' = Md*Mpi-Jd*Mdi'*Mp'*Jp, Sd =
%                            u*sqrt(s), Sp=v*sqrt(s) 
%                          - Let T1 and T2 denote the collection of
%                            eigenvectors corresponding to the positive and
%                            negative eigenvalues of 
%
%                            eig(-(Sd\Jd*Md)*Mp'*Sp - ... 
%                                            Sp'*Mp*J'/(J'*Pp*J)*J'*Mp'*Sp)
%
%                            and
%
%                            eig(-(Sd\Jd*Md)*Mp'*Sp - ... 
%                                         Sp'*Mp*Jt/(Jt'*Pp*Jt)*Jt'*Mp'*Sp)
%
%                            respectively.
%                          - Let T = [T1,T2]
%                          - Then Pe = [ Qe , Se ] = ...
%                                      [ Se', Re ] 
%
%                                 = [ Mpi    , 0       ]^-1[ Jp*Mp , Sp*T ]
%                                   [ Jd*Mdi', Sd/(T') ]   [ Md    , 0    ]
%
%--------------------------------------------------------------------------
npos            = sum(eig(Pp) > 0);
nneg            = sum(eig(Pp) < 0);
J               = fJ(npos,nneg);
Jt              = fJt(nneg,npos);
Je              = kron(eye(2),J);
Jte             = kron(eye(2),Jt);
switch type
    case 1
        Pi      = fInv(Pp-fInv(Pd));
        [V1,D1] = eig(Pi-J/(J'*Pp*J)*J');
        [V2,D2] = eig(Pi-Jt/(Jt'*Pp*Jt)*Jt');
        T1      = V1*orth(double(diag(diag(D1)>0)));
        T2      = V2*orth(double(diag(diag(D2)<0)));
        T       = [T1,T2];
        Qe      = Je'*[Pp,T;T',T'*Pi*T]*Je;
        Se      = Je'*[Pp,T;T',T'*Pi*T]*Jte;
        Re      = Jte'*[Pp,T;T',T'*Pi*T]*Jte;        
    case 2
        Pi      = Pp-fInv(Pd);
        [V1,D1] = eig(Pi-Pi*J/(J'*Pp*J)*J'*Pi);
        [V2,D2] = eig(Pi-Pi*Jt/(Jt'*Pp*Jt)*Jt'*Pi);
        T1      = V1*orth(double(diag(diag(D1)>0)));
        T2      = V2*orth(double(diag(diag(D2)<0)));
        T       = [T1,T2];
        Qe      = Je'*[Pp,Pi*T;T'*Pi,T'*Pi*T]*Je;
        Se      = Je'*[Pp,Pi*T;T'*Pi,T'*Pi*T]*Jte;
        Re      = Jte'*[Pp,Pi*T;T'*Pi,T'*Pi*T]*Jte;
    case 3
        ncDe    = npos;
        nrDe    = nneg;
        n       = size(Pp,1);
        neg     = nrDe+1:nrDe+ncDe;
        pos     = setdiff(1:n,neg);
        [Mp,Jp] = fJfac(Pp,neg);
        Mpi     = fInv(Mp);    
        [Md,Jd] = fJfac(-Pd,pos);
        Jd      = -Jd;
        Mdi     = fInv(Md);
        [u,s,v] = svd(Mdi'*Mpi-Jd*Md*Mp'*Jp);
        Sd      = u*sqrt(s);
        Sp      = v*sqrt(s);
        P11     = Pp(pos,pos);
        P21     = Sp'*Mp(:,pos);
        H1      = 0.5*fHe(-(Sd\Jd*Md)*Mp'*Sp-P21/P11*P21');
        [~,j,T] = fEigfac(H1);
        Tpos    = T*null(double(j<0));
        P11     = Pp(neg,neg);
        P21     = Sp'*Mp(:,neg);
        H2      = 0.5*fHe(-(Sd\Jd*Md)*Mp'*Sp-P21/P11*P21');
        [~,j,T] = fEigfac(H2);
        Tneg    = T*null(double(j>0));
        T       = [Tpos,Tneg];
        nrDc    = size(Tneg,2);
        ncDc    = size(Tpos,2);
        Pe1     = [Mpi',zeros(size(Mpi,2),size(Sd,2));Jd*Md,Sd/(T')];
        Pe2     = [Jp*Mp,Sp*T;Mdi',zeros(size(Mdi,2),size(Sp,2))];
        Pe      = 0.5*fHe(Pe1\Pe2);
        Qind    = [1:ncDe,ncDe+nrDe+1:ncDe+nrDe+ncDc];
        Rind    = [ncDe+1:ncDe+nrDe,ncDe+nrDe+ncDc+1:ncDe+nrDe+ncDc+nrDc];
        Qe      = Pe(Qind,Qind);
        Se      = Pe(Qind,Rind);
        Re      = Pe(Rind,Rind);

%         ncDe    = npos;
%         nrDe    = nneg;
%         n       = size(Pp,1);
%         neg     = nrDe+1:nrDe+ncDe;
%         pos     = setdiff(1:n,neg);
%         Pp11    = Pp(pos,pos);
%         Pp12    = Pp(pos,neg);
%         Pp22    = Pp(neg,neg);
%         J       = fJ(npos,nneg);
%         Jt      = [-Pp11\Pp12;eye(nneg)];
% 
%         To      = [zeros(npos,nneg),eye(npos);eye(nneg),zeros(nneg,npos)];
%         [Mp,Jp] = fJfac(-To'*Pp*To,neg);
%         Mp      = To'*Mp*To;
%         Jp      = -To'*Jp*To;
%         Mpi     = fInv(Mp);
%         
%         To      = [zeros(nneg,npos),eye(nneg);eye(npos),zeros(npos,nneg)];
%         [Md,Jd] = fJfac(To'*Pd*To,pos);
%         Md      = To'*Md*To;
%         Jd      = To'*Jd*To;
%         Mdi     = fInv(Md);
%         
%         [u,s,v] = svd(Mdi'*Mpi-Jd*Md*Mp'*Jp);
%         Sd      = u*sqrt(s);
%         Sp      = v*sqrt(s);
% 
%         P11     = J'*Pp*J;
%         P21     = Sp'*Mp*J;
%         H1      = 0.5*fHe(-(Sd\Jd*Md)*Mp'*Sp-P21/P11*P21');
%         [~,j,T] = fEigfac(H1);
%         Tpos    = T*null(double(j<0));
% 
%         P11     = Jt'*Pp*Jt;
%         P21     = Sp'*Mp*Jt;
%         H2      = 0.5*fHe(-(Sd\Jd*Md)*Mp'*Sp-P21/P11*P21');
%         [~,j,T] = fEigfac(H2);
%         Tneg    = T*null(double(j>0));
% 
%         T       = [Tpos,Tneg];
%         nrDc    = size(Tneg,2);
%         ncDc    = size(Tpos,2);
%         Pe1     = [Mpi',zeros(size(Mpi,2),size(Sd,2));Jd*Md,Sd/(T')];
%         Pe2     = [Jp*Mp,Sp*T;Mdi',zeros(size(Mdi,2),size(Sp,2))];
%         Qind    = [1:nrDe,nrDe+ncDe+1:nrDe+ncDe+nrDc];
%         Rind    = [nrDe+1:nrDe+ncDe,nrDe+ncDe+nrDc+1:nrDe+ncDe+nrDc+ncDc];
%         Pe      = 0.5*fHe(Pe1\Pe2);
%         Qe      = Pe(Qind,Qind);
%         Se      = Pe(Qind,Rind);
%         Re      = Pe(Rind,Rind);
%         Pet = Pe^-1;
%         Pd(1:7,1:7)-Pet(1:7,1:7)
%         Pd(1:7,8:14)-Pet(1:7,15:21)
%         Pd(8:14,8:14)-Pet(15:21,15:21)
    case 4
        Pi      = 0.5*fHe(fInv(Pp-Pd));
        [V1,D1] = eig(Pi-J/(J'*Pp*J)*J');
        [V2,D2] = eig(Pi-Jt/(Jt'*Pp*Jt)*Jt');
        T1      = V1*orth(double(diag(diag(D1)>0)));
        T2      = V2*orth(double(diag(diag(D2)<0)));
        T       = [T1,T2];
        Qe      = Je'*[Pp,T;T',T'/Pi*T]*Je;
        Se      = Je'*[Pp,T;T',T'/Pi*T]*Jte;
        Re      = Jte'*[Pp,T;T',T'/Pi*T]*Jte;
    case 5
        Pi      = Pp-Pd;
        [V1,D1] = eig(Pi-Pi*J/(J'*Pp*J)*J'*Pi);
        [V2,D2] = eig(Pi-Pi*Jt/(Jt'*Pp*Jt)*Jt'*Pi);
        T1      = V1*orth(double(diag(diag(D1)>0)));
        T2      = V2*orth(double(diag(diag(D2)<0)));
        T       = [T1,T2];
        Qe      = Je'*[Pp,Pi*T;T'*Pi,T'*Pi*T]*Je;
        Se      = Je'*[Pp,Pi*T;T'*Pi,T'*Pi*T]*Jte;
        Re      = Jte'*[Pp,Pi*T;T'*Pi,T'*Pi*T]*Jte;
    case 6
        ncDe    = npos;
        nrDe    = nneg;
        n       = size(Pp,1);
        neg     = nrDe+1:nrDe+ncDe;
        pos     = setdiff(1:n,neg);        
        [Mp,Jp] = fJfac(Pp,neg);
        Mpi     = fInv(Mp);
        [Md,Jd] = fJfac(-Pd,pos);
        Jd      = -Jd;
        Mdi     = fInv(Md);
        [u,s,v] = svd(Md*Mpi-Jd*Mdi'*Mp'*Jp);
        Sd      = u*sqrt(s);
        Sp      = v*sqrt(s);
        P11     = Pp(pos,pos);
        P21     = Sp'*Mp(:,pos);
        H1      = 0.5*fHe(-Sd\Jd*Mdi'*Mp'*Sp-P21/P11*P21');
        [~,j,T] = fEigfac(H1);
        Tpos    = T*null(double(j<0));
        P11     = Pp(neg,neg);
        P21     = Sp'*Mp(:,neg);
        H2      = 0.5*fHe(-Sd\Jd*Mdi'*Mp'*Sp-P21/P11*P21');
        [~,j,T] = fEigfac(H2);
        Tneg    = T*null(double(j>0));
        T       = [Tpos,Tneg];
        nrDc    = size(Tneg,2);
        ncDc    = size(Tpos,2);
        Pe1     =[Mpi',zeros(size(Mpi,2),size(Sd,2));Jd*Mdi',Sd/(T')];
        Pe2     = [Jp*Mp,Sp*T;Md,zeros(size(Mdi,2),size(Sp,2))];
        Qind    = [1:ncDe,ncDe+nrDe+1:ncDe+nrDe+ncDc];
        Rind    = [ncDe+1:ncDe+nrDe,ncDe+nrDe+ncDc+1:ncDe+nrDe+ncDc+nrDc];
        Pe      = 0.5*fHe(Pe1\Pe2);
        Qe      = Pe(Qind,Qind);
        Se      = Pe(Qind,Rind);
        Re      = Pe(Rind,Rind);
end
Pe              = [Qe,Se;Se',Re];
end