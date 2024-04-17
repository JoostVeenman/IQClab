function Phi=fBasis(l,pl,nr,type,Ts)

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
% Description: This function creates state-space realizations for various
%              basis function for the IQC multipliers:
%
%              type = 1:
%
%              col(1, (iw+pl)/(iw-pl), (iw+pl)^2/(iw-pl)^2, ... ,
%                                                   (iw+pl)^l/(iw-pl)^l)
%
%              type = 2:
%
%              col(1, (iw)/(iw-pl)^(l-1), (iw)^2/(iw-pl)^(l-1), ... ,
%                                                 (iw)^(l-1)/(iw-pl)^(l-1))
%
%              type = 3:
%
%              col(1, 1/(iw-pl), 1/(iw-pl)^2, ... , 1/(iw-pl)^l)
%
%              type = 4: 
%
%              (-1, -1/(iw+pl), 1/(iw+pl)^2, ... ,
%                                    (-1)^(l+1)/(iw-pl)^l) = -(type 3)^*
%
%              type = 5
%
%              this is the discrete time version of type = 1
%
%              type = 7
%
%              this is the discrete time version of type = 3
%
%              type = 101
%
%              This is a discrete time FIR basis function
%
%              col(1, 1/z, 1/z^2, ... , 1/z^l)
%
% Syntax:      Phi = fBasis(l,pl,nr,type)
%
% Usage:       Phi = fBasis(l,pl,nr,type) creates a state space realization
%              of a basis function of typeX with length l > 1, pole
%              location pl < 0, and number of repetitions nr > 0.
% -------------------------------------------------------------------------
if nargin < 5
    Ts                  = 0;
end
if Ts > 0 || Ts == -1
    if type == 1
        type            = 5;
    elseif type == 3
        type            = 7;
    end
end
switch type
     case 1
        A1              = diag(ones(l-1,1));A2=[];
        for i           = 1:l-1
            A2(i,1:l-1) = [zeros(1,i),ones(1,l-i-1)];    
        end
        pl              = -pl;
        Aphi            = kron(-pl*A1-2*pl*A2,eye(nr));
        Bphi            = kron(-sqrt(2*pl)*ones(l-1,1),eye(nr));
        Cphi            = kron([zeros(1,l-1);sqrt(2*pl)*rot90(A1+A2,3)],eye(nr));
        Dphi            = kron(ones(l,1),eye(nr));
        Phi             = ss(Aphi,Bphi,Cphi,Dphi,Ts);
    case 2
        if l == 1
            Phi         = ss([],Ts);
        else
            Aphi        = diag(pl*ones(l-2,1),0)+diag(ones(l-3,1),-1);
            Bphi        = eye(l-2,1);
            T           = zeros(0,l-1); 
            p           = 1;
            for i       = 1:l-1
                T       = [T;zeros(1,l-1-i) p];
                p       = conv(p,[1 pl]);
            end
            Cphi        = T(:,2:l-1);
            Dphi        = T(:,1);
            Phi         = ss(Aphi,Bphi,Cphi,Dphi,Ts);
        end        
    case 3
        Aphi            = kron(diag(pl*ones(l-1,1),0)+diag(ones(l-2,1),-1),eye(nr));
        Bphi            = kron(eye(l-1,1),eye(nr));
        Cphi            = kron([zeros(1,l-1);eye(l-1)],eye(nr));
        Dphi            = kron(eye(l,1),eye(nr));
        Phi             = ss(Aphi,Bphi,Cphi,Dphi,Ts);
    case 4
        Aphi            = kron(diag(pl*ones(l-1,1),0)+diag(ones(l-2,1),-1),eye(nr));
        Bphi            = kron(eye(l-1,1),eye(nr));
        Cphi            = kron([zeros(1,l-1);eye(l-1)],eye(nr));
        Dphi            = kron(eye(l,1),eye(nr));
        Phi             = ss(-Aphi,-Bphi,-Cphi,-Dphi,Ts);
    case 5 % discrete time version of type 1
        if pl == 0
           error('For this basis function it is not possible to have a pole at zero'); 
        end
        if Ts == 0
           error('For discrete time basis functions the sampling time should be strictly larger than 0');
        end
        A1              = diag(ones(l-1,1));A2=[];
        for i           = 1:l-1
            A2(i,1:l-1) = [zeros(1,i),ones(1,l-i-1)];    
        end
        Aphi            = kron(pl*A1+2*pl*A2,eye(nr));
        Bphi            = kron(-sign(pl)*sqrt(2*abs(pl))*ones(l-1,1),eye(nr));
        Cphi            = kron([zeros(1,l-1);-sqrt(2*abs(pl))*rot90(A1+A2,3)],eye(nr));
        Dphi            = kron(ones(l,1),eye(nr));
        Phi             = ss(Aphi,Bphi,Cphi,Dphi,Ts);
    case 6 % discrete time version of type 2
        
    case 7 % discrete time version of type 3
        if Ts == 0
           error('For discrete time basis functions the sampling time should be strictly larger than 0');
        end
        Aphi            = kron(pl*eye(l-1)+[zeros(1,l-1);eye(l-2,l-1)],eye(nr));
        Bphi            = kron(eye(l-1,1),eye(nr));
        Cphi            = kron([zeros(1,l-1);eye(l-1)],eye(nr));
        Dphi            = kron([eye(l,1)],eye(nr));
        Phi             = ss(Aphi,Bphi,Cphi,Dphi,Ts);
    case 101 % discrete time FIR basis function
        if Ts == 0
           error('For discrete time basis functions the sampling time should be strictly larger than 0');
        end
        if l == 1
            Phi         = ss([],[],[],kron(1,eye(nr)),Ts);
        else
            Aphi        = kron([zeros(l-1,1),eye(l-1,l-2)],eye(nr));
            Bphi        = kron([zeros(l-2,1);1],eye(nr));
            Cphi        = kron([zeros(1,l-1);rot90(eye(l-1),1)],eye(nr));
            Dphi        = kron([1;zeros(l-1,1)],eye(nr));
            Phi         = ss(Aphi,Bphi,Cphi,Dphi,Ts);
        end
end
end