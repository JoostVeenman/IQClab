function Psi = fAugss(Psi,Phi,nr)
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
% Date:        17-04-2012
% 
% -------------------------------------------------------------------------
%
% Description: Augment the state space object Phi to the state space object
%              Psi (This can be done nr times).
%
% Syntax:      Psi = fAugss(Psi,Phi,iter)
%
% Usage:       Provide Psi and Phi as inputs as well as the number nr that
%              Phi should be repeated. The function generates the following
%              augmented state space as output:
%
%                    [ Apsi            | Bpsi            ]
%                    [     Aphi        |     Bphi        ]
%                    [         ...     |         ...     ]
%                    [            Aphi |            Bphi ]
%              Psi = [-----------------------------------]
%                    [ Cpsi            | Dpsi            ]
%                    [     Cphi        |     Dphi        ]
%                    [         ...     |         ...     ]
%                    [            Cphi |            Dphi ]
%
% -------------------------------------------------------------------------
for i=1:nr
    Ts   = Phi.Ts;
    Apsi = blkdiag(Psi.a,Phi.a);
    Bpsi = blkdiag(Psi.b,Phi.b);
    Cpsi = blkdiag(Psi.c,Phi.c);
    Dpsi = blkdiag(Psi.d,Phi.d);
    Psi  = ss(Apsi,Bpsi,Cpsi,Dpsi,Ts);
end
end