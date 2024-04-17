function iqcprob = iqclmi(iqcprob,P,J,A,B,C)

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
% Date:        24-11-2019
%
% -------------------------------------------------------------------------
%
% Description: fLMI creates new LMI constraints and adds them to the
%              considered set of LMI.
%
% Syntax:      iqcprob = iqclmi(iqcprob,P,J)
%              iqcprob = iqclmi(iqcprob,P,J,A)
%              iqcprob = iqclmi(iqcprob,P,J,A,B)
%              iqcprob = iqclmi(iqcprob,P,J,A,B,C)
%
% Usage:       This function creates LMIs of the form
%
%              newlmi = C^T*(J*B'*P*B >= J*A)*C
%
%              As input one can specify:
%
%                1.) "iqcprob" is a structure containing information about
%                    the LMI problem (should be defined by the iqc-class
%                    "iqc")
%                2.) P = symmetric matrix variable (should be defined by
%                    the iqc-class "iqc"
%                3.) J = Inertia of LMI (pos|neg definite = 1|-1)
%                4.) A = Constant symmetric matrix
%                5.) B = Constant matrix
%                6.) C = Constant matrix
%
%              As output one obtains:
%              - iqcprob with an additional LMI constraint
%
% -------------------------------------------------------------------------

if strcmp(iqcprob.Parser,'Yalmip')
    Q = P.var;
    if nargin == 3
        Aeps = iqcprob.eps*eye(size(Q,1));
        iqcprob.lmi = iqcprob.lmi + [J*Q >= Aeps];         % J*P >= eps*I
    elseif nargin == 4
        Aeps = iqcprob.eps*eye(size(A,1));
        iqcprob.lmi = iqcprob.lmi + [J*Q >= J*A+Aeps];      % J*P >= J*A+eps*I
    elseif nargin == 5
        Aeps = iqcprob.eps*eye(size(A,1));
        iqcprob.lmi = iqcprob.lmi + [J*B'*Q*B >= J*A+Aeps]; % J*B'*P*B >= J*A+eps*I
    elseif nargin == 6
        T = B*C;
        D = C'*A*C;
        Deps = iqcprob.eps*eye(size(D,1));
        iqcprob.lmi = iqcprob.lmi + [J*T'*Q*T >= J*D+Deps]; % C^T*(J*B'*P*B)*C >= C^T*(J*A+eps*I)*C
    end
elseif strcmp(iqcprob.Parser,'LMIlab')
    ni = newlmi;
    Q  = lmivar(3,P.svar);
    if nargin == 3
        Aeps = J*iqcprob.eps*eye(size(Q,1));
        lmiterm([-J*ni 1 1 Q],1,1);                       % J*P >= eps*I
        lmiterm([J*ni 1 1 0],Aeps);
    elseif nargin == 4
        Aeps = A+J*iqcprob.eps*eye(size(A,1));
        lmiterm([-J*ni 1 1 Q],1,1);                       % J*P >= J*A+eps*I
        lmiterm([J*ni 1 1 0],Aeps);
    elseif nargin == 5
        Aeps = A+J*iqcprob.eps*eye(size(A,1));
        lmiterm([-J*ni 1 1 Q],B',B);                      % J*B'*P*B >= J*A+eps*I
        lmiterm([J*ni 1 1 0],Aeps);
    elseif nargin == 6
        T = B*C;
        D = C'*A*C;
        Deps = D+J*iqcprob.eps*eye(size(D,1));
        lmiterm([-J*ni 1 1 Q],T',T);                      % C^T*(J*B'*P*B)*C >= C^T*(J*A+eps*I)*C
        lmiterm([J*ni 1 1 0],Deps);
    end
end
end