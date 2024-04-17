function X = fInvproj(A,B,P)
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
% Date:        20-04-2011
% 
% -------------------------------------------------------------------------
%
% Description: Under the hypothesis that Aperp'*P*Aperp < 0 and
%              Bperp'*P*Bperp < 0, where Aperp and Bperp denote arbitrary
%              matrices whose columns form a basis of ker(A) and ker(B)
%              respectively, this function returns the unstructured matrix
%              X such that A'*X*B + B'*X'*A + P < 0 
%
% Syntax:      X = fInvproj(A,B,P)
%
% Usage:       Let S=[S1,S2,S3,S4] be a nonsingular matrix such that
%                - the columns of S3 span ker(A)\cap ker(B)
%                - the columns of [S1,S3] span ker(A)
%                - the columns of [S2,S3] span ker(B)
% 
%              Then one can consider the inequality:
%
%              S'*(A'*X*B+B'*X'*A+P)*S<0
%
%              It follows that:
%                                                          [0  , 0, 0, 0  ]
%              (AS)'*X*(BS) = [0,A2,0,A4]'*X*[B1,0,0,B4] = [Z21, 0, 0, Z24]
%                                                          [0  , 0, 0, 0  ]
%                                                          [Z41, 0, 0, Z44]
%              and 
%                                       
%              S'*(A'*X*B+B'*X*A+P)*S = 
%
%                 [Q11      , Q12 + Z21', Q13, Q14 + Z14'      ]
%               = [Q21 + Z21, Q22       , Q23, Q24 + Z24       ] < 0
%                 [Q31      , Q32       , Q33, Q34             ]
%                 [Q41 + Z41, Q42 + Z24', Q43, Q44 + Z44 + Z44']
%
%              This allows to directly obtain X.
%
% -------------------------------------------------------------------------
[S1,S2,S3,S4]=fIntkerext(A,B);

Q11         = S1'*P*S1;
Q12         = S1'*P*S2;
Q13         = S1'*P*S3;
Q14         = S1'*P*S4;
Q22         = S2'*P*S2;
Q23         = S2'*P*S3;
Q24         = S2'*P*S4;
Q33         = S3'*P*S3;
Q34         = S3'*P*S4;
Q44         = S4'*P*S4;

% check if: [Q11,Q13]<0 and  [Q22,Q23]<0
%           [Q31,Q33]        [Q32,Q33]

I1          = [Q11,Q13;Q13',Q33];
I2          = [Q22,Q23;Q23',Q33];
if max(eig(I1)) > 0 || max(eig(I2)) > 0
    error('Check1: The main LMI is not satisfied');
end

T33         = chol(-Q33);
T33i        = fInv(T33);
Q33i        = -T33i*T33i';
Z21         = Q23*Q33i*Q13'-Q12';
Z41         = Q34'*Q33i*Q13'-Q14';
Z24         = Q23*Q33i*Q34-Q24;
Q44t        = Q44-Q34'*Q33i*Q34;
Z44         = -0.5*Q44t-1e-4*norm(I2)*eye(size(Q44,1));

L           = [A*S2,A*S4]';
R           = [B*S1,B*S4];
X           = L\[Z21,Z24;Z41,Z44]/R;

if max(eig(fHe(A'*X*B) + P)) > 0
    alp     = 2;
    i       = 1;
    while max(eig(fHe(A'*X*B) + P)) > 0 && i < 2
        Z44 = -0.5*Q44t-0.01*alp*eye(size(Q44,1));
        L   = [A*S2,A*S4]';
        R   = [B*S1,B*S4];
        X   = L\[Z21,Z24;Z41,Z44]/R;
        i   = i+1;
        alp = 2*alp;
    end
end
        
% check if A'*X*B + B'*X*A + P < 0
if max(eig(fHe(A'*X*B)+P)) > 0
    error('Check2: The main LMI is not satisfied');
end
end