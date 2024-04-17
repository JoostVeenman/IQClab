function Mex=fplantex(M,IO,options)
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
% Description: This function creates the extended plant Mex from M.
%
% Syntax:      Mex = fPlantex(prob,M,options)
%
% Usage:       This function creates extended plants as follows:
%
%              1.) Given the original plant "M", with realization
%
%                  [A  | B1  B2  ... Bn ]
%                  [---|--------------- ]
%                  [C1 | D11 D12 ... B1n]
%                  [C2 | D21 D22 ... B2n]
%                  [...| ... ... ... ...]
%                  [Cm | Dm1 Dm2 ... Bmn]
%
%                  the function creates the extended realization "Mex":
%
%                  [A  | B1  B2  ... Bn ]
%                  [---|----------------]
%                  [C1 | D11 D12 ... B1n]
%                  [ 0 |  I   0  ...  0 ]
%                  [C2 | D21 D22 ... B2n]
%                  [ 0 |  0   I  ...  0 ]
%                  [...| ... ... ... ...]
%                  [...| ... ... ... ...]
%                  [Cm | Dm1 Dm2 ... Bmn]
%                  [ 0 |  0   0  ...  I ]
%
%              2.) Or from the transposed plant "-M^*", with realization
%
%                  [-A^T  | -C1^T  -C2^T  ... -Cn^T ]
%                  [------|-------------------------]
%                  [-B1^T | -D11^T -D21^T ... -Bn1^T]
%                  [-B2^T | -D12^T -D22^T ... -Bn2^T]
%                  [ ...  |   ...    ...  ...   ... ]
%                  [-Bn^T | -D1n^T -D2n^T ... -Bnn^T]
%
%                  the function creates the extended realization "Mex":
%
%                  [-A^T  | -C1^T  -C2^T  ... -Cn^T ]
%                  [------|-------------------------]
%                  [  0   |   I       0   ...    0  ]
%                  [-B1^T | -D11^T -D21^T ... -Bn1^T]
%                  [  0   |   0       I   ...    0  ]
%                  [-B2^T | -D12^T -D22^T ... -Bn2^T]
%                  [ ...  |   ...    ...  ...   ... ]
%                  [  0   |   0       0   ...    I  ]
%                  [-Bn^T | -D1n^T -D2n^T ... -Bnn^T]
%
%
%             These augmented plants are compatible with the block
%             structure of the IQC-multiplier: 
% 
%             Pi = blkdiag(Pi_1,Pi_2, ... ,Pi_n).
%
%             The following inputs are to be provided:
%             1.) M:  The plant
%             2.) IO: The input/output data of the IQC-multipliers
%             3.) options (optional): Structure with
%                 - options.n  : With n = [n1,n2] this defines the start
%                                and stop criteria such that the code
%                                starts at IO(n1,:) and stops at IO(n2,:) 
%                 - options.pd : Use options.pd='primal_augm' for primal
%                                augmentation or options.pd = 'dual_augm',
%                                for dual augmentation (Default: options.pd
%                                = 'primal_augm')
%
% -------------------------------------------------------------------------
Mex                     = M;
if nargin == 3
    if isfield(options,'n')
        n1              = options.n(1);
        n2              = options.n(2);
    else
        n1              = 1;
        n2              = size(IO,1);
    end
    if isfield(options,'pd')
        options.pd      = options.pd;
    else
        options.pd      = 'primal_augm';
    end
elseif nargin < 3
    n1                  = 1;
    n2                  = size(IO,1);
    options.pd          = 'primal_augm';
end
switch options.pd
    case 'primal_augm'
        counter         = 0;
        for i = n1:n2
            if IO(i,3) == 0
                Mex     = Mex;
                counter = counter+IO(i,2);
            elseif IO(i,3)>0
                Mb      = [Mex.b(:,1:counter+IO(i,3)),zeros(size(Mex.b,1),IO(i,4)),Mex.b(:,counter+IO(i,3)+1:end)];
                Md      = [Mex.d(:,1:counter+IO(i,3)),zeros(size(Mex.d,1),IO(i,4)),Mex.d(:,counter+IO(i,3)+1:end)];
                counter = counter + IO(i,4)+1;
                Mex     = ss(Mex.a,Mb,Mex.c,Md,Mex.Ts);
            end
        end

        T1              = [];
        T2              = [];
        for i = n1:n2
            if IO(i,3) == 0
                T1      = blkdiag(T1,[eye(IO(i,1));zeros(IO(i,2),IO(i,1))]);
                T2      = blkdiag(T2,[zeros(IO(i,1),IO(i,2));eye(IO(i,2))]);
            elseif IO(i,3)>0
                T1      = blkdiag(T1,[eye(IO(i,1));zeros(IO(i,3)+IO(i,4),IO(i,1))]);
                T2      = blkdiag(T2,[zeros(IO(i,1),IO(i,3)+IO(i,4));eye(IO(i,3)+IO(i,4))]);
            end
        end
    case 'dual_augm'
        counter         = 0;
        for i = n1:n2
            if IO(i,3) == 0
                Mex     = Mex;
                counter = counter + IO(i,2);
            elseif IO(i,3) > 0
                Mb      = [Mex.b(:,1:counter+IO(i,3)),zeros(size(Mex.b,1),IO(i,4)),Mex.b(:,counter+IO(i,3)+1:end)];
                Md      = [Mex.d(:,1:counter+IO(i,3)),zeros(size(Mex.d,1),IO(i,4)),Mex.d(:,counter+IO(i,3)+1:end)];
                counter = counter + IO(i,4)+1;
                Mex     = ss(Mex.a,Mb,Mex.c,Md,Mex.Ts);
            end
        end

        T1              = [];
        T2              = [];
        for i = n1:n2
            if IO(i,3) == 0
                T1      = blkdiag(T1,[zeros(IO(i,2),IO(i,1));eye(IO(i,1))]);
                T2      = blkdiag(T2,[eye(IO(i,2));zeros(IO(i,1),IO(i,2))]);
            elseif IO(i,3)>0
                T1      = blkdiag(T1,[zeros(IO(i,3)+IO(i,4),IO(i,1));eye(IO(i,1))]);
                T2      = blkdiag(T2,[eye(IO(i,3)+IO(i,4));zeros(IO(i,1),IO(i,3)+IO(i,4))]);
            end
        end
end

T11                     = blkdiag(T1,eye(size(Mex.d,1)-size(T1,2)));
T22                     = blkdiag(T2,zeros(size(T11*Mex.d,1)-size(T2,1),size(T11*Mex.d,2)-size(T2,2)));
Mex                     = ss(Mex.a,Mex.b,T11*Mex.c,T11*Mex.d+T22,Mex.Ts);
end