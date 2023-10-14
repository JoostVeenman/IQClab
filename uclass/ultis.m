classdef ultis < handle & matlab.mixin.SetGetExactNames

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
% Date:        30-10-2019
% 
% -------------------------------------------------------------------------
%
% Description: Class definition for LTI static diagonally repeated
%              parametric uncertainties of the form:
%
%                      [ delta_1 * I_nr1,     ,        0        ]
%              delta = [                , ... ,                 ]
%                      [        0       ,     , delta_N * I_nrN ]
%
%              where:
%              - delta_i \in S \subset R, i \in {1,...,N}
%              - S is star convex [0,1]S\subset S
%              - nri is the number of repetitions of delta_i
%
% Syntax:      delta = ultis('name')
%              delta = ultis('name',varargin)
%
% Usage:       "delta = ultis('name')" defines an LTI parametric
%              uncertainty on the interval [-1,1], which is repeated once.
%
%              For "delta = ultis('name',varargin)", the varargin inputs
%              come in pairs and can be defined as: 
%
%              delta = ultis('name','prop1','value1','prop2','value2',...)
%
%              Alternatively, the properties can be set by defining the
%              structure:
%
%              ultisOpt.prop1 = value1
%              ultisOpt.prop2 = value2
%                       ...              
%              ultisOpt.propN = valueN
% 
%              delta = ultis('name',ultisOpt)
%
%              Finally, the properties can be set and retrieved by using:
%
%              set(delta,'propX','valueX') and get(delta,'propX').
%
%              The properties that can be specified are:
%            
%               1.) 'NumberOfRepetitions' specifies the number of
%                   repetitions of delta_i (default =  1).
%
%                   Since it is possible to define multiple diagonally
%                   repeated uncertainties one should specify:
%
%                   delta = ultv_rb('name','NumberOfRepetitions',[nr1;nr2;...])
%
%               2.) 'Bounds' specifies the domain on which the uncertainty
%                   is defined (default = [-1;1]). It is assumed that the
%                   value zero is contained in the uncertainty set:
%
%                          Bounds = [-delta1, ... ,-deltan] 
%                                   [+delta1, ... ,+deltan]
%
%               3.) Alternatively, instead of using the option 'Bounds',
%                   one can specify the option 'Polytope'
%
%                                     [delta1p1, ... ,deltamp1]
%                      b.) Polytope = [   ...  , ... ,   ...  ]
%                                     [delta1pn, ... ,deltampn]
%
%               4.) 'InputChannel' specifies which input-channels of the
%                   uncertain plant are affected by delta. For each deltai,
%                   the channels should be specified as: 
%
%                         row_i = [Chnl_deltai_x,...,Chnl_deltai_y]
%
%                    Here the order of the channels is not relevant, while
%                    the row-length equals the number of repetitions of
%                    delta_i. 
%
%                   'InputChannel' should then be specified as a cell
%
%                           InputChannel = {row_1,...,row_N}
%
%               5.) 'OutputChannel' specifies which output-channels of the
%                   uncertain plant are affected by delta. For each deltai,
%                   the channels should be specified as: 
%
%                         row_i = [Chnl_deltai_x,...,Chnl_deltai_y]
%
%                    Here the order of the channels is not relevant, while
%                    the row-length equals the number of repetitions of
%                    delta_i.
%
%                   'OutputChannel' should then be specified as a cell
%
%                           OutputChannel = {row_1,...,row_N}
%
%
%              In addition, various properties of the corresponding
%              IQC-multiplier can be specified: 
%
%               6.) 'BasisFunctionType' specifies the type of basis
%                   function (default = 1).
%
%               7.) 'Length' specifies the length of the basis function
%                   (default =  1). In case of multiple diagonally repeated
%                   uncertainties, once can specify one common length as
%
%                              delta = ultis('Length',l) 
%
%                   or a different length for each delta_i
%
%                           delta = ultis('Length',[l1;l2;l3;...]).
%                                  
%               8.) 'PoleLocation' specifies the pole location of the basis
%                   function (default = -1). In case of multiple diagonally
%                   repeated uncertainties, once cam specify one common
%                   pole-location 
%
%                             delta = ultis('PoleLocation',pl)
%
%                   or a different pole-location for each delta_i
%
%                   delta = ultis('PoleLocation',[pl1;pl2;pl3;...])
%
%               9.) 'SampleTime' specifies the sampling time of the plant
%                   (default = 0, is continuous time)
%
%              10.) 'RelaxationType' specifies the multiplier relaxation
%                   type, options are (default = 'DG'):                     
%
%                    a.) 'DG'   DG-scalings (only applicable for 'Bounds')
%                    b.) 'PR'   Positive real (equivalent to DG) (only
%                               applicable for 'Bounds')
%                    c.) 'CH'   Convex hull relaxation
%                    d.) 'PC'   Partial Convexity
%                    e.) 'ZP'   Zeroth order Polya relaxation
%
%              11.) 'TerminalCost' defines the terminal cost constraint, Z,
%                   which is used in invariance analysis. Options are
%
%                    a.) 'off' (default)
%                    a.) 'on'
%
%              12.) 'PrimalDual' specifies whether the multiplier should be
%                   a primal/dual parametrization (default = 'Primal').
%                   
%                   Note: For a standard IQC-analysis, all multipliers are
%                   Primal ones.
% -------------------------------------------------------------------------

properties (SetAccess = public)
    % Uncertainty name
    Name string = 'delta';
    
    % Basis function type
    BasisFunctionType double {mustBeReal,mustBeFinite} = 1;

    % Length of the muliplier basis function
    Length double {mustBeReal,mustBeFinite} = 1;

    % Pole location of the multiplier basis function
    PoleLocation double {mustBeReal,mustBeFinite} = -1;
    
    % Sampling time of the multiplier basis function
    SampleTime double {mustBeReal,mustBeFinite} = 0;

    % Number of repetitions of the uncertain parameter
    NumberOfRepetitions double {mustBeReal,mustBeFinite} = 1;

    % Relaxation type of the uncertainty (DG scaling (DG), Convex hull (CH) Partial Convexity (PC), First order Polya (ZP), Positive Real (PR))
    RelaxationType string {mustBeMember(RelaxationType,{'DG','CH','PC','ZP','PR'})} = 'DG';

    % Terminal cost certificate Z
    TerminalCost string {mustBeMember(TerminalCost,{'on','off'})} = 'off';

    % Bounds of the uncertainty
    Bounds double {mustBeReal} = [-1;1];
    
    % Bounds of the uncertainty region defined by a polytope
    Polytope double {mustBeReal} = [];

    % Primal or Dual multiplier
    PrimalDual string {mustBeMember(PrimalDual,{'Primal','Dual'})} = 'Primal';
    
    % InputChannel
    InputChannel = {};
    
    % OutputChannel
    OutputChannel = {};
end
methods
    function obj = ultis(varargin)
        if nargin == 0
            error('Error: To utilize "ultis" one should specify at least one input argument (see "help ultis" for further details)');
        elseif nargin == 1
            if ischar(varargin{1})
                obj.Name = varargin{1};
            else
                error('Error: The name of the uncertainty should be defined by a string (see "help ultis" for further details)');
            end
        elseif nargin > 1
            if nargin == 2
                if ischar(varargin{1})
                    obj.Name = varargin{1};
                else
                    error('Error: The name of the uncertainty should be defined by a string (see "help ultis" for further details)');
                end
                if isstruct(varargin{2})
                else
                    error('Error: The second input argument should be defined by a structure (see "help ultis" for further details)');
                end
                if isfield(varargin{2},'BasisFunctionType')
                    a = 1:1:10;
                    if isreal(varargin{2}.BasisFunctionType) && ismatrix(varargin{2}.BasisFunctionType) && min(ismember(varargin{2}.BasisFunctionType,a)) > 0
                        obj.BasisFunctionType = varargin{2}.BasisFunctionType;
                    else
                        error('Error: The type of basis function should be defined as Natural numbers in N1');
                    end
                end
                if isfield(varargin{2},'Length')
                    a = 1:1:10;
                    if isreal(varargin{2}.Length) && ismatrix(varargin{2}.Length) && min(ismember(varargin{2}.Length,a)) > 0
                        obj.Length = varargin{2}.Length;
                    else
                        error('Error: The length of a basis function should be defined as a natural number in N1 (see "help ultis" for further details).');
                    end
                end
                if isfield(varargin{2},'PoleLocation')
                    if isreal(varargin{2}.PoleLocation) && ismatrix(varargin{2}.PoleLocation)
                        obj.PoleLocation = varargin{2}.PoleLocation;
                    else
                        error('Error: The pole location of a basis function should be defined as a real rational number. For continuous and distrete time systems, it is not allowed to have poles on the imaginary axis and unit circle respectively (see "help ultis" for further details).');
                    end
                end
                if isfield(varargin{2},'SampleTime')
                    if isreal(varargin{2}.SampleTime) && ismatrix(varargin{2}.SampleTime)
                        obj.SampleTime = varargin{2}.SampleTime;
                    else
                        error('Error: The Sampling time of a basis function should be defined as a nonnegative real rational number (see "help ultis" for further details).');
                    end
                end
                if isfield(varargin{2},'NumberOfRepetitions')
                    if isreal(varargin{2}.NumberOfRepetitions) && ismatrix(varargin{2}.NumberOfRepetitions) && min(varargin{2}.NumberOfRepetitions > 0) > 0 
                        obj.NumberOfRepetitions = varargin{2}.NumberOfRepetitions;
                    else
                        error('Error: The number of repetitions of the uncertainty should be defined as a vector with positive natural numbers in N1 (see "help ultis" for further details).');
                    end
                end
                if isfield(varargin{2},'RelaxationType')
                    obj.RelaxationType = varargin{2}.RelaxationType;
                end
                if isfield(varargin{2},'TerminalCost')
                    obj.TerminalCost   = varargin{2}.TerminalCost;
                end
                if isfield(varargin{2},'Bounds')
                    if isreal(varargin{2}.Bounds) && ismatrix(varargin{2}.Bounds)
                        obj.Bounds     = varargin{2}.Bounds;
                        obj.Polytope   = [];
                    else
                        error('Error: The bounds of the uncertainty should be defined as a real matrix (see "help ultis" for further details).');
                    end
                end
                if isfield(varargin{2},'Polytope')
                    if isreal(varargin{2}.Polytope) && ismatrix(varargin{2}.Polytope)
                        obj.Polytope   = varargin{2}.Polytope;
                        obj.Bounds     = [];
                    else
                        error('Error: The Polytope should be defined as a real matrix (see "help ultis" for further details).');
                    end
                end
                if isfield(varargin{2},'PrimalDual')
                    obj.PrimalDual     = varargin{2}.PrimalDual;
                end
                if isfield(varargin{2},'InputChannel')
                    obj.InputChannel   = varargin{2}.InputChannel;
                end
                if isfield(varargin{2},'OutputChannel')
                    obj.OutputChannel  = varargin{2}.OutputChannel;
                end
            elseif nargin > 2
                if ischar(varargin{1})
                    obj.Name = varargin{1};
                else
                    error('Error: The name of the uncertainty should be defined by a string (see "help ultis" for further details)');
                end 
                n = length(varargin);
                if mod(n-1,2) ~= 0
                    error('Error: The input arguments should be defined in pairs')
                else
                    j = linspace(3,n,n/2);
                end
                for i = 1:length(j)
                    prop = varargin{j(i)-1};
                    switch prop
                        case 'BasisFunctionType'
                            a = 1:1:10;
                            if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)}) && min(ismember(varargin{j(i)},a)) > 0
                                obj.BasisFunctionType = varargin{j(i)};
                            else
                                error('Error: The type of basis function should be defined as natural numbers in N1 (see "help ultis" for further details)');
                            end
                        case 'Length'
                            a = 1:1:10;
                            if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)}) && min(ismember(varargin{j(i)},a)) > 0
                                obj.Length = varargin{j(i)};
                            else
                                error('Error: The length of a basis function should be defined as a natural number in N1 (see "help ultis" for further details).');
                            end
                        case 'PoleLocation'
                            if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)})
                                obj.PoleLocation = varargin{j(i)};
                            else
                                error('Error: The pole location of a basis function should be defined as a real rational number. For continuous and distrete time systems, it is not allowed to have poles on the imaginary axis and unit circle respectively (see "help ultis" for further details).');
                            end
                        case 'SampleTime'
                            if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)})
                                obj.SampleTime = varargin{j(i)};
                            else
                                error('Error: The sampling time of a basis function should be defined as a nonnegative real rational number (see "help ultis" for further details).');
                            end
                        case 'NumberOfRepetitions'
                            if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)}) && min(varargin{j(i)} > 0) > 0 
                                obj.NumberOfRepetitions = varargin{j(i)};
                            else
                                error('Error: The number of repetitions of the uncertainty should be defined as a vector with positive natural numbers in N1 (see "help ultis" for further details).');
                            end
                        case 'RelaxationType'
                            obj.RelaxationType = varargin{j(i)};
                        case 'Bounds'
                            if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)})
                                obj.Bounds = varargin{j(i)};
                                obj.Polytope = [];
                            else
                                error('Error: The uncertainty domain should be defined as a real matrix (see "help ultis" for further details).');
                            end
                        case 'Polytope'
                            if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)})
                                obj.Polytope = varargin{j(i)};
                                obj.Bounds = [];
                            else
                                error('Error: The uncertainty domain should be defined as a real matrix (see "help ultis" for further details).');
                            end
                        case 'TerminalCost'
                            obj.TerminalCost = varargin{j(i)};
                        case 'PrimalDual'
                            obj.PrimalDual = varargin{j(i)};
                        case 'InputChannel'
                            obj.InputChannel = varargin{j(i)};
                        case 'OutputChannel'
                            obj.OutputChannel = varargin{j(i)};
                    end
                end    
            end
        end
    end
    function prob = iqcltis_p(Delta,prob)
    % ---------------------------------------------------------------------
    % This function generates two possible IQC-multipliers for the class
    % ultis:
    %
    %   1.) Dynamic DG-scalings
    %   2.) Dynamic full-block multipliers
    %
    % These are precisely implemented as discussed in Section 5.3 of: 
    %
    % [1] J. Veenman, C.W. Scherer, H. Koroglu, Robust stability and
    %     performance analysis based on integral quadratic constraints,
    %     European Journal of Control, Vol. 31, pp. 1-32, 2016.
    % ---------------------------------------------------------------------
    % Dynamic DG-scalings
    % ---------------------------------------------------------------------
    %
    % DG scalings are defined as
    %
    % Pi_DG = [D  , G]
    %         [G^*,-D]
    %
    % where:
    % - D = blkdiag(D_1,...,D_n)\in\RLinf
    % - G = blkdiag(G_1,...,G_n)\in\RLinf
    %
    % We have
    %
    % [I    ]^T[D  , G][I    ] >= 0 for all omega\in[0,\infty]
    % [Delta]  [G^*,-D][Delta]
    %
    % and thus
    %
    % D + Delta^T G + G^* Delta - Delta^T D Delta\geq0
    %
    % and thus (I-Delta^TDelta) D\geq0 if:
    % - D(i\omega)>0 
    % - G(i\omega)+G(i\omega)^*=0
    %
    % This function generates the LMI constraints
    %
    %      [ 0, X, 0   ][  I  ,  0   ]
    % [x]^T[ X, 0, 0   ][ Apsi, Bpsi ] >0 and P12 = -P12^T.             (1)
    %      [ 0, 0, P11 ][ Cpsi, Dpsi ]
    %
    % Which are the KYP certificates for:
    % 1.) D = Psi^*P11 Psi > 0
    % 2.) G = Psi^*P12 Psi, P12 = -P12^T
    % 
    % with the basis realization Psi = (Apsi,Bpsi,Cpsi,Dpsi)\in\RHinf (see
    % the function "fBasis").
    %
    % ---------------------------------------------------------------------
    % Dynamic PR (Positive real) multipliers
    % ---------------------------------------------------------------------
    % Positive real multipliers are defined as
    %
    % Pi_PR = [x][0   Phi^*][bet  I, -I]
    %         [x][Phi,     ][-alp I,  I]
    %
    % with alp < 0 < beta and for delta\in[alp,bet]
    %
    % where Phi\in RL_\infty and satisfies the constraint Phi^* + Phi > 0
    %
    % With the realization 
    %
    %  Psi = [psi1,0   ][bnd I, I] = [Apsi | Bpsi ]
    %        [0,   psi2][bnd I, I]   [-----|------]
    %                                [Cpsi | Dpsi ] 
    %
    % This function generates the LMI constraint
    %
    %      [ 0, Y, 0 ][  I  ,  0   ]
    % [x]^T[ Y, 0, 0 ][ Apsi, Bpsi ]E >0
    %      [ 0, 0, P ][ Cpsi, Dpsi ]
    %
    % where P = [0,  P12^T] and E = [I]
    %           [P12,0    ]         [0]
    %
    % ---------------------------------------------------------------------
    % Dynamic Full-block multipliers
    % ---------------------------------------------------------------------
    % With
    %
    % Pi_FB = [Q11  , Q12]\in\RLinf
    %         [Q12^*, Q22]
    %
    % we have:
    %
    % Fij(Q) = [I      ]^T[Q11  , Q12][I      ]\geq0,
    %          [Delta^i]  [Q12^*, Q22][Delta^j]
    %
    % where i,j\in 1,...n and Delta^n denote the vertices of the convex
    % hull of Delta (Property "Polytope"), defined as:
    %
    % Polytope = [delta1_p1,...,deltam_p1]
    %            [  ...    ,...,   ...   ]
    %            [delta1_pn,...,deltam_pn]
    %
    % We thus have
    %
    % Fij(Q) = Q11 + (Delta^i)^T Q12 + Q12^* Delta^j + ...
    %            (Delta^i)^T Q22 Delta^j >= 0, i,j\in 1,...n.
    %
    % In order to relax this semi-infinite dimensional matrix inequality,
    % this function considers two optional relaxation schemes:
    %
    % 1.)  Partial Convexity (PC): Fij(Q)\geq0 if
    %      a.) Q22 < 0
    %      b.) Fjj(Q) > 0
    %
    %      This function generates the LMI constraints
    %
    %             [ 0, X, 0   ][  I  ,  0   ]
    %        [x]^T[ X, 0, 0   ][ Apsi, Bpsi ] < 0                       (2)
    %             [ 0, 0, P22 ][ Cpsi, Dpsi ]
    %
    %      and
    %
    %             [ 0, X, 0   ][  I  ,  0   ]
    %        [x]^T[ X, 0, 0   ][ Apsi, Bpsi ]>0, j=1,...,n              (3)
    %             [ 0, 0, Pjj ][ Cpsi, Dpsi ]
    %
    %      where Pjj = P11+He{(Delta^j)^TP12}+(Delta^j)^T P22 Delta^j
    %
    %      based on the parameterizations
    %      a.) Q22    = Psi^* P22 Psi < 0
    %      b.) Fjj(Q) = Psi^* Pjj Psi > 0, j = 1,...,n
    %
    %      with the basis realization (see the function "fBasis)"
    %
    %            [ Apsi | Bpsi ]
    %      Psi = [------|------]\in\RHinf,
    %            [ Cpsi | Dpsi ]
    %
    %      as well as the factorization
    %
    %      Pi_PC = [Psi, 0 ]^*[P11  ,P12][Psi, 0 ]
    %              [0  ,Psi]  [P12^T,P22][0  ,Psi]
    %
    % II.) Zeroth order Polya relaxation (ZP): equation (3) is satisfied if
    %      a.) Fjj(Q) >= 0,          j = 1,...,n
    %      b.) Fij(Q) + Fji(Q) >= 0, j = 1,...,n, i = j + 1,...,n.
    %
    %      Then this function generates the LMI constraints
    %
    %             [ 0, X, 0   ][  I  ,  0   ]
    %        [x]^T[ X, 0, 0   ][ Apsi, Bpsi 0 > 0, j = 1,...,n,         (4)
    %             [ 0, 0, Pjj ][ Cpsi, Dpsi ]
    %
    %      and
    %
    %            [ 0, X, 0       ][  I  ,  0   ]
    %       [x]^T[ X, 0, 0       ][ Apsi, Bpsi ]>0, j = 1,...,n,        (5)
    %            [ 0, 0, Pij+Pji ][ Cpsi, Dpsi ]    i = j + 1,...,n.
    %
    %      where 
    %      a.) Pjj = P11+He{(Delta^j)^TP12}+(Delta^j)^T P22 Delta^j
    %      b.) Pij = P11+(Delta^i)^TP12+P12^TDelta^j+(Delta^i)^T P22 Delta^j,
    %
    %      based on the parameterizations
    %      a.) Fjj(Q) = Psi^* Pjj Psi < 0,     j = 1,...,n
    %      b.) Fij(Q) = Psi^*(Pij+Pji)Psi > 0, j = 1,...,n, i = j + 1,...,n
    % 
    %      with the basis realization Psi = (Apsi,Bpsi,Cpsi,Dpsi)\in\RHinf
    %      (see the function "fBasis").
    %
    %      Note: The first order Polya relaxation is not more
    %      conservative than the convex Hull relaxation.
    %
    %----------------------------------------------------------------------
    
    % Define outer-factor (basis function) of IQC-multiplier
    bft                            = get(Delta,'BasisFunctionType');
    nr                             = get(Delta,'NumberOfRepetitions');
    l                              = get(Delta,'Length');
    pl                             = get(Delta,'PoleLocation');
    Ts                             = get(Delta,'SampleTime');
    tc                             = get(Delta,'TerminalCost');
    lnr                            = length(nr);
    ll                             = length(l);
    lpl                            = length(pl);
    if ll < lnr
        l                          = l(1)*ones(lnr,1);
    end
    if lpl < lnr
        pl                         = pl(1)*ones(lnr,1);
    end
    Phi                            = ss([],[],[],[],Ts);
    if Ts == 0 && sum(pl == 0) > 0
        disp('Warning: In continuous time it is not allowed to have poles on the imaginary axis. The "PoleLocation" was set to its nominal value (i.e. -1).');
        pl                         = -ones(size(pl));
    elseif Ts > 0 && sum(abs(pl) == 1) > 0
        disp('Warning: In discrete time it is not allowed to have poles on the unit circle. The "PoleLocation" was set to its nominal value (i.e. 0).');
        pl                         = zeros(size(pl));
    elseif Ts == -1 && sum(abs(pl) == 1) > 0
        disp('Warning: In discrete time it is not allowed to have poles on the unit circle. The "PoleLocation" was set to its nominal value (i.e. 0).');
        pl                         = zeros(size(pl));
    end
    for i = 1:length(nr)
        Phi                        = fAugss(Phi,fBasis(l(i),pl(i),nr(i),bft,Ts),1);
    end
    prob                           = fPsi(prob,Phi,Phi);

    % Define multiplier variables plus lmi constraints
    rt                             = get(Delta,'RelaxationType');
    
    % Define KYP certificate structure (for continuous or discrete time analysis)
    if Ts == 0
        kyp                        = [0,1;1,0];
    else
        kyp                        = [-1,0;0,1];
    end
    switch rt
        case 'DG'
            if isempty(Delta.Bounds)
                error('Error: To apply DG-scalings one should specify the property "Bounds" of the uncertainty.');
            end
            P11                    = iqcvar(prob,[0,0],'symmetric');
            P12                    = iqcvar(prob,[0,0],'skew');
            sc11                   = [];
            for i = 1:length(nr)
                P11                = blkdiag(P11,iqcvar(prob,[l(i)*nr(i),l(i)*nr(i)],'symmetric'));
                P12                = blkdiag(P12,iqcvar(prob,[l(i)*nr(i),l(i)*nr(i)],'skew'));
                sc11               = blkdiag(sc11,Delta.Bounds(2,i)*eye(l(i)*nr(i)));
            end
            P                      = [P11,P12;P12.',-P11];
            sc                     = blkdiag(sc11,eye(size(sc11,1)));
            prob                   = fP(prob,P);
            prob                   = fsc(prob,sc);
            switch tc
                case 'on'
                    Z              = iqcvar(prob,[0,0],'symmetric');
                    Tz11           = [];
                    Tz22           = [];
                    for i = 1:length(nr)
                        if (l(i)-1)*nr(i) ~= 0
                            Z      = blkdiag(Z,iqcvar(prob,[(l(i)-1)*nr(i),(l(i)-1)*nr(i)],'symmetric'));
                            Tz11   = blkdiag(Tz11,Delta.Bounds(2,i)*eye((l(i)-1)*nr(i)));
                            Tz22   = blkdiag(Tz22,eye((l(i)-1)*nr(i)));
                        end
                    end
                    prob           = fZ(prob,blkdiag(Z,-Z),blkdiag(Tz11,Tz22));
                    A              = fss2m(Phi);
                    V              = blkdiag(kron(kyp,Z),P11);
                    prob           = iqclmi(prob,V,1,0,A);           
                case 'off'
                    prob           = iqclmi(prob,P11,1);        
            end
        case 'PR'
            if isempty(Delta.Bounds)
                error('Error: To apply PR-multipliers one should specify the property "Bounds" of the uncertainty.');
            end
            R                      = iqcvar(prob,[0,0],'symmetric');
            P12                    = iqcvar(prob,[0,0],'full');
            sc11                   = [];
            sc12                   = [];
            sc21                   = [];
            sc22                   = [];
            for i = 1:length(nr)
                P12                = blkdiag(P12,iqcvar(prob,[l(i)*nr(i),l(i)*nr(i)],'full'));
                sc11               = blkdiag(sc11,Delta.Bounds(2,i)*eye(l(i)*nr(i)));
                sc12               = blkdiag(sc12,-eye(l(i)*nr(i)));
                sc21               = blkdiag(sc21,-Delta.Bounds(1,i)*eye(l(i)*nr(i)));
                sc22               = blkdiag(sc22,eye(l(i)*nr(i)));
                if (l(i)-1)*nr(i) ~= 0
                    R              = blkdiag(R,iqcvar(prob,2*[(l(i)-1)*nr(i),(l(i)-1)*nr(i)],'symmetric'));
                end
            end
            P                      = oblkdiag(P12.');
            sc                     = [sc11,sc12;sc21,sc22];
            prob                   = fP(prob,P);
            prob                   = fsc(prob,sc);
            
            nPhi2                  = size(Phi,2);
            B                      = fss2m(sc*blkdiag(Phi,Phi)*[eye(nPhi2);zeros(nPhi2)]);
            V                      = blkdiag(kron(kyp,R),P);
            prob                   = iqclmi(prob,V,1,0,B);
            switch tc
                case 'on'
                    Z12            = iqcvar(prob,[0,0],'symmetric');
                    for i = 1:length(nr)
                        if (l(i)-1)*nr(i) ~= 0
                            Z12    = blkdiag(Z12,iqcvar(prob,[(l(i)-1)*nr(i),(l(i)-1)*nr(i)],'symmetric'));
                        end
                    end
                    T              = [eye(R.Dim(1));eye(R.Dim(1))];
                    Z              = [-Z12,Z12;Z12,Z12];
                    prob           = fZ(prob,Z,eye(Z.Dim(1)));
                    W              = blkdiag(R,-Z);
                    if W.Dim(1) > 0 && W.Dim(2) > 0
                        prob       = iqclmi(prob,Z12,-1);
                        prob       = iqclmi(prob,W,-1,0,T);
                    end
            end
        case {'CH','PC','ZP'}
            P11                    = iqcvar(prob,[size(Phi.c,1),size(Phi.c,1)],'symmetric');
            P12                    = iqcvar(prob,[size(Phi.c,1),size(Phi.c,1)],'full');
            P22                    = iqcvar(prob,[size(Phi.c,1),size(Phi.c,1)],'symmetric');
            P                      = [P11,P12;P12.',P22];
            sc                     = blkdiag(eye(size(Phi.c,1)),eye(size(Phi.c,1)));
            prob                   = fP(prob,P);
            prob                   = fsc(prob,sc);
            switch tc
                case 'on'
                    Z              = iqcvar(prob,2*[size(Phi.a,1),size(Phi.a,1)],'symmetric');
                    prob           = fZ(prob,Z,eye(Z.Dim(1)));
                    B              = fss2m(blkdiag(Phi,Phi)*[zeros(size(Phi,2));eye(size(Phi,2))]);
                    W              = blkdiag(kron(kyp,Z),P);
                    prob           = iqclmi(prob,W,-1,0,B);
            end
            if ~isempty(Delta.Bounds)
                La                 = polydec(pvec('box',Delta.Bounds'))';
            elseif ~isempty(Delta.Polytope)
                La                 = Delta.Polytope;
            end
            for i = 1:size(La,1)
                A                  = [];
                for j = 1:length(nr)
                    A              = blkdiag(A,La(i,j)*eye(l(j)*nr(j)));
                end
                A                  = [eye(size(A,1));A];
                prob               = iqclmi(prob,P,1,0,A);
            end
            switch rt
                case 'CH' % Convex Hull Relaxation
                    prob           = iqclmi(prob,P22,-1);
                case 'PC'
                    sLa1           = size(La,1);
                    if ll-sum(l)/length(l) < 0
                       error('Error: This relaxation can only be applied for "Length" is 1.') 
                    end
                    if sLa1 ~= 2^lnr
                        error('Error: For the relaxation option "PC" the individual uncertainties should be defined on the interval ci*[-1,1].');
                    end
                    for i = 1:length(nr)
                        A          = fJ(l(i)*nr(i),l(i)*sum(nr(i+1:end)),l(i)*sum(nr(1:i-1)));
                        prob       = iqclmi(prob,P22,-1,0,A);
                    end                    
                case 'ZP' % First Order Polya Relaxation
                    for i = 1:length(La)
                        for j = i+1:length(La)
                            A1     = [];
                            A2     = [];
                            for k = 1:length(nr)
                                A1 = blkdiag(A1,La(i,k)*eye(l(k)*nr(k)));
                                A2 = blkdiag(A2,La(j,k)*eye(l(k)*nr(k)));
                            end
                            prob   = iqclmi(prob,oblkdiag(P),1,0,[eye(size(A1,2));A1;eye(size(A2,2));A2]);
                        end
                    end
            end
    end
    IO                             = [sum(nr),sum(nr),0,0,size(Phi.c,1),size(Phi.c,1),size(Phi.c,1),size(Phi.c,1),size(Phi.a,1),size(Phi.a,1)];
    prob                           = fIO(prob,IO);
    end
    function prob = iqcltis_d(Delta,prob)
        error('Error: Dual IQC-multipliers are not supported for this class of uncertainties.');
    end
    function display(obj)
        a1 = [];b1 = [];
        for i = 1:length(obj.OutputChannel)
            if i == 1
                a2 = ['[',num2str(obj.InputChannel{i}),']'];
                b2 = ['[',num2str(obj.OutputChannel{i}),']'];
            else
                a2 = [', [',num2str(obj.InputChannel{i}),']'];
                b2 = [', [',num2str(obj.OutputChannel{i}),']'];
            end
                a1 = [a1,a2];
                b1 = [b1,b2];
        end
        disp('-------------------------------------------------------------');
        disp(' LTI static diagonally repeated parametric uncertainty block');
        disp('-------------------------------------------------------------');
        disp(' Uncertainty details:');
        disp(['  - Name:                           ',num2str(obj.Name)]);
        disp(['  - Number of uncertain parameters: ',num2str(length(obj.NumberOfRepetitions))]);
        disp(['  - Number of repetitions:          ',num2str(fDim(obj.NumberOfRepetitions))]);
        disp(['  - Input channels plant:           ',a1]);
        disp(['  - Output channels plant:          ',b1]);
        if isempty(obj.Bounds)
            disp('  - Convex Hull generator points:  ');
            disp('           [delta1^1, ... ,deltaN^1]');
            disp('           [   ...  , ... ,   ...  ]');
            disp('           [delta1^M, ... ,deltaN^M]');
            disp('La =');
            disp(' ')
            disp(obj.Polytope);
        else
            disp('  - Bounds:');
            disp(' ')
            disp(obj.Bounds);
        end
        disp(' Multiplier details:');
        disp(['  - Basis function type:            ',num2str(fDim(obj.BasisFunctionType))]);
        disp(['  - Basis function length:          ',num2str(fDim(obj.Length))]);
        disp(['  - Basis function pole location:   ',num2str(fDim(obj.PoleLocation))]);
        disp(['  - Relaxation type:                ',num2str(obj.RelaxationType)]);
        disp(['  - Multiplier type:                ',num2str(obj.PrimalDual)]);
    end
end
end
function y = fDim(x)
    [a1,a2] = size(x);
    if a1 > a2
        y = x';
    else
        y = x;
    end
end