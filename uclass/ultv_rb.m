classdef ultv_rb < handle & matlab.mixin.SetGetExactNames

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
% Date:        06-10-2019
% 
% -------------------------------------------------------------------------
%
% Description: Class definition for LTV rate bounded diagonally repeated
%              parametric uncertainties of the form:
%
%                      [ delta_1 * I_nr1,     ,        0        ]
%              Delta = [                , ... ,                 ]
%                      [        0       ,     , delta_N * I_nrN ]
%
%              where:
%              - delta_i\in S
%              - S = {delta_i(.)\in C^1[[0,\infty),R): ...
%                              (delta(t),deltadot(t))\in La for all t\geq0}
%              - 0\in S
%              - La is a compact polytope:
%                       La = co{delta_i^1, ... ,delta_i^M} = ...
%                                 = {sum_{a=1}^M b_a*zeta^a: ...
%                               b_a\leq0, sum_{a=1}^M b_a = 1}
%                                 
%                with zeta^j = (delta_i^j,deltadot_i^j), j\in{1,...,M},
%                as generator points and with 0\in La.
%              - nri is the number of repetitions of delta_i.
%
% Syntax:      delta = ultv_rb('name')
%              delta = ultv_rb('name',varargin)
%
% Usage:       "delta = ultv_rb('name')" defines an LTV parametric
%              uncertainty on the interval [-1,1], which is repeated once
%              and which can vary arbitrarily fast.
%
%              For "delta = ultv_rb('name',varargin)", the varargin inputs
%              come in pairs and can be defined as: 
%
%              delta = ultv_rb('name','prop1','value1','prop2','value2',..)
%
%              Alternatively, the properties can be set by defining the
%              structure:
%
%              ultv_rbOpt.prop1 = value1
%              ultv_rbOpt.prop2 = value2
%                       ...              
%              ultv_rbOpt.propN = valueN
% 
%              delta = ultv_rb('name',ultv_rbOpt)
%
%              Finally, the properties can be set and retrieved by using:
%
%              set(delta,'propX','valueX') and get(delta,'propX')
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
%                          Bounds = [-delta1, ... ,-deltaN] 
%                                   [+delta1, ... ,+deltaN]
%
%               3.) 'RateBounds' specifies the rate-bounds for delta_i
%                   (default = [-1;1]). These should be specified in the
%                   same fashion as 'Bounds':
%
%                      RateBounds = [-delta1dot, ... ,-deltaNdot] 
%                                   [+delta1dot, ... ,+deltaNdot]
%
%               4.) Alternatively, instead of using the option 'Bounds/
%                   RateBounds', one can specify the option 'Polytope':
%
%                                   [d1p1,d1dotp1, ... ,dNp1,dNdotp1]
%                        Polytope = [... ,  ...  , ... ,... ,  ...  ]
%                                   [d1pM,d1dotpM, ... ,dNpM,dNdotpM]
%
%               5.) 'InputChannel' specifies which input-channels of the
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
%               6.) 'OutputChannel' specifies which output-channels of the
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
%               7.) 'BasisFunctionType' specifies the type of basis
%                   function (default = 1).
%
%               8.) 'Length' specifies the length of the basis function
%                   (default =  1). In case of multiple diagonally repeated
%                   uncertainties, once can specify one common length as
%
%                              delta = ultis('Length',l) 
%
%                   or a different length for each delta_i
%
%                           delta = ultv_rb('Length',[l1;l2;l3;...]).
%                                  
%               9.) 'PoleLocation' specifies the pole location of the basis
%                   function (default = -1). In case of multiple diagonally
%                   repeated uncertainties, once cam specify one common
%                   pole-location 
%
%                             delta = ultv_rb('PoleLocation',pl)
%
%                   or a different pole-location for each delta_i
%
%                   delta = ultv_rb('PoleLocation',[pl1;pl2;pl3;...])
%
%              10.) 'SampleTime' specifies the sampling time of the plant
%                   (default = 0, is continuous time)
%
%              11.) 'RelaxationType' specifies the multiplier relaxation
%                   type, options are (default = 'DG'):                     
%
%                     a.) 'DG'  DG-scalings (only applicable for 'Bounds')
%                     b.) 'CH'  Convex Hull
%                     c.) 'PC'  Partial Convexity
%                     d.) 'ZP'  Zeroth order Polya relaxation
%
%              12.) 'TerminalCost' defines the terminal cost constraint, Z,
%                   which is used in invariance analysis. Options are
%
%                    a.) 'off' (default)
%                    a.) 'on'
%
%              13.) 'PrimalDual' specifies whether the multiplier should be
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

    % Relaxation type of the uncertainty (DG scaling (DG), Zeroth order Polya (ZP), Partial Convexity (PC))
    RelaxationType string {mustBeMember(RelaxationType,{'DG','CH','PC','ZP'})} = 'DG';

    % Bounds of the uncertainty
    Bounds double {mustBeReal} = [-1;1];
    
    % Rate bounds of the uncertainty
    RateBounds double {mustBeReal} = [-1;1];
    
    % Bounds of the uncertainty region defined by a polytope
    Polytope double {mustBeReal} = [];

    % Terminal cost certificate Z
    TerminalCost string {mustBeMember(TerminalCost,{'on','off'})} = 'off';
    
    % Primal or Dual multiplier
    PrimalDual string {mustBeMember(PrimalDual,{'Primal','Dual'})} = 'Primal';
    
    % InputChannel
    InputChannel = {};
    
    % OutputChannel
    OutputChannel = {};
end
methods
    function obj = ultv_rb(varargin)
        if nargin == 0
            error('Error: To utilize "ultv_rb" one should specify at least one input argument (see "help ultv_rb" for further details)');
        elseif nargin == 1
            if ischar(varargin{1})
                obj.Name = varargin{1};
            else
                error('Error: The name of the uncertainty should be defined by a string (see "help ultv_rb" for further details)');
            end
        elseif nargin > 1
            if nargin == 2
                if ischar(varargin{1})
                    obj.Name = varargin{1};
                else
                    error('Error: The name of the uncertainty should be defined by a string (see "help ultv_rb" for further details)');
                end
                if isstruct(varargin{2})
                else
                    error('Error: The second input argument should be defined by a structure (see "help ultv_rb" for further details)');
                end
                if isfield(varargin{2},'BasisFunctionType')
                    a = 1:1:5;
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
                        error('Error: The length of a basis function should be defined as a natural number in N1 (see "help ultv_rb" for further details).');
                    end
                end
                if isfield(varargin{2},'PoleLocation')
                    if isreal(varargin{2}.PoleLocation) && ismatrix(varargin{2}.PoleLocation)
                        obj.PoleLocation = varargin{2}.PoleLocation;
                    else
                        error('Error: The pole location of a basis function should be defined as a real rational number. For continuous and distrete time systems, it is not allowed to have poles on the imaginary axis and unit circle respectively (see "help ultv_rb" for further details).');
                    end
                end
                if isfield(varargin{2},'SampleTime')
                    if isreal(varargin{2}.SampleTime) && ismatrix(varargin{2}.SampleTime)
                        obj.SampleTime = varargin{2}.SampleTime;
                    else
                        error('Error: The Sampling time of a basis function should be defined as a nonnegative real rational number (see "help ultv_rb" for further details).');
                    end
                end
                if isfield(varargin{2},'NumberOfRepetitions')
                    if isreal(varargin{2}.NumberOfRepetitions) && ismatrix(varargin{2}.NumberOfRepetitions) && min(varargin{2}.NumberOfRepetitions > 0) > 0 
                        obj.NumberOfRepetitions = varargin{2}.NumberOfRepetitions;
                    else
                        error('Error: The number of repetitions of the uncertainty should be defined as a vector with positive natural numbers in N1 (see "help ultv_rb" for further details).');
                    end
                end
                if isfield(varargin{2},'RelaxationType')
                    obj.RelaxationType = varargin{2}.RelaxationType;
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
                if isfield(varargin{2},'RateBounds')
                    if isreal(varargin{2}.RateBounds) && ismatrix(varargin{2}.RateBounds)
                        obj.RateBounds     = varargin{2}.RateBounds;
                    else
                        error('Error: The rate bounds of the uncertainty should be defined as a real matrix (see "help ultis" for further details).');
                    end
                end
                if isfield(varargin{2},'TerminalCost')
                    obj.TerminalCost   = varargin{2}.TerminalCost;
                end
                if isfield(varargin{2},'PrimalDual')
                    obj.PrimalDual = varargin{2}.PrimalDual;
                end
                if isfield(varargin{2},'InputChannel')
                    obj.InputChannel = varargin{2}.InputChannel;
                end
                if isfield(varargin{2},'OutputChannel')
                    obj.OutputChannel = varargin{2}.OutputChannel;
                end
            elseif nargin > 2
                if ischar(varargin{1})
                    obj.Name = varargin{1};
                else
                    error('Error: The name of the uncertainty should be defined by a string (see "help ultv_rb" for further details)');
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
                            a = 1:1:5;
                            if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)}) && min(ismember(varargin{j(i)},a)) > 0
                                obj.BasisFunctionType = varargin{j(i)};
                            else
                                error('Error: The type of basis function should be defined as natural numbers in N1 (see "help ultv_rb" for further details)');
                            end
                        case 'Length'
                            a = 1:1:10;
                            if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)}) && min(ismember(varargin{j(i)},a)) > 0
                                obj.Length = varargin{j(i)};
                            else
                                error('Error: The length of a basis function should be defined as a natural number in N1 (see "help ultv_rb" for further details).');
                            end
                        case 'PoleLocation'
                            if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)})
                                obj.PoleLocation = varargin{j(i)};
                            else
                                error('Error: The pole location of a basis function should be defined as a real rational number. For continuous and distrete time systems, it is not allowed to have poles on the imaginary axis and unit circle respectively (see "help ultv_rb" for further details).');
                            end
                        case 'SampleTime'
                            if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)})
                                obj.SampleTime = varargin{j(i)};
                            else
                                error('Error: The sampling time of a basis function should be defined as a nonnegative real rational number (see "help ultv_rb" for further details).');
                            end
                        case 'NumberOfRepetitions'
                            if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)}) && min(varargin{j(i)} > 0) > 0 
                                obj.NumberOfRepetitions = varargin{j(i)};
                            else
                                error('Error: The number of repetitions of the uncertainty should be defined as a vector with positive natural numbers in N1 (see "help ultv_rb" for further details).');
                            end
                        case 'RelaxationType'
                            obj.RelaxationType = varargin{j(i)};
                        case 'Bounds'
                            if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)})
                                obj.Bounds = varargin{j(i)};
                                obj.Polytope = [];
                            else
                                error('Error: The bounds of the uncertainty should be defined as a real matrix (see "help ultv_rb" for further details).');
                            end
                        case 'Polytope'
                            if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)})
                                obj.Polytope = varargin{j(i)};
                                obj.Bounds = [];
                            else
                                error('Error: The Polytope should be defined as a real matrix (see "help ultis" for further details).');
                            end
                        case 'RateBounds'
                            if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)})
                                obj.RateBounds = varargin{j(i)};
                            else
                                error('Error: The rate bounds of the uncertainty should be defined as a real matrix (see "help ultis" for further details).');
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
    function prob = iqcltv_rb_p(Delta,prob)
    % ---------------------------------------------------------------------
    % This function generates IQC-multipliers for the class ultv_rb. This
    % multiplier is precisely implemented as discussed in Section 5.4 of:
    %
    % [1] J. Veenman, C.W. Scherer, H. Koroglu, Robust stability and
    %     performance analysis based on integral quadratic constraints,
    %     European Journal of Control, Vol. 31, pp. 1-32, 2016.
    %
    % With Delta(delta) = blkdiag(delta1 In1 ,..., deltaN*IN), we have
    %
    % we have the multiplier
    %
    % Pi = [Psi1, 0   ]^*[P11  , P12][Psi1, 0  ]
    %      [0   , Psi2]  [P12^T, P22][0   ,Psi2]
    %
    % Here P = [P11,P12;P12.',P22] is a symmetric matrix variable which
    % must satisfy the constraints:
    %
    % Fij(Delta) = [I    ]^T[P11  , P12][I    ] >= 0,
    %              [Delta]  [P12^T, P22][Delta]
    %
    % for all Delta\in Polytope (Property "Polytope"), which can be
    % defined, indirectly, by
    %
    % Bounds = [-delta1, ... ,-deltaN] 
    %          [+delta1, ... ,+deltaN]
    %
    % and
    %
    % RateBounds = [-delta1dot, ... ,-deltaNdot] 
    %              [+delta1dot, ... ,+deltaNdot]
    %
    % or directly by:
    %
    %            [d1p1,d1dotp1, ... ,dNp1,dNdotp1]
    % Polytope = [... ,  ...  , ... ,... ,  ...  ]
    %            [d1pM,d1dotpM, ... ,dNpM,dNdotpM]
    %
    % In order to relax the semi-infinite dimensional matrix inequality,
    % this function considers 4 optional relaxation schemes:
    %
    %       1.) DG-scalings:
    %
    %           With Delta(delta) = blkdiag(delta1*In1,...,deltaN*InN) one
    %           can apply DG-scalings. For each deltai, i\in{1,...N} we
    %           then have Pi_i=[Di,Gi;Gi^T,-Di] with Di > 0 and Gi = -Gi^T.
    %
    % 2.) - 4.) Convex hull / Partial convexity / 1st order Polya:
    %
    %           These 3 relaxation schemes assume the Fij(Delta) is
    %           satisfied at the generator points of the polytope:
    %
    %           Fij(Delta^i) = [I      ]^T[P11  , P12][I      ] >= 0    (1)
    %                          [Delta^i]  [P12^T, P22][Delta^i]
    %
    %       2.) Convex hull relaxation. In addition to (1), the following
    %           LMI constraint is added to render the map Fij(Delta^i)
    %           concave:
    %
    %           P22 < 0
    %
    %       3.) Partial convexity: Just like with the DG-scalings, if
    %           Delta(delta) = blkdiag(delta1*In1,..., deltaN*InN) one can
    %           relax the previous constraint by only rendering P22
    %           partially negative definite:
    %
    %           Ji^T*P22*Ji < 0, i\in{1,...,N} 
    %
    %           with 
    %
    %           Ji = blkdiag(0,...,Ini,...,0)
    %
    %       4.) First order Polya relaxation: In addition to (1), one can
    %           also render the problem tractable if adding the additional
    %           constraints:
    %
    %           Fij(Delta^i,Delta^j) =[I      ]^T[P11  , P12][I      ] >= 0
    %                                 [Delta^i]  [P12^T, P22][Delta^j]
    %
    %           for i = 1,...,N, and j = i,...,N
    %
    %----------------------------------------------------------------------
    
    % Define outer-factor (basis function) of IQC-multiplier
    bft                             = get(Delta,'BasisFunctionType');
    nr                              = get(Delta,'NumberOfRepetitions');
    l                               = get(Delta,'Length');
    pl                              = get(Delta,'PoleLocation');
    Ts                              = get(Delta,'SampleTime');
    tc                              = get(Delta,'TerminalCost');
    lnr                             = length(nr);
    ll                              = length(l);
    lpl                             = length(pl);
    if ll < lnr
        l                           = l(1)*ones(lnr,1);
    end
    if lpl < lnr
        pl                          = pl(1)*ones(lnr,1);
    end
    if Ts == 0 && sum(pl == 0) > 0
        disp('Warning: In continuous time it is not allowed to have poles on the imaginary axis. The "PoleLocation" was set to its nominal value (i.e. -1).');
        pl                          = -ones(size(pl));
    elseif Ts > 0 && sum(abs(pl) == 1) > 0
        disp('Warning: In discrete time it is not allowed to have poles on the unit circle. The "PoleLocation" was set to its nominal value (i.e. 0).');
        pl                          = zeros(size(pl));
    elseif Ts == -1 && sum(abs(pl) == 1) > 0
        disp('Warning: In discrete time it is not allowed to have poles on the unit circle. The "PoleLocation" was set to its nominal value (i.e. 0).');
        pl                          = zeros(size(pl));
    end
    H                               = ss([],[],[],[],Ts);
    for i = 1:length(nr)
        H                           = fAugss(H,fBasis(l(i),pl(i),nr(i),bft,Ts),1);
        q(i,1)                      = l(i)*nr(i);
        k(i,1)                      = (l(i)-1)*nr(i);
    end
    if Ts == 0
        Phi1                        = ss(H.a,H.b,[H.c;eye(size(H.c,2))],[H.d;zeros(size(H.c,2),size(H.d,2))],Ts);
    else
        Phi1                        = ss(H.a,H.b,[H.c;H.a],[H.d;H.b],Ts);
    end
    Phi2                            = ss(H.a,[H.b,eye(size(H.b,1))],[H.c;zeros(size(H.a))],blkdiag(H.d,eye(size(H.a,2))),Ts);
    prob                            = fPsi(prob,Phi1,Phi2);
    
    % Define KYP certificate structure (for continuous or discrete time analysis)
    if Ts == 0
        kyp                         = [0,1;1,0];
    else
        kyp                         = [-1,0;0,1];
    end
    
    % Define multiplier variables plus lmi constraints
    rt                              = get(Delta,'RelaxationType');
    
    switch rt
        case 'DG'
            if isempty(Delta.Bounds) || isempty(Delta.RateBounds)
                error('Error: To apply DG-scalings one should specify the properties "Bounds" and "RateBounds" of the uncertainty.');
            end
            P11a                    = iqcvar(prob,[0,0],'symmetric');
            P11b                    = iqcvar(prob,[0,0],'symmetric');
            P12a                    = iqcvar(prob,[0,0],'skew');
            P12b                    = iqcvar(prob,[0,0],'skew');
            sc11a                   = [];
            sc11b                   = [];
            for i = 1:length(nr)
                P11a                = blkdiag(P11a,iqcvar(prob,[q(i),q(i)],'symmetric'));
                P11b                = blkdiag(P11b,iqcvar(prob,[k(i),k(i)],'symmetric'));
                P12a                = blkdiag(P12a,iqcvar(prob,[q(i),q(i)],'skew'));
                P12b                = blkdiag(P12b,iqcvar(prob,[k(i),k(i)],'skew'));
                sc11a               = blkdiag(sc11a,Delta.Bounds(2,i)*eye(q(i)));
                sc11b               = blkdiag(sc11b,Delta.RateBounds(2,i)*eye(k(i)));
            end
            P11                     = blkdiag(P11a,P11b);
            P12                     = blkdiag(P12a,P12b);
            sc11                    = blkdiag(sc11a,sc11b);
            P22                     = -P11;
            P                       = [P11,P12;P12.',P22];
            sc                      = blkdiag(sc11,eye(size(sc11,1)));
            prob                    = fP(prob,P);
            prob                    = fsc(prob,sc);
            prob                    = iqclmi(prob,P11,1);
            switch tc
                case 'on'
                    xPhi            = size(Phi1.a,1) + size(Phi2.a,1);
                    Z               = iqcvar(prob,[xPhi,xPhi],'symmetric');
                    prob            = fZ(prob,Z,eye(Z.Dim(1)));
                    nPhi1           = size(Phi1,2);
                    nPhi2           = size(Phi2,2);
                    B               = fss2m(blkdiag(Phi1,Phi2)*[zeros(nPhi1,nPhi2);eye(nPhi2)]);
                    W               = blkdiag(kron(kyp,Z),P);
                    prob            = iqclmi(prob,W,-1,0,B);
            end
        case {'CH','PC','ZP'}
            P11                     = iqcvar(prob,[sum([q;k]),sum([q;k])],'symmetric');
            P12                     = iqcvar(prob,[sum([q;k]),sum([q;k])],'full');
            P22                     = iqcvar(prob,[sum([q;k]),sum([q;k])],'symmetric');
            P                       = [P11,P12;P12.',P22];
            sc                      = blkdiag(eye(P11.Dim(1)),eye(P22.Dim(1)));
            prob                    = fP(prob,P);
            prob                    = fsc(prob,sc);
            if ~isempty(Delta.Bounds) && ~isempty(Delta.RateBounds)
                La                  = polydec(pvec('box',[Delta.Bounds';Delta.RateBounds']))';
            elseif ~isempty(Delta.Polytope)
                La                  = Delta.Polytope;
                m1                  = 1:2:size(La,2);
                m2                  = 2:2:size(La,2);
                La                  = [La(:,m1),La(:,m2)];
            end
            ijk                     = length(nr);
            for i = 1:size(La,1)
                A1a                 = [];
                A1b                 = [];    
                for j = 1:ijk
                    A1a             = blkdiag(A1a,La(i,j)*eye(q(j)));
                    A1b             = blkdiag(A1b,La(i,j+ijk)*eye(k(j)));
                end
                A1                  = blkdiag(A1a,A1b);
                A                   = [eye(size(A1,2));A1];
                prob                = iqclmi(prob,P,1,0,A);
            end
            switch rt
                case 'CH'
                    prob            = iqclmi(prob,P22,-1);
                case 'PC'
                    sLa1            = size(La,1);
                    if sLa1 ~= 2^(2*lnr)
                        error('Error: For the relaxation option "PC" the individual uncertainties should be defined on the interval ci*[-1,1].');
                    end
                    for i = 1:length(nr)
                        A           = fJ(q(i),sum(q(i+1:end))+sum(k),sum(q(1:i-1)));
                        prob        = iqclmi(prob,P22,-1,0,A);
                        B           = fJ(k(i),sum(k(i+1:end)),sum(k(1:i-1))+sum(q));
                        prob        = iqclmi(prob,P22,-1,0,B);
                    end
                case 'ZP'
                    ijk             = length(nr);
                    for i = 1:length(La)
                        for j = i+1:length(La)
                            A1a     = [];
                            A1b     = [];
                            A2a     = [];
                            A2b     = [];
                            for r = 1:ijk
                                A1a = blkdiag(A1a,La(i,r)*eye(q(r)));
                                A1b = blkdiag(A1b,La(i,r+ijk)*eye(k(r)));
                                A2a = blkdiag(A2a,La(j,r)*eye(q(r)));
                                A2b = blkdiag(A2b,La(j,r+ijk)*eye(k(r)));
                            end
                            A1      = blkdiag(A1a,A1b);
                            A2      = blkdiag(A2a,A2b);
                            A       = [eye(size(A1,2));A1;eye(size(A2,2));A2];
                            V       = oblkdiag(P);
                            prob    = iqclmi(prob,V,1,0,A);
                        end
                    end
            end
            switch tc
                case 'on'
                    xPhi            = size(Phi1.a,1) + size(Phi2.a,1);
                    Z               = iqcvar(prob,[xPhi,xPhi],'symmetric');
                    prob            = fZ(prob,Z,eye(Z.Dim(1)));
                    nPhi1           = size(Phi1,2);
                    nPhi2           = size(Phi2,2);
                    B               = fss2m(blkdiag(Phi1,Phi2)*[zeros(nPhi1,nPhi2);eye(nPhi2)]);
                    W               = blkdiag(kron(kyp,Z),P);
                    prob            = iqclmi(prob,W,-1,0,B);
            end
    end
    IO                              = [sum(nr),sum(nr),sum(nr),sum(k),P11.Dim(1),P22.Dim(1),P11.Dim(1),P22.Dim(1),size(Phi1.a,1),size(Phi2.a,1)];
    prob                            = fIO(prob,IO);
    end
    function prob = iqcltv_rb_d(Delta,prob)
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
       disp('------------------------------------------------------------------');
       disp(' LTV rate bounded diagonally repeated parametric uncertainty block');
       disp('------------------------------------------------------------------');
       disp(' Uncertainty details:');
       disp(['  - Name:                           ',num2str(obj.Name)]);
       disp(['  - Number of uncertain parameters: ',num2str(length(obj.NumberOfRepetitions))]);
       disp(['  - Number of repetitions:          ',num2str(fDim(obj.NumberOfRepetitions))]);
       disp(['  - Input channels plant:           ',a1]);
       disp(['  - Output channels plant:          ',b1]);
       if isempty(obj.Bounds)
           disp('  - Polytope:');
           disp(' ')
           disp(obj.Polytope);
       else
           disp('  - Bounds:');
           disp(' ')
           disp(obj.Bounds);
           disp('  - Rate bounds:');
           disp(' ')
           disp(obj.RateBounds);
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