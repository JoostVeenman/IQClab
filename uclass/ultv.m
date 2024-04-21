classdef ultv < handle & matlab.mixin.SetGetExactNames

% -------------------------------------------------------------------------
%
% IQClab:      Version 3.4.0
% Copyright:   This is copyrighted material owned by Novantec B.V.
% Terms:       IQClab is available for non-commercial usage under a
%              Creative Commons (Attribution-NoDerivatives 4.0
%              International (CC BY-ND 4.0)) license:  
%              https://creativecommons.org/licenses/by-nd/4.0/
%              For further information please visit iqclab.eu
% Author:      J.Veenman
% Date:        06-11-2019
% 
% -------------------------------------------------------------------------
%
% Description: Class definition for LTV parametric uncertainties of the
%              form: 
%
%              Delta(delta) = sum_{i=1}^N delta_i*Ti = ...
%                                 = delta1*T1 + ... + deltaN*TN,
%
%              where
%              - Ti are some fixed matrices Ti \in R^{n x m}
%              - delta:[0,\infty)->La is a piecewise continuous
%                time-varying parameter vector that takes its values from
%                the compact polytope 
%                      La = co{delta^1, ... ,delta^M} = ...
%                                = {sum_{a=1}^M b_a*delta^a: ...
%                              b_a\leq0, sum_{a=1}^M b_a = 1}
%                                 
%                with delta^j = (delta_1^j, ... ,delta_N^j), j\in{1,...,M},
%                as generator point and with 0\in La.
%              - La is assumed to be star convex: [0,1]La\subset La.
%
% Syntax:      delta = ultv('name')
%              delta = ultv('name',varargin)
%
% Usage:       "delta = ultv('name')" defines an LTV parametric
%              uncertainty on the interval [-1,1], which is repeated once:
%              I.e. delta1\in[-1,1] and T1=1.
%
%              For "delta = ultv('name',varargin)", the varargin inputs
%              come in pairs and can be defined as: 
%
%              delta = ultv('name','prop1','value1','prop2','value2',...)
%
%              Alternatively, the properties can be set by defining the
%              structure:
%
%              ultvOpt.prop1 = value1
%              ultvOpt.prop2 = value2
%                       ...              
%              ultvOpt.propN = valueN
% 
%              delta = ultv('name',ultvOpt)
%
%              Finally, the properties can be set and retrieved by using:
%
%              set(delta,'propX','valueX') and get(delta,'propX').
%
%              The properties that can be specified are:
%            
%               1.) 'Polytope' specifies the generator points of delta by:
%
%                                 [delta1^1, ... ,deltaN^1]
%                      Polytope = [   ...  , ... ,   ...  ]
%                                 [delta1^M, ... ,deltaN^M]
%
%               2.) 'UncertaintyMap' specifies the matrices T1,...,TN as
%                   a cell array T = {T1,...,TN}. This defines the
%                   uncertainty mapping: 
%
%                      Delta(delta) = sum_{i=1}^N delta_i*Ti = ...
%                                       = delta1*T1 + ... + deltaN*TN.
%
%               3.) 'InputChannel' specifies which input-channels of the
%                   uncertain plant are affected by Delta. The channels
%                   should be specified as a cell:
%
%                           InputChannel = {[Chnl_x,...,Chnl_y]}
%
%                    Here the order of the channels is not relevant, while
%                    the row-length equals the number of outputs of Delta.
%
%               4.) 'OutputChannel' specifies which input-channels of the
%                   uncertain plant are affected by Delta. The channels
%                   should be specified as a cell:
%
%                           OutputChannel = {[Chnl_x,...,Chnl_y]}
%
%                    Here the order of the channels is not relevant, while
%                    the row-length equals the number of inputs of Delta.
%
%
%              In addition, various properties of the corresponding
%              IQC-multiplier can be specified:
%
%               5.) 'RelaxationType' specifies the multiplier relaxation
%                   type, options are (default = 'DG'):                     
%
%                      a.) 'DG'        DG scalings
%                      b.) 'CH'        Convex Hull relaxation
%                      c.) 'PC'        Partial Convexity*
%                      d.) 'ZP'        Zeroth order Polya relaxation
%
%                   % Note: The partial convexity relaxation can only be
%                           applied in the case that Delta(delta) is
%                           defined by diagonally repeated parameters: 
%
%                           Delta(delta) = diag(delta1,...,deltaN)
%
%               6.) 'PrimalDual' specifies whether the multiplier should be
%                   a primal/dual parametrization (default = 'Primal').
%                   
%                   Note: For a standard IQC-analysis, all multipliers are
%                   Primal ones.
% -------------------------------------------------------------------------

properties (SetAccess = public)
    % Uncertainty name
    Name string = 'delta';
    
    % Generator points of the uncertainty
    Polytope double {mustBeReal,mustBeFinite} = [-1,1];
    
    % Uncertainty matrices defining the full-block uncertainty mapping
    UncertaintyMap cell = {1};

    % Relaxation type of the uncertainty (convex hull(CH), partial convexity (PC), Zeroth order Polya (ZP))
    RelaxationType string {mustBeMember(RelaxationType,{'DG','CH','PC','ZP'})} = 'DG';

    % Sampling time of the multiplier basis function
    SampleTime double {mustBeReal,mustBeFinite} = 0;
    
    % Primal or Dual multiplier
    PrimalDual string {mustBeMember(PrimalDual,{'Primal','Dual'})} = 'Primal';
    
    % InputChannel
    InputChannel = {};
    
    % OutputChannel
    OutputChannel = {};
end
methods
    function obj = ultv(varargin)
        if nargin == 0
            error('Error: To utilize "ultv" one should specify at least one input argument (see "help ultv" for further details)');
        elseif nargin == 1
            if ischar(varargin{1})
                obj.Name = varargin{1};
            else
                error('Error: The name of the uncertainty should be defined by a string (see "help ultv" for further details)');
            end
        elseif nargin > 1
            if nargin == 2
                if ischar(varargin{1})
                    obj.Name = varargin{1};
                else
                    error('Error: The name of the uncertainty should be defined by a string (see "help ultv" for further details)');
                end
                if isstruct(varargin{2})
                else
                    error('Error: The second input argument should be defined by a structure (see "help ultv" for further details)');
                end
                if isfield(varargin{2},'Polytope')
                    if isreal(varargin{2}.Polytope) && ismatrix(varargin{2}.Polytope)
                        obj.Polytope = varargin{2}.Polytope;
                    else
                        error('Error: The generator points of the uncertainty should be defined as a real valued matrix (see "help ultv" for further details).');
                    end
                end
                if isfield(varargin{2},'UncertaintyMap')
                    if isreal(varargin{2}.UncertaintyMap) && iscell(varargin{2}.UncertaintyMap)
                        obj.UncertaintyMap = varargin{2}.UncertaintyMap;
                    else
                        error('Error: "UncertaintyMap" should be defined as cell arry of real valued matrices all of the same dimension (see "help ultv" for further details).');
                    end
                end
                if isfield(varargin{2},'SampleTime')
                    if isreal(varargin{2}.SampleTime) && ismatrix(varargin{2}.SampleTime)
                        obj.SampleTime = varargin{2}.SampleTime;
                    else
                        error('Error: The Sampling time of a basis function should be defined as a nonnegative real rational number (see "help ultv" for further details).');
                    end
                end
                if isfield(varargin{2},'RelaxationType')
                    obj.RelaxationType = varargin{2}.RelaxationType;
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
                    error('Error: The name of the uncertainty should be defined by a string (see "help ultv" for further details)');
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
                        case 'Polytope'
                            if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)})
                                obj.Polytope = varargin{j(i)};
                            else
                                error('Error: The generator points of the uncertainty should be defined as a real valued matrix (see "help ultv" for further details).');
                            end
                        case 'UncertaintyMap'
                        if iscell(varargin{j(i)})
                            obj.UncertaintyMap = varargin{j(i)};
                        else
                            error('Error: "UncertaintyMap" should be defined as cell arry of real valued matrices all of the same dimension (see "help ultv" for further details).');
                        end
                        case 'RelaxationType'
                            obj.RelaxationType = varargin{j(i)};
                        case 'RelaxProp'
                            obj.RelaxProp = varargin{j(i)};
                        case 'UncertaintyDomain'
                            if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)})
                                obj.UncertaintyDomain = varargin{j(i)};
                            else
                                error('Error: The uncertainty domain should be defined as a real matrix (see "help ultv" for further details).');
                            end
                        case 'SampleTime'
                            if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)})
                                obj.SampleTime = varargin{j(i)};
                            else
                                error('Error: The sampling time of a basis function should be defined as a nonnegative real rational number (see "help ultv" for further details).');
                            end
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
    function prob = iqcltv_p(Delta,prob)
    % ---------------------------------------------------------------------
    % This function generates IQC-multipliers for the class ultv. These are
    % precisely implemented as discussed in Section 5.2 of:  
    %
    % [1] J. Veenman, C.W. Scherer, H. Koroglu, Robust stability and
    %     performance analysis based on integral quadratic constraints,
    %     European Journal of Control, Vol. 31, pp. 1-32, 2016.
    %
    % With
    %
    % Delta(delta) = H1*delta1 + ,..., + HN*deltaN, 
    % delta        = col(delta1,...,deltaN) \in Polytope
    %
    % we have:
    %
    % Fij(Delta) = [I    ]^T[P11  , P12][I    ] >= 0,
    %              [Delta]  [P12^T, P22][Delta]
    %
    % for all Delta\in Polytope (Property "Polytope"), defined as:
    %
    % Polytope = [delta1_p1,...,deltam_p1]
    %            [  ...    ,...,   ...   ]
    %            [delta1_pn,...,deltam_pn]
    %
    % In order to relax this semi-infinite dimensional matrix inequality,
    % this function considers 4 optional relaxation schemes:
    %
    %       1.) DG-scalings:
    %
    %           If Delta(delta) = blkdiag(delta1*In1,...,deltaN*InN) one
    %           can apply DG-scalings. For each deltai, i\in{1,...N} we
    %           have Pi_i=[Di,Gi;Gi^T,-Di] with Di > 0 and Gi = -Gi^T.
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

    % Get uncertaity info
    La                             = get(Delta,'Polytope');
    H                              = get(Delta,'UncertaintyMap');
    Ts                             = get(Delta,'SampleTime');
    nou                            = length(H);
    [suo,sui]                      = size(H{1});
    
    % Create outer factor (basis function) of the uncertainty
    Phi1                           = fAugss(ss([],[],[],[],Ts),fBasis(1,-1,suo,1,Ts),1);
    Phi2                           = fAugss(ss([],[],[],[],Ts),fBasis(1,-1,sui,1,Ts),1);
    prob                           = fPsi(prob,Phi1,Phi2);
    
    % Define multiplier variables plus lmi constraints
    rt                             = get(Delta,'RelaxationType');
    switch rt
        case 'DG'
            J                      = zeros(suo,sui);
            for i = 1:nou
               J                   = J + H{i};
            end
            if suo == sui && nou <= suo
                if norm(J-eye(suo)) ~= 0 || trace(J) ~= suo
                    error('Error: It is only possible to apply DG scalings if the uncertainty block has all uncertainties repeated on the diagonal.');
                end
            else
                error('Error: It is only possible to apply DG scalings if the uncertainty block is square.');
            end
            if size(La,1) ~= 2^nou
                error('Error: DG-scaling can only be applied for Hyper rectangular and symmetric polytopes');
            end
            for i = 1:nou
                if sum(La(:,i)) == 0
                    Bounds(:,i)    = abs(La(1,i))*[-1;1];
                else
                    error('Error: DG-scaling can only be applied for Hyper rectangular and symmetric polytopes');
                end
            end
            P11                    = iqcvar(prob,[0,0],'symmetric');
            P12                    = iqcvar(prob,[0,0],'skew');
            sc11                   = zeros(size(H{1}));
            for i = 1:nou
                k                  = 0;
                n                  = 0;
                for j = 1:suo
                    if H{i}(j,j) == 0
                        k          = k + 1;
                    elseif H{i}(j,j) == 1
                        n          = n + 1;
                    else
                        error('Error: It is only possible to apply DG scalings if the uncertainty block has all uncertainties repeated on the diagonal.');
                    end
                end
                if k + n ~= suo
                    error('Error: It is only possible to apply DG scalings if the uncertainty block has all uncertainties repeated on the diagonal.');
                end
                P11                = blkdiag(P11,iqcvar(prob,[n,n],'symmetric'));
                P12                = blkdiag(P12,iqcvar(prob,[n,n],'skew'));
                sc11               = sc11 + Bounds(2,i)*H{i};
            end
            P22                    = -P11;
            P                      = [P11,P12;P12.',P22];
            sc                     = blkdiag(sc11,eye(size(sc11,1)));
            prob                   = fP(prob,P);
            prob                   = fsc(prob,sc);
            prob                   = iqclmi(prob,P11,1);
        case {'CH','PC','ZP'}
            P11                    = iqcvar(prob,[sui,sui],'symmetric');
            P12                    = iqcvar(prob,[sui,suo],'full');
            P22                    = iqcvar(prob,[suo,suo],'symmetric');
            P                      = [P11,P12;P12.',P22];
            sc                     = blkdiag(eye(sui+suo));
            prob                   = fP(prob,P);
            prob                   = fsc(prob,sc);
            for i = 1:size(La,1)
                De{i}              = zeros(suo,sui);
                for j = 1:nou
                    De{i}          = De{i} + La(i,j)*H{j};
                end
                A                  = [eye(sui);De{i}];
                prob               = iqclmi(prob,P,1,0,A);
            end
            switch rt
                case 'CH' % Convex Hull Relaxation
                    prob           = iqclmi(prob,P22,-1);
                case 'PC'
                    J              = zeros(suo,sui);
                    for i = 1:nou
                       J           = J + H{i};
                    end
                    if suo == sui && nou <= suo
                        if norm(J-eye(suo)) ~= 0 || trace(J) ~= suo
                            error('Error: It is only possible to apply the PC relaxation scheme if the uncertainty block has all uncertainties repeated on the diagonal.');
                        end
                    else
                        error('Error: It is only possible to apply  the PC relaxation scheme if the uncertainty block is square.');
                    end
                    if size(La,1) ~= 2^nou
                        error('Error:  The PC relaxation scheme can only be applied for Hyper rectangular and symmetric polytopes');
                    end
                    for i = 1:nou
                        if sum(La(:,i)) ~= 0
                            error('Error:  The PC relaxation scheme can only be applied for Hyper rectangular and symmetric polytopes');
                        end
                    end
                    for i = 1:nou
                        k          = 0;
                        n          = 0;
                        for j = 1:suo
                            if H{i}(j,j) == 0
                                k  = k + 1;
                            elseif H{i}(j,j) == 1
                                n  = n + 1;
                            else
                                error('Error: It is only possible to apply the PC relaxation scheme if the uncertainty block has all uncertainties repeated on the diagonal.');
                            end
                        end
                        if k + n ~= suo
                            error('Error: It is only possible to the PC relaxation scheme if the uncertainty block has all uncertainties repeated on the diagonal.');
                        end
                        m          = sum(H{i});
                        n1         = 0;
                        n2         = 0;
                        n3         = 0;
                        for ijk = 1:length(m)
                            if m(ijk) == 0
                                n1 = n1 + 1;
                            else
                                break
                            end
                        end
                        for ijk = n1+1:length(m)
                            if m(ijk) == 1
                                n2 = n2 + 1;
                            else
                                break
                            end
                        end
                        for ijk = n1+n2+1:length(m)
                            if m(ijk) == 0
                                n3 = n3 + 1;
                            else
                                break
                            end
                        end
                        prob       = iqclmi(prob,P22,-1,0,fJ(n2,n3,n1));
                    end
                case 'ZP'
                    for i = 1:length(De)
                        for j = i + 1:length(De)
                            prob   = iqclmi(prob,oblkdiag(P),1,0,[eye(sui);De{i};eye(sui);De{j}]);
                        end
                    end
            end
    end
    IO                             = [suo,sui,0,0,P11.Dim(2),P22.Dim(2),P11.Dim(2),P22.Dim(2),0,0];
    prob                           = fIO(prob,IO);
    end
    function prob = iqcltv_d(Delta,prob)
    % ---------------------------------------------------------------------
    % This function generates dual IQC-multipliers for the class ultv.
    % These are implemented as discussed in Section 5.2 of (but in the dual
    % form thereof).
    %
    % [1] J. Veenman, C.W. Scherer, H. Koroglu, Robust stability and
    %     performance analysis based on integral quadratic constraints,
    %     European Journal of Control, Vol. 31, pp. 1-32, 2016.
    %
    % With
    %
    % Delta(delta) = H1*delta1 + ,..., + HN*deltaN, 
    % delta        = col(delta1,...,deltaN) \in Polytope
    %
    % we have:
    %
    % Fij(Delta) = [-Delta^T]^T[P11 , P12][-Delta^T] <= 0
    %              [I       ] [P12^T, P22][I       ]
    %
    % for all Delta\in Polytope (Property "Polytope"), defined as:
    %
    % Polytope = [delta1_p1,...,deltam_p1]
    %            [  ...    ,...,   ...   ]
    %            [delta1_pn,...,deltam_pn]
    %
    % In order to relax this semi-infinite dimensional matrix inequality,
    % this function considers 4 optional relaxation schemes:
    %
    %       1.) DG-scalings:
    %
    %           If Delta(delta) = blkdiag(delta1*In1,...,deltaN*InN) one
    %           can apply DG-scalings. For each deltai, i\in{1,...N} we
    %           have Pi_i=[Di,Gi;Gi^T,-Di] with Di > 0 and Gi = -Gi^T.
    %
    % 2.) - 4.) Convex hull / Partial convexity / 1st order Polya:
    %
    %           These 3 relaxation schemes assume the Fij(Delta) is
    %           satisfied at the generator points of the polytope:
    %
    %           Fij(Delta^i) = [*]^T[P11  ,P12][-(Delta^i)^T] <= 0      (1)
    %                          [*]  [P12^T,P22][I           ]
    %
    %       2.) Convex hull relaxation. In addition to (1), the following
    %           LMI constraint is added to render the map Fij(Delta^i)
    %           concave:
    %
    %           P11 > 0
    %
    %       3.) Partial convexity: Just like with the DG-scalings, if
    %           Delta(delta) = blkdiag(delta1*In1,..., deltaN*InN) one can
    %           relax the previous constraint by only rendering P11
    %           partially positive definite:
    %
    %           Ji^T*P11*Ji > 0, i\in{1,...,N} 
    %
    %           with 
    %
    %           Ji = blkdiag(0,...,Ini,...,0)
    %
    %       4.) First order Polya relaxation: In addition to (1), one can
    %           also render the problem tractable if adding the additional
    %           constraints:
    %
    %           Fij(Delta^i,Delta^j) = ...
    %              =[-(Delta^j)^T]^T[P11  , P12][-(Delta^j)^T] <= 0
    %               [I           ]  [P12^T, P22][I           ]
    %
    %           for i = 1,...,N, and j = i,...,N
    %
    %----------------------------------------------------------------------

    % Get uncertaity info
    La                             = get(Delta,'Polytope');
    H                              = get(Delta,'UncertaintyMap');
    Ts                             = get(Delta,'SampleTime');
    nou                            = length(H);
    [suo,sui]                      = size(H{1});
    
    % Create outer factor (basis function) of the uncertainty
    Phi1                           = fAugss(ss([],[],[],[],Ts),fBasis(1,-1,suo,1,Ts),1);
    Phi2                           = fAugss(ss([],[],[],[],Ts),fBasis(1,-1,sui,1,Ts),1);
    prob                           = fPsi(prob,Phi1,Phi2);
    
    % Define multiplier variables plus lmi constraints
    rt                             = get(Delta,'RelaxationType');
    switch rt
        case 'DG'
            J                      = zeros(suo,sui);
            for i = 1:nou
               J                   = J + H{i};
            end
            if suo == sui && nou <= suo
                if norm(J-eye(suo)) ~= 0 || trace(J) ~= suo
                    error('Error: It is only possible to apply DG scalings if the uncertainty block has all uncertainties repeated on the diagonal.');
                end
            else
                error('Error: It is only possible to apply DG scalings if the uncertainty block is square.');
            end
            if size(La,1) ~= 2^nou
                error('Error: DG-scaling can only be applied for Hyper rectangular and symmetric polytopes');
            end
            for i = 1:nou
                if sum(La(:,i)) == 0
                    Bounds(:,i)    = abs(La(1,i))*[-1;1];
                else
                    error('Error: DG-scaling can only be applied for Hyper rectangular and symmetric polytopes');
                end
            end
            P11                    = iqcvar(prob,[0,0],'symmetric');
            P12                    = iqcvar(prob,[0,0],'skew');
            sc11                   = zeros(size(H{1}));
            for i = 1:nou
                k                  = 0;
                n                  = 0;
                for j = 1:suo
                    if H{i}(j,j) == 0
                        k          = k + 1;
                    elseif H{i}(j,j) == 1
                        n          = n + 1;
                    else
                        error('Error: It is only possible to apply DG scalings if the uncertainty block has all uncertainties repeated on the diagonal.');
                    end
                end
                if k + n ~= suo
                    error('Error: It is only possible to apply DG scalings if the uncertainty block has all uncertainties repeated on the diagonal.');
                end
                P11                = blkdiag(P11,iqcvar(prob,[n,n],'symmetric'));
                P12                = blkdiag(P12,iqcvar(prob,[n,n],'skew'));
                sc11               = sc11 + Bounds(2,i)*H{i};
            end
            P22                    = -P11;
            P                      = [P11,P12;P12.',P22];
            sc                     = blkdiag(sc11,eye(size(sc11,1)));
            prob                   = fP(prob,P);
            prob                   = fsc(prob,sc);
            prob                   = iqclmi(prob,P11,1);
        case {'CH','PC','ZP'}
            P11                    = iqcvar(prob,[sui,sui],'symmetric');
            P12                    = iqcvar(prob,[sui,suo],'full');
            P22                    = iqcvar(prob,[suo,suo],'symmetric');
            P                      = [P11,P12;P12.',P22];
            sc                     = blkdiag(eye(sui+suo));
            prob                   = fP(prob,P);
            prob                   = fsc(prob,sc);
            for i = 1:size(La,1)
                De{i}              = zeros(suo,sui);
                for j = 1:nou
                    De{i}          = De{i} + La(i,j)*H{j};
                end
                A                  = [-De{i}';eye(suo)];
                prob               = iqclmi(prob,P,-1,0,A);
            end
            switch rt
                case 'CH' % Convex Hull Relaxation
                    prob           = iqclmi(prob,P11,1);
                case 'PC'
                    J              = zeros(suo,sui);
                    for i = 1:nou
                       J           = J + H{i};
                    end
                    if suo == sui && nou <= suo
                        if norm(J-eye(suo)) ~= 0 || trace(J) ~= suo
                            error('Error: It is only possible to apply the PC relaxation scheme if the uncertainty block has all uncertainties repeated on the diagonal.');
                        end
                    else
                        error('Error: It is only possible to apply  the PC relaxation scheme if the uncertainty block is square.');
                    end
                    if size(La,1) ~= 2^nou
                        error('Error:  The PC relaxation scheme can only be applied for Hyper rectangular and symmetric polytopes');
                    end
                    for i = 1:nou
                        if sum(La(:,i)) ~= 0
                            error('Error:  The PC relaxation scheme can only be applied for Hyper rectangular and symmetric polytopes');
                        end
                    end
                    for i = 1:nou
                        k          = 0;
                        n          = 0;
                        for j = 1:suo
                            if H{i}(j,j) == 0
                                k  = k + 1;
                            elseif H{i}(j,j) == 1
                                n  = n + 1;
                            else
                                error('Error: It is only possible to apply the PC relaxation scheme if the uncertainty block has all uncertainties repeated on the diagonal.');
                            end
                        end
                        if k + n ~= suo
                            error('Error: It is only possible to the PC relaxation scheme if the uncertainty block has all uncertainties repeated on the diagonal.');
                        end
                        m          = sum(H{i});
                        n1         = 0;
                        n2         = 0;
                        n3         = 0;
                        for ijk = 1:length(m)
                            if m(ijk) == 0
                                n1 = n1 + 1;
                            else
                                break
                            end
                        end
                        for ijk = n1+1:length(m)
                            if m(ijk) == 1
                                n2 = n2 + 1;
                            else
                                break
                            end
                        end
                        for ijk = n1+n2+1:length(m)
                            if m(ijk) == 0
                                n3 = n3 + 1;
                            else
                                break
                            end
                        end
                        prob       = iqclmi(prob,P11,1,0,fJ(n2,n3,n1));%%
                    end
                case 'ZP'
                    for i = 1:length(De)
                        for j = i + 1:length(De)
                            prob   = iqclmi(prob,oblkdiag(P),-1,0,[-De{i}';eye(suo);-De{j}';eye(suo)]);
                        end
                    end
            end
    end
    IO                             = [suo,sui,0,0,P11.Dim(2),P22.Dim(2),P11.Dim(2),P22.Dim(2),0,0];
    prob                           = fIO(prob,IO);
    end
    function display(obj)
        s1 = size(obj.UncertaintyMap{1},1);
        s2 = size(obj.UncertaintyMap{1},2);
        um = length(obj.UncertaintyMap);
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
        disp('----------------------------------------------');
        disp(' LTV parametric (full-block) uncertainty');
        disp('----------------------------------------------');
        disp(' Uncertainty details:');
        disp(['  - Name:                               ',num2str(obj.Name)]);
        disp(['  - Number of uncertain parameters:     ',num2str(um)]);
        disp(['  - Size of uncertainty block Delta:    ',num2str(s1),' x ',num2str(s2)]);
        disp(['  - Input channels plant:               ',a1]);
        disp(['  - Output channels plant:              ',b1]);
        disp('  - Convex Hull generator points:  ');
        disp('           [delta1^1, ... ,deltaN^1]');
        disp('           [   ...  , ... ,   ...  ]');
        disp('           [delta1^M, ... ,deltaN^M]');
        disp('    La =');
        disp(' ');
        disp(obj.Polytope);
        disp('  - Matrices Ti defining Delta(delta)=');
        disp('           =delta1*T1 + ... + deltaN*TN:');
        disp(' ');
        for ijk = 1:um
            disp(['T',num2str(ijk),' =']);
            disp(obj.UncertaintyMap{ijk})
        end
        disp(' Multiplier details:');
        disp(['  - Relaxation type:                    ',num2str(obj.RelaxationType)]);
        disp(['  - Multiplier type:                    ',num2str(obj.PrimalDual)]);
    end 
end
end
