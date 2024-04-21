classdef udel < handle & matlab.mixin.SetGetExactNames

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
% Date:        04-11-2019
% 
% -------------------------------------------------------------------------
%
% Description: Class definition for LTI diagonally repeated delay
%              uncertainties of the following forms: 
%
%              Type 1:     delta =     (e^(-s*tau)-1)*I_nr
%              Type 2:     delta = 1/s*(e^(-s*tau)-1)*I_nr
%
%              where:
%              - s is the Laplace operator
%              - tau\in[0,tau_max] and tau_max>0 is the maximum delay.
%              - nr is the number of repetitions
%
% Syntax:      delta = udel('name')
%              delta = udel('name',varargin)
%
% Usage:       "delta = udel('name')" defines an LTI diagonally repeated
%              delay uncertainty of Type 1, which has a maximum delay of
%              tau = 1 s, and is repeated once.
%
%              For "delta = udel('name',varargin)", the varargin inputs
%              come in pairs and can be defined as: 
%
%              delta = udel('name','prop1','value1','prop2','value2',...)
%
%              Alternatively, the properties can be set by defining the
%              structure:
%
%              udelOpt.prop1 = value1
%              udelOpt.prop2 = value2
%                       ...              
%              udelOpt.propN = valueN
% 
%              delta = udel('name',udelOpt)
%
%              Finally, the properties can be set and retrieved by using:
%
%                   set(delta,'propX','valueX') and get(delta,'propX')
%
%              The properties that can be specified are:
%            
%               1.) 'DelayType' specifies the type of delay operator (i.e
%                   Type 1 or Type 2) (default = 1).
%
%               2.) 'NumberOfRepetitions' specifies the number of
%                   repetitions of delta (default =  1).
%
%               3.) 'DelayTime' specifies the maximum time delay tau
%                   (default = 1). It is assumed that the value zero is
%                   contained in the uncertainty set.
%
%               4.) 'InputChannel' specifies which input-channels of the
%                   uncertain plant are affected by Delta. The channels
%                   should be specified as a cell:
%
%                           InputChannel = {[Chnl_x,...,Chnl_y]}
%
%                    Here the order of the channels is not relevant, while
%                    the row-length equals the number of outputs of Delta.
%
%               5.) 'OutputChannel' specifies which input-channels of the
%                   uncertain plant are affected by Delta. The channels
%                   should be specified as a cell:
%
%                           OutputChannel = {[Chnl_x,...,Chnl_y]}
%
%
%              In addition, various properties of the corresponding
%              IQC-multiplier can be specified: 
%
%               6.) 'BasisFunctionType' specifies the type of basis
%                   function (default = 1).
%
%               7.) 'Length' specifies the length of the basis function
%                   (default =  1).
%                                  
%               8.) 'PoleLocation' specifies the pole location of the basis
%                   function (default = -1).
%
%               9.) 'SampleTime' specifies the sampling time of the plant
%                   (default = 0, is continuous time)
%
%              10.) 'AddIQC' specifies whether one would like to add an
%                   additional LMI constraint at the cost of more
%                   computational complexity (default = 'yes').
%
%              11.) 'MstrictlyProp' specifies whether the part of the plant
%                   that is seen by the uncertainty (i.e. M) satisfies
%                   M(\infty) == 0. If so, the additional LMI constraint
%                   following from 9.) can be dropped allowing more freedom
%                   in the optimization problem  (default = 'no').
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
    NumberOfRepetitions double {mustBeReal,mustBeFinite,mustBePositive} = 1;

    % Delay operator type
    DelayType double {mustBeReal,mustBeFinite} = 1;
    
    % Add additional IQC constraint
    AddIQC string {mustBeMember(AddIQC,{'yes','no'})} = 'yes';
    
    % If the part of the plant seen by the uncertainty (i.e. M) satisfies
    % M(\infty) == 0 the one can set the property
    MstrictlyProp string {mustBeMember(MstrictlyProp,{'yes','no'})} = 'no';
    
    % Maximum time delay tau
    DelayTime double {mustBeReal,mustBeFinite} = 1;
    
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
    function obj = udel(varargin)
        if nargin == 0
            error('Error: To utilize "udel" one should specify at least one input argument (see "help udel" for further details)');
        elseif nargin == 1
            if ischar(varargin{1})
                obj.Name = varargin{1};
            else
                error('Error: The name of the uncertainty should be defined by a string (see "help udel" for further details)');
            end
        elseif nargin > 1
            if nargin == 2
                if ischar(varargin{1})
                    obj.Name = varargin{1};
                else
                    error('Error: The name of the uncertainty should be defined by a string (see "help udel" for further details)');
                end
                if isstruct(varargin{2})
                else
                    error('Error: The second input argument should be defined by a structure (see "help udel" for further details)');
                end
                if isfield(varargin{2},'BasisFunctionType')
                    a = 1:1:5;
                    if isreal(varargin{2}.BasisFunctionType) && isscalar(varargin{2}.BasisFunctionType) && ismember(varargin{2}.BasisFunctionType,a)
                        obj.BasisFunctionType = varargin{2}.BasisFunctionType;
                    else
                        error('Error: The type of basis function should be defined as Natural numbers in N1');
                    end
                end
                if isfield(varargin{2},'Length')
                    a = 1:1:10;
                    if isreal(varargin{2}.Length) && isscalar(varargin{2}.Length) && ismember(varargin{2}.Length,a)
                        obj.Length = varargin{2}.Length;
                    else
                        error('Error: The length of a basis function should be defined as a natural number in N1 (see "help udel" for further details).');
                    end
                end
                if isfield(varargin{2},'PoleLocation')
                    if isreal(varargin{2}.PoleLocation) && isscalar(varargin{2}.PoleLocation)
                        obj.PoleLocation = varargin{2}.PoleLocation;
                    else
                        error('Error: The pole location of a basis function should be defined as a real rational number. For continuous and distrete time systems, it is not allowed to have poles on the imaginary axis and unit circle respectively (see "help udel" for further details).');
                    end
                end
                if isfield(varargin{2},'SampleTime')
                    if isreal(varargin{2}.SampleTime) && ismatrix(varargin{2}.SampleTime)
                        obj.SampleTime = varargin{2}.SampleTime;
                    else
                        error('Error: The Sampling time of a basis function should be defined as a nonnegative real rational number (see "help udel" for further details).');
                    end
                end
                if isfield(varargin{2},'NumberOfRepetitions')
                    if isreal(varargin{2}.NumberOfRepetitions) && isscalar(varargin{2}.NumberOfRepetitions) && varargin{2}.NumberOfRepetitions > 0 
                        obj.NumberOfRepetitions = varargin{2}.NumberOfRepetitions;
                    else
                        error('Error: The number of repetitions of the uncertainty should be defined as a vector with positive natural numbers in N1 (see "help udel" for further details).');
                    end
                end
                if isfield(varargin{2},'DelayType')
                    a = [1,2];
                    if isreal(varargin{2}.DelayType) && isscalar(varargin{2}.DelayType) && ismember(varargin{2}.DelayType,a)
                        obj.DelayType = varargin{2}.DelayType;
                    else
                        error('Error: It is only possible to consider two different types of delay operators (i.e. Type 1 or Type 2) (see "help udel" for further details).');
                    end
                end
                if isfield(varargin{2},'DelayTime')
                    if isreal(varargin{2}.DelayTime) && isscalar(varargin{2}.DelayTime) && varargin{2}.DelayTime > 0
                        obj.DelayTime = varargin{2}.DelayTime;
                    else
                        error('Error: The maximum time delay should be defined by a positive scalar (see "help udel" for further details).');
                    end
                end
                if isfield(varargin{2},'AddIQC')
                    obj.AddIQC = varargin{2}.AddIQC;
                end
                if isfield(varargin{2},'MstrictlyProp')
                    obj.MstrictlyProp = varargin{2}.MstrictlyProp;
                end
                if isfield(varargin{2},'TerminalCost')
                    obj.TerminalCost   = varargin{2}.TerminalCost;
                end
                if isfield(varargin{2},'PrimalDual')
                    obj.PrimalDual = varargin{2}.PrimalDual;
                end
            elseif nargin > 2
                if ischar(varargin{1})
                    obj.Name = varargin{1};
                else
                    error('Error: The name of the uncertainty should be defined by a string (see "help udel" for further details)');
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
                            if isreal(varargin{j(i)}) && isscalar(varargin{j(i)}) && ismember(varargin{j(i)},a)
                                obj.BasisFunctionType = varargin{j(i)};
                            else
                                error('Error: The type of basis function should be defined as natural numbers in N1 (see "help udel" for further details)');
                            end
                        case 'Length'
                            a = 1:1:10;
                            if isreal(varargin{j(i)}) && isscalar(varargin{j(i)}) && ismember(varargin{j(i)},a)
                                obj.Length = varargin{j(i)};
                            else
                                error('Error: The length of a basis function should be defined as a natural number in N1 (see "help udel" for further details).');
                            end
                        case 'PoleLocation'
                            if isreal(varargin{j(i)}) && isscalar(varargin{j(i)}) && varargin{j(i)} ~= 0
                                obj.PoleLocation = varargin{j(i)};
                            else
                                error('Error: The pole location of a basis function should be defined as a real rational number. For continuous and distrete time systems, it is not allowed to have poles on the imaginary axis and unit circle respectively (see "help udel" for further details).');
                            end
                        case 'SampleTime'
                            if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)})
                                obj.SampleTime = varargin{j(i)};
                            else
                                error('Error: The sampling time of a basis function should be defined as a nonnegative real rational number (see "help udel" for further details).');
                            end
                        case 'NumberOfRepetitions'
                            if isreal(varargin{j(i)}) && isscalar(varargin{j(i)}) && varargin{j(i)} > 0
                                obj.NumberOfRepetitions = varargin{j(i)};
                            else
                                error('Error: The number of repetitions of the uncertainty should be defined as a vector with positive natural numbers in N1 (see "help udel" for further details).');
                            end
                        case 'DelayType'
                            a = [1,2];
                            if isreal(varargin{j(i)}) && isscalar(varargin{j(i)}) && ismember(varargin{j(i)},a)
                                obj.DelayType = varargin{j(i)};
                            else
                                error('Error: It is only possible to consider two different types of delay operators (i.e. Type 1 or Type 2) (see "help udel" for further details).');
                            end
                        case 'DelayTime'
                            if isreal(varargin{j(i)}) && isscalar(varargin{j(i)}) && varargin{j(i)} > 0
                                obj.DelayTime = varargin{j(i)};
                            else
                                error('Error: The maximum time delay should be defined by a positive scalar (see "help udel" for further details).');
                            end
                        case 'TerminalCost'
                            obj.TerminalCost = varargin{j(i)};
                        case 'AddIQC'
                            obj.AddIQC = varargin{j(i)};
                        case 'MstrictlyProp'
                            obj.MstrictlyProp = varargin{j(i)};
                        case 'PrimalDual'
                            obj.PrimalDual = varargin{j(i)};
                    end
                end    
            end
        end
    end
    function prob = iqcdel_p(Delta,prob)
    % ---------------------------------------------------------------------
    % This function generates two types of multipliers for the class udel.
    % These are precisely implemented as discussed in Section 5.5 of:
    %
    % [1] J. Veenman, C.W. Scherer, H. Koroglu, Robust stability and
    %     performance analysis based on integral quadratic constraints,
    %     European Journal of Control, Vol. 31, pp. 1-32, 2016.
    % ---------------------------------------------------------------------
    % Type 1: delta =     (e^(-s*tau)-1)*I_nr
    % ---------------------------------------------------------------------
    %
    % This function generates the multiplier factorization
    %
    %       [x]^*[ D, 0 ][ kron(Phi*W,I_nr), 0              ]
    %       [x]  [ 0,-D ][ 0               , kron(Phi,I_nr) ]
    %
    % with:
    % - (kron(Phi,I_nr))^* D(kron(Phi,I_nr)) > 0
    %
    % Here Phi\in RHinf is a basis function and W a weighting function that
    % ensures that ||Delta(i\omega)||<=|W(i\omega)| for all omega\in R and
    % for all time delays in [0,bet], with "bet" being the maximum
    % time-delay.
    %
    % If the property "AddIQC" is set to "yes" an additional IQC constraint
    % is added at the cost of more computational complexity. This yields
    % the multiplier factorization: 
    %
    %       [x]^*[ D, 0, 0,  0 ][ kron(Phi*W,I_nr), 0              ]
    %       [x]  [ 0, 0, 0,  R ][ kron(Phi,I_nr)  , 0              ]
    %       [x]  [ 0, 0, R,  0 ][ 0               , kron(Phi,I_nr) ]
    %       [x]  [ 0, R, 0, -D ][ 0               , kron(Phi,I_nr) ]
    %
    % with:
    % - (kron(Phi,I_nr))^* D(kron(Phi,I_nr)) > 0
    % - (kron(Phi,I_nr))^*(R-D)(kron(Phi,I_nr)) < 0
    %
    % The last constraint can be dropped if the part of the plant that is
    % seen by the uncertainty is strictly proper.
    %
    % ---------------------------------------------------------------------
    % Type 2: delta = 1/s*(e^(-s*tau)-1)*I_nr
    % ---------------------------------------------------------------------
    % This function generates the following multiplier factorization, which
    % forms a valid multiplier class for the delay operator of Type 2:
    %
    %       [x]^*[ D,  0 ][ bet*kron(Phi,I_nr), 0              ]
    %       [x]  [ 0, -D ][ 0                 , kron(Phi,I_nr) ]
    %
    % Here "bet" is the maximum time-delay, while Phi\in RHinf is the basis
    % function and D a matrix variable. The this multiplier parametrization
    % is valid for the uncertainty class if the following constraint holds:
    %
    % - (kron(Phi,I_nr))^* D(kron(Phi,I_nr)) > 0
    %
    % ---------------------------------------------------------------------
    
    bft           = get(Delta,'BasisFunctionType');
    nr            = get(Delta,'NumberOfRepetitions');
    l             = get(Delta,'Length');
    pl            = get(Delta,'PoleLocation');
    Ts            = get(Delta,'SampleTime');
    bet           = get(Delta,'DelayTime');
    xtrIQC        = get(Delta,'AddIQC');
    MstrP         = get(Delta,'MstrictlyProp');
    tc            = get(Delta,'TerminalCost');
    lnr           = length(nr);
    ll            = length(l);
    lpl           = length(pl);
    
    if lnr > 1 || ll > 1 || lpl > 1
        error('Error: The properties "Length", "PoleLocation", and "NumberOfRepetitions" should be defined as scalar inputs.');
    end
    if Ts == 0 && pl == 0
        disp('Warning: In continuous time it is not allowed to have poles on the imaginary axis. The "PoleLocation" was set to its nominal value (i.e. -1).');
        pl                = -1;
    elseif Ts > 0 && abs(pl) == 1
        disp('Warning: In discrete time it is not allowed to have poles on the unit circle. The "PoleLocation" was set to its nominal value (i.e. 0).');
        pl                = 0;
    elseif Ts == -1 && abs(pl) == 1
        disp('Warning: In discrete time it is not allowed to have poles on the unit circle. The "PoleLocation" was set to its nominal value (i.e. 0).');
        pl                = 0;
    end
    if Ts > 0
        q                 = ceil(bet/Ts);
        if q-bet/Ts > 0
            disp('The "(DelayTime/SampleTime)" was rounded off to the nearest integer');
            bet           = q;
        else
            bet           = bet/Ts;
        end
    end
    
    % Define KYP certificate structure (for continuous or discrete time analysis)
    if Ts == 0
        kyp               = [0,1;1,0];
    else
        kyp               = [-1,0;0,1];
    end
    if Delta.DelayType == 1
        Phi               = fBasis(l,pl,nr,bft,Ts);
        W                 = fW(bet,1e-6,3,Ts);
        W                 = fAugss(ss([],[],[],[],Ts),W,nr);
        if strcmp(xtrIQC,'no')
            Phi1          = Phi*W;
            Phi2          = Phi;
            prob          = fPsi(prob,Phi1,Phi2);

            % Define LMI variables
            a             = size(Phi.a,1);
            c             = size(Phi.c,1);
            D             = iqcvar(prob,[c,c],'symmetric');
            X             = iqcvar(prob,[a,a],'symmetric');
            P             = [D,zeros(c);zeros(c),-D];
            sc            = eye(2*c);
            IO            = [nr,nr,0,0,c,c,c,c,a,a];

            % Define LMI constraints
            V             = blkdiag(kron(kyp,X),D);
            A             = fss2m(Phi);
            prob          = iqclmi(prob,V,1,0,A);
            switch tc
                case 'on'
                    Z     = blkdiag(X,zeros(size(W.a,1)),-X);
                    Tz    = eye(Z.Dim(1));
                    prob  = fZ(prob,Z,Tz);
            end
        elseif strcmp(xtrIQC,'yes')
            Phi1          = blkdiag(Phi,Phi)*[W;eye(size(W,2))];
            Phi2          = [eye(size(Phi,1));eye(size(Phi,1))]*Phi;
            prob          = fPsi(prob,Phi1,Phi2);
        
            % Define LMI variables
            a             = size(Phi.a,1);
            c             = size(Phi.c,1);
            D             = iqcvar(prob,[c,c],'symmetric');
            R             = iqcvar(prob,[c,c],'symmetric');
            X             = iqcvar(prob,[a,a],'symmetric');
            P11           = blkdiag(D,zeros(c));
            P12           = blkdiag(zeros(c),R);
            P22           = blkdiag(R,-D);
            P             = [P11,P12;P12.',P22];
            sc            = eye(4*c);
            IO            = [nr,nr,0,0,2*c,2*c,2*c,2*c,size(Phi1.a,1),size(Phi2.a,1)];
        
            % Define LMI constraints
            V             = blkdiag(kron(kyp,X),D);
            A             = fss2m(Phi);
            prob          = iqclmi(prob,V,1,0,A);
        
            if strcmp(MstrP,'no')
                Y         = iqcvar(prob,[a,a],'symmetric');
                Q         = blkdiag(kron(kyp,Y),-D,R);
                B         = fss2m(Phi2);
                prob      = iqclmi(prob,Q,-1,0,B);
            end
            switch tc
                case 'on'
                    xPhi  = size(Phi1.a,1) + size(Phi2.a,1);
                    Z     = iqcvar(prob,[xPhi,xPhi],'symmetric');
                    prob  = fZ(prob,Z,eye(Z.Dim(1)));
                    nPhi1 = size(Phi1,2);
                    nPhi2 = size(Phi2,2);
                    B     = fss2m(blkdiag(Phi1,Phi2)*[zeros(nPhi1,nPhi2);eye(nPhi2)]);
                    W     = blkdiag(kron(kyp,Z),P);
                    prob  = iqclmi(prob,W,-1,0,B);
            end
        end
    elseif Delta.DelayType == 2
        % Define basis functions
        Phi               = fBasis(l,pl,nr,bft,Ts);
        prob              = fPsi(prob,Phi,Phi);
        
        % Define LMI variables
        a                 = size(Phi.a,1);
        c                 = size(Phi.c,1);
        D                 = iqcvar(prob,[c,c],'symmetric');
        X                 = iqcvar(prob,[a,a],'symmetric');        
        P                 = blkdiag(D,-D);
        sc                = blkdiag(bet*eye(c),eye(c));
        IO                = [nr,nr,0,0,c,c,c,c,a,a];
        
        % Define LMI constraint
        V                 = blkdiag(kron(kyp,X),D);
        A                 = fss2m(Phi);
        prob              = iqclmi(prob,V,1,0,A);
        switch tc
            case 'on'
                Z         = blkdiag(X,-X);
                Tz        = blkdiag(bet*eye(X.Dim(1)),eye(X.Dim(1)));
                prob      = fZ(prob,Z,Tz);
        end
    else
        error('Error: The property "DelayType" is not assigned.');
    end
    prob                  = fP(prob,P);
    prob                  = fsc(prob,sc);
    prob                  = fIO(prob,IO);
    end
    function prob = iqcdel_d(Delta,prob)
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
        disp('------------------------------------------------');
        disp(' LTI diagonally repeated delay uncertainty block');
        disp('------------------------------------------------');
        disp(' Uncertainty details:');
        disp(['  - Name:                           ',num2str(obj.Name)]);
        disp(['  - Delay operator type:            ',num2str(obj.DelayType)]);
        disp(['  - Maximum time delay:             ',num2str(obj.DelayTime)]);
        disp(['  - Number of repetitions:          ',num2str(obj.NumberOfRepetitions)]);
        disp(['  - Input channels plant:           ',a1]);
        disp(['  - Output channels plant:          ',b1]);
        disp(['  - Plant strictly proper:          ',num2str(obj.MstrictlyProp)]);
        disp(' ');
        disp(' Multiplier details:');
        disp(['  - Basis function type:            ',num2str(obj.BasisFunctionType)]);
        disp(['  - Basis function length:          ',num2str(obj.Length)]);
        disp(['  - Basis function pole location:   ',num2str(obj.PoleLocation)]);
        disp(['  - Add Extra Constraint:           ',num2str(obj.AddIQC)]);
        disp(['  - Multiplier type:                ',num2str(obj.PrimalDual)]);
    end
end
end