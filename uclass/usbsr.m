classdef usbsr < handle & matlab.mixin.SetGetExactNames
    
% -------------------------------------------------------------------------
%
% Copyright:   This is copyrighted material owned by Novantec B.V.
% Terms:       IQClab is available under a Creative Commons
%              (Attribution-NoDerivatives 4.0 International (CC BY-ND 4.0))
%              license: https://creativecommons.org/licenses/by-nd/4.0/
%              For further information please visit iqclab.eu
%
% Author:      J.Veenman
% Date:        02-12-2019
% 
% -------------------------------------------------------------------------
%
% Description: Class definition for diagonally repeated sector bounded and
%              slope restricted nonlinearities of the form:
%
%              p = Delta(q)=col(delta(q1) , ... , delta(qN)), 
%
%              which satisfy:
%
%               1.) The sector constraint:
%
%                 (delta(qi)-a*qi)*(b*qi-delta(qi)) >= 0, i\in{1,...,N}
%
%                   for all qi\in L2 and with fixed constants a <= 0 <= b
%
%               2.) And optionally, the incremental sector constraint
%                   (slope restriction): 
%
%                   c <= (delta(x1)-delta(x2))/(x1-x2) <= d
%
%                   for all x1,x2\in R, x1\neq x2 and with fixed constants
%                   c <= 0 <= d.
%
%              Note: delta\in slope(c,d) -> delta\in sector(c,d).
%              Thus, if desired, one can consider tighter sector bounds: 
%              b < d and c < a (but not the other way around.
%
% Syntax:      Delta = usbsr('name')
%              Delta = usbsr('name',varargin)
%
% Usage:       "Delta = usbsr('name')" defines a sector bounded
%              nonlinearity, which is repeated once and which satisfies the
%              sector constraint with [a,b] = [0,1];
%
%              For "delta = usbsr('name',varargin)", the varargin inputs
%              come in pairs and can be defined as: 
%
%              delta = usbsr('name','prop1','value1','prop2','value2',...)
%
%              Alternatively, the properties can be set by defining the
%              structure:
%
%              usbsrOpt.prop1 = value1
%              usbsrOpt.prop2 = value2
%                       ...              
%              usbsrOpt.propN = valueN
% 
%              Delta = usbsr('name',usbsrOpt)
%
%              Finally, the properties can be set and retrieved by using:
%
%              set(Delta,'propX','valueX') and get(Delta,'propX').
%
%              The properties that can be specified are:
%
%               1.) 'NumberOfRepetitions' specifies the number of
%                   repetitions of delta (default =  1).
%            
%               2.) 'SectorBounds' specifies the sector constraint as
%                   (default: [0,1]):
%
%                             SectorBounds = [a,b]
%
%               3.) 'SlopeBounds' specifies the slope constraints as
%                   (default: [0,0]):
%
%                              SlopeBounds = [c,d]
%
%               3.) 'Odd' specifies if delta is odd, i.e. delta(-qi) =
%                   -delta(qi) ('yes' or 'no') (default: 'no');
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
%              10.) 'TerminalCost' defines the terminal cost constraint, Z,
%                   which is used in invariance analysis. Options are
%
%                    a.) 'off' (default)
%                    a.) 'on'
%
%              11.) 'PrimalDual' specifies whether the multiplier should be
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
   
    % Norm alpha of the uncertainty delta
    SectorBounds double {mustBeReal,mustBeFinite} = [0,1];
    
    % Norm alpha of the uncertainty delta
    SlopeBounds double {mustBeReal,mustBeFinite} = [0,0];
    
    %  Is the nonlinearity odd or even
    Odd string {mustBeMember(Odd,{'yes','no',''})} = 'no';

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
    function obj = usbsr(varargin)
        if nargin == 0
            error('Error: To utilize "usbsr" one should specify at least one input argument (see "help usbsr" for further details)');
        elseif nargin == 1
            if ischar(varargin{1})
                obj.Name = varargin{1};
            else
                error('Error: The name of the uncertainty should be defined by a string (see "help usbsr" for further details)');
            end
        elseif nargin > 1
            if nargin == 2
                if ischar(varargin{1})
                    obj.Name = varargin{1};
                else
                    error('Error: The name of the uncertainty should be defined by a string (see "help usbsr" for further details)');
                end
                if isstruct(varargin{2})
                else
                    error('Error: The second input argument should be defined by a structure (see "help usbsr" for further details)');
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
                        error('Error: The length of a basis function should be defined as a natural number in N1 (see "help usbsr" for further details).');
                    end
                end
                if isfield(varargin{2},'PoleLocation')
                    if isreal(varargin{2}.PoleLocation) && isscalar(varargin{2}.PoleLocation)
                        obj.PoleLocation = varargin{2}.PoleLocation;
                    else
                        error('Error: The pole location of a basis function should be defined as a real rational number. For continuous and distrete time systems, it is not allowed to have poles on the imaginary axis and unit circle respectively (see "help usbsr" for further details).');
                    end
                end
                if isfield(varargin{2},'SampleTime')
                    if isreal(varargin{2}.SampleTime) && ismatrix(varargin{2}.SampleTime)
                        obj.SampleTime = varargin{2}.SampleTime;
                    else
                        error('Error: The Sampling time of a basis function should be defined as a nonnegative real rational number (see "help usbsr" for further details).');
                    end
                end
                if isfield(varargin{2},'NumberOfRepetitions')
                    if isreal(varargin{2}.NumberOfRepetitions) && isscalar(varargin{2}.NumberOfRepetitions) && varargin{2}.NumberOfRepetitions > 0 
                        obj.NumberOfRepetitions = varargin{2}.NumberOfRepetitions;
                    else
                        error('Error: The number of repetitions of the uncertainty should be defined as a vector with positive natural numbers in N1 (see "help usbsr" for further details).');
                    end
                end
                if isfield(varargin{2},'SectorBounds')
                    a = varargin{2}.SectorBounds;
                    if isreal(a) && ismatrix(a)
                        [s1,s2] = size(a);
                        if s2 ~= 2
                           error('The option SectorBounds should be specified as a n by 2 matrix.');
                        end
                        for ijk = 1:s1
                            if a(ijk,1) > 0 || a(ijk,2) < 0
                                error('Sector bounds should satisfy ai <= 0 <= bi.');
                            end
                        end
                        obj.varargin{2}.SectorBounds;
                    else
                        error('Error: The size of "SectorBounds" should be defined as [n_out,n_in], with n_out and n_in being positive natural numbers in N1 (see "help usbsr" for further details).');
                    end
                end
                if isfield(varargin{2},'SlopeBounds')
                    a = varargin{2}.SlopeBounds;
                    if isreal(a) && ismatrix(a)
                        [s1,s2] = size(a);
                        if s2 ~= 2
                           error('The option SlopeBounds should be specified as a n by 2 matrix.');
                        end
                        for ijk = 1:s1
                            if a(ijk,1) > 0 || a(ijk,2) < 0
                                error('Slope bounds should satisfy ai <= 0 <= bi.');
                            end
                        end
                        obj.varargin{2}.SlopeBounds;
                    else
                        error('Error: The size of "SlopeBounds" should be defined as [n_out,n_in], with n_out and n_in being positive natural numbers in N1 (see "help usbsr" for further details).');
                    end
                end
                if isfield(varargin{2},'Odd')
                    obj.Odd = varargin{2}.Odd;
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
                    error('Error: The name of the uncertainty should be defined by a string (see "help usbsr" for further details)');
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
                                error('Error: The type of basis function should be defined as natural numbers in N1 (see "help usbsr" for further details)');
                            end
                        case 'Length'
                            a = 1:1:10;
                            if isreal(varargin{j(i)}) && isscalar(varargin{j(i)}) && ismember(varargin{j(i)},a)
                                obj.Length = varargin{j(i)};
                            else
                                error('Error: The length of a basis function should be defined as a natural number in N1 (see "help usbsr" for further details).');
                            end
                        case 'PoleLocation'
                            if isreal(varargin{j(i)}) && isscalar(varargin{j(i)}) && varargin{j(i)} ~= 0
                                obj.PoleLocation = varargin{j(i)};
                            else
                                error('Error: The pole location of a basis function should be defined as a real rational number. For continuous and distrete time systems, it is not allowed to have poles on the imaginary axis and unit circle respectively (see "help usbsr" for further details).');
                            end
                        case 'SampleTime'
                            if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)})
                                obj.SampleTime = varargin{j(i)};
                            else
                                error('Error: The sampling time of a basis function should be defined as a nonnegative real rational number (see "help usbsr" for further details).');
                            end
                        case 'NumberOfRepetitions'
                            if isreal(varargin{j(i)}) && isscalar(varargin{j(i)}) && varargin{j(i)} > 0
                                obj.NumberOfRepetitions = varargin{j(i)};
                            else
                                error('Error: The number of repetitions of the uncertainty should be defined as a vector with positive natural numbers in N1 (see "help usbsr" for further details).');
                            end
                        case 'SectorBounds'
                            if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)})
                                a = varargin{j(i)};
                                [s1,s2] = size(a);
                                if s2 ~= 2
                                   error('The option "SectorBounds" should be specified as a n by 2 matrix.');
                                end
                                for ijk = 1:s1
                                    if a(ijk,1) > 0 || a(ijk,2) < 0
                                        error('Sector bounds should satisfy ai <= 0 <= bi.');
                                    end
                                end
                                obj.SectorBounds = varargin{j(i)};
                            else
                                error('Error: The size of "SectorBounds" should be defined as [n_out,n_in], with n_out and n_in being positive natural numbers in N1 (see "help usbsr" for further details).');
                            end
                        case 'SlopeBounds'
                            if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)})
                                a = varargin{j(i)};
                                [s1,s2] = size(a);
                                if s2 ~= 2
                                   error('The option "SlopeBounds" should be specified as a n by 2 matrix.');
                                end
                                for ijk = 1:s1
                                    if a(ijk,1) > 0 || a(ijk,2) < 0
                                        error('Slope bounds should satisfy ai <= 0 <= bi.');
                                    end
                                end
                                obj.SlopeBounds = varargin{j(i)};
                            else
                                error('Error: The size of "Slope" should be defined as [n_out,n_in], with n_out and n_in being positive natural numbers in N1 (see "help usbsr" for further details).');
                            end
                        case 'TerminalCost'
                            obj.TerminalCost = varargin{j(i)};
                        case 'PrimalDual'
                            obj.PrimalDual = varargin{j(i)};
                        case 'Odd'
                            obj.Odd = varargin{j(i)};
                        case 'InputChannel'
                            obj.InputChannel = varargin{j(i)};
                        case 'OutputChannel'
                            obj.OutputChannel = varargin{j(i)};
                    end
                end    
            end
        end
    end
    function prob = iqcsbsr_p(Delta,prob)
    % ---------------------------------------------------------------------
    % This function generates 2 IQC-multipliers for the class usbsr:
    %
    % 1.) Full-block sector criterion for sector-bounded nonlinearities
    % 2.) Full-block sector criterion + full-block Zames-Falb multipliers
    %     for sector-bounded and slope restricted nonlinearities.
    %
    % These multiplieres are precisely implemented as described in Section
    % 5.8 of [1] and Section 4.3 of [2] respectively. The discrete time
    % version is addressed in [3].
    %
    % [1] J. Veenman, C.W. Scherer, H. Koroglu, Robust stability and
    %     performance analysis based on integral quadratic constraints,
    %     European Journal of Control, Vol. 31, pp. 1-32, 2016.
    % [2] M. Fetzer, C.W. Scherer, Full-block multipliers for repeated
    %     slope restricted scalar nonlinearities, International Journal of
    %     Robust and Nonlinear Control, Vol. 27, pp. 3376-3411, 2017.
    % [3] M. Fetzer, C. W. Scherer, ?Absolute stability analysis of
    %     discrete time feedback interconnections?, Proceedings of the 20th
    %     IFAC World Congress, pp. 8447-8453, 2017.
    %
    % ---------------------------------------------------------------------

    % Retrieve uncertaity info
    nr                                       = get(Delta,'NumberOfRepetitions');
    OE                                       = get(Delta,'Odd');
    l                                        = get(Delta,'Length');
    pl                                       = get(Delta,'PoleLocation');
    Ts                                       = get(Delta,'SampleTime');
    sect                                     = get(Delta,'SectorBounds');
    slope                                    = get(Delta,'SlopeBounds');
    tc                                       = get(Delta,'TerminalCost');
    lnr                                      = length(nr);
    ll                                       = length(l);
    lpl                                      = length(pl);
    
    if lnr > 1 || ll > 1 || lpl > 1
        error('Error: The properties "Length", "PoleLocation", and "NumberOfRepetitions" should be defined as scalar inputs.');
    end
    if sect(1) > 0 || sect(1) > sect(2)
        error('Error: The sector should satisfy: SectorBounds(1) <= 0 <= SectorBounds(2)');
    end
    if slope(1) > 0 || slope(1) > slope(2)
        error('Error: The sector should satisfy: SlopeBounds(1) <= 0 <= SlopeBounds(2)');
    end
    if Ts == 0 && pl == 0
        disp('Warning: In continuous time it is not allowed to have poles on the imaginary axis. The "PoleLocation" was set to its nominal value (i.e. -1).');
        pl                                   = -1;
    end
    if slope(1) == 0 && slope(2) == 0
        sw                                   = 'sector';
    elseif slope(1) ~= 0 || slope(2) ~= 0
        sw                                   = 'slope';
        if sect(1) == 0 && sect(2) == 0
            sect                             = slope;
            disp('Warning: If the sector upper-bound is not specified, it is set to be equal to the slope upper-bound.');
        end
    end
    % Define KYP certificate structure (for continuous or discrete time analysis)
    if Ts == 0
        kyp                                  = [0,1;1,0];
    else
        kyp                                  = [-1,0;0,1];
    end
    switch sw
        case 'sector'
            % Create outer factor (basis function) of the uncertainty
            prob                             = fPsi(prob,ss([],[],[],eye(nr),Ts),ss([],[],[],eye(nr),Ts));
            
            % Define LMI Variables
            P11                              = iqcvar(prob,[nr,nr],'symmetric');
            P12                              = iqcvar(prob,[nr,nr],'full');
            P22                              = iqcvar(prob,[nr,nr],'symmetric');
            P                                = [P11,P12;P12.',P22];
            sc                               = eye(2*nr);
            prob                             = fP(prob,P);
            prob                             = fsc(prob,sc);
            
            % Define LMI constraints
%             prob                             = iqclmi(prob,P22,-1);           % Convex hull relaxation
            prob                             = iqclmi(prob,diag(diag(P22)),-1); % Partial convexity relaxation
            
            po                               = polydec(pvec('box',ones(nr,1)*sect))';
            for i = 1:length(po)
                A                            = [eye(nr);diag(po(i,:))];
                prob                         = iqclmi(prob,P,1,0,A);
            end

            % Provide input output data
            IO                               = [nr,nr,0,0,nr,nr,nr,nr,0,0];
            prob                             = fIO(prob,IO);
        case 'slope'
            if Ts == 0
                % Create outer factor (basis function) of the uncertainty
                Phi                          = [ss([],[],[],eye(nr),Ts);fBasis(l,pl,nr,3,Ts)];
                prob                         = fPsi(prob,Phi,Phi);

                % Define LMI Variables
                Q                            = iqcvar(prob,[nr,nr],'symmetric');
                S                            = iqcvar(prob,[nr,nr],'full');
                R                            = iqcvar(prob,[nr,nr],'symmetric');
                C0                           = iqcvar(prob,[nr,nr],'full');
                if l == 1
                    P11                      = blkdiag(Q,zeros(nr));
                    P12                      = blkdiag(S,C0.');
                    P22                      = blkdiag(R,zeros(nr));
                    P                        = [P11,P12;P12.',P22];
                    sc                       = [eye(nr),zeros(nr,3*nr);zeros(nr),slope(2)*eye(nr),zeros(nr),-eye(nr);...
                                                zeros(nr,2*nr),eye(nr),zeros(nr);zeros(nr),-slope(1)*eye(nr),zeros(nr),eye(nr)];
                    prob                     = fP(prob,P);
                    prob                     = fsc(prob,sc);

                    % Define LMI sector constraints
%                     prob                     = iqclmi(prob,R,-1);             % Convex hull relaxation
                    prob                     = iqclmi(prob,diag(diag(R)),-1); % Partial convexity relaxation
                    Psect                    = [Q,S;S.',R];

                    po                       = polydec(pvec('box',ones(nr,1)*sect))';
                    for i = 1:length(po)
                        A                    = [eye(nr);diag(po(i,:))];
                        prob                 = iqclmi(prob,Psect,1,0,A);
                    end
                    if strcmp(OE,'no')
                        for i = 1:nr
                            for j = 1:nr
                                if i ~= j
                                    prob     = iqclmi(prob,C0(i,j),-1);
                                end
                            end
                        end
                    end
                elseif l > 1
                    % Define basis function slope constraint
                    Phit                     = fBasis(l,pl,1,2);

                    % Define common LMI variables
                    C1                       = iqcvar(prob,[nr,nr*(l-1)],'full');
                    C3                       = iqcvar(prob,[nr,nr*(l-1)],'full');

                    if strcmp(OE,'yes')
                        % Define LMI variables
                        C2                   = iqcvar(prob,[nr,nr*(l-1)],'full');
                        C4                   = iqcvar(prob,[nr,nr*(l-1)],'full');
                        P11                  = blkdiag(Q,zeros(nr+2*nr*(l-1)));
                        P12                  = blkdiag(S,[C0.',C2,C1;C4.',zeros(nr*(l-1),2*nr*(l-1));C3.',zeros(nr*(l-1),2*nr*(l-1))]);
                        P22                  = blkdiag(R,zeros(nr+2*nr*(l-1)));
                        sca                  = kron(eye(2),blkdiag(eye(2*nr),[eye(nr*(l-1));-eye(nr*(l-1))]));

                        % Define LMI constraints: Ci*phi(t) > 0, i = 1,...,4
                        Rv                   = sqrt(diag(factorial(1:1:l-1)))^-1;
                        J                    = kron(ones(l-1,1),eye(nr))';
                        A                    = fss2m(Rv*Phit);
                        for i = 1:nr
                            for j = 1:nr
                                X1           = iqcvar(prob,[l-2,l-2],'symmetric');
                                X2           = iqcvar(prob,[l-2,l-2],'symmetric');
                                X3           = iqcvar(prob,[l-2,l-2],'symmetric');
                                X4           = iqcvar(prob,[l-2,l-2],'symmetric');
                                V1           = blkdiag(oblkdiag(X1),diag(C1(i,J(j,:) > 0)));
                                V2           = blkdiag(oblkdiag(X2),diag(C2(i,J(j,:) > 0)));
                                V3           = blkdiag(oblkdiag(X3),diag(C3(i,J(j,:) > 0)));
                                V4           = blkdiag(oblkdiag(X4),diag(C4(i,J(j,:) > 0)));
                                prob         = iqclmi(prob,V1,1,0,A);
                                prob         = iqclmi(prob,V2,1,0,A);
                                prob         = iqclmi(prob,V3,1,0,A);
                                prob         = iqclmi(prob,V4,1,0,A);
                            end
                        end

                        % Define LMI constraints: ||H||_1 < C0
                        A                    = [0.5;-ones(4*nr*(l-1),1);ones(nr,1)];
                        if nr == 1
                            V                = blkdiag(oblkdiag(-[C1,C2,C3,C4]),-C0);
                            prob             = iqclmi(prob,V,-1,0,A);
                        else
                            for i = 1:nr
                                L1           = [];
                                L2           = [];
                                M1           = [];
                                M2           = [];
                                for j = 1:nr
                                    L1       = [L1,[C1(i,J(j,:) > 0),C2(i,J(j,:) > 0),C3(i,J(j,:) > 0),C4(i,J(j,:) > 0)]];
                                    M1       = [M1,[C1(j,J(i,:) > 0),C2(j,J(i,:) > 0),C3(j,J(i,:) > 0),C4(j,J(i,:) > 0)]];
                                    if i ~= j
                                        SL2  = iqcvar(prob,[1,1],'symmetric');
                                        SM2  = iqcvar(prob,[1,1],'symmetric');
                                        L2   = blkdiag(L2,SL2);
                                        M2   = blkdiag(M2,SM2);
                                        prob = iqclmi(prob,[SL2,C0(i,j);C0(i,j).',SL2],1);
                                        prob = iqclmi(prob,[SM2,C0(j,i);C0(j,i).',SM2],1);
                                    end
                                end
                                L2           = blkdiag(L2,-C0(i,i));
                                M2           = blkdiag(M2,-C0(i,i));
                                V            = blkdiag(oblkdiag(-L1),L2);
                                W            = blkdiag(oblkdiag(-M1),M2);
                                prob         = iqclmi(prob,V,-1,0,A);
                                prob         = iqclmi(prob,W,-1,0,A);
                            end
                        end
                    elseif strcmp(OE,'no')
                        P11                  = blkdiag(Q,zeros(nr*l));
                        P12                  = blkdiag(S,[C0.',-C1;-C3.',zeros(nr*(l-1))]);
                        P22                  = blkdiag(R,zeros(nr*l));
                        sca                  = eye(2*(nr+nr*l));

                        % Define LMI constraints: Ci*phi(t) > 0, i = 1,...,4
                        Rv                   = sqrt(diag(factorial(1:1:l-1)))^-1;
                        J                    = kron(ones(l-1,1),eye(nr))';
                        A                    = fss2m(Rv*Phit);
                        for i = 1:nr
                            for j = 1:nr
                                X1           = iqcvar(prob,[l-2,l-2],'symmetric');
                                X3           = iqcvar(prob,[l-2,l-2],'symmetric');
                                V1           = blkdiag(oblkdiag(X1),diag(C1(i,J(j,:) > 0)));
                                V3           = blkdiag(oblkdiag(X3),diag(C3(i,J(j,:) > 0)));
                                prob         = iqclmi(prob,V1,1,0,A);
                                prob         = iqclmi(prob,V3,1,0,A);
                            end
                        end

                        % Define LMI constraints: Gij < 0 for all i ~= j
                        for i = 1:nr
                            for j = 1:nr
                                if i ~= j
                                    prob     = iqclmi(prob,C0(i,j),-1);
                                end
                            end
                        end

                        % Define LMI constraints: ||H||_1 < C0
                        A                    = [0.5;-ones(2*nr*(l-1),1);ones(nr,1)];
                        if nr == 1
                            V                = blkdiag(oblkdiag(-[C1,C3]),-C0);
                            prob             = iqclmi(prob,V,-1,0,A);
                        else
                            for i = 1:nr
                                L1           = [];
                                L2           = [];
                                M1           = [];
                                M2           = [];
                                for j = 1:nr
                                    L1       = [L1,[C1(i,J(j,:) > 0),C3(i,J(j,:) > 0)]];
                                    M1       = [M1,[C1(j,J(i,:) > 0),C3(j,J(i,:) > 0)]];
                                    if i ~= j
                                        SL2  = iqcvar(prob,[1,1],'symmetric');
                                        SM2  = iqcvar(prob,[1,1],'symmetric');
                                        L2   = blkdiag(L2,SL2);
                                        M2   = blkdiag(M2,SM2);
                                        prob = iqclmi(prob,[SL2,C0(i,j);C0(i,j).',SL2],1);
                                        prob = iqclmi(prob,[SM2,C0(j,i);C0(j,i).',SM2],1);
                                    end
                                end
                                L2           = blkdiag(L2,-C0(i,i));
                                M2           = blkdiag(M2,-C0(i,i));
                                V            = blkdiag(oblkdiag(-L1),L2);
                                W            = blkdiag(oblkdiag(-M1),M2);
                                prob         = iqclmi(prob,V,-1,0,A);
                                prob         = iqclmi(prob,W,-1,0,A);
                            end
                        end
                    end
                    P                        = [P11,P12;P12.',P22];
                    prob                     = fP(prob,P);
                    sc11b                    = blkdiag(eye(nr),slope(2)*eye(l*nr));
                    sc12b                    = blkdiag(zeros(nr),-eye(l*nr));
                    sc21b                    = blkdiag(zeros(nr),-slope(1)*eye(l*nr));
                    sc22b                    = eye(nr*(1+l));
                    scb                      = [sc11b,sc12b;sc21b,sc22b];
                    sc                       = sca*scb;
                    prob                     = fsc(prob,sc);
                    
                    % Define LMI sector constraints
%                     prob                     = iqclmi(prob,R,-1);             % Convex hull relaxation
                    prob                     = iqclmi(prob,diag(diag(R)),-1); % Partial convexity relaxation
                    Psect                    = [Q,S;S.',R];
                    
                    po                       = polydec(pvec('box',ones(nr,1)*sect))';
                    for i = 1:length(po)
                        A                    = [eye(nr);diag(po(i,:))];
                        prob                 = iqclmi(prob,Psect,1,0,A);
                    end
                    switch tc
                        case 'on'
                            Z                = iqcvar(prob,2*[size(Phi.a,1),size(Phi.a,1)],'symmetric');
                            prob             = fZ(prob,Z,eye(Z.Dim(1)));
                            B                = fss2m(sc*blkdiag(Phi,Phi)*[zeros(size(Phi,2));eye(size(Phi,2))]);
                            W                = blkdiag(kron(kyp,Z),P);
                            prob             = iqclmi(prob,W,-1,0,B);
                    end
                end
            else
                % Create outer factor (basis function) of the uncertainty
                Phi                          = [ss([],[],[],eye(nr),Ts);fBasis(l,0,nr,101,Ts)];
                prob                         = fPsi(prob,Phi,Phi);

                % Define LMI Variables
                Q                            = iqcvar(prob,[nr,nr],'symmetric');
                S                            = iqcvar(prob,[nr,nr],'full');
                R                            = iqcvar(prob,[nr,nr],'symmetric');
                M0                           = iqcvar(prob,[nr,nr],'full');
                if l == 1
                    P11                      = blkdiag(Q,zeros(nr));
                    P12                      = blkdiag(S,M0.');
                    P22                      = blkdiag(R,zeros(nr));
                    P                        = [P11,P12;P12.',P22];
                    sc                       = [eye(nr),zeros(nr,3*nr);zeros(nr),slope(2)*eye(nr),zeros(nr),-eye(nr);...
                                                zeros(nr,2*nr),eye(nr),zeros(nr);zeros(nr),-slope(1)*eye(nr),zeros(nr),eye(nr)];
                    prob                     = fP(prob,P);
                    prob                     = fsc(prob,sc);

                    % Define LMI sector constraints
%                     prob                     = iqclmi(prob,R,-1);             % Convex hull relaxation
                    prob                     = iqclmi(prob,diag(diag(R)),-1); % Partial convexity relaxation
                    Psect                    = [Q,S;S.',R];

                    po                       = polydec(pvec('box',ones(nr,1)*sect))';
                    for i = 1:length(po)
                        A                    = [eye(nr);diag(po(i,:))];
                        prob                 = iqclmi(prob,Psect,1,0,A);
                    end

                    e                        = ones(nr,1);
                    if strcmp(OE,'no')
                        for i = 1:nr
                            % Mij\leq0 for i\neq j
                            for j = 1:nr
                                if i ~= j
                                    prob     = iqclmi(prob,M0(i,j),-1);
                                end
                            end
                            % Me\geq0, e^TM\geq0
                            prob             = iqclmi(prob,diag(M0(i,:)),1,0,e);
                            prob             = iqclmi(prob,diag(M0(:,i)),1,0,e);
                        end
                    elseif strcmp(OE,'yes')
                        V                    = iqcvar(prob,[nr*(nr-1),1],'full');
                        W                    = diag(diag(M0));
                        cnt                  = 1;
                        for i = 1:nr
                            for j = 1:nr
                                if i ~= j
                                    W(i,j)   = V(cnt,1);
                                    cnt      = cnt + 1;
                                end
                            end    
                        end
                        for i = 1:nr 
                            T                = -eye(nr);
                            T(i,i)           = -T(i,i);
                            prob             = iqclmi(prob,T*diag(W(i,:)),1,0,e);
                            prob             = iqclmi(prob,T*diag(W(:,i)),1,0,e);
                        end
                        for i = 1:nr
                            for j = 1:nr
                                if i~=j
                                   prob      = iqclmi(prob,blkdiag(-M0(i,j),-W(i,j)),-1,0,[1;1]);
                                   prob      = iqclmi(prob,blkdiag( M0(i,j),-W(i,j)),-1,0,[1;1]);
                                end
                            end
                        end
                    end
                elseif l > 1
                    P11                      = blkdiag(Q,zeros(nr*l));
                    Mp                       = [];
                    Mm                       = [];
                    for i = 1:l-1
                        Mpi{i}               = iqcvar(prob,[nr,nr],'full');
                        Mmi{i}               = iqcvar(prob,[nr,nr],'full');
                        Mp                   = [Mp',Mpi{i}']';
                        Mm                   = [Mm,Mmi{i}];
                    end
                    P12                      = blkdiag(S,[M0,Mm;Mp,zeros((l-1)*nr)].');
                    P22                      = blkdiag(R,zeros(nr*l));
                    P                        = [P11,P12;P12.',P22];
                    sc11                     = blkdiag(eye(nr),slope(2)*eye(l*nr));
                    sc12                     = blkdiag(zeros(nr),-eye(l*nr));
                    sc21                     = blkdiag(zeros(nr),-slope(1)*eye(l*nr));
                    sc22                     = eye(nr*(1+l));
                    sc                       = [sc11,sc12;sc21,sc22];
                    prob                     = fP(prob,P);
                    prob                     = fsc(prob,sc);

                    % Define LMI sector constraints
%                     prob                     = iqclmi(prob,R,-1);             % Convex hull relaxation
                    prob                     = iqclmi(prob,diag(diag(R)),-1); % Partial convexity relaxation
                    Psect                    = [Q,S;S.',R];
                    
                    po                       = polydec(pvec('box',ones(nr,1)*sect))';
                    for i = 1:length(po)
                        A                    = [eye(nr);diag(po(i,:))];
                        prob                 = iqclmi(prob,Psect,1,0,A);
                    end

                    Mv                       = M0;
                    Mh                       = M0;
                    for k = 1:l-1
                        Mv                   = [Mpi{k}',Mv',Mmi{k}']';
                        Mh                   = [Mpi{k},Mh,Mmi{k}];
                    end
                    e                        = ones(2*nr*(l-1)+nr,1);
                    if strcmp(OE,'no')
                        % M0ij\leq0 for i\neq j
                        for k = 1:nr
                            for j = 1:nr
                                if j ~= k
                                    prob     = iqclmi(prob,M0(k,j),-1);
                                end
                            end
                        end
                        for i = 1:l-1
                            for j = 1:nr
                                for k = 1:nr
                                    prob     = iqclmi(prob,Mpi{i}(k,j),-1);
                                    prob     = iqclmi(prob,Mmi{i}(k,j),-1);
                                end
                            end
                        end
                        % Me\geq0, e^TM\geq0
                        for k = 1:nr
                            prob             = iqclmi(prob,diag(Mv(:,k)),1,0,e);
                            prob             = iqclmi(prob,diag(Mh(k,:)),1,0,e);
                        end
                    elseif strcmp(OE,'yes')
                        V                    = iqcvar(prob,[nr*(nr-1),1],'full');
                        W                    = diag(diag(M0));
                        cnt                  = 1;
                        for i = 1:nr
                            for j = 1:nr
                                if i ~= j
                                    W(i,j)   = V(cnt,1);
                                    cnt      = cnt + 1;
                                end
                            end
                        end
                        Wv                   = W;
                        Wh                   = W;
                        for i = 1:l-1
                            Wpi{i}           = iqcvar(prob,[nr,nr],'full');
                            Wmi{i}           = iqcvar(prob,[nr,nr],'full');
                            Wv               = [Wpi{i}',Wv',Wmi{i}']';
                            Wh               = [Wpi{i},Wh,Wmi{i}];
                        end
                        for i = 1:nr
                            T{i}             = -eye(Wv.Dim(1));
                            j                = (l-1)*nr+i;
                            T{i}(j,j)        = -T{i}(j,j);
                        end
                        for i = 1:nr
                            prob             = iqclmi(prob,T{i}*diag(Wv(:,i)),1,0,e);
                            prob             = iqclmi(prob,T{i}*diag(Wh(i,:)),1,0,e);
                        end
                        for i = 1:nr
                            k                = (l-1)*nr+i;
                            for j = 1:Wv.Dim(1)
                                if j ~= k
                                    prob     = iqclmi(prob,blkdiag(-Mv(j,i),-Wv(j,i)),-1,0,[1;1]);
                                    prob     = iqclmi(prob,blkdiag( Mv(j,i),-Wv(j,i)),-1,0,[1;1]);
                                    prob     = iqclmi(prob,blkdiag(-Mh(i,j),-Wh(i,j)),-1,0,[1;1]);
                                    prob     = iqclmi(prob,blkdiag( Mh(i,j),-Wh(i,j)),-1,0,[1;1]);
                                end
                            end
                        end
                    end
                    switch tc
                        case 'on'
                            Z                = iqcvar(prob,2*[size(Phi.a,1),size(Phi.a,1)],'symmetric');
                            prob             = fZ(prob,Z,eye(Z.Dim(1)));
                            B                = fss2m(sc*blkdiag(Phi,Phi)*[zeros(size(Phi,2));eye(size(Phi,2))]);
                            W                = blkdiag(kron(kyp,Z),P);
                            prob             = iqclmi(prob,W,-1,0,B);
                    end
                end
            end
            % Provide input output data
            IO                               = [nr,nr,0,0,size(Phi.c,1),size(Phi.c,1),P11.Dim(1),P22.Dim(1),size(Phi.a,1),size(Phi.a,1)];
            prob                             = fIO(prob,IO);
    end    
    end
    function prob = iqcsbsr_d(Delta,prob)
        error('Error: Dual IQC-multipliers are not supported for this class of uncertainties.');
    end
    function display(obj)
        a1 = [];b1 = [];
        lu = length(obj.OutputChannel);
        for i = 1:lu
            if i == 1
                a2 = ['[',num2str(obj.InputChannel{i}),']'];
                b2 = ['[',num2str(obj.OutputChannel{i}),']'];
            else
                a2 = [', [',num2str(obj.InputChannel{i}),']'];
                b2 = [', [',num2str(obj.OutputChannel{i}),']'];
            end
            a1 = [a1,a2];
            b1 = [b1,b2];
            nr(1,i) = length(obj.OutputChannel{i});
        end
        disp('--------------------------------------------------------------');
        disp(' Sector bounded and (optionally) slope restricted nonlinearity');
        disp('--------------------------------------------------------------');
        disp(' Uncertainty details:');
        disp(['  - Name:                         ',num2str(obj.Name)]);
        disp(['  - Number of repetitions:        ',num2str(obj.NumberOfRepetitions)]);
        disp(['  - Input channels plant:         ',a1]);
        disp(['  - Output channels plant:        ',b1]);
        disp(['  - Odd:                          ',num2str(obj.Odd)]);
        disp(['  - SectorBounds:                 [',num2str(obj.SectorBounds(1)),'  ',num2str(obj.SectorBounds(2)),']']);
        disp(['  - SlopeBounds:                  [',num2str(obj.SlopeBounds(1)),'  ',num2str(obj.SlopeBounds(2)),']']);
        disp(' ');
        disp(' Multiplier details:');
        disp(['  - Basis function type:          ',num2str(obj.BasisFunctionType)]);
        disp(['  - Basis function length:        ',num2str(obj.Length)]);
        disp(['  - Basis function pole location: ',num2str(obj.PoleLocation)]);
        disp(['  - Multiplier type:              ',num2str(obj.PrimalDual)]);
    end
end
end

