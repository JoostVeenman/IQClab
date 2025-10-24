function obj = iqcassign(Delta,Multiplier,varargin)
% -------------------------------------------------------------------------
%
% Copyright:   This is copyrighted material owned by Novantec B.V.
% Terms:       IQClab is available under a Creative Commons
%              (Attribution-NoDerivatives 4.0 International (CC BY-ND 4.0))
%              license: https://creativecommons.org/licenses/by-nd/4.0/
%              For further information please visit iqclab.eu
%
% Author:      J.Veenman
% Date:        29-11-2019
%              20-10-2025: Added error message in case of wrong inputs
% 
% -------------------------------------------------------------------------
%
% Description: Function to assign IQC-multipliers to uncertainty objects
%              that belong to the class iqcdelta 
%
% Syntax:      obj = iqcassign(Delta,Multiplier)
%              obj = iqcassign(Delta,Multiplier,varargin)
%
% Usage:       obj = iqcassign(Delta,Multiplier) assigns an IQC-multiplier
%              'Multiplier' to the uncertainty object Delta from the class
%              iqcdelta.
%
%              The following multiplier classes are supported:
%
%              1.) 'ultid':   Uncertain LTI dynamics
%              2.) 'ultis':   LTI static diagonally repeated parametric
%                             uncertainties
%              3.) 'ultv':    LTV parametric uncertainties
%              4.) 'ultv_rb': LTV rate bounded parametric uncertainties
%              5.) 'udel':    Uncertain delay
%              6.) 'usbsr':   Sector bounded and slope restricted
%                             nonlinearities
%              7.) 'usg':     Norm-bounded uncertainties/nonlinearities
%                            (small-gain)
%              8.) 'ups':     Passive uncertainties/nonlinearities
%
%              See help u* for further information on each of these
%              classes.
%
%              For obj = iqcassign(Delta,Multiplier,varargin) the input
%              arguments come in pairs. Depending on the particular
%              mutliplier, the properties that can be specified are:
%
%              1.) 'BasisFunctionType' specifies the type of basis
%                   function (1,2,...) (default: 1).
%              2.) 'Length' specifies the length of the basis function > 1.
%                   For multiple repeated uncertainties one can specify
%                   [l1;l2;...lN] (default: 1).
%              3.) 'PoleLocation' specifies the pole locationof the basis
%                   function < 0. For multiple repeated uncertainties one
%                   can specify [pl1;pl2;...plN] (default: -1).
%              4.) 'RelaxationType' specifies the relaxation type of the
%                   multiplier ('DG', 'CH', 'PC', 'ZP') (default: 'DG').
%              5.) 'RelaxProp' specifies an additional relaxation property
%                   for dynamic uncertainties in the class ultis ('S', 'D')
%                   (default: 'S')
%              6.) 'PrimalDual' Specify whether the multiplier should be a
%                   primal or dual parametrization (default: 'Primal').
%                   Note: For standard IQC-analysis, all multipliers are
%                   primal ones.
%
%              Remark: The remaining properties, like the bounds on the
%              uncertainties are automatically assigned to the multiplier
%              class. 
% -------------------------------------------------------------------------
switch Multiplier
    case 'ultis'
        nD   = length(Delta.InputChannel);
        nA   = length(Delta.LinNonlin);
        for i = 1:nA
            if strcmp(Delta.LinNonlin{i},'NL')
                error('Error: This multiplier is not applicable for nonlinear uncertainties.');
            end
            if strcmp(Delta.TimeInvTimeVar{i},'TV')
                disp('Warning: This multiplier is only applicable for time-invariant parametric uncertainties, unless');
                disp('the option "Length" of the basisfunction is set to 1. Then also the class of arbitrarily fast');
                disp('varying uncertainties is covered. The option "RateBounds" is ignored.');
            end
            if strcmp(Delta.StaticDynamic{i},'D')
                error('Error: This multiplier is not applicable for dynamic uncertainties.');
            end
            if strcmp(Delta.Structure{i},'FB')
                error('Error: This multiplier is not applicable for FB uncertainties. Consider "ultv" instead.');
            end
            if ~isempty(Delta.NormBounds{i}) || ~isempty(Delta.SectorBounds{i}) || ~isempty(Delta.SlopeBounds{i}) || ~isempty(Delta.DelayTime{i}) || ~isempty(Delta.DelayType{i})
                error('Error: For this multiplier one should specify either the property "Bounds" or "Polytope".  The others should be kept empty.');
            end
        end
        obj = ultis(Delta.Name{1});
        for i = 1:nD
            n(i,1) = length(Delta.InputChannel{i});
            m(i,1) = length(Delta.OutputChannel{i});
            if n(i,1) ~= m(i,1)
                error('Error: This multiplier requires the uncertainty to have the same number of in- and outputs.');
            end
        end
        set(obj,'NumberOfRepetitions',n);
        set(obj,'InputChannel',Delta.InputChannel);
        set(obj,'OutputChannel',Delta.OutputChannel);
        for i = 1:nD
            if length(Delta.Bounds) == 1 && nD > 1
                i = nD;
            else
                if ~isempty(Delta.Bounds{i})
                    B(:,i) = Delta.Bounds{i}';
                end 
            end
        end
        if exist('B')
            if ~isempty(B)
                [Bs1,Bs2] = size(B);
                if Bs1 ~= 2 && Bs2 ~= nD
                    text = ['Error: The property "Bounds" should have dimension [2,',num2str(nD),'].'];
                    error(text);
                end
                set(obj,'Bounds',B);
            end
        end
        if ~isempty(Delta.Polytope)
            B = Delta.Polytope;
            if ~isempty(Delta.Bounds{1})
                disp('Warning: The option "Bounds" is ignored if the option "Polytope" is specified.');
            end
            set(obj,'Polytope',B);
            set(obj,'Bounds',[]);
        end
        if nargin > 2
            n = length(varargin);
            if mod(n,2) ~= 0
                error('Error: The input arguments should be defined in pairs')
            else
                j = linspace(2,n,n/2);
            end
            for i = 1:length(j)
                prop = varargin{j(i)-1};
                switch prop
                    case 'BasisFunctionType'
                        a = 1:1:10;
                        if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)}) && min(ismember(varargin{j(i)},a)) > 0
                            set(obj,'BasisFunctionType',varargin{j(i)});
                        else
                            error('Error: The type of basis function should be defined as natural numbers in N1 (see "help ultis" for further details)');
                        end
                     case 'Length'
                        a = 1:1:15;
                        if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)}) && min(ismember(varargin{j(i)},a)) > 0
                            set(obj,'Length',varargin{j(i)});
                        else
                            error('Error: The length of a basis function should be defined as a natural number in N1 (see "help ultis" for further details).');
                        end
                    case 'PoleLocation'
                        if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)})
                            set(obj,'PoleLocation',varargin{j(i)});
                        else
                            error('Error: The pole location of a basis function should be defined as a real rational number. For continuous and distrete time systems, it is not allowed to have poles on the imaginary axis and unit circle respectively (see "help ultis" for further details).');
                        end
                    case 'RelaxationType'
                        set(obj,'RelaxationType',varargin{j(i)});
%                     case 'RelaxProp'
%                         set(obj,'RelaxProp',varargin{j(i)});
                    case 'TerminalCost'
                        set(obj,'TerminalCost',varargin{j(i)});
                    case 'PrimalDual'
                        set(obj,'PrimalDual',varargin{j(i)});
                    otherwise
                        error('Error: You can only specify the options "BasisFunctionType", "Length", "PoleLocation", "RelaxationType", "TerminalCost", and "PrimalDual" for this multiplier.');
                end
            end
        end
    case 'ultid'
        nD   = length(Delta.InputChannel);
        nA   = length(Delta.LinNonlin);
        if nA > 1
           error('Error: This multiplier is only applicable for LTI dynamic uncertainties of the form kron(I_nr,Delta)');
        end
        if strcmp(Delta.LinNonlin{1},'NL')
            error('Error: This multiplier is not applicable for nonlinear uncertainties');
        end
        if strcmp(Delta.TimeInvTimeVar{1},'TV')
            error('Error: This multiplier is not applicable for time-varying uncertainties');
        end
        if strcmp(Delta.StaticDynamic{1},'S')
            error('Error: This multiplier is not applicable for static uncertainties ');
        end
        if ~isempty(Delta.Bounds{1}) || ~isempty(Delta.SectorBounds{1}) || ~isempty(Delta.SlopeBounds{1}) || ~isempty(Delta.Polytope)  || ~isempty(Delta.DelayTime{1}) || ~isempty(Delta.DelayType{1})
            error('Error: For this multiplier one should specify the propery "NormBound". The others should be kept empty.');
        end
        obj       = ultid(Delta.Name{1});
        [li1,li2] = size(Delta.InputChannel{1});
        [lo1,lo2] = size(Delta.OutputChannel{1});
        if strcmp(Delta.Structure{1},'D')
            d     = [1,1];
        else
            d     = [li2,lo2];
        end
        if lo1 ~= li1
            error('Error: The in- outputs channels should have the same number of repetitions.')
        else
            n     = li1;
        end
        set(obj,'NumberOfRepetitions',length(Delta.InputChannel));
        set(obj,'Dimensions',d);
        set(obj,'InputChannel',Delta.InputChannel);
        set(obj,'OutputChannel',Delta.OutputChannel);
        if ~isempty(Delta.NormBounds{1})
            B = Delta.NormBounds{1};
        end
        set(obj,'NormBounds',B);
        if nargin > 2
            n = length(varargin);
            if mod(n,2) ~= 0
                error('Error: The input arguments should be defined in pairs')
            else
                j = linspace(2,n,n/2);
            end
            for i = 1:length(j)
                prop = varargin{j(i)-1};
                switch prop
                    case 'BasisFunctionType'
                        a = 1:1:10;
                        if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)}) && min(ismember(varargin{j(i)},a)) > 0
                            set(obj,'BasisFunctionType',varargin{j(i)});
                        else
                            error('Error: The type of basis function should be defined as natural numbers in N1 (see "help ultid" for further details)');
                        end
                     case 'Length'
                        a = 1:1:10;
                        if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)}) && min(ismember(varargin{j(i)},a)) > 0
                            set(obj,'Length',varargin{j(i)});
                        else
                            error('Error: The length of a basis function should be defined as a natural number in N1 (see "help ultid" for further details).');
                        end
                    case 'PoleLocation'
                        if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)})
                            set(obj,'PoleLocation',varargin{j(i)});
                        else
                            error('Error: The pole location of a basis function should be defined as a real rational number. For continuous and distrete time systems, it is not allowed to have poles on the imaginary axis and unit circle respectively (see "help ultid" for further details).');
                        end
                    case 'TerminalCost'
                        set(obj,'TerminalCost',varargin{j(i)});
                    case 'PrimalDual'
                        set(obj,'PrimalDual',varargin{j(i)});
                    otherwise
                        error('Error: You can only specify the options "BasisFunctionType", "Length", "PoleLocation", "TerminalCost", and "PrimalDual" for this multiplier.');
                end
            end
        end
    case 'ultv'
        nD   = length(Delta.InputChannel);
        nA   = length(Delta.LinNonlin);
        for i = 1:nA
            if strcmp(Delta.LinNonlin{i},'NL')
                error('Error: This multiplier is not applicable for nonlinear uncertainties.');
            end
            if strcmp(Delta.TimeInvTimeVar{i},'TI')
                disp('Warning: This multiplier is primarily applicable for the class of time-varying parametric uncertainties.');
                disp('It can be applied to the class of time-invariant parametric uncertainties as well, but this might yield');
                disp('conservative results.');
            end
            if strcmp(Delta.StaticDynamic{i},'D')
                error('Error: This multiplier is not applicable for dynamic uncertainties.');
            end
            if strcmp(Delta.Structure{i},'D')
                error('Error: This multiplier is only applicable for full-block uncertainties. Note, however, that one can consider diagonal uncertainty as well by properly defining "UncertaintyMap".');
            end
            if ~isempty(Delta.Bounds{i}) || ~isempty(Delta.NormBounds{i}) || ~isempty(Delta.SectorBounds{i}) || ~isempty(Delta.SlopeBounds{i}) || ~isempty(Delta.DelayTime{i}) || ~isempty(Delta.DelayType{i})
                error('Error: For this multiplier one should specify the propery "Polytope". The others should be kept empty.');
            end
        end
        obj = ultv(Delta.Name{1});
        set(obj,'InputChannel',Delta.InputChannel);
        set(obj,'OutputChannel',Delta.OutputChannel);
        if ~isempty(Delta.Polytope)
            B = Delta.Polytope;
        end
        if ~isempty(Delta.UncertaintyMap)
            C = Delta.UncertaintyMap;
        end
        if ~exist('B') || ~exist('C')
            error('Error: For this multiplier one should specify the properties: "Polytope" and "UncertaintyMap".');
        end
        set(obj,'Polytope',B);
        set(obj,'UncertaintyMap',C);
        if nargin > 2
            n = length(varargin);
            if mod(n,2) ~= 0
                error('Error: The input arguments should be defined in pairs')
            else
                j = linspace(2,n,n/2);
            end
            for i = 1:length(j)
                prop = varargin{j(i)-1};
                switch prop
                    case 'RelaxationType'
                        set(obj,'RelaxationType',varargin{j(i)});
                    case 'PrimalDual'
                        set(obj,'PrimalDual',varargin{j(i)});
                    otherwise
                        error('Error: You can only specify the options "RelaxationType" and "PrimalDual" for this multiplier.');
                end
            end
        end
    case 'udel'
        nD   = length(Delta.InputChannel);
        if nD > 1
           error('Error: This multiplier is only applicable for a single repeated time-delay uncertainty.');
        end
        if strcmp(Delta.LinNonlin{1},'NL')
            error('Error: This multiplier is not applicable for nonlinear uncertainties.');
        end
        if strcmp(Delta.TimeInvTimeVar{1},'TV')
            disp('Warning: This multiplier is only applicable for time-invariant delay uncertainties, unless');
            disp('the option "Length" of the basisfunction is set to 1. Then also the class of arbitrarily fast');
            disp('varying delay uncertainties is covered.');
        end
        if strcmp(Delta.StaticDynamic{1},'S')
            error('Error: This multiplier is not applicable for static uncertainties.');
        end
        if strcmp(Delta.Structure{1},'FB')
            error('Error: This multiplier not applicable for full-block uncertainties.');
        end
        if ~isempty(Delta.Bounds{1}) || ~isempty(Delta.NormBounds{1}) || ~isempty(Delta.SectorBounds{1}) || ~isempty(Delta.SlopeBounds{1}) || ~isempty(Delta.Polytope)
            error('Error: For this multiplier one should specify the propery "DelayTime" and "DelayType". The others should be kept empty.');
        end
        obj = udel(Delta.Name{1});
        set(obj,'InputChannel',Delta.InputChannel);
        set(obj,'OutputChannel',Delta.OutputChannel);
        set(obj,'NumberOfRepetitions',length(Delta.OutputChannel{1}));
        set(obj,'DelayTime',Delta.DelayTime{1});
        set(obj,'DelayType',Delta.DelayType{1});
        if nargin > 2
            n = length(varargin);
            if mod(n,2) ~= 0
                error('Error: The input arguments should be defined in pairs')
            else
                j = linspace(2,n,n/2);
            end
            for i = 1:length(j)
                prop = varargin{j(i)-1};
                switch prop
                    case 'BasisFunctionType'
                        a = 1:1:10;
                        if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)}) && min(ismember(varargin{j(i)},a)) > 0
                            set(obj,'BasisFunctionType',varargin{j(i)});
                        else
                            error('Error: The type of basis function should be defined as natural numbers in N1 (see "help udel" for further details)');
                        end
                     case 'Length'
                        a = 1:1:15;
                        if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)}) && min(ismember(varargin{j(i)},a)) > 0
                            set(obj,'Length',varargin{j(i)});
                        else
                            error('Error: The length of a basis function should be defined as a natural number in N1 (see "help udel" for further details).');
                        end
                    case 'PoleLocation'
                        if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)})
                            set(obj,'PoleLocation',varargin{j(i)});
                        else
                            error('Error: The pole location of a basis function should be defined as a real rational number. For continuous and distrete time systems, it is not allowed to have poles on the imaginary axis and unit circle respectively (see "help udel" for further details).');
                        end
                    case 'AddIQC'
                        set(obj,'AddIQC',varargin{j(i)});
                    case 'MstrictlyProp'
                        set(obj,'MstrictlyProp',varargin{j(i)});
                    case 'TerminalCost'
                        set(obj,'TerminalCost',varargin{j(i)});
                    case 'PrimalDual'
                        set(obj,'PrimalDual',varargin{j(i)});
                    otherwise
                        error('Error: You can only specify the options "BasisFunctionType", "Length", "PoleLocation", "AddIQC", "MstrictlyProper", "TerminalCost", and "PrimalDual" for this multiplier.');
                end
            end
        end
    case 'ultv_rb'
        nD   = length(Delta.InputChannel);
        nA   = length(Delta.LinNonlin);
        for i = 1:nA
            if strcmp(Delta.LinNonlin{i},'NL')
                error('Error: This multiplier is not applicable for nonlinear uncertainties.');
            end
             if strcmp(Delta.TimeInvTimeVar{i},'TI')
                disp('Error: This multiplier is only applicable for the class of time-varying rate-bounded parametric uncertainties.');
            end
            if strcmp(Delta.StaticDynamic{i},'D')
                error('Error: This multiplier is not applicable for dynamic uncertainties.');
            end
            if strcmp(Delta.Structure{i},'FB')
                error('Error: This multiplier is not applicable for FB uncertainties.');
            end
            if ~isempty(Delta.NormBounds{i}) || ~isempty(Delta.SectorBounds{i}) || ~isempty(Delta.SlopeBounds{i}) || ~isempty(Delta.DelayTime{i}) || ~isempty(Delta.DelayType{i})
                error('Error: For this multiplier one should specify either the property "Bounds/RateBounds" or "Polytope".  The others should be kept empty.');
            end
        end
        obj = ultv_rb(Delta.Name{1});
        for i = 1:nD
            n(i,1) = length(Delta.InputChannel{i});
            m(i,1) = length(Delta.OutputChannel{i});
            if n(i,1) ~= m(i,1)
                error('Error: This multiplier requires the uncertainty to have the same number of in- and outputs.');
            end
        end
        set(obj,'NumberOfRepetitions',n);
        set(obj,'InputChannel',Delta.InputChannel);
        set(obj,'OutputChannel',Delta.OutputChannel);
        for i = 1:nD
            if length(Delta.Bounds) == 1 && nD > 1
                i = nD;
            else
               if ~isempty(Delta.Bounds{i}) && ~isempty(Delta.RateBounds{i})
                    B(:,i) = Delta.Bounds{i}';
                    Br(:,i) = Delta.RateBounds{i}';
                else
                    error('Error: For this multiplier, if specifying the property "Bounds", one should also specify the property "RateBounds".');
                end 
            end
        end
        if exist('B') && exist('Br')
            if ~isempty(B) && ~isempty(Br)
                [Bs1,Bs2] = size(B);
                if Bs1 ~= 2 && Bs2 ~= nD
                    text = ['Error: The property "Bounds" should have dimension [2,',num2str(nD),'].'];
                    error(text);
                end
                [Brs1,Brs2] = size(Br);
                if Brs1 ~= 2 && Brs2 ~= nD
                    text = ['Error: The property "Bounds" should have dimension [2,',num2str(nD),'].'];
                    error(text);
                end
                set(obj,'Bounds',B);
                set(obj,'RateBounds',Br);
            end
        end
        if ~isempty(Delta.Polytope)
            B  = Delta.Polytope;
            if ~isempty(Delta.Bounds{1}) || ~isempty(Delta.RateBounds{1})
                disp('Warning: The option "Bounds/RateBounds" is ignored if the option "Polytope" is specified.');
            end
            set(obj,'Polytope',B);
            set(obj,'Bounds',[]);
            set(obj,'RateBounds',[]);
        end
        if nargin > 2
            n = length(varargin);
            if mod(n,2) ~= 0
                error('Error: The input arguments should be defined in pairs')
            else
                j = linspace(2,n,n/2);
            end
            for i = 1:length(j)
                prop = varargin{j(i)-1};
                switch prop
                    case 'BasisFunctionType'
                        a = 1:1:10;
                        if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)}) && min(ismember(varargin{j(i)},a)) > 0
                            set(obj,'BasisFunctionType',varargin{j(i)});
                        else
                            error('Error: The type of basis function should be defined as natural numbers in N1 (see "help ultv_rb" for further details)');
                        end
                     case 'Length'
                        a = 1:1:10;
                        if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)}) && min(ismember(varargin{j(i)},a)) > 0
                            set(obj,'Length',varargin{j(i)});
                        else
                            error('Error: The length of a basis function should be defined as a natural number in N1 (see "help ultv_rb" for further details).');
                        end
                    case 'PoleLocation'
                        if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)})
                            set(obj,'PoleLocation',varargin{j(i)});
                        else
                            error('Error: The pole location of a basis function should be defined as a real rational number. For continuous and distrete time systems, it is not allowed to have poles on the imaginary axis and unit circle respectively (see "help ultv_rb" for further details).');
                        end
                    case 'RelaxationType'
                        set(obj,'RelaxationType',varargin{j(i)});
                    case 'TerminalCost'
                        set(obj,'TerminalCost',varargin{j(i)});
                    case 'PrimalDual'
                        set(obj,'PrimalDual',varargin{j(i)});
                    otherwise
                        error('Error: You can only specify the options "BasisFunctionType", "Length", "PoleLocation", "RelaxationType", "TerminalCost", and "PrimalDual" for this multiplier.');
                end
            end
        end
    case 'usbsr'
        nD   = length(Delta.InputChannel);
        if nD > 1
           error('Error: This multiplier is only applicable for a single repeated sector bounded and slope restricted nonlinearity.');
        end
        if strcmp(Delta.LinNonlin{1},'L')
            error('Error: This multiplier is not applicable for linear uncertainties.');
        end
        if strcmp(Delta.TimeInvTimeVar{1},'TV')
            disp('Error: This multiplier is only applicable for the class of time-invariant uncertainties.');
        end
        if strcmp(Delta.StaticDynamic{1},'D')
            error('Error: This multiplier is not applicable for dynamic uncertainties.');
        end
        if strcmp(Delta.Structure{1},'FB')
            error('Error: This multiplier is not applicable for FB uncertainties.');
        end
        if ~isempty(Delta.Bounds{1}) || ~isempty(Delta.RateBounds{1}) || ~isempty(Delta.DelayTime{1}) || ~isempty(Delta.DelayType{1}) || ~isempty(Delta.Polytope)
            error('Error: For this multiplier one should specify the property "SectorBounds" and optionally "RateBounds".  The others should be kept empty.');
        end
        obj = usbsr(Delta.Name{1});
        n(1,1) = length(Delta.InputChannel{1});
        m(1,1) = length(Delta.OutputChannel{1});
        if n(1,1) ~= m(1,1)
            error('Error: This multiplier requires the uncertainty to have the same number of in- and outputs.');
        end
        set(obj,'NumberOfRepetitions',n);
        set(obj,'InputChannel',Delta.InputChannel);
        set(obj,'OutputChannel',Delta.OutputChannel);
        if ~isempty(Delta.Odd)
            set(obj,'Odd',Delta.Odd);
        end
        if ~isempty(Delta.SectorBounds{1})
            B = Delta.SectorBounds{1};
            set(obj,'SectorBounds',B);
        else
            error('Error: For this multiplier one should specify the property "SectorBounds".')
        end
        if ~isempty(Delta.SlopeBounds{1})
            C = Delta.SlopeBounds{1};
            set(obj,'SlopeBounds',C);
        end
        if nargin > 2
            n = length(varargin);
            if mod(n,2) ~= 0
                error('Error: The input arguments should be defined in pairs')
            else
                j = linspace(2,n,n/2);
            end
            for i = 1:length(j)
                prop = varargin{j(i)-1};
                switch prop
                    case 'BasisFunctionType'
                        a = 1:1:10;
                        if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)}) && min(ismember(varargin{j(i)},a)) > 0
                            set(obj,'BasisFunctionType',varargin{j(i)});
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
                        if isreal(varargin{j(i)}) && isscalar(varargin{j(i)})
                            obj.PoleLocation = varargin{j(i)};
                        else
                            error('Error: The pole location of a basis function should be defined as a real rational number. For continuous and distrete time systems, it is not allowed to have poles on the imaginary axis and unit circle respectively (see "help usbsr" for further details).');
                        end
                    case 'PrimalDual'
                        set(obj,'PrimalDual',varargin{j(i)});
                    case 'Odd'
                        set(obj,'Odd',varargin{j(i)});
                    case 'TerminalCost'
                        set(obj,'TerminalCost',varargin{j(i)});
                    otherwise
                        error('Error: You can only specify the options "BasisFunctionType", "Length", "PoleLocation", "Odd", "TerminalCost", and "PrimalDual" for this multiplier.');
                end
            end
        end
    case 'usg'
        nD   = length(Delta.InputChannel);
        if nD > 1
           error('Error: This multiplier is only applicable for a single full-block uncertainty.');
        end
        if ~isempty(Delta.Bounds{1}) || ~isempty(Delta.RateBounds{1}) || ~isempty(Delta.SectorBounds{1}) || ~isempty(Delta.SlopeBounds{1}) || ~isempty(Delta.Polytope) || ~isempty(Delta.DelayTime{1}) || ~isempty(Delta.DelayType{1})
            error('Error: For this multiplier one should specify the propery "NormBounds". The others should be kept empty.');
        end
        obj = usg(Delta.Name{1});
        set(obj,'InputChannel',Delta.InputChannel);
        set(obj,'OutputChannel',Delta.OutputChannel);
        set(obj,'NormBounds',Delta.NormBounds{1});
        set(obj,'Dimensions',[size(Delta.InputChannel{1},2),size(Delta.OutputChannel{1},2)]);
        if nargin > 2
            n = length(varargin);
            if mod(n,2) ~= 0
                error('Error: The input arguments should be defined in pairs')
            else
                j = linspace(2,n,n/2);
            end
            for i = 1:length(j)
                prop = varargin{j(i)-1};
                switch prop
                    case 'BasisFunctionType'
                        a = 1:1:10;
                        if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)}) && min(ismember(varargin{j(i)},a)) > 0
                            set(obj,'BasisFunctionType',varargin{j(i)});
                        else
                            error('Error: The type of basis function should be defined as natural numbers in N1 (see "help usg" for further details)');
                        end
                     case 'Length'
                        a = 1:1:10;
                        if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)}) && min(ismember(varargin{j(i)},a)) > 0
                            set(obj,'Length',varargin{j(i)});
                            if varargin{j(i)} > 1
                                disp('Warning: "Length" is only allowed to be larger than 0 if the integrant of the IQC is positive for all omega in [0,\infty].');
                            end
                        else
                            error('Error: The length of a basis function should be defined as a natural number in N1 (see "help usg" for further details).');
                        end
                    case 'PoleLocation'
                        if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)})
                            set(obj,'PoleLocation',varargin{j(i)});
                        else
                            error('Error: The pole location of a basis function should be defined as a real rational number. For continuous and distrete time systems, it is not allowed to have poles on the imaginary axis and unit circle respectively (see "help usg" for further details).');
                        end
                    case 'TerminalCost'
                        set(obj,'TerminalCost',varargin{j(i)});
                    case 'PrimalDual'
                        set(obj,'PrimalDual',varargin{j(i)});
                    otherwise
                        error('Error: You can only specify the options "BasisFunctionType", "Length", "PoleLocation", "TerminalCost", and "PrimalDual" for this multiplier.');
                end
            end
        end
    case 'ups'
        nD   = length(Delta.InputChannel);
        if nD > 1
           error('Error: This multiplier is only applicable for a single full-block uncertainty.');
        end
        if ~isempty(Delta.Bounds{1}) || ~isempty(Delta.RateBounds{1}) || ~isempty(Delta.NormBounds{1}) || ~isempty(Delta.SectorBounds{1}) || ~isempty(Delta.SlopeBounds{1}) || ~isempty(Delta.Polytope) || ~isempty(Delta.DelayTime{1}) || ~isempty(Delta.DelayType{1})
            error('Error: For this multiplier one should specify the propery "Passive". The others should be kept empty.');
        end
        obj = ups(Delta.Name{1});
        set(obj,'InputChannel',Delta.InputChannel);
        set(obj,'OutputChannel',Delta.OutputChannel);
        set(obj,'Dimensions',[size(Delta.InputChannel{1},2),size(Delta.OutputChannel{1},2)]);
        if nargin > 2
            n = length(varargin);
            if mod(n,2) ~= 0
                error('Error: The input arguments should be defined in pairs')
            else
                j = linspace(2,n,n/2);
            end
            for i = 1:length(j)
                prop = varargin{j(i)-1};
                switch prop
                    case 'BasisFunctionType'
                        a = 1:1:10;
                        if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)}) && min(ismember(varargin{j(i)},a)) > 0
                            set(obj,'BasisFunctionType',varargin{j(i)});
                        else
                            error('Error: The type of basis function should be defined as natural numbers in N1 (see "help ups" for further details)');
                        end
                     case 'Length'
                        a = 1:1:10;
                        if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)}) && min(ismember(varargin{j(i)},a)) > 0
                            set(obj,'Length',varargin{j(i)});
                            if varargin{j(i)} > 1
                                disp('Warning: "Length" is only allowed to be larger than 0 if the integrant of the IQC is positive for all omega in [0,\infty].');
                            end
                        else
                            error('Error: The length of a basis function should be defined as a natural number in N1 (see "help ups" for further details).');
                        end
                    case 'PoleLocation'
                        if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)})
                            set(obj,'PoleLocation',varargin{j(i)});
                        else
                            error('Error: The pole location of a basis function should be defined as a real rational number. For continuous and distrete time systems, it is not allowed to have poles on the imaginary axis and unit circle respectively (see "help ups" for further details).');
                        end
                    case 'TerminalCost'
                        set(obj,'TerminalCost',varargin{j(i)});
                    case 'PrimalDual'
                        set(obj,'PrimalDual',varargin{j(i)});
                    otherwise
                        error('Error: You can only specify the options "BasisFunctionType", "Length", "PoleLocation", "TerminalCost", and "PrimalDual" for this multiplier.');
                end
            end
        end
    otherwise
        error(['You cannot specify this multiplier: ',Multiplier]);
end
end