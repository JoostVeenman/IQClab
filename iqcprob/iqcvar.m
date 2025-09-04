classdef iqcvar < matlab.mixin.SetGetExactNames

%--------------------------------------------------------------------------
%
% Copyright:   This is copyrighted material owned by Novantec B.V.
% Terms:       IQClab is available under a Creative Commons
%              (Attribution-NoDerivatives 4.0 International (CC BY-ND 4.0))
%              license: https://creativecommons.org/licenses/by-nd/4.0/
%              For further information please visit iqclab.eu
%
% Author:      J.Veenman
% Date:        24-11-2019
% 
%--------------------------------------------------------------------------
%
% Description: Class definition to setup an IQC problem ...
%
% Syntax:      iqcprob = iqcvar(iqcprob,dim,type)
%
% Usage:       iqcprob = iqcvar(iqcprob,dim,type) defines an iqc (matrix)
%              variable. The following inputs should be specified:
%
%                1.) iqcprob - iqcprob object defining the IQC-LMI problem
%                2.) dim:    - Dimension of the variable
%                3.) type:   - Type of variable ('symmetric','full','skew')
%
%              Each variable has 4 properties (var, nvar, svar, Parser).
%
%              If the considered parser is LMIlab, then:
%
%                1.) var  - is the variable identifier (is a number)
%                2.) nvar - is the number of scalar variable in the var
%                3.) svar - is the structure of the variable
%
%              If the considered parser is Yalmip, then:
%
%                1.) var  - is the sdp variable
%                2.) nvar - is zero
%                3.) svar - is the 0-matrix with the same dimension as var
%
%              The following operations are supported for the class iqcvar:
%
%                1.) Horizontal concatenation:   Xnew = [X1,X2,...,XN]
%                2.) Vertical concatenation:     Xnew = [X1;X2;...;XN]
%                3.) Transposition (real):       Xnew = X^T
%                4.) Negation:                   Xnew = -X
%                5.) Multiplication by identity: Xnew = X*I, X = scalar var
%                6.) Blkdiagonal concatenation:  Xnew = blkdiag(X1,...,XN)
%                7.) Skew-diagonal block matrix: Xnew = oblkdiag(X) =
%                                                           = [0,X;X^T,0] 
%
%                8.) Kronecker product:          Xnew = kron(x,y) =
%                                                 = x\otimes y   or
%                                                 = x\otimes y 
%                                                where resp. x and y are
%                                                matrices whos entries are
%                                                limited to 1,0,-1.
%
%                9.) diag:                       Xnew = diag(X)
%               10.) subindexing                 Xnew = X(a:b,c:d)
%
%--------------------------------------------------------------------------
    
properties (SetAccess = public, GetAccess = public)
    var     % Name of the variables
    nvar    % Number of the variables
    svar    % Structure of variable
    Dim     % Dimension of variable
end
properties (SetAccess = private, GetAccess = public)
    Parser  % Parser
end
methods
    function obj = iqcvar(iqcprob,dim,type)
        if isprop(iqcprob,'Parser')
            ptype                                       = iqcprob.Parser;
        elseif ischar(iqcprob)
            ptype                                       = iqcprob;
        end
        if nargin == 1
            obj.var                                     = [];
            obj.nvar                                    = [];
            obj.svar                                    = [];
            obj.Parser                                  = ptype;
            obj.Dim                                     = [];
        elseif nargin == 3
            if strcmp(ptype,'Yalmip')
                switch type
                    case {'symmetric','full','skew'}
                        obj.var                         = sdpvar(dim(1),dim(2),type,'real');
                        obj.nvar                        = 0;
                        obj.svar                        = zeros(size(obj.var));
                        obj.Parser                      = 'Yalmip';
                        obj.Dim                         = [dim(1),dim(2)];
                end            
            elseif strcmp(ptype,'LMIlab')
                switch type
                    case {'symmetric'}
                        [obj.var,obj.nvar,obj.svar]     = lmivar(1,[dim(1) 1]);
                        obj.Parser                      = 'LMIlab';
                        obj.Dim                         = [dim(1),dim(2)];
                    case {'full'}
                        [obj.var,obj.nvar,obj.svar]     = lmivar(2,[dim(1),dim(2)]);
                        obj.Parser                      = 'LMIlab';
                        obj.Dim                         = [dim(1),dim(2)];
                    case {'skew'}
                        if dim(1)<1
                            obj.var                     = [];
                            obj.nvar                    = [];
                            obj.svar                    = [];
                            obj.Parser                  = 'LMIlab';
                            obj.Dim                     = [];
                        elseif dim(1) == 1
                            obj.var                     = 0;
                            obj.nvar                    = 0;
                            obj.svar                    = 0;
                            obj.Parser                  = 'LMIlab';
                            obj.Dim                     = [1,1];
                        elseif dim(1) > 1
                            [obj.var,obj.nvar,obj.svar] = lmivar(1,[0 1]);
                            [obj.var,obj.nvar,obj.svar] = lmivar(3,skewdec(dim(1),obj.nvar));
                            obj.Parser                  = 'LMIlab';
                            obj.Dim                     = [dim(1),dim(2)];
                        end
                end
            end
        end
    end
    function obj = horzcat(varargin)
        % Define the horizontal concatenation of N matrices variables that
        % have the same number of rows: Xnew = [X1,X2,...,XN].
        n = length(varargin);
        if n < 2
            error('Error using iqcvar/horzcat. The function requires at least 2 input arguments.');
        end
        for i = 1:n
            if i == 1
                if ismatrix(varargin{1}) && ~isobject(varargin{1})
                    if norm(varargin{1}) ~= 0
                        error('Error: It is only possible to concatenate IQC-variables and zero matrices.');
                    end
                    for j = 2:n
                        if isobject(varargin{j}) && ~isobject(varargin{1})
                            m               = varargin{1};
                            if strcmp(varargin{j}.Parser,'Yalmip')
                                varargin{1} = iqcvar('Yalmip');
                                set(varargin{1},'var',zeros(size(m)));
                            elseif strcmp(varargin{j}.Parser,'LMIlab')
                                varargin{1} = iqcvar('LMIlab');
                            end
                            set(varargin{1},'svar',zeros(size(m)));
                            set(varargin{1},'Dim',size(m));
                        end
                    end
                end
                if ismatrix(varargin{1}) && ~isobject(varargin{1})
                   error('Error: At least one input should be defined as an object from the iqcvar class.');
                end
            else
                if ismatrix(varargin{i}) && ~isobject(varargin{i})
                    if norm(varargin{i}) ~= 0
                        error('Error: It is only possible to concatenate IQC-variables and zero matrices.');
                    end
                    m = varargin{i};
                    if strcmp(varargin{1}.Parser,'Yalmip')
                        varargin{i} = iqcvar('Yalmip');
                        set(varargin{i},'var',zeros(size(m)));
                    elseif strcmp(varargin{1}.Parser,'LMIlab')
                        varargin{i} = iqcvar('LMIlab');
                        set(varargin{i},'var',[]);
                    end
                    set(varargin{i},'svar',zeros(size(m)));
                    set(varargin{i},'Dim',size(m));
                end
            end
        end
        if strcmp(varargin{1}.Parser,'Yalmip')
            x = 1;
            for i = 1:n
                x      = x*strcmp(varargin{i}.Parser,'Yalmip');
                if i == n && x ~= 1
                    error('Error using iqcvar/horzcat. Each variable should be defined with the same Parser.');
                end
            end
            m  = horzcat(varargin{1}.var,varargin{2}.var);
            if n > 2
                for i = 3:n
                    m = horzcat(m,varargin{i}.var);
                end
            end
            obj = iqcvar('Yalmip');
            set(obj,'var',m);
            set(obj,'nvar',0);
            set(obj,'svar',zeros(size(m)));
            set(obj,'Dim',size(m));
        elseif strcmp(varargin{1}.Parser,'LMIlab')
            x = 1;
            for i = 1:n
                x = x*strcmp(varargin{i}.Parser,'LMIlab');
                if i == n && x ~= 1
                    error('Error using iqcvar/horzcat. Each variable should be defined with the same Parser.');
                end
            end
            m = horzcat(varargin{1}.svar,varargin{2}.svar);
            if n > 2
                for i = 3:n
                    m  = horzcat(m,varargin{i}.svar);
                end
            end
            obj = iqcvar('LMIlab');
            if isempty(m)
                set(obj,'var',[]);
                set(obj,'nvar',[]);
                set(obj,'svar',[]);
                set(obj,'Dim',size(m));
            elseif m == 0
                set(obj,'var',[]);
                set(obj,'nvar',[]);
                set(obj,'svar',zeros(size(m)));
                set(obj,'Dim',size(m));
            else
                [m1,m2,m3] = lmivar(3,m);
                set(obj,'var',m1);
                set(obj,'nvar',m2);
                set(obj,'svar',m3);
                set(obj,'Dim',size(m));
            end
        end
    end
    function obj = vertcat(varargin)
        % Define the vertical concatenation of N matrices variables that
        % have the same number of column: Xnew = [X1;X2;...;XN].
        n              = length(varargin);
        if n < 2
            error('Error using iqcvar/vertcat. The function requires at least 2 input arguments.');
        end
        for i = 1:n
            if i == 1
                if ismatrix(varargin{1}) && ~isobject(varargin{1})
                    if norm(varargin{1}) ~= 0
                        error('Error: It is only possible to concatenate IQC-variables and zero matrices.');
                    end
                    for j = 2:n
                        if isobject(varargin{j}) && ~isobject(varargin{1})
                            m               = varargin{1};
                            if strcmp(varargin{j}.Parser,'Yalmip')
                                varargin{1} = iqcvar('Yalmip');
                                set(varargin{1},'var',zeros(size(m)));
                            elseif strcmp(varargin{j}.Parser,'LMIlab')
                                varargin{1} = iqcvar('LMIlab');
                            end
                            set(varargin{1},'svar',zeros(size(m)));
                            set(varargin{1},'Dim',size(m));
                        end
                    end
                end
                if ismatrix(varargin{1}) && ~isobject(varargin{1})
                   error('Error: At least one input should be defined as an object from the iqcvar class.');
                end
            else
                if ismatrix(varargin{i}) && ~isobject(varargin{i})
                    if norm(varargin{i}) ~= 0
                        error('Error: It is only possible to concatenate IQC-variables and zero matrices.');
                    end
                    m = varargin{i};
                    if strcmp(varargin{1}.Parser,'Yalmip')
                        varargin{i} = iqcvar('Yalmip');
                        set(varargin{i},'var',zeros(size(m)));
                    elseif strcmp(varargin{1}.Parser,'LMIlab')
                        varargin{i} = iqcvar('LMIlab');
                        set(varargin{i},'var',[]);
                    end
                    set(varargin{i},'svar',zeros(size(m)));
                    set(varargin{i},'Dim',size(m));
                end
            end
        end
        if strcmp(varargin{1}.Parser,'Yalmip')
            x          = 1;
            for i = 1:n
                x      = x*strcmp(varargin{i}.Parser,'Yalmip');
                if i == n && x ~= 1
                    error('Error using iqcvar/vertcat. Each variable should be defined with the same Parser.');
                end
            end
            m          = vertcat(varargin{1}.var,varargin{2}.var);
            if n > 2
                for i = 3:n
                    m  = vertcat(m,varargin{i}.var);
                end
            end
            obj        = iqcvar('Yalmip');
            set(obj,'var',m);
            set(obj,'nvar',0);
            set(obj,'svar',zeros(size(m)));
            set(obj,'Dim',size(m));
        elseif strcmp(varargin{1}.Parser,'LMIlab')
            for i = 1:n
                x      = 1;
                x      = x*strcmp(varargin{i}.Parser,'LMIlab');
                if i == n && x ~= 1
                    error('Error using iqcvar/vertcat. Each variable should be defined with the same Parser.');
                end
            end
            m          = vertcat(varargin{1}.svar,varargin{2}.svar);
            if n > 2
                for i = 3:n
                    m  = vertcat(m,varargin{i}.svar);
                end
            end
            obj        = iqcvar('LMIlab');
            if isempty(m)
                set(obj,'var',[]);
                set(obj,'nvar',[]);
                set(obj,'svar',[]);
                set(obj,'Dim',size(m));
            elseif m == 0
                set(obj,'var',[]);
                set(obj,'nvar',[]);
                set(obj,'svar',zeros(size(m)));
                set(obj,'Dim',size(m));
            else
                [m1,m2,m3] = lmivar(3,m);
                set(obj,'var',m1);
                set(obj,'nvar',m2);
                set(obj,'svar',m3);
                set(obj,'Dim',size(m));
            end
        end
    end
    function obj = transpose(X)
        % Define the transpose of the matrix variable X: Xnew = X^T.
        if ismatrix(X) && ~isobject(X)
           obj         = X.';
        end
        if strcmp(X.Parser,'Yalmip')
            obj        = iqcvar('Yalmip');
            m          = transpose(X.var);
            set(obj,'var',m);
            set(obj,'nvar',0);
            set(obj,'svar',zeros(size(m)));
            set(obj,'Dim',size(m));
        elseif strcmp(X.Parser,'LMIlab')
            obj        = iqcvar('LMIlab');
            m          = get(X,'svar');
            m          = m.';
            if isempty(m)
                set(obj,'var',[]);
                set(obj,'nvar',[]);
                set(obj,'svar',[]);
                set(obj,'Dim',size(m));
            elseif m == 0
                set(obj,'var',[]);
                set(obj,'nvar',[]);
                set(obj,'svar',zeros(size(m)));
                set(obj,'Dim',size(m));
            else
                [m1,m2,m3] = lmivar(3,m);
                set(obj,'var',m1);
                set(obj,'nvar',m2);
                set(obj,'svar',m3);
                set(obj,'Dim',size(m));
            end
        end
    end
    function obj = ctranspose(X)
        % Define the complex transpose of the matrix variable X: Xnew = X^*.
        if ismatrix(X) && ~isobject(X)
           obj         = X';
        end
        if strcmp(X.Parser,'Yalmip')
            obj        = iqcvar('Yalmip');
            m          = ctranspose(X.var);
            set(obj,'var',m);
            set(obj,'nvar',0);
            set(obj,'svar',zeros(size(m)));
            set(obj,'Dim',size(m));
        elseif strcmp(X.Parser,'LMIlab')
            obj        = iqcvar('LMIlab');
            m          = get(X,'svar');
            m          = m';
            if isempty(m)
                set(obj,'var',[]);
                set(obj,'nvar',[]);
                set(obj,'svar',[]);
                set(obj,'Dim',size(m));
            elseif m == 0
                set(obj,'var',[]);
                set(obj,'nvar',[]);
                set(obj,'svar',zeros(size(m)));
                set(obj,'Dim',size(m));
            else
                [m1,m2,m3] = lmivar(3,m);
                set(obj,'var',m1);
                set(obj,'nvar',m2);
                set(obj,'svar',m3);
                set(obj,'Dim',size(m));
            end
        end
    end
    function obj = uminus(X)
        % Define the negative of the matrix variable X: Xnew = -X.
        if  ismatrix(X) && ~isobject(X)
           obj         = -X;
        end
        if strcmp(X.Parser,'Yalmip')
            obj        = iqcvar('Yalmip');
            m          = -X.var;
            set(obj,'var',m);
            set(obj,'nvar',0);
            set(obj,'svar',zeros(size(m)));
            set(obj,'Dim',size(m));
        elseif strcmp(X.Parser,'LMIlab')
            obj        = iqcvar('LMIlab');
            m          = get(X,'svar');
            if m == 0 
                set(obj,'var',[]);
                set(obj,'nvar',[]);
                set(obj,'svar',zeros(size(m)));
                set(obj,'Dim',size(m));
            elseif isempty(m)
                set(obj,'var',[]);
                set(obj,'nvar',[]);
                set(obj,'svar',zeros(size(m)));
                set(obj,'Dim',size(m));
            else
                [m1,m2,m3] = lmivar(3,-m);
                set(obj,'var',m1);
                set(obj,'nvar',m2);
                set(obj,'svar',m3);
                set(obj,'Dim',size(m));
            end
        end
    end
    function obj = mtimes(X1,X2)
        % Define the multiplication of a scalar variable X with the
        % identity matrix: Xnew = X*I. Both inputs X1 and X2 can be the
        % scalar variable or the identity matrix.
        if isobject(X1)
            if ismatrix(X2)
                [d1,d2] = size(X2);
                n       = norm(X2);
                trc     = sum(abs(diag(X2)));
                if d1 ~= d2 || n ~= 1 || trc ~= d1
                    error('If x1 is an iqc-variable, then x2 should be defined as a diagonal matrix with ones and minus-ones on the diagonal.');
                end
            else
                error('If x1 is an iqc-variable, then x2 should be defined as a diagonal matrix with ones and minus-ones on the diagonal.');
            end
            if strcmp(X1.Parser,'Yalmip')
                obj     = iqcvar('Yalmip');
                m       = X1.var*X2;
                set(obj,'var',m);
                set(obj,'nvar',0);
                set(obj,'svar',zeros(size(m)));
                set(obj,'Dim',size(m));
            elseif strcmp(X1.Parser,'LMIlab')
                obj     = iqcvar('LMIlab');
                m       = get(X1,'svar')*X2;
                if m == 0
                    set(obj,'var',[]);
                    set(obj,'nvar',[]);
                    set(obj,'svar',zeros(size(m)));
                    set(obj,'Dim',size(m));
                else
                    [m1,m2,m3] = lmivar(3,m);
                    set(obj,'var',m1);
                    set(obj,'nvar',m2);
                    set(obj,'svar',m3);
                    set(obj,'Dim',size(m));
                end
            end
        elseif isobject(X2)
            if ismatrix(X1)
                [d1,d2] = size(X1);
                n       = norm(X1);
                trc     = sum(abs(diag(X1)));
                if d1 ~= d2 || n ~= 1 || trc ~= d1
                    error('If x1 is an iqc-variable, then x2 should be defined as a diagonal matrix with ones and minus-ones on the diagonal.');
                end
            else
                error('If x1 is an iqc-variable, then x2 should be defined as a diagonal matrix with ones and minus-ones on the diagonal.');
            end
            if strcmp(X2.Parser,'Yalmip')
                obj     = iqcvar('Yalmip');
                m       = X1*X2.var;
                set(obj,'var',m);
                set(obj,'nvar',0);
                set(obj,'svar',zeros(size(m)));
                set(obj,'Dim',size(m));
            elseif strcmp(X2.Parser,'LMIlab')
                obj     = iqcvar('LMIlab');
                m       = X1*get(X2,'svar');
                if m == 0
                    set(obj,'var',[]);
                    set(obj,'nvar',[]);
                    set(obj,'svar',zeros(size(m)));
                    set(obj,'Dim',size(m));
                else
                    [m1,m2,m3] = lmivar(3,m);
                    set(obj,'var',m1);
                    set(obj,'nvar',m2);
                    set(obj,'svar',m3);
                    set(obj,'Dim',size(m));
                end
            end
        end
    end
    function obj = blkdiag(varargin)
        % Define the block diagonal concatenation of N matrices variables:
        % Xnew = diag(X1,X2,...,XN).
        n = length(varargin);
        if n < 2
            error('Error using iqcvar/vertcat. The function requires at least 2 input arguments.');
        end
        for i = 1:n
            if i == 1
                if ismatrix(varargin{1}) && ~isobject(varargin{1})
                    if norm(varargin{1}) ~= 0
                        error('Error: It is only possible to concatenate IQC-variables and zero matrices.');
                    end
                    for j = 2:n
                        if isobject(varargin{j}) && ~isobject(varargin{1})
                            m               = varargin{1};
                            if strcmp(varargin{j}.Parser,'Yalmip')
                                varargin{1} = iqcvar('Yalmip');
                                set(varargin{1},'var',zeros(size(m)));
                            elseif strcmp(varargin{j}.Parser,'LMIlab')
                                varargin{1} = iqcvar('LMIlab');
                            end
                            set(varargin{1},'svar',zeros(size(m)));
                            set(varargin{1},'Dim',size(m));
                        end
                    end
                end
                if ismatrix(varargin{1}) && ~isobject(varargin{1})
                   error('Error: At least one input should be defined as an object from the iqcvar class.');
                end
            else
                if ismatrix(varargin{i}) && ~isobject(varargin{i})
                    if norm(varargin{i}) ~= 0
                        error('Error: It is only possible to concatenate IQC-variables and zero matrices.');
                    end
                    m = varargin{i};
                    if strcmp(varargin{1}.Parser,'Yalmip')
                        varargin{i} = iqcvar('Yalmip');
                        set(varargin{i},'var',zeros(size(m)));
                    elseif strcmp(varargin{1}.Parser,'LMIlab')
                        varargin{i} = iqcvar('LMIlab');
                        set(varargin{i},'var',[]);
                    end
                    set(varargin{i},'svar',zeros(size(m)));
                    set(varargin{i},'Dim',size(m));
                end
            end
        end
        if strcmp(varargin{1}.Parser,'Yalmip')
            x              = 1;
            for i = 1:n
                x = x*strcmp(varargin{i}.Parser,'Yalmip');
                if i == n && x ~= 1
                    error('Error using iqcvar/vertcat. Each variable should be defined with the same Parser.');
                end
            end
            m              = varargin{1}.var;
            for i = 2:n
                [r1,r2]    = size(m);
                [c1,c2]    = size(varargin{i}.var);
                m          = [m,zeros(r1,c2);zeros(c1,r2),varargin{i}.var];
            end
            obj = iqcvar('Yalmip');
            set(obj,'var',m);
            set(obj,'nvar',0);
            set(obj,'svar',zeros(size(m)));
            set(obj,'Dim',size(m));
        elseif strcmp(varargin{1}.Parser,'LMIlab')
            for i = 1:n
                x          = 1;
                x          = x*strcmp(varargin{i}.Parser,'LMIlab');
                if i == n && x ~= 1
                    error('Error using iqcvar/vertcat. Each variable should be defined with the same Parser.');
                end
            end
            x          = 1;
            for i = 1:n
                x          = x*strcmp(varargin{i}.Parser,'LMIlab');
                if i == n && x ~= 1
                    error('Error using iqcvar/vertcat. Each variable should be defined with the same Parser');
                end
            end
            m              = varargin{1}.svar;
            for i = 2:n
                [r1,r2]    = size(m);
                [c1,c2]    = size(varargin{i}.svar);
                m          = [m,zeros(r1,c2);zeros(c1,r2),varargin{i}.svar];
            end
            obj            = iqcvar('LMIlab');
            if isempty(m)
                set(obj,'var',[]);
                set(obj,'nvar',[]);
                set(obj,'svar',[]);
                set(obj,'Dim',size(m));
            elseif m == 0
                set(obj,'var',[]);
                set(obj,'nvar',[]);
                set(obj,'svar',zeros(size(m)));
                set(obj,'Dim',size(m));
            else
                [m1,m2,m3] = lmivar(3,m);
                set(obj,'var',m1);
                set(obj,'nvar',m2);
                set(obj,'svar',m3);
                set(obj,'Dim',size(m));
            end
        end
    end
    function obj = oblkdiag(x)
        % Construct the matrix Xnew = [0,X;X^T,0].
        if strcmp(x.Parser,'Yalmip')
            [s1,s2]        = size(x.svar);
            m              = [zeros(s1),x.var;x.var.',zeros(s2)];
            obj            = iqcvar('Yalmip');
            set(obj,'var',m);
            set(obj,'nvar',0);
            set(obj,'svar',zeros(size(m)));
            set(obj,'Dim',size(m));
        elseif strcmp(x.Parser,'LMIlab')
            [s1,s2]        = size(x.svar);
            m              = [zeros(s1),x.svar;x.svar.',zeros(s2)];
            obj            = iqcvar('LMIlab');
            if isempty(m)
                set(obj,'var',[]);
                set(obj,'nvar',[]);
                set(obj,'svar',[]);
                set(obj,'Dim',size(m));
            elseif m == 0
                set(obj,'var',[]);
                set(obj,'nvar',[]);
                set(obj,'svar',zeros(size(m)));
                set(obj,'Dim',size(m));
            else
                [m1,m2,m3] = lmivar(3,m);
                set(obj,'var',m1);
                set(obj,'nvar',m2);
                set(obj,'svar',m3);
                set(obj,'Dim',size(m));
            end
        end
    end
    function obj = kron(x,y)
        % Contstruct the matrix Xnes = kron(x,y) with x being a matrix
        % consisting of ones and minus ones, while y is an LMI variable, or
        % where x is an LMI variable and y is a matrix consisting of ones
        % and minus ones.
        if ismatrix(x) && isreal(x)
            if norm(x-sign(x)) ~= 0
                error('The entries of the matrix x can only consist of ones and minus ones.');
            end
            if strcmp(y.Parser,'Yalmip')
                m          = kron(x,y.var);
                obj        = iqcvar('Yalmip');
                set(obj,'var',m);
                set(obj,'nvar',0);
                set(obj,'svar',zeros(size(m)));
                set(obj,'Dim',size(m));
            elseif strcmp(y.Parser,'LMIlab')
                m          = kron(x,y.svar);
                obj        = iqcvar('LMIlab');
                if isempty(m)
                    set(obj,'var',[]);
                    set(obj,'nvar',[]);
                    set(obj,'svar',[]);
                    set(obj,'Dim',size(m));
                elseif m == 0
                    set(obj,'var',[]);
                    set(obj,'nvar',[]);
                    set(obj,'svar',zeros(size(m)));
                    set(obj,'Dim',size(m));
                else
                    [m1,m2,m3] = lmivar(3,m);
                    set(obj,'var',m1);
                    set(obj,'nvar',m2);
                    set(obj,'svar',m3);
                    set(obj,'Dim',size(m));
                end
            end
        elseif ismatrix(y) && isreal(y)
            if norm(y-sign(y)) ~= 0
                error('The entries of the matrix y can only consist of ones and minus ones.');
            end
            if strcmp(x.Parser,'Yalmip')
                m          = kron(x.var,y);
                obj        = iqcvar('Yalmip');
                set(obj,'var',m);
                set(obj,'nvar',0);
                set(obj,'svar',zeros(size(m)));
                set(obj,'Dim',size(m));
            elseif strcmp(x.Parser,'LMIlab')
                m          = kron(x.svar,y);
                obj        = iqcvar('LMIlab');
                if isempty(m)
                    set(obj,'var',[]);
                    set(obj,'nvar',[]);
                    set(obj,'svar',[]);
                    set(obj,'Dim',size(m));
                elseif m == 0
                    set(obj,'var',[]);
                    set(obj,'nvar',[]);
                    set(obj,'svar',zeros(size(m)));
                    set(obj,'Dim',size(m));
                else
                    [m1,m2,m3] = lmivar(3,m);
                    set(obj,'var',m1);
                    set(obj,'nvar',m2);
                    set(obj,'svar',m3);
                    set(obj,'Dim',size(m));
                end
            end
        else
            error('Error: One of the inputs should be defined as a positive scalar.');
        end
    end
    function obj = diag(x)
        % Pick the diagonal elements of a full square matrix and put them
        % on diagonal of a new matrix: Xnew = diag(X).
        if  ismatrix(x) && ~isobject(x)
           obj         = diag(x);
        end
        if strcmp(x.Parser,'Yalmip')
            [s1,s2]        = size(x.var);
            if s1 ~= s2
                if s1 > 1 && s2 > 1
                    error('Error: This command only applies to square matrices and 1 by n or n by 1 matrices .');
                end
            end
            m              = diag(x.var);
            obj            = iqcvar('Yalmip');
            set(obj,'var',m);
            set(obj,'nvar',0);
            set(obj,'svar',zeros(size(m)));
            set(obj,'Dim',size(m));
        elseif strcmp(x.Parser,'LMIlab')
            [s1,s2]        = size(x.svar);
            if s1 ~= s2
                if s1 > 1 && s2 > 1
                    error('Error: This command only applies to square matrices and 1 by n or n by 1 matrices .');
                end
            end
            m              = diag(x.svar);
            obj            = iqcvar('LMIlab');
            if isempty(m)
                set(obj,'var',[]);
                set(obj,'nvar',[]);
                set(obj,'svar',[]);
                set(obj,'Dim',size(m));
            elseif m == 0
                set(obj,'var',[]);
                set(obj,'nvar',[]);
                set(obj,'svar',zeros(size(m)));
                set(obj,'Dim',size(m));
            else
                [m1,m2,m3] = lmivar(3,m);
                set(obj,'var',m1);
                set(obj,'nvar',m2);
                set(obj,'svar',m3);
                set(obj,'Dim',size(m));
            end
        end
    end
    function obj = subsref(x,s)   
        switch s(1).type
            case '()'
                s1       = s.subs{1};
                s2       = s.subs{2};
                if  ismatrix(x) && ~isobject(x)
                    obj  = x(s.subs{1},s.subs{2});
                end
                if strcmp(x.Parser,'Yalmip')
                    m    = x.var(s1,s2);
                    obj  = iqcvar('Yalmip');
                    set(obj,'var',m);
                    set(obj,'nvar',0);
                    set(obj,'svar',zeros(size(m)));
                    set(obj,'Dim',size(m));
                elseif strcmp(x.Parser,'LMIlab')
                    m    = x.svar(s1,s2);
                    obj  = iqcvar('LMIlab');
                    if isempty(m)
                        set(obj,'var',[]);
                        set(obj,'nvar',[]);
                        set(obj,'svar',[]);
                        set(obj,'Dim',size(m));
                    elseif m == 0
                        set(obj,'var',[]);
                        set(obj,'nvar',[]);
                        set(obj,'svar',zeros(size(m)));
                        set(obj,'Dim',size(m));
                    else
                        [m1,m2,m3] = lmivar(3,m);
                        set(obj,'var',m1);
                        set(obj,'nvar',m2);
                        set(obj,'svar',m3);
                        set(obj,'Dim',size(m));
                    end
                end
            case '.'
                l = size(s,2);
                if l == 1
                    obj = get(x,s(1).subs);
                elseif l > 1
                    y   = get(x,s(1).subs);
                    obj = y(s(2).subs{1});
                end
        end
    end
    function obj = subsasgn(x,s,y)
        if ~isobject(y) && ismatrix(y)
            if norm(y) ~= 0
                error('Error: It is only possible to subassign IQC-variables or zero matrices.');
            end
        end
        switch s(1).type
            case '()'
                s1       = s.subs{1};
                s2       = s.subs{2};
                if isobject(x)
                    if isobject(y)
                        if ~strcmp(x.Parser,y.Parser)
                            error('Error: Each variable should be defined with the same Parser.');
                        end
                        if strcmp(x.Parser,'Yalmip')
                            m        = x.var;
                            m(s1,s2) = y.var;
                            obj      = iqcvar('Yalmip');
                            set(obj,'var',m);
                            set(obj,'nvar',0);
                            set(obj,'svar',zeros(size(m)));
                            set(obj,'Dim',size(m));
                        elseif strcmp(x.Parser,'LMIlab')
                            m        = x.svar;
                            m(s1,s2) = y.svar;
                            obj      = iqcvar('LMIlab');
                            if isempty(m)
                                set(obj,'var',[]);
                                set(obj,'nvar',[]);
                                set(obj,'svar',[]);
                                set(obj,'Dim',size(m));
                            elseif m == 0
                                set(obj,'var',[]);
                                set(obj,'nvar',[]);
                                set(obj,'svar',zeros(size(m)));
                                set(obj,'Dim',size(m));
                            else
                                [m1,m2,m3] = lmivar(3,m);
                                set(obj,'var',m1);
                                set(obj,'nvar',m2);
                                set(obj,'svar',m3);
                                set(obj,'Dim',size(m));
                            end
                        end
                    elseif ismatrix(y)
                        if strcmp(x.Parser,'Yalmip')
                            m        = x.var;
                            m(s1,s2) = 0*y;
                            obj      = iqcvar('Yalmip');
                            set(obj,'var',m);
                            set(obj,'nvar',0);
                            set(obj,'svar',zeros(size(m)));
                            set(obj,'Dim',size(m));
                        elseif strcmp(x.Parser,'LMIlab')
                            m        = x.svar;
                            m(s1,s2) = y;
                            obj      = iqcvar('LMIlab');
                            if isempty(m)
                                set(obj,'var',[]);
                                set(obj,'nvar',[]);
                                set(obj,'svar',[]);
                                set(obj,'Dim',size(m));
                            elseif m == 0
                                set(obj,'var',[]);
                                set(obj,'nvar',[]);
                                set(obj,'svar',zeros(size(m)));
                                set(obj,'Dim',size(m));
                            else
                                [m1,m2,m3] = lmivar(3,m);
                                set(obj,'var',m1);
                                set(obj,'nvar',m2);
                                set(obj,'svar',m3);
                                set(obj,'Dim',size(m));
                            end
                        end
                    end
                end
        end
    end
end
end
