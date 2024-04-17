function out = sensanalysis(P,meth,N,plotonoff)
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
% Date:        06-04-2024
% 
% -------------------------------------------------------------------------
%
% Description: This function performs sensitivity analyses for LTI systems
%              that are affected by LTI parametric uncertain with the aim
%              to identify driving uncertainties (i.e., identify which
%              uncertainties affect the system's performance the most and
%              quantify the effect). 
%
% Syntax:      out = sensanalysis(P,meth,N)
%              out = sensanalysis(P,meth,N,plotonoff)
%
% Usage:       out = sensanalysis(P,meth,N) performs a sensitivity analysis
%              for the uncertain plant P. Here:
%                1.) P = lft(Delta,M) is an LFT, which is composed of 
%                        a.) the stable LTI plant M (continuous or discrete
%                            time)
%                        b.) the uncertainty block 
%
%                            Delta = diag(delta_1*I_r1, ... delta_k*I_rk),
%
%                            where the uncertain parameters are objects
%                            from the class ureal.
%
%                    Note: In case (some of) the uncertainties are not
%                    normalized, then they are normalized as
%                    delta_i\in[-1,1] in the analysis.
%
%                2.) 'meth' specifies the sensitivity analysis method. 
%
%                    All methods compute the effect of each of the
%                    uncertainties in the output for the function y =
%                    f(Delta), which is the H_\infty norm of P.
%
%                    The function provides the following options for
%                    computing sensitivity indices:
%
%                        a.) 'spea': This option computes the Spearman rank
%                            correlation 'rho_{x,y}'. This is a
%                            "one-at-the-time" (OAT) analysis approach that
%                            requires the least amount of computational
%                            effort.
%
%                            As output one obtains the structure "out" with
%                            the fields: 
%
%                               - out.rho: This is the sensitivity index
%                                 rho in descending order.
%                               - out.unc: This is the uncertainty vector
%                                 in descending order corresponding to rho.
%
%                            The larger the value of rho, the more the
%                            corresponding uncertainty affects the output
%                            of the system.
%
%                        b.) 'src': This option computes the standard
%                            regression coefficients 
% 
%                                      src_i = b_i*sdelta_i/sy,
%
%                            where, respectively, sdelta_i and sy are the
%                            standard deviations of y and delta_i, i =
%                            1,...,k, while b_i are the regression
%                            coefficients of the linear regression model
% 
%                               yhat_j = b0 + sum_{i=1}^k b_i delta_i^j
%
%                            such that yhat_j equals approximately y_j =
%                            f(Delta_j) for some $Delta_j).
%
%                            As output one obtains the structure "out" with
%                            the fields: 
%
%                               - out.src: This is the sensitivity index
%                                 src_i, i = 1,...,k, in descending order.
%                               - out.unc: This is the uncertainty vector
%                                 in descending order corresponding to the
%                                 src indices.
%                               - out.Rs: This is a measure for the
%                                 approximation quality of the linear model
%                                 with Rs = sum_{i=1}^N src_i^2. If Rs >
%                                 0.7, then the approximation can be
%                                 considered (relatively) good. Otherwise
%                                 it is recommended to apply either the
%                                 morris or the variance based method.
%
%                        c.) 'morris': This option computes the elementary
%                            effects 
% 
%                            di(delta) = [f(delta_1,...,delta_{i-1}, ...
%                                           delta_i + De, delta_{i+1}, ...
%                                              delta_k) - f(delta)]/De
%
%                            even though this is an OAT effect, the
%                            sampling method explores the entire parameter
%                            space (see iqclab.eu for further details).
%
%                            Based on the elementary effects, two
%                            sensitivity metrics are computed:
%
%                               - mu_i = mean(abs(di(delta_j)))
%                               - std_i = std(di(delta_j)
%
%                            These are, respectively, the mean (of the
%                            absolute) and the standard deviation of
%                            elementary effects. The mean indicates
%
%                            As output one obtains the structure "out" with
%                            the fields: 
%
%                               - out.mu: This is the sensitivity index
%                                 that reflects the total effect the
%                                 uncertainty has on the output
%                               - out.std: This is a second sensitivity
%                                 index which provides an indication on
%                                 the correlation with other uncertainties
%                                 or the nonlinear influence of the
%                                 parameter.
%                               - out.unc: This is the uncertainty vector
%                                 in descending order corresponding to mu.
%
%                        d.) 'var': This approach quantifies the effect of
%                            each uncertainty with respect to its variance
%                            it generates in the output.
%
%                            This is done by computing two metrics:
%
%                               - Si, i = 1,...,k, which are called the
%                                 first order effects.
%                               - STi, i = 1,...,k, which are called the
%                                 total effects.
%
%                            See iqclab.eu for further details).
%
%                            As output one obtains the structure "out" with
%                            the fields: 
%
%                               - out.STi: This is the sensitivity index
%                                 that reflects the total effect the
%                                 uncertainty has on the output
%                               - out.Si: This is the sensitivity index
%                                 that reflects the first order effects.
%                               - out.unc: This is the uncertainty vector
%                                 in descending order corresponding to STi.
%
%                            Note: Since this method is computational the
%                            most expensive, one could consider first
%                            running Morris method (i.e., 'morris'), and
%                            run the variance based approach for a subset
%                            of paramters that were identified by the
%                            Morris method to be the most sensitive).
%
%                3.) N is the number of samples
%
%                4.) plotonoff = 'on'/'off' is an optional input that
%                    allows to turn 'on' or 'off' the plot generation
%                    within the function.
%
% -------------------------------------------------------------------------
if nargin == 3
    plotonoff = 'on';
elseif nargin == 4
    if ~ischar(plotonoff)
        plotonoff = 'on';
    end
    if ~strcmp(plotonoff,'on') && ~strcmp(plotonoff,'off')
        plotonoff = 'on';
    end
end
methset = {'spea','src','morris','var'};
if ~ismember(meth,methset)
    error('Please specify either of the following methods: spea, src, morris, or var.');
end

[M,d]                       = fLFTdata(P);
k                           = length(d.unc);

for i = 1:k
    if strcmp(d.type{i},'ureal') ~= 1
        error('Error: This function can only consider uncertain systems that are affected by real parametric uncertainties.');
    end
end

un = [];
for i = 1:k
    if d.bounds{i}(1) ~= -1 || d.bounds{i}(2) ~= 1
        un                  = [un,i];
    end
end
if ~isempty(un)
    disp('The following uncertainties are not normalized and will be normalized in the analysis:');
    disp(' ');
    for i = 1:length(un)
        disp(['   - ',d.name{un(i)}]),
    end
    del                     = [];
    for i = 1:k
        if d.bounds{i}(1) ~= -1 || d.bounds{i}(2) ~= 1
            da              = 0.5*(d.bounds{i}(2)+d.bounds{i}(1));
            dd              = 0.5*(d.bounds{i}(2)-d.bounds{i}(1));
            un              = da + dd*ureal(d.name{i},0);
        else
            un              = ureal(d.name{i},0);
        end
        del                 = blkdiag(del,un*eye(length(d.in{i})));
    end
    Pn                      = lft(del,M);
    [M,d]                   = fLFTdata(Pn);
    k                       = length(d.unc);
    eM                      = max(real(eig(M)));
    if eM > -1e-16
        error('The normalized plant is not nominally stable.');
    end
end

switch meth
    case 'morris'
        p                   = 2000;
        DEL                 = p/(2*(p-1));
        v                   = 0:1/(p-1):1-DEL;
        nv                  = length(v);
        B                   = tril(ones(k+1,k),-1);
        Jmk                 = ones(k+1,k);
        Jm1                 = ones(k+1,1);
        for j = 1:N    
            Ds              = diag(sign(randn(k,1)));
            Ps              = eye(k);
            Ps              = Ps(randperm(k),:);
            xs              = randi([0,nv-1],1,k)/(p-1);
            Bs              = (Jm1*xs + DEL/2*((2*B-Jmk)*Ds+Jmk))*Ps;
            Bsc             = 2*Bs-Jm1;
            for w = 1:k
                J(j,w)      = find(abs(Bsc(w+1,:)-Bsc(w,:)));
            end
            for n = 1:k+1
                del         = [];
                for m = 1:k
                    del     = [del,Bsc(n,m)*ones(1,length(d.in{m}))];
                end
                Del{n}      = diag(del);
            end
            for i = 1:k
                di(j,i)     = abs((norm(lft(Del{i+1},M),inf)-norm(lft(Del{i},M),inf))/DEL);
            end
        end
        for n = 1:k
            mu(n)           = mean(di(J == n));
            si(n)           = std(di(J == n));
        end

        [muo,muI]           = sort(mu,'descend');
        out.Unc             = d.name(muI);
        out.mu              = muo;
        out.std             = si(muI);

        if strcmp(plotonoff,'on')
            set(0,'DefaultAxesFontSize',14);
            clrs            = jet(k);
    
            figure;
            for i = 1:k
                plot(out.mu(i),out.std(i),'Marker','o','MarkerEdgeColor',clrs(i,:),'MarkerFaceColor',clrs(i,:),'Color',clrs(i,:));hold on
            end
    
            grid on
            axis([0,1.1*max(out.mu),0,1.1*max(out.std),]);
            title('Sensitivity analysis: Morris method');
            xlabel('Mean absolute elementary effects (\mu_i^*)');
            ylabel('Std elementary effects (\sigma_i)');
            legend(out.Unc,'Location','northwest');
            set(gcf,'Position',[23,44,1321,640]);
            fCutFig(1,1);
        end
    case 'var'
        A                   = 2*rand(N,k)-1;
        B                   = 2*rand(N,k)-1;
        C                   = [A;B];
        for n = 1:2*N
            del             = [];
            for m = 1:k
                del         = [del,C(n,m)*ones(1,length(d.in{m}))];
            end
            ev(n)           = norm(lft(diag(del),M),inf);
        end
        varY                = var(ev);
        for m = 1:k
            ABi             = A;
            ABi(:,m)        = B(:,m);
            for n = 1:N
                delA        = [];
                delB        = [];
                delABi      = [];
                for i = 1:k
                    delA    = [delA,  A(n,i)*ones(1,length(d.in{i}))];
                    delB    = [delB,  B(n,i)*ones(1,length(d.in{i}))];
                    delABi  = [delABi,ABi(n,i)*ones(1,length(d.in{i}))];
                end
                nrmA        = norm(lft(diag(delA),M),inf);
                nrmB        = norm(lft(diag(delB),M),inf);
                nrmABi      = norm(lft(diag(delABi),M),inf);
        
                evSi(n)     = (nrmB-nrmABi)^2;
                evSTi(n)    = (nrmA-nrmABi)^2;
            end
            Si(m)           = abs(varY - sum(evSi)/2/N)/varY;
            STi(m)          = sum(evSTi)/2/N/varY;
        end
        [STio,STiI]         = sort(STi,'descend');
        out.Unc             = d.name(STiI);
        out.STi             = STio;
        out.Si              = Si(STiI);
        

        if strcmp(plotonoff,'on')
            set(0,'DefaultAxesFontSize',14);
            figure;
            for i = 1:k
                bar(i,out.STi(i));hold on
            end
            grid on
            title('Sensitivity analysis: Variance based method');
            xlabel('Uncertain parameters');
            ylabel('Total effect (STi)');
            xticks(1:1:k);
            xticklabels(d.name(STiI));
            xtickangle(45);
            set(gcf,'Position',[23,44,1321,640]);
            fCutFig(1,1);
            
            figure;
            for i = 1:k
                bar(i,out.Si(i));hold on
            end
            grid on
            title('Sensitivity analysis: Variance based method');
            xlabel('Uncertain parameters');
            ylabel('First order effect (Si)');
            xticks(1:1:k);
            xticklabels(d.name(STiI));
            xtickangle(45);
            set(gcf,'Position',[23,44,1321,640]);
            fCutFig(1,1);
        end
    case 'spea'
        for j = 1:k
            xi              = 2*rand(1,N)-1;
            for i = 1:N
                del         = [];
                for m = 1:k
                    if m ~= j
                        del = [del,0*ones(1,length(d.in{m}))];
                    elseif m == j
                        del = [del,xi(i)*ones(1,length(d.in{m}))];
                    end
                end
                yi(i)       = norm(lft(diag(del),M),inf);
            end
            vxi             = xi - mean(xi);
            vyi             = yi - mean(yi);
            rho(j)          = abs(vxi*vyi'/sqrt(sum(vxi.^2)*sum(vxi.^2)));
        end
        [rhoo,rhoI]         = sort(rho,'descend');
        out.Unc             = d.name(rhoI);
        out.rho             = rhoo;

        if strcmp(plotonoff,'on')
            set(0,'DefaultAxesFontSize',14);
            figure;
            for i = 1:k
                bar(i,out.rho(i));hold on
            end
            grid on
            title('Sensitivity analysis: Spearman rank correlation coefficient (SPEA)');
            xlabel('Uncertain parameters');
            ylabel('SPEA (\rho_i)');
            xticks(1:1:k);
            xticklabels(d.name(rhoI));
            xtickangle(45);
            set(gcf,'Position',[23,44,1321,640]);
            fCutFig(1,1);
        end
    case 'src'
        A                   = 2*rand(N,k)-1;
        for i = 1:N
            del             = [];
            for m = 1:k
                del         = [del,A(i,m)*ones(1,length(d.in{m}))];
            end
            y(i)            = norm(lft(diag(del),M),inf);
        end
        
        yb                  = mean(y);
        ys                  = std(y);
        yt                  = (y-yb)/ys;
        
        for j = 1:k
            xb(j)           = mean(A(:,j));
            xs(j)           = std(A(:,j));
            Mt(:,j)         = (A(:,j)-xb(j))/xs(j);
        end
        
        bet                 = abs(pinv(Mt)*yt');
        [beto,betI]         = sort(abs(bet),'descend');
        out.Unc             = d.name(betI);
        out.src             = beto;
        out.Rs              = out.src'*out.src;
        
        if strcmp(plotonoff,'on')
            set(0,'DefaultAxesFontSize',14);
            figure;
            for i = 1:k
                bar(i,abs(out.src(i)));hold on
            end
            grid on
            title('Sensitivity analysis: Standard regression quotients (SRC)');
            xlabel('Uncertain parameters');
            ylabel('SRC (src_i)');
            xticks(1:1:k);
            xticklabels(d.name(betI));
            xtickangle(45);
            set(gcf,'Position',[23,44,1321,640]);
            fCutFig(1,1);
        end
end