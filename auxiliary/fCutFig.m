function fCutFig(n,m)
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
% Date:        23-12-2019
% 
% -------------------------------------------------------------------------
%
% Description: This function cuts away the empty spaces around figures.
%
% Syntax:      fCutFig(n,m)
%
% Usage:       The inputs correspond to the subplots of the figure:
%
%               1.) n = nr of rows
%               2.) m = nr of columns
%
% -------------------------------------------------------------------------
set(gcf,'PaperPositionMode','auto');
if n == 1 && m == 1
    hold on
    h1                = get(gca,'OuterPosition');
    h2                = get(gca,'TightInset');
    set(gca,'Position',[h1(1)+h2(1),h1(2)+h2(2),h1(3)-h2(1)-h2(3),h1(4)-h2(2)-h2(4)]);
else
    q                 = 0;
    x                 = 1;
    a                 = x/n;
    b                 = 0.5*(1-x);
    y                 = 1;
    c                 = y/m;
    d                 = 0.5*(1-y);
    for i = 1:n
        for j = 1:m
            q         = q + 1;
            subplot(n,m,q)
            hold on
            h1(q,1:4) = get(gca,'OuterPosition');
            h2(q,1:4) = get(gca,'TightInset');
        end
    end
    h2                = max(h2);
    q                 = 0;
    for i = 1:n
        for j = 1:m
            q         = q + 1;
            p(q,1:4)  = [d+(j-1)*c+h2(1),b+(n-i)*a+h2(2),c-h2(1)-h2(3),a-h2(2)-h2(4)];
        end
    end
    b                 = reshape(1:1:n*m,m,n)';
    r1                = 1:1:m;
    r2                = m:-1:1;
    v                 = [];
    for j = 1:m
        v             = [v,r1(j),r2(j)];
    end
    v                 = v(1:m);
    for j = 1:m
        e(1:n,j)      = b(:,v(j));
    end
    q                 = reshape(e',1,n*m); 
    for i = 1:n*m
        subplot(n,m,q(i))
        hold on
        set(gca,'Position',p(q(i),:));
    end
end
end