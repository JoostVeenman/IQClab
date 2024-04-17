function Gred = fMOR(G,w,T,n,disp)
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
% Date:        13-03-2020
% 
% -------------------------------------------------------------------------
%
% Description:  This function performs a model (i.e. state dimension)
%               reduction on the stable multi-variable system G.
%
% Syntax:       Gred = fMOR(G,w)
%               Gred = fMOR(G,w,T)
%               Gred = fMOR(G,w,T,n)
%               Gred = fMOR(G,w,T,n,'off')
%
% Usage:       This function performs a model (i.e. state dimension)
%              reduction on the stable multi-variable system G. The
%              algorithm exploits balanced truncation techniques and Hankel
%              singular values.
%              
%              As input one should provide:
%
%              if nargin == 2
%                - The system G
%                - The frequency vector w
%
%              if nargin == 3
%                - The system G
%                - The frequency vector w
%                - The title of the plot T = '...'
%
%              if nargin == 4
%                - The system G
%                - The frequency vector w
%                - The title of the plot T = '...' (if non-empty)
%                - The desired state dimension k of the reduced system Gred
%
%              if nargin == 5
%                - The system G
%                - The frequency vector w
%                - The title of the plot T = '...' (if non-empty)
%                - The desired state dimension k of the reduced system Gred
%                - disp = 'off' turns of the visual output.
%   
%              Note: if k is chosen as 0 < k < 1, the system will only
%              remove the uncontrollable and unobservable modes with a
%              tolerance of 1. This, hence, is an alternative approach for
%              the applying the command "minreal".
%
%             As output:
%               
%             The reduced system Gred
% -------------------------------------------------------------------------

G       = ss(G);
A       = G.a;
B       = G.b;
C       = G.c;
D       = G.d;

U       = lyapchol(A,B);
L       = lyapchol(A',C');
L       = L';

[Y,S,X] = svd(U*L);
hs      = diag(S);

if nargin < 5
    disp = 1;
elseif nargin == 5
    if strcmp('off',disp) == 1
        disp = 0;
    else
        disp = 1;
    end
end

if disp == 1
    h       = get(0,'CurrentFigure');
    if isempty(h)
        h = 1;
    else
        h = h.Number;
    end
    
    [n1,n2] = size(G);
    m       = reshape(1:1:2*n1*(n2+1),n2+1,2*n1)';
    
    figure(h+1)
    set(gcf,'Position', [23,300,1321,640]);
    set(gcf,'DefaultAxesFontSize',10);
    for i = 1:n1
        for j = 1:n2        
            [mag,phase]   = bode(G(i,j),w);
            mag           = 20*log10(squeeze(mag));
            phase         = squeeze(phase);
            maxmag(i,j)   = max(mag);
            minmag(i,j)   = min(mag);
            maxphase(i,j) = max(phase);
            minphase(i,j) = min(phase);
            
            subplot(2*n1,n2+1,m(2*i-1,j))
            semilogx(w,mag,'b-');
            hold on;
            grid on;
            axis tight
            set(gca,'XMinorGrid','off');
            title(['In(',num2str(j),') -> Out(',num2str(i),')']);
            xlabel('Frequency (rad/s)');
            ylabel('Magnitude (dB)');
            
            subplot(2*n1,n2+1,m(2*i,j))
            semilogx(w,phase,'b-');
            hold on;
            grid on;
            axis tight
            set(gca,'XMinorGrid','off');
            title(['In(',num2str(j),') -> Out(',num2str(i),')']);
            xlabel('Frequency (rad/s)');
            ylabel('Phase (deg)');
        end
    end
    mama = max(max(maxmag));
    mima = min(min(minmag));
    maph = max(max(maxphase));
    miph = min(min(minphase));
    
    magmax   = mama + 0.15*abs(mama);
    magmin   = mima - 0.15*abs(mima);
    phasemax = maph + 0.15*abs(maph);
    phasemin = miph - 0.15*abs(miph);
    
    intvm    = abs(magmax)+abs(magmin);
    
    if intvm > 0.1 && intvm < 1
        com = 0.1;
    elseif intvm > 1 && intvm < 10
        com = 1;
    elseif intvm > 10 && intvm < 80
        com = 20;
    elseif intvm > 80 && intvm < 140
        com = 40;
    elseif intvm > 140
        com = 50;
    end

    intvp     = abs(phasemax)+abs(phasemin);
    if intvp > 0.1 && intvp < 100
        cop = 45;
    elseif intvp > 100
        cop = 90;
    end
    
    am       = -520:com:520;
    ap       = -720:cop:720;
    am1 = am(am>magmin);
    am2 = am1(am1<magmax);
    ap1 = ap(ap>phasemin);
    ap2 = ap1(ap1<phasemax);

    for i = 1:n1
        for j = 1:n2
            subplot(2*n1,n2+1,m(2*i-1,j));
            axis([w(1),w(end),magmin,magmax]);
            set(gca,'YTick',am2);
            subplot(2*n1,n2+1,m(2*i,j))
            axis([w(1),w(end),phasemin,phasemax]);
            set(gca,'YTick',ap2);
        end
    end
    
    legend('Original Model','Location','SouthEast');
    subplot(2*n1,n2+1,m(:,end));
    bar(hs(hs>0),'r');
    title('Hankel singular values');
    pause(1)
end

if nargin < 4
    k  = input('Enter the desired model order: ');
    if nargin == 3
        set(gcf,'Name',T);
    end
elseif nargin == 4
    k = n;
    if ischar(T)
        set(gcf,'Name',T);
    end
elseif nargin > 4
    if isempty(n)
        k  = input('Enter the desired model order: ');
    else
        k = n;
    end
end
if k < 1e-2
    k = length(hs(hs>n));
end

hsk     = hs(1:k);
Yk      = Y(:,1:k);
Xk      = X(:,1:k);
Wk      = L*Xk*diag(1./sqrt(hsk));
Vk      = U'*Yk*diag(1./sqrt(hsk));
Ak      = Wk'*A*Vk;
Bk      = Wk'*B;
Ck      = C*Vk;
Dk      = D;
Gred    = ss(Ak,Bk,Ck,Dk);

if disp == 1
    figure(gcf)
    [n1,n2] = size(G);
    m1      = n1;
    m2      = n2;
    for i = 1:n1
        for j = 1:n2
            [magr,phaser] = bode(Gred(i,j),w);
            magr          = 20*log10(squeeze(magr));
            phaser        = squeeze(phaser);
            maxmagr(i,j)   = max(magr);
            minmagr(i,j)   = min(magr);
            maxphaser(i,j) = max(phaser);
            minphaser(i,j) = min(phaser);
            
            figure(h+1)
            subplot(2*n1,n2+1,m(2*i-1,j));
            semilogx(w,magr,'g--');
            hold on;
            grid on;
            axis tight
            set(gca,'XMinorGrid','off');
            title(['In(',num2str(j),') -> Out(',num2str(i),')']);
            xlabel('Frequency (rad/s)');
            ylabel('Magnitude (dB)');

            subplot(2*n1,n2+1,m(2*i,j))
            semilogx(w,phaser,'g--');
            hold on;
            grid on;
            axis tight
            set(gca,'XMinorGrid','off');
            title(['In(',num2str(j),') -> Out(',num2str(i),')']);
            xlabel('Frequency (rad/s)');
            ylabel('Phase (deg)');
        end
    end
    mama = max(max(maxmagr));
    mima = min(min(minmagr));
    maph = max(max(maxphaser));
    miph = min(min(minphaser));
    
    magmax   = mama + 0.15*abs(mama);
    magmin   = mima - 0.15*abs(mima);
    phasemax = maph + 0.15*abs(maph);
    phasemin = miph - 0.15*abs(miph);
    
    intvm    = abs(magmax)+abs(magmin);
    
    if intvm > 0.1 && intvm < 1
        com = 0.1;
    elseif intvm > 1 && intvm < 10
        com = 1;
    elseif intvm > 10 && intvm < 80
        com = 20;
    elseif intvm > 80 && intvm < 140
        com = 40;
    elseif intvm > 140
        com = 50;
    end

    intvp   = abs(phasemax)+abs(phasemin);
    if intvp > 0.1 && intvp < 100
        cop = 45;
    elseif intvp > 100
        cop = 90;
    end
    
    am  = -520:com:520;
    ap  = -720:cop:720;
    am1 = am(am>magmin);
    am2 = am1(am1<magmax);
    ap1 = ap(ap>phasemin);
    ap2 = ap1(ap1<phasemax);

    for i = 1:n1
        for j = 1:n2
            subplot(2*n1,n2+1,m(2*i-1,j));
            axis([w(1),w(end),magmin,magmax]);
            set(gca,'YTick',am2);
            subplot(2*n1,n2+1,m(2*i,j))
            axis([w(1),w(end),phasemin,phasemax]);
            set(gca,'YTick',ap2);
        end
    end
    
    legend('Original Model','Reduced Model','Location','SouthEast')
    pause(1)
end
end