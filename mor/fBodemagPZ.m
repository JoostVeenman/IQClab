function fBodemagPZ(w,fignr,G1,sty1,leg1,G2,sty2,leg2)
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
% Date:        13-03-2020
% 
% -------------------------------------------------------------------------
%
% Description:  This function generates a bode magnitude plot of the
%               (possibly multi-variable) system G together with the
%               location of the zeros and poles indicated by circles and
%               crosses along the magnitude and phase plots. 
%
% Syntax:       fBodemagPZ(w,fignr,G1,sty1,leg1)
%               fBodemagPZ(w,fignr,G1,sty1,leg1,G2,sty2,leg2)
%
% Usage:        With nargin == 5, 
%
%               fBodemagPZ(w,fignr,G1,sty1,leg1) 
%
%               generates a bodeplot for the system G1 together with the
%               location of the zeros and poles indicated by circles and
%               crosses respectively along the magnitude and phase plots,
%               where: 
%
%                 - w = [w1,w2,...,wN] : is the frequency vector 
%                 - fignr              : is the figure number
%                 - G1                 : is the system 
%                 - sty1 = '...'       : is the plot style  
%                 - leg1 = '...'       : is the plot legend 
%
%               In addition, if nargin == 8,
%
%               fBodemagPZ(w,fignr,G1,sty1,leg1,G2,sty2,leg2)
%
%               generates the bodeplots for two systems (i.e. G1 and G2),
%               where:
%
%                 - w = [w1,w2,...,wN] : is the frequency vector 
%                 - fignr              : is the figure number
%                 - G1                 : is the first system 
%                 - sty1 = '...'       : is the first plot style  
%                 - leg1 = '...'       : is the first plot legend 
%                 - G2                 : is the second system 
%                 - sty2 = '...'       : is the second plot style  
%                 - leg2 = '...'       : is the second plot legend 
%
%               As output one obtaines:
%
%               - with nargin == 5: The Bode magnitude plot of G1,
%                                   toghether with its poles and zeros.
%
%               - with nargin == 8: The Bode magnitude plots of G1 and G2,
%                                   toghether with the poles and zeros of
%                                   G1.
%
%--------------------------------------------------------------------------

[n1,n2] = size(G1);
% [z,p,k] = zpkdata(G1,'v');

for i=1:n1
    for j = 1:n2
        G1ij = fMOR(G1(i,j),w,[],1e-6,'off');
        [z{i,j},p{i,j},k(i,j)] = zpkdata(G1ij,'v');
    end
end

if nargin > 5
    [m1,m2] = size(G2);
    if m1 ~= n1 && m2 ~= n2
        error('G1 and G2 must have the same I/O dimensions.')
    end
end

% if n1 == 1 && n2 == 1
%     zstr{1,1} = z;
%     pstr{1,1} = p;
%     z         = zstr;
%     p         = pstr;
% end
figure(fignr)
set(gcf,'Position', [23,300,1321,640]);

for i = 1:n1    
    for j = 1:n2
        w1           = w;
        [mag1,~]     = bode(G1(i,j),w1);
        mag1         = 20*log10(squeeze(mag1));
        
        if nargin > 5
            w2       = w;
            [mag2,~] = bode(G2(i,j),w2);
            mag2     = 20*log10(squeeze(mag2));
        end

        w3           = sort(abs(z{i,j}));
        if isempty(w3)
            mag3     = [];
        else
            [mag3,~] = bode(G1(i,j),w3);
            mag3     = 20*log10(squeeze(mag3));
        end
        
        w4           = sort(abs(p{i,j}));
        if isempty(w4)
            mag4     = [];
        else
            [mag4,~] = bode(G1(i,j),w4);
            mag4     = 20*log10(squeeze(mag4));
        end

        subplot(n1,n2,j+(i-1)*n1)
        if nargin > 5
            semilogx(w2,mag2,sty2);hold on;
        end
        semilogx(w1,mag1,sty1);hold on;
        semilogx(w3,mag3,'ko','MarkerEdgeColor','k','MarkerSize',10,'LineWidth',2);hold on
        semilogx(w4,mag4,'rx','MarkerSize',10,'LineWidth',2);hold on;
        grid on;
        set(gca,'XMinorGrid','off');
        title('Bode diagram with pole/zero indication');
        xlabel('Frequency (rad/s)');
        ylabel('Magnitude (dB)');
        if nargin <= 5
            legend(leg1);
        elseif nargin > 5
            legend(leg2,leg1);
        end
    end
end
pause(.1)
end