function [Deff, Dph] = ESLeffect(D,HD)
% (Safaeian 2011)
% Effects of endothelial surface layer (ESL) on flow resistance are
% involved with Dph, the physical diameter of vessel reduced by ESL
% thickness, and Deff, the effective diameter.
% input: vessel diameter, hematocrit discharge [as vectors]
% output: effective and physical diameters

%% Experimental parameters for ESL
Wmax = 2.6;    % maximum thickness
Doff = 2.4;    % lower threshold diameter
Dcrit = 10.5;  % critical diameter
D50 = 100;     %
Eamp = 1.1;    %
Ewidth = 0.03; %
Epeak = 0.6;   %
Ehd = 1.18;    % coefficient of sensitivity to HD

%% Solve
nseg = length(D);
Was = zeros(nseg,1);
Wpeak = zeros(nseg,1);
Wph = zeros(nseg,1);
Weff = zeros(nseg,1);
Deff = zeros(nseg,1);
Dph = zeros(nseg,1);

for iseg = 1:nseg
    if D(iseg) <= Doff
        Was(iseg) = 0;
    elseif Doff <= D(iseg)
        Was(iseg) = Wmax*(D(iseg) - Doff)/(D(iseg) + D50 - 2*Doff);
    end
    
    if D(iseg) <= Doff
        Wpeak(iseg) = 0;
    elseif (D(iseg) > Doff) && (D(iseg) <= Dcrit)
        Wpeak(iseg) = Eamp*(D(iseg) - Doff)/(Dcrit - Doff);
    elseif Dcrit < D(iseg)
        Wpeak(iseg) = Eamp*exp(-Ewidth*(D(iseg) - Dcrit));
    end
    
    Wph(iseg) = Was(iseg) + Wpeak(iseg)*Epeak;
    Weff(iseg) = Was(iseg) + Wpeak(iseg)*(1 + HD(iseg)*Ehd);
    
    Deff(iseg) = D(iseg) - 2*Weff(iseg);
    Dph(iseg) = D(iseg) - 2*Wph(iseg);
end

