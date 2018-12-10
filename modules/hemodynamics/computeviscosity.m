function u = computeviscosity(D,Deff,Dph,HD)
% (Safaeian 2011)
% calculate apparent in vivo viscosity
% input: diameter, hematocrit
% output: in vivo viscosity

nseg = length(D);

C = zeros(nseg,1);
u45 = zeros(nseg,1);
uvitro = zeros(nseg,1);
u = zeros(nseg,1);

for iseg = 1:nseg
    C(iseg) = (0.8 + exp(-0.075*Dph(iseg)))*(-1 + 1/(1 + 10*(Dph(iseg)/10)^12))...
        + (1/(10*(Dph(iseg)/10)^12));
    
    u45(iseg) = 220*exp(-1.3*Dph(iseg)) + 3.2 - 2.44*exp(-0.06*Dph(iseg)^0.645);
    
    uvitro(iseg) = 1 + (u45(iseg) - 1)*((1 - HD(iseg))^C(iseg) - 1)/((1-0.45)^C(iseg) - 1);
    
    u(iseg) = uvitro(iseg)*(D(iseg)/Deff(iseg))^4; % in vivo viscosity
end



