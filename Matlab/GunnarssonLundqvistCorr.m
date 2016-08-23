function [Vc, ec] = GunnarssonLundqvistCorr(rs, polarization)
% Gunnarsson-Lundqvist correlation

if polarization == 0
    Ap = 11.4;
    Cp = 0.0666;
else
    Ap = 15.9;
    Cp = 0.0406;
end

xp = rs./Ap;

ec = -Cp*( (1+xp.^3).*log(1+1./xp)  + 1/2*xp - xp.^2 - 1/3);
Vc = -Cp*log(1+1./xp);

I_inf = isinf(rs);
Vc(I_inf) = 0;
ec(I_inf) = 0;

end