function [Vc, ec] = PerdewZungerCorr(rs, polarization)
% Perdew-Zunger correlation
% rs - radius
% polarization = 0 (non-polarized) or 1 (polarized)

if polarization == 0
    gamma = -0.1423; beta1 = 1.0529; beta2 = 0.3334;
    A = 0.0311; B = -0.048; C = 0.0020; D = -0.0116;
else
    gamma = -0.0843; beta1 = 1.3981; beta2 = 0.2611;
    A = 0.01555; B = -0.0269; C = 0.0007; D = -0.0048;
end

I_big = rs >= 1;
I_small = rs < 1;
I_inf = isinf(rs);
rs_big = rs(I_big);
rs_small = rs(I_small);

ec(I_big) = gamma./(1 + beta1*sqrt(rs_big) + beta2*rs_big);
Vc(I_big) = ec(I_big).*(1+(7/6)*beta1*sqrt(rs_big)+(4/3)*beta2*rs_big)./(1+beta1*sqrt(rs_big)+beta2*rs_big);

ec(I_small) = A*log(rs_small) + B + C*rs_small.*log(rs_small) + D*rs_small;
Vc(I_small) = A*log(rs_small) + B - A/3 + (2/3)*C*rs_small.*log(rs_small) + (2*D-C)*rs_small./3;

Vc(I_inf) = 0;
ec(I_inf) = 0;

end