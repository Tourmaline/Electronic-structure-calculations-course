function [Vc, ec] = VoskoWilkNusairCorr(rs, polarization)
% Vosko-Wilk-Nusair correlation

if polarization == 0
    A = 0.0621814; b = 3.72744;
    c = 12.9352; x0 = -0.10498;
    b1 = 9.81379; b2 = 2.82224;
    b3 = 0.736412;
else
    A = 0.0310907; b = 7.06042;
    c = 18.0578; x0 = -0.32500;
    b1 = 3.46791; b2 = 1.25842;
    b3 = 0.170393;
end

Xrs = rs + b*sqrt(rs)+c;
Xx0 = x0^2 + b*sqrt(x0^2)+c;
Q = sqrt(4*c-b^2);

temp2 =Q./(2*sqrt(rs)+b);
temp1 = log((sqrt(rs) - x0).^2./Xrs) + 2*(b+2*x0)/Q*atan(temp2);

temp = log(rs./Xrs)+2*b*atan(temp2)/Q - b*x0*temp1./Xx0;

ec = A*temp;
Vc = ec - A/3*(1+b1*sqrt(rs))./(1+b1*sqrt(rs) + b2*rs + b3*rs.^(3/2));

I_inf = isinf(rs);
Vc(I_inf) = 0;
ec(I_inf) = 0;

end