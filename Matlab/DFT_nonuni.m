function [r, r2Density, Etotal, Eigenvalue] = DFT_nonuni(corr, Nmax, delta)
% Self-consistent DFT program for the helium atom
% Non-uniform grid
% INPUT:
% corr:   approximation of exchange-correlation energy
%       1 = Gunnarsson-Lundqvist
%       2 = Perdew-Zunger
%       3 = Vosko-Wilk-Nusair
%
% Written by Anastasia Kruchinina
% anastasia.kruchinina@it.uu.se
% using matherials and scripts provided on the course
% Electronic Structure Theory and Calculations


Ex = [];
Ec = [];
polarization = 0;

% Construct grid with Nmax points.
% rmin and rmax are the end points of the mesh.
rmin = 0.00001;
rmax = 10;

if nargin == 1
    Nmax = 300;
    delta = 0.05;
end
if nargin == 2
    delta = 0.05;
end

i = 0:Nmax-1;
r = rmin + (rmax-rmin)*(exp(i*delta) - 1)/(exp(Nmax*delta) - 1);
h = r(2:Nmax) - r(1:Nmax-1);

% figure(33)
% plot(r, ones(1, max(size(r))), '*r');

% Number of electrons 
q = 2;

% Convergence criteria
convcrit = 1e-7;

% Construct first guess of charge density * r
r2DensityInit = exp(-2*r)/pi;
r2DensityInit(1) = 0; r2DensityInit(Nmax) = 0;

r2DensityOld = r2DensityInit;


%%%%%%% Hartree potential %%%%%%%%%%%
% Form of problem: Ay = b.
% Construct trigonal square matrix AH. AH has dimension (Nmax-2).
% AH is the same for every iteration, different just b, which depends on
% density

AH = sparse(diag(h(1:Nmax-3), 1) - ...
    diag(h(2:Nmax-1)+h(1:Nmax-2)) + ...
    diag(h(3:Nmax-1), -1) );

d = zeros(1, Nmax-2);
d(1)=1;
b = diag(AH, 1);
c = diag(AH, -1);
for i=1:Nmax-3
    d(i+1) = d(i)*sqrt(b(i)/c(i));
end
Dinv_AH = diag(1./d);
D_AH = diag(d);
AH = D_AH*AH*Dinv_AH;


% Boundary conditions
BC0 = 0; % any constant
BCNmax = q; % charge contained inside a sphere

ediff = inf; % change in energy between iterations
niter = 0;

while ediff >= convcrit
    
    DensityOld = r2DensityOld./(r.^2);
    
    % Calculate b on the mesh.
    bH = -4/2*pi*(h(2:Nmax-1)+h(1:Nmax-2)).*h(2:Nmax-1).*h(1:Nmax-2).*r2DensityOld(2:Nmax-1)./r(2:Nmax-1);
    
    % Add the boundary conditions.
    bH(1) = bH(1) - BC0*h(2);
    bH(Nmax-2) = bH(Nmax-2) - BCNmax*h(Nmax-2);
    bH = D_AH*bH';
    
    % Solve Ay=b.
    yH = AH\bH;
    yH = Dinv_AH*yH;
    
    % Finalizing the result.
    rVH = [BC0; yH; BCNmax]';
    
    %%%%%%% Exchange term %%%%%%%%%%%
    ex = -(3/4)*(3/pi)^(1/3)*DensityOld.^(1/3);
    Vx = -(3/pi)^(1/3)*DensityOld.^(1/3);
    
    Ex = [Ex, trapz(r,(ex).*DensityOld)];
        
    %%%%%%% Correlation term %%%%%%%%%%%
    
    rs = (3./(DensityOld*4*pi)).^(1/3);
    
    if corr == 0
        Vc(1:max(size(rs))) = 0;
        ec(1:max(size(rs))) = 0;
    elseif corr == 1
        [Vc, ec] = GunnarssonLundqvistCorr(rs, polarization);
    elseif corr == 2
        [Vc, ec] = PerdewZungerCorr(rs, polarization);
    elseif corr == 3
        [Vc, ec] = VoskoWilkNusairCorr(rs, polarization);
    else
        display('Wrong correlation parameter!');
        return;
    end
    
    Ec = [Ec, trapz(r,(ec).*DensityOld)];
           
           
    %%%%%%% Schroedinger equation %%%%%%%%%%%
    % Construct tridiagonal matrix A. This square matrix has dimension Nmax-2.
    AS = diag(-1./(h(2:Nmax-2).*(h(2:Nmax-2)+h(1:Nmax-3))), 1) + ...
        diag( 1./(h(2:Nmax-1).*h(1:Nmax-2)) + ...
        (rVH(2:Nmax-1)./r(2:Nmax-1) + Vx(2:Nmax-1) + Vc(2:Nmax-1) - q./r(2:Nmax-1)) ) + ...
        diag(-1./(h(2:Nmax-2).*(h(2:Nmax-2)+h(3:Nmax-1))), -1);
    
    
    d(1)=1;
    b = diag(AS, 1);
    c = diag(AS, -1);
    for i=1:Nmax-3
        d(i+1) = d(i)*sqrt(b(i)/c(i));
    end
    Dinv = diag(1./d);
    AS = diag(d)*AS*Dinv;


    % Boundary conditions. Wave function is zero at first and last point, y(1) = y(Nmax) = 0,
    % enabling the simple matrix eigenvalue formulation of the equation.
    
    % Solve (A + eE)y = 0. E = unity matrix.
    % Find the ground state.
    
    [V,L] = eig(AS);
    [lambdaMin, ilambdaMin] = min(diag(L));
    ymin = V(:,ilambdaMin); % ground-state 
    ymin = (Dinv*ymin)';
    
    
    % Set the end point values to the boundary condition values.
    rphi = [0 ymin 0];
    
    % Normalize the calculated (rphi)^2 solution to q/(4*pi)
    % This is the radial part of the density normalization.
    r2Density = rphi.^2.*q/(4*pi)./trapz(r,rphi.^2);
    
    % Copy over the new density. Mix.
    mix = 0.5;
    r2DensityOld = mix*r2Density + (1-mix)*r2DensityOld;
    
    % Calculate the total energy (Eq.34)
    Etotal(niter + 1) = q*lambdaMin - 8*pi*trapz(r,(1/q)*r2Density.*(0.5*rVH./r - ex - ec + Vx + Vc));
    Eigenvalue(niter + 1) = lambdaMin;
    
    if niter > 1
        ediff = abs(Etotal(niter) - Etotal(niter-1));
    end
    
    niter = niter + 1;
    
end

Ex(end)
Ec(end)
Etotal(end)


end


