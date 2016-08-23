function [r, r2Density, Etotal, Eigenvalue] = DFT(corr, Nmax)
% Self-consistent DFT program for the helium atom
% Uniform grid
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

polarization = 0;

% Construct grid with Nmax points.
% rmin and rmax are the end points of the mesh.
rmin = 0.00001;
rmax = 10;

if nargin < 2
    Nmax = 200;
end

r = linspace(rmin,rmax, Nmax);

% figure(34)
% plot(r, ones(1, max(size(r))), '*r');

h = r(2)-r(1);


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

AH = sparse(diag(ones(1,Nmax-3),1) - ...
    2*diag(ones(1,Nmax-2)) + ...
    diag(ones(1,Nmax-3),-1));


% Boundary conditions
BC0 = 0; % any constant
BCNmax = q; % charge contained inside a sphere

ediff = inf; % change in energy between iterations
niter = 0;

while ediff >= convcrit
    
    DensityOld = r2DensityOld./(r.^2);
    
    % Calculate b on the mesh.
    bH = -h^2*4*pi*r2DensityOld(2:Nmax-1)./r(2:Nmax-1);
    
    % Add the boundary conditions.
    bH(1) = bH(1) - BC0;
    bH(Nmax-2) = bH(Nmax-2) - BCNmax;
    
    % Solve Ay=b.
    yH = AH\bH';
    
    % Finalizing the result.
    rVH = [BC0; yH; BCNmax]';
    
    %%%%%%% Exchange term %%%%%%%%%%%
    ex = -(3/4)*(3/pi)^(1/3)*DensityOld.^(1/3);
    Vx = -(3/pi)^(1/3)*DensityOld.^(1/3);
   
    % Vx(1:max(size(r))) = 0;   ex(1:max(size(r))) = 0;
    
    
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
    
    %%%%%%% Schroedinger equation %%%%%%%%%%%
    % Construct tridiagonal matrix A. This square matrix has dimension Nmax-2.
    AS = diag(ones(1,Nmax-3),1) + ...
        diag((-2 - 2*h^2*(rVH(2:Nmax-1)./r(2:Nmax-1) + Vx(2:Nmax-1) + Vc(2:Nmax-1) - q./r(2:Nmax-1)))) + ...
        diag(ones(1,Nmax-3),-1);
    
    
    % Rescale A so that we get an ordinary eigenvalue problem.
    AS = AS./(-2*h^2);
    
    % Boundary conditions. Wave function is zero at first and last point, y(1) = y(Nmax) = 0,
    % enabling the simple matrix eigenvalue formulation of the equation.
    
    % Solve (A + eE)y = 0. E = unity matrix.
    % Find the ground state.
    
    [V,L] = eig(AS);
    [lambdaMin, ilambdaMin] = min(diag(L));
    ymin = V(:,ilambdaMin)'; % ground-state eigenvector
    
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

end


