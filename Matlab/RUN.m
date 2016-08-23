% Written by Anastasia Kruchinina
% anastasia.kruchinina@it.uu.se
% using matherials and scripts provided on the course
% Electronic Structure Theory and Calculations

uniform = 0;
N = 200;
delta = 0.01;

if uniform == 1
    [r1, r2Density1, Etotal1, Eigenvalue1] = DFT(1);
    [r2, r2Density2, Etotal2, Eigenvalue2] = DFT(2);
    [r3, r2Density3, Etotal3, Eigenvalue3] = DFT(3);
else
    [r1, r2Density1, Etotal1, Eigenvalue1] = DFT_nonuni(1);
    [r2, r2Density2, Etotal2, Eigenvalue2] = DFT_nonuni(2);
    [r3, r2Density3, Etotal3, Eigenvalue3] = DFT_nonuni(3);
end

figure
plot(r1,r2Density1,'r-')
hold on
plot(r2,r2Density2,'b-')
plot(r3,r2Density3,'g-')
xlabel('r(au)');
ylabel('r^2n(r)(au)');
legend('Gunnarsson-Lundqvist', 'Perdew-Zunger', 'Vosko-Wilk-Nusair');
%title('Density*r^2');
hold off

exact=-2.9037; %(au)

figure
plot(Etotal1,'r*-')
hold on
plot(Etotal2,'b*-')
plot(Etotal3,'g*-')
plot([0, 25], exact*ones(1,2), 'm-')
xlabel('Iteration number')
ylabel('Etot (Ha)')
legend('Gunnarsson-Lundqvist', 'Perdew-Zunger', 'Vosko-Wilk-Nusair');
%title('Etotal');
hold off

%exact_eig = -0.903724;


figure
plot(Eigenvalue1,'r*-')
hold on
plot(Eigenvalue2,'b*-')
plot(Eigenvalue3,'g*-')
%plot([0, 25], exact_eig*ones(1,2), 'm-')
xlabel('Iteration number')
ylabel('\lambda _{min}')
legend('Gunnarsson-Lundqvist', 'Perdew-Zunger', 'Vosko-Wilk-Nusair');
%title('Eigenvalue');
hold off