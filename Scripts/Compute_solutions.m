clear all
close all

%%%The script generates stresses for a TPE inclusion 
%%% different p and T
%%% Input parameters are in S.I.
%%% Massimo Nespoli 01/03/2023

%%%%%%INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H=10*10^9;               % Constant of Biot (Pa)
alfa=3*10^(-5);          % thermal expansion (1/K)
dp=10e6;                 % pore pressure change of the inclusion (Pa)
dT=100;                  % Temperature change of the inclusion  (°C)
a=2500;                  % disk radius (m)
db=500;                  % disk height  (m)
mu=6*10^9;               % Shear modulus (Pa)
lambda=4*10^9;           % Lamè constant (Pa)
MedianPlane=2000;        % TPE inclusion, depth   of median plane  (m)
limiteplot=8000;         % Limit in plot (max(x)) (m)
k=25;                    % step for plot in x (m)
Zmin=0;                  % min z for computation (m)
Zmax=0;                  % max z for computation (m)
Zstep=25;                % step for plot in z (m)
eta=10^16;               % Viscosity (Pa*s)
TV=2;                    % t1/tau
TV2=0.5;                 % t2/tau
plotshow=1;              % Depth to plot x(plotshow,:)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ni=lambda/(2*(lambda+mu)); %Poisson's ratio

Zetav=Zmin:Zstep:Zmax;

for i=1:length(Zetav)
    
    disp(i)
    zlm=Zetav(i);
[xA(i,:),e11A(i,:),e22A(i,:),e33A(i,:),u1A(i,:),u3A(i,:),e111veA(i,:),e221veA(i,:),e331veA(i,:),ux1veA(i,:),uz1veA(i,:),e112veA(i,:),e222veA(i,:),e332veA(i,:),ux2veA(i,:),uz2veA(i,:),tau11ve(i,:),tau22ve(i,:),tau33ve(i,:),tau111ve(i,:),tau221ve(i,:),tau331ve(i,:),tau112ve(i,:),tau222ve(i,:),tau332ve(i,:),e13A(i,:),e131veA(i,:),e132veA(i,:),tau13veA(i,:),tau131veA(i,:),tau132veA(i,:)]=TPVE_SOLUTIONS_DISK(H,alfa,dp,dT,a,db,ni,mu,lambda,MedianPlane,limiteplot,zlm,k,eta,TV,TV2);

x(i,:)=xA(i,:);
z(i,1:length(x(1,:)))=zlm;
end

figure()
plot(x(plotshow,:),-u3A(plotshow,:),x,u1A(plotshow,:));
hold on
plot(x(plotshow,:),-uz2veA(plotshow,:),x(plotshow,:),ux2veA(plotshow,:));
plot(x(plotshow,:),-uz1veA(plotshow,:),x(plotshow,:),ux1veA(plotshow,:));
legend('-u_3','u_1','-u_3ve (\tau =0.5)','u_1ve (\tau =0.5)','-u_3ve (\tau =2)','u_1ve (\tau =2)');
title('Displacement');
xlabel(' x (m)');
ylabel(' Displacement (m)');
%xlim([0 8000]);
set(gca,'FontSize',12)
set(gca, 'FontName', 'Times');

figure()
plot(x(plotshow,:),e11A(plotshow,:),x(plotshow,:),e22A(plotshow,:),x(plotshow,:),e33A(plotshow,:),x(plotshow,:),e13A(plotshow,:));
hold on
plot(x(plotshow,:),e112veA(plotshow,:),x(plotshow,:),e222veA(plotshow,:),x(plotshow,:),e332veA(plotshow,:),x(plotshow,:),e132veA(plotshow,:));
plot(x(plotshow,:),e111veA(plotshow,:),x(plotshow,:),e221veA(plotshow,:),x(plotshow,:),e331veA(plotshow,:),x(plotshow,:),e131veA(plotshow,:));
legend('e_{11}','e_{22}','e_{33}','e_{13}','e_{11}ve (\tau =0.5)','e_{22}ve (\tau =0.5)','e_{33}ve (\tau =0.5)','e_{13}ve (\tau =0.5)','e_{11}ve (\tau =2)','e_{22}ve (\tau =2)','e_{33}ve (\tau =2)','e_{13}ve (\tau =2)');
title(' Strain');
xlabel(' x (m)');
ylabel(' e_{ij}');
%ylim([-4e-3 4e-3]);
%xlim([0 8000]);
set(gca,'FontSize',12)
set(gca, 'FontName', 'Times')

figure()
plot(x(plotshow,:),tau11ve(plotshow,:),x(plotshow,:),tau22ve(plotshow,:),x(plotshow,:),tau33ve(plotshow,:),x(plotshow,:),tau13veA(plotshow,:));
hold on
plot(x(plotshow,:),tau112ve(plotshow,:),x(plotshow,:),tau222ve(plotshow,:),x(plotshow,:),tau332ve(plotshow,:),x(plotshow,:),tau132veA(plotshow,:));
plot(x(plotshow,:),tau111ve(plotshow,:),x(plotshow,:),tau221ve(plotshow,:),x(plotshow,:),tau331ve(plotshow,:),x(plotshow,:),tau131veA(plotshow,:));
legend('\tau_{11}','\tau_{22}','\tau_{33}','\tau_{13}','\tau_{11}ve (\tau =0.5)','\tau_{22}ve (\tau =0.5)','\tau_{33}ve (\tau =0.5)','\tau_{13}ve (\tau =0.5)','\tau_{11}ve (\tau =2)','\tau_{22}ve (\tau =2)','\tau_{33}ve (\tau =2)','\tau_{13}ve (\tau =2)');
title(' Stress');
xlabel(' x (m)');
ylabel(' \sigma_{ij}');
%ylim([-3e7 2e7]);
%xlim([0 8000]);
set(gca,'FontSize',12)
set(gca, 'FontName', 'Times')

