function [x,e11,e22,e33,u1,u3,e111ve,e221ve,e331ve,ux1ve,uz1ve,e112ve,e222ve,e332ve,ux2ve,uz2ve,tau11ve,tau22ve,tau33ve,tau111ve,tau221ve,tau331ve,tau112ve,tau222ve,tau332ve,e13,e131ve,e132ve,tau13ve,tau131ve,tau132ve]=TPVE_SOLUTIONS_DISK(H,alfa,dp,dT,a,db,ni,mu,lambda,MedianPlane,limiteplot,zlm,k,eta,TV,TV2)


%%%The function generates stresses for a 1 disk TPE inclusion 
%%% Input parameters are in S.I.
%%% Massimo Nespoli 01/03/2023

%%%%EXAMPLES OF PARAMETERS
% %%%%%INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H=10*10^9;               % Constant of Biot (Pa)
% alfa=3*10^(-5);          % thermal expansion (1/K)
% dp=10e6;                 % pore pressure change of disk 1 (Pa)
% dT=100;                  % Temperature change of disk 1 (Pa)
% a=2500;                  % disk radius (m)
% db=500;                  % disk height  (m)
% mu=6*10^9;               % Shear modulus (Pa)
% lambda=4*10^9;           % Lamè constant (Pa)
% MedianPlane=2000;        % TPE inclusion, depth   of median plane  (m)
% limiteplot=8000;         % Limit in plot (max(x)) (m)
% k=25;                    % step for plot in x (m)
% Zmin=0;                  % min z for computation (m)
% Zmax=0;                  % max z for computation (m)
% Zstep=25;                % step for plot in z (m)
% eta=10^16;               % Viscosity (Pa*s)
% TV=2;                    % t1/tau
% TV2=0.5;                 % t2/tau
% plotshow=1;              % Depth to plot x(plotshow,:)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K=lambda+2/3*mu;
c=MedianPlane; 

C=(1-2*ni)/(8*pi*mu*(1-ni));
high=db; 
e0=(1/(3*H))*dp+(1/3)*alfa*dT;
e1=e0*(1+ni)/(1-ni);
A=(e1*high)/(2*a);
%DEFINING RECEIVER
x=0.1:k:limiteplot;
y=0;
z=c-zlm;
r=sqrt(x.^2+y.^2+z.^2);
ct=z./r;
st=sqrt(x.^2+y.^2)./r;
sp=y./sqrt(x.^2+y.^2);
cp=x./sqrt(x.^2+y.^2);
errl(length(r))=0;
err(length(r))=0;
ettl(length(r))=0;
ett(length(r))=0;
ertl(length(r))=0;
ert(length(r))=0;
effl(length(r))=0;
eff(length(r))=0;
url(length(r))=0;
utl(length(r))=0;
ur(length(r))=0;
ut(length(r))=0;
M=1000;
L=2*M;
cc(L)=0;
cr(L)=0; 
tau11s(length(r))=0;
tau22s(length(r))=0;
tau33s(length(r))=0;
tau13s(length(r))=0;
distancefroma=0.01*a;
P(L)=0;
pp(L)=0;
ppc(L)=0;
ppp(L)=0;

% %%%%%%%%non_singular part computation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a1=-a;
f2=@(w) sqrt(a^2-w.^2);
Rpp= @(w,x) sqrt((y-w).^2+(zlm+c+high/2)^2+(x-f2(w)).^2);
Rpm= @(w,x) sqrt((y-w).^2+(zlm+c-high/2)^2+(x-f2(w)).^2);
Rmp= @(w,x) sqrt((y-w).^2+(zlm+c+high/2)^2+(x+f2(w)).^2);
Rmm= @(w,x) sqrt((y-w).^2+(zlm+c-high/2)^2+(x+f2(w)).^2);
Rpp2= @(w,x) sqrt((x-w).^2+(zlm+c+high/2)^2+(y-f2(w)).^2);
Rpm2= @(w,x) sqrt((x-w).^2+(zlm+c-high/2)^2+(y-f2(w)).^2);
Rmp2= @(w,x) sqrt((x-w).^2+(zlm+c+high/2)^2+(y+f2(w)).^2);
Rmm2= @(w,x) sqrt((x-w).^2+(zlm+c-high/2)^2+(y+f2(w)).^2);

fun1= @(w,x)  log((Rpp(w,x)+zlm+c+high/2)./(Rpm(w,x)+zlm+c-high/2))-log((Rmp(w,x)+zlm+c+high/2)./(Rmm(w,x)+zlm+c-high/2));
fun2= @(w,x) zlm*((1./Rpp(w,x))-(1./Rmp(w,x))-(1./Rpm(w,x))+(1./Rmm(w,x)));
fun3= @(w,x)  log((Rpp2(w,x)+zlm+c+high/2)./(Rpm2(w,x)+zlm+c-high/2))-log((Rmp2(w,x)+zlm+c+high/2)./(Rmm2(w,x)+zlm+c-high/2));
fun4= @(w,x) zlm*((1./Rpp2(w,x))-(1./Rmp2(w,x))-(1./Rpm2(w,x))+(1./Rmm2(w,x)));
fun5=@(w,x) (log((Rpp(w,x)-x+f2(w))./(Rmp(w,x)-x-f2(w)))-log((Rpm(w,x)-x+f2(w))./(Rmm(w,x)-x-f2(w))));
%fun6old=@(w,x) zlm*((zlm+c+high/2)*((x-f2(w))./(Rpp(w,x).*((y-w).^2+(zlm+c+high/2)^2))-(x+f2(w))./(Rmp(w,x).*((y-w).^2+(zlm+c+high/2)^2)))-(zlm+c-high/2).*((x-f2(w))./(Rpm(w,x).*((y-w).^2+(zlm+c-high/2)^2))-(x+f2(w))./(Rmm(w,x).*((y-w).^2+(zlm+c-high/2)^2))));
fun6=@(w,x) zlm*(zlm+c+high/2)*((x-f2(w))./(Rpp(w,x))-(x+f2(w))./(Rmp(w,x)))./((y-w).^2+(zlm+c+high/2)^2)-zlm*(zlm+c-high/2)*((x-f2(w))./(Rpm(w,x))-(x+f2(w))./(Rmm(w,x)))./((y-w).^2+(zlm+c-high/2)^2);
  
  %  ((x-f2(w))./(Rpp(w,x).*((y-w).^2+(zlm+c+high/2)^2))-(x+f2(w))./(Rmp(w,x).*((y-w).^2+(zlm+c+high/2)^2)))-(zlm+c-high/2).*((x-f2(w))./(Rpm(w,x).*((y-w).^2+(zlm+c-high/2)^2))-(x+f2(w))./(Rmm(w,x).*((y-w).^2+(zlm+c-high/2)^2))));

fun7=@(w,x) ((1./Rpp(w,x))-(1./Rmp(w,x))-(1./Rpm(w,x))+(1./Rmm(w,x)));
fun8=@(w,x) 2*((1./Rpp(w,x))-(1./Rmp(w,x))-(1./Rpm(w,x))+(1./Rmm(w,x)))-zlm*((zlm+c+high/2)./Rpp(w,x).^3-(zlm+c-high/2)./Rpm(w,x).^3-(zlm+c+high/2)./Rmp(w,x).^3+(zlm+c-high/2)./Rmm(w,x).^3);
fun9=@(w,x) -((1./Rpp(w,x))-(1./Rmp(w,x))-(1./Rpm(w,x))+(1./Rmm(w,x)));
fun10=@(w,x) -2*zlm*((zlm+c+high/2)./Rpp(w,x).^3-(zlm+c-high/2)./Rpm(w,x).^3-(zlm+c+high/2)./Rmp(w,x).^3+(zlm+c-high/2)./Rmm(w,x).^3);


fun11=@(w,x) (x-f2(w))./((y-w).^2+(x-f2(w)).^2).*((zlm+c+high/2)./Rpp(w,x)-(zlm+c-high/2)./Rpm(w,x));
fun12=@(w,x) (x+f2(w))./((y-w).^2+(x+f2(w)).^2).*((zlm+c+high/2)./Rmp(w,x)-(zlm+c-high/2)./Rmm(w,x));
fun21=@(w,x) (x-f2(w))./(Rpp(w,x)).^3-(x-f2(w))./(Rpm(w,x)).^3-(x+f2(w))./(Rmp(w,x)).^3+(x+f2(w))./(Rmm(w,x)).^3;

fun31=@(w,x) (y-f2(w))./((x-w).^2+(y-f2(w)).^2).*((zlm+c+high/2)./Rpp2(w,x)-(zlm+c-high/2)./Rpm2(w,x)); %alternative 
fun32=@(w,x) (y+f2(w))./((x-w).^2+(y+f2(w)).^2).*((zlm+c+high/2)./Rmp2(w,x)-(zlm+c-high/2)./Rmm2(w,x)); %alternative
%fun31bis = @(w,x) -2*f2(w)./((x-w).^2+(y+f2(w)).^2).*((zlm+c+high/2)./Rpp2(w,x)-(zlm+c-high/2)./Rpm2(w,x)); %alternative to fun31 and fun 32


fun41=@(w,x) zlm*((y-f2(w))./(Rpp2(w,x)).^3-(y-f2(w))./(Rpm2(w,x)).^3-(y+f2(w))./(Rmp2(w,x)).^3+(y+f2(w))./(Rmm2(w,x)).^3); %alternative
%fun41bis= @(w,x) -4*zlm*f2(w).*(1./(Rpp2(w,x).^3)-1./(Rpm2(w,x).^3)); %alternative to fun41


fun51=@(w,x) -(zlm+c+high/2).*(-(x-f2(w))./(((y-w).^2+(zlm+c+high/2)^2).*Rpp(w,x))+(x+f2(w))./(((y-w).^2+(zlm+c+high/2)^2).*Rmp(w,x)))-(zlm+c-high/2).*((x-f2(w))./(((y-w).^2+(zlm+c-high/2)^2).*Rpm(w,x))-(x+f2(w))./(((y-w).^2+(zlm+c-high/2)^2).*Rmm(w,x)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i=1:length(r)
for m=2:L
P(1,i)=1;
P(2,i)=ct(i);
P(m+1,i)=((2*m-1)*ct(i)*P(m,i)-(m-1)*P(m-1,i))/m; % polinomio di Legendre
pp_zero(i)=0;
ppc(m,i)=(m*ct(i)*P(m+1,i)-(m)*P(m,i))/(1-ct(i)^2);
pp(m,i)=ppc(m,i)*st(i); %First derivative of  Legendre p
ppp(m,i)=ct(i)*pp(m,i)/(-st(i))-m*(m+1)*P(m+1,i);  %Second derivative of  Legendre p
end
for m=1:M
l=2*m;
if m<80
            cc(l)=(((-1)^m)/4^m)*factorial(2*m)/(factorial(m))^2;
         else
             cc(l)=((-1)^m)/(sqrt(pi*m));
 end
cr(l)=cc(l);
cr(1)=1;
%INTERNAL DOMAIN
if abs(r(i))<=a
url(1,i)= cr(1).*P(1,i)*((-1/2));
    
url(l,i)= cc(l).*P(l+1,i).*((1/(l+2)+1/(l-1))-l/(l-1).*(r(i)./a)^(l-1));
utl(l,i)= cc(l).*pp(l,i).*(1/(l+2)+1/(l-1)-1/(l-1).*(r(i)./a)^(l-1));
    
errl(l,i)=l*cc(l).*P(l+1,i).*(abs(r(i))./a).^(l-2);
ettl(l,i)=(cc(l)/(l-1)).*(ppp(l,i)+l.*P(l+1,i)).*(abs(r(i))./a).^(l-2);
effl(l,i)=(cc(l)/(l-1)).*((ct(i)./st(i)).*pp(l,i)+l.*P(l+1,i)).*(abs(r(i))./a).^(l-2);
ertl(l,i)=cc(l).*pp(l,i).*(abs(r(i))./a).^(l-2);
ert(i)=A*sum(ertl(:,i));
% if abs(z)==high/2
% ett(i)=A*(sum(ettl(:,i)))+2*A*a/high;
% else
ett(i)=A*(sum(ettl(:,i)));
% end
err(i)=A*sum(errl(:,i));
eff(i)=A*sum(effl(:,i));

ur(i)=-a*sum(url(:,i));
ut(i)=-a*sum(utl(:,i));

% %CONFINE
% elseif abs(r(i))==a
% err(i)=err(i-1);
% ett(i)=ett(i-1);
% eff(i)=eff(i-1);
% ert(i)=ert(i-1);
%EXTERNAL DOMAIN
elseif abs(r(i))>a


url(1,i)=   cr(1).*P(1,i)*((1/2).*(a./abs(r(i))).^(2));
    
url(l,i)=  cc(l).*P(l+1,i).*(l+1)./(l+2).*(a./abs(r(i))).^(l+2);
utl(l,i)=  cc(l).*pp(l,i).*(1/(l+2).*(a./abs(r(i))).^(l+2));
    
errl(1,i)=(1)*cr(1).*P(1,i).*(a./abs(r(i))).^(3);
effl(1,i)=(cr(1)./(2)).*((-ct(i)./st(i)).*pp_zero(i)+(1).*P(1,i)).*(a./abs(r(i))).^(3); %il meno va nel coseno
errl(l,i)=(l+1)*cr(l).*P(l+1,i).*(a./abs(r(i))).^(l+3);
ertl(l,i)=cc(l).*pp(l,i).*(a./abs(r(i))).^(l+3);
ettl(l,i)=(cc(l)/(m+1)).*(ppp(l,i)-(l+1).*P(l+1,i)).*(a./abs(r(i))).^l;
effl(l,i)=(cc(l)./(l+2)).*((-ct(i)./st(i)).*pp(l,i)+(l+1).*P(l+1,i)).*(a./abs(r(i))).^(l+3); %il meno va nel coseno
ett(i)=(A/2)*(a./abs(r(i))).^3.*(1-sum(ettl(:,i)));
ert(i)=A*sum(ertl(:,i));
err(i)=-A*sum(errl(:,i));
eff(i)=A*sum(effl(:,i));

ur(i)=a*sum(url(:,i));
ut(i)=-a*sum(utl(:,i));   
end
end


Ia1(i)=(3-4*ni).*integral(@(w) fun1(w,x(i)),a1,a);
Ib1(i)=2.*integral(@(w) fun2(w,x(i)),a1,a);
% Ia2(i)=(3-4*ni).*integral(@(w) fun3(w,x(i)),a1,a);
% Ib2(i)=2.*integral(@(w) fun4(w,x(i)),a1,a);
Ia3(i)=-(3-4*ni).*integral(@(w) fun5(w,x(i)),a1,a);
Ib3(i)=2.*integral(@(w) fun6(w,x(i)),a1,a);
%Ib3old(i)=2.*integral(@(w) fun6old(w,x(i)),a1,a)

Ia1ve(i)=integral(@(w) fun1(w,x(i)),a1,a);
Ib1ve(i)=Ib1(i);
% Ia2(i)=integral(@(w) fun3(w,x(i)),a1,a);
% Ib2(i)=2.*integral(@(w) fun4(w,x(i)),a1,a);
Ia3ve(i)=-integral(@(w) fun5(w,x(i)),a1,a);
Ib3ve(i)=Ib3(i);

R2=@(w,j,k) sqrt((x(i)-j).^2+(y-w).^2+(zlm+k).^2); 
fun61=@(w,j) (((zlm+c+high/2)./R2(w,j,c+high/2).^3-(zlm +c- high/2)./R2(w,j,c-high/2).^3)+zlm./R2(w,j,c+high/2).^3-zlm./R2(w,j,c-high/2).^3-3*zlm*(zlm +c+high/2)^2./R2(w,j,c+high/2).^5+3*zlm*(zlm +c-high/2)^2./R2(w,j,c-high/2).^5);

dIa1(i)=-(3-4*ni)*(integral(@(w) fun11(w,x(i)),a1,a))+(3-4*ni)*(integral(@(w) fun12(w,x(i)),a1,a));
dIb1(i)=-2*zlm*integral(@(w) fun21(w,x(i)),a1,a);
dIa2(i)=(3-4*ni)*(integral(@(w) fun31(w,x(i)),a1,a)-integral(@(w) fun32(w,x(i)),a1,a)); %alternative dIa2 integral
%dIa2_alternative(i)=(3-4*ni)*(integral(@(w) fun31bis(w,x(i)),a1,a)); %alternative dIa2 integral

dIb2(i)=2*integral(@(w) fun41(w,x(i)),a1,a); %alternative dIb2 integral
%dIb2_alternative(i)=integral(@(w) fun41bis(w,x(i)),a1,a); %alternative dIb2 integral

dIa3(i)=-(3-4*ni)*integral(@(w) fun51(w,x(i)),a1,a);
%dIb3(i)=((4*zlm+2*c+high)*integral(@(w) fun61(w,x(i)),a1,a)-4*zlm*(zlm+c+high/2)^2*integral(@(w) fun62(w,x(i)),a1,a)-2*zlm*(zlm+c+high/2)^2*integral(@(w) fun63(w,x(i)),a1,a)-(4*zlm+2*c-high)*integral(@(w) fun64(w,x(i)),a1,a)+4*zlm*(zlm+c-high/2)^2*integral(@(w) fun65(w,x(i)),a1,a)+2*zlm*(zlm+c-high/2)^2*integral(@(w) fun66(w,x(i)),a1,a));
dIb3_1=@(w) integral(@(j) fun61(w,j),-sqrt(a^2-w.^2), sqrt(a^2-w.^2));
dIb3(i)=-2*integral(@(w) dIb3_1(w),a1,a, 'ArrayValued', true);

dIa13(i)=  (3-4*ni)*(integral(@(w) fun7(w,x(i)),a1,a));
dIb13(i)=  integral(@(w) fun8(w,x(i)),a1,a);
dIa31(i)=  (3-4*ni)*(integral(@(w) fun9(w,x(i)),a1,a));
dIb31(i)=  integral(@(w) fun10(w,x(i)),a1,a);

dIa1ve(i)=-(integral(@(w) fun11(w,x(i)),a1,a))+(integral(@(w) fun12(w,x(i)),a1,a));
dIb1ve(i)=dIb1(i);
dIa2ve(i)=(integral(@(w) fun31(w,x(i)),a1,a)-integral(@(w) fun32(w,x(i)),a1,a));
dIb2ve(i)=dIb2(i);
dIa3ve(i)=-integral(@(w) fun51(w,x(i)),a1,a);
%dIb3(i)=((4*zlm+2*c+high)*integral(@(w) fun61(w,x(i)),a1,a)-4*zlm*(zlm+c+high/2)^2*integral(@(w) fun62(w,x(i)),a1,a)-2*zlm*(zlm+c+high/2)^2*integral(@(w) fun63(w,x(i)),a1,a)-(4*zlm+2*c-high)*integral(@(w) fun64(w,x(i)),a1,a)+4*zlm*(zlm+c-high/2)^2*integral(@(w) fun65(w,x(i)),a1,a)+2*zlm*(zlm+c-high/2)^2*integral(@(w) fun66(w,x(i)),a1,a));
dIb3_1ve=dIb3_1;
dIb3ve(i)=dIb3(i);

dIa13ve(i)=  dIa13(i);
dIb13ve(i)=  integral(@(w) fun8(w,x(i)),a1,a);
dIa31ve(i)=  dIa31(i);
dIb31ve(i)=  integral(@(w) fun10(w,x(i)),a1,a);
end

ett = interp1(x(r<a-distancefroma | r>a+distancefroma ),ett(r<a-distancefroma | r>a+distancefroma ),x);
ert = interp1(x(r<a-distancefroma | r>a+distancefroma ),ert(r<a-distancefroma | r>a+distancefroma ),x);
err = interp1(x(r<a-distancefroma | r>a+distancefroma ),err(r<a-distancefroma | r>a+distancefroma ),x);
eff = interp1(x(r<a-distancefroma | r>a+distancefroma ),eff(r<a-distancefroma | r>a+distancefroma ),x);
ur = interp1(x(r<a-distancefroma | r>a+distancefroma ),ur(r<a-distancefroma | r>a+distancefroma ),x);
ut = interp1(x(r<a-distancefroma | r>a+distancefroma ),ut(r<a-distancefroma | r>a+distancefroma ),x);


err=err/A;
ett=ett/A;
ert=ert/A;
eff=eff/A;
% 
% SPOSTAMENTO NON SINGOLARE
u1ns=3*K*C*e0*(Ia1+Ib1);
%u2ns=3*K*C*e0*(Ia2+Ib2);
u3ns=-3*K*C*e0*(Ia3+Ib3);
% %STRAIN NON SINGOLARE
du1=3*K*C*e0*(dIa1+dIb1);
du2=-3*K*C*e0*(dIa2+dIb2);
du3=-3*K*C*e0*(dIa3+dIb3);
du13=3*K*C*e0/2*(dIa13+dIb13+dIa31+dIb31);


e11ns=du1;
e22ns=du2;
e33ns=du3;
e13ns=du13;


syms s
lambdas=(s*lambda+K*mu/eta)/(s+mu/eta);
mus=(s/(s/mu+1/eta));
Ks=K;
nus=lambdas/(2*(lambdas+mus));


At=e0/s*(3*lambdas+2*mus)/(lambdas+2*mus)*db/(2*a);
IL=ilaplace(At);
gf=matlabFunction(IL);

Ic1=3*Ks*(1-2*nus)*(3-4*nus)/(8*pi*mus*(1-nus))*e0/s;
Ic2=3*Ks*(1-2*nus)/(8*pi*mus*(1-nus))*e0/s;

IL1=ilaplace(Ic1);
gf1=matlabFunction(IL1);
IL2=ilaplace(Ic2);
gf2=matlabFunction(IL2);

%%%%%%%%%%%%%%%%%for stress%%%%%%%%%%%

Atmu=e0*mus/s*(3*lambdas+2*mus)/(lambdas+2*mus)*db/(2*a);

ILAmu=ilaplace(Atmu);
gfAmu=matlabFunction(ILAmu);

Atlambda=e0*lambdas/s*(3*lambdas+2*mus)/(lambdas+2*mus)*db/(2*a);
ILAlambda=ilaplace(Atlambda);
gfAlambda=matlabFunction(ILAlambda);

Y1mu=mus*3*Ks*(1-2*nus)*(3-4*nus)/(8*pi*mus*(1-nus))*e0/s;
Y1lambda=lambdas*3*Ks*(1-2*nus)*(3-4*nus)/(8*pi*mus*(1-nus))*e0/s;

ILY1mu=ilaplace(Y1mu);
ILY1lambda=ilaplace(Y1lambda);

gfY1mu=matlabFunction(ILY1mu);
gfY1lambda=matlabFunction(ILY1lambda);

Y2mu=mus*3*Ks*(1-2*nus)/(8*pi*mus*(1-nus))*e0/s;
Y2lambda=lambdas*3*Ks*(1-2*nus)/(8*pi*mus*(1-nus))*e0/s;

ILY2mu=ilaplace(Y2mu);
ILY2lambda=ilaplace(Y2lambda);

gfY2mu=matlabFunction(ILY2mu);
gfY2lambda=matlabFunction(ILY2lambda);


urE=A*ur;
utE=A*ut;
tau=eta/mu;

t=TV*tau;
urt1=ur*gf(t);
utt1=ut*gf(t);
ttau=t/tau;

t2=TV2*tau;
urt2=ur*gf(t2);
utt2=ut*gf(t2);
ttau2=t2/tau;

errE=err*A;
ertE=ert*A;
effE=eff*A;
ettE=ett*A;

err1=err*gf(t);
ert1=ert*gf(t);
eff1=eff*gf(t);
ett1=ett*gf(t);

err2=err*gf(t2);
ert2=ert*gf(t2);
eff2=eff*gf(t2);
ett2=ett*gf(t2);

taurr=2*gfAmu(0)*err+gfAlambda(0)*(err+eff+ett);
tauff=2*gfAmu(0)*eff+gfAlambda(0)*(err+eff+ett);
tautt=2*gfAmu(0)*ett+gfAlambda(0)*(err+eff+ett);
taurt=2*gfAmu(0)*ert;


taurr1=2*gfAmu(t)*err+gfAlambda(t)*(err+eff+ett);
tauff1=2*gfAmu(t)*eff+gfAlambda(t)*(err+eff+ett);
tautt1=2*gfAmu(t)*ett+gfAlambda(t)*(err+eff+ett);
taurt1=2*gfAmu(t)*ert;

taurr2=2*gfAmu(t2)*err+gfAlambda(t2)*(err+eff+ett);
tauff2=2*gfAmu(t2)*eff+gfAlambda(t2)*(err+eff+ett);
tautt2=2*gfAmu(t2)*ett+gfAlambda(t2)*(err+eff+ett);
taurt2=2*gfAmu(t2)*ert;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
  uxns1=gf1(t)*Ia1ve+gf2(t)*Ib1ve;
  uzns1=-gf1(t)*Ia3ve-gf2(t)*Ib3ve;
%  
   uxns2=gf1(t2)*Ia1ve+gf2(t2)*Ib1ve;
   uzns2=-gf1(t2)*Ia3ve-gf2(t2)*Ib3ve;
  
  e11ns1=gf1(t)*dIa1ve+gf2(t)*dIb1ve;
  e22ns1=-gf1(t)*dIa2ve-gf2(t)*dIb2ve;
  e33ns1=-gf1(t)*dIa3ve-gf2(t)*dIb3ve;
  e13ns1=gf1(t)/2*(dIa13ve+dIa31ve)+gf2(t)/2*(dIb13ve+dIb31ve);
 
   e11ns2=gf1(t2)*dIa1ve+gf2(t2)*dIb1ve;
   e22ns2=-gf1(t2)*dIa2ve-gf2(t2)*dIb2ve;
   e33ns2=-gf1(t2)*dIa3ve-gf2(t2)*dIb3ve;
   e13ns2=gf1(t2)/2*(dIa13ve+dIa31ve)+gf2(t2)/2*(dIb13ve+dIb31ve);
   
   tau11ns=2*gfY1mu(0)*dIa1ve+2*gfY2mu(0)*dIb1ve+gfY1lambda(0)*(dIa1ve-dIa2ve-dIa3ve)+gfY2lambda(0)*(dIb1ve-dIb2ve-dIb3ve);
   tau22ns=-2*gfY1mu(0)*dIa2ve-2*gfY2mu(0)*dIb2ve+gfY1lambda(0)*(dIa1ve-dIa2ve-dIa3ve)+gfY2lambda(0)*(dIb1ve-dIb2ve-dIb3ve);
   tau33ns=-2*gfY1mu(0)*dIa3ve-2*gfY2mu(0)*dIb3ve+gfY1lambda(0)*(dIa1ve-dIa2ve-dIa3ve)+gfY2lambda(0)*(dIb1ve-dIb2ve-dIb3ve);
   tau13ns=gfY1mu(0)*(dIa13ve+dIa31ve)+gfY2mu(0)*(dIb13ve+dIb31ve);
  
  
   tau11ns1=2*gfY1mu(t)*dIa1ve+2*gfY2mu(t)*dIb1ve+gfY1lambda(t)*(dIa1ve-dIa2ve-dIa3ve)+gfY2lambda(t)*(dIb1ve-dIb2ve-dIb3ve);
   tau22ns1=-2*gfY1mu(t)*dIa2ve-2*gfY2mu(t)*dIb2ve+gfY1lambda(t)*(dIa1ve-dIa2ve-dIa3ve)+gfY2lambda(t)*(dIb1ve-dIb2ve-dIb3ve);
   tau33ns1=-2*gfY1mu(t)*dIa3ve-2*gfY2mu(t)*dIb3ve+gfY1lambda(t)*(dIa1ve-dIa2ve-dIa3ve)+gfY2lambda(t)*(dIb1ve-dIb2ve-dIb3ve);
   tau13ns1=gfY1mu(t)*(dIa13ve+dIa31ve)+gfY2mu(t)*(dIb13ve+dIb31ve);
   
   tau11ns2=2*gfY1mu(t2)*dIa1ve+2*gfY2mu(t2)*dIb1ve+gfY1lambda(t2)*(dIa1ve-dIa2ve-dIa3ve)+gfY2lambda(t2)*(dIb1ve-dIb2ve-dIb3ve);
   tau22ns2=-2*gfY1mu(t2)*dIa2ve-2*gfY2mu(t2)*dIb2ve+gfY1lambda(t2)*(dIa1ve-dIa2ve-dIa3ve)+gfY2lambda(t2)*(dIb1ve-dIb2ve-dIb3ve);
   tau33ns2=-2*gfY1mu(t2)*dIa3ve-2*gfY2mu(t2)*dIb3ve+gfY1lambda(t2)*(dIa1ve-dIa2ve-dIa3ve)+gfY2lambda(t2)*(dIb1ve-dIb2ve-dIb3ve);
   tau13ns2=gfY1mu(t2)*(dIa13ve+dIa31ve)+gfY2mu(t2)*(dIb13ve+dIb31ve);


for ju=1:length(errE)
   matrixR(1,1)=st(ju)*cp(ju);
   matrixR(1,2)=ct(ju)*cp(ju);
   matrixR(1,3)=-sp(ju);
   matrixR(2,1)=st(ju)*sp(ju);
   matrixR(2,2)=ct(ju)*sp(ju);
   matrixR(2,3)=cp(ju);
   matrixR(3,1)=+ct(ju);
   matrixR(3,2)=-st(ju);
   matrixR(3,3)=0;
   
   matrixR2(1,1)=1;
   matrixR2(1,2)=0;
   matrixR2(1,3)=0;
   matrixR2(2,1)=0;
   matrixR2(2,2)=1;
   matrixR2(2,3)=0;
   matrixR2(3,1)=0;
   matrixR2(3,2)=0;
   matrixR2(3,3)=-1;
   

  e_polar=[errE(ju),ertE(ju),0;ertE(ju),ettE(ju),0;0,0,effE(ju)];
  e_cartesianA=matrixR*e_polar*matrixR';
  e_cartesian=matrixR2*e_cartesianA*matrixR2';
  
  e11smio(ju)=e_cartesian(1,1);
  e12smio(ju)=e_cartesian(1,2);
  e13smio(ju)=e_cartesian(1,3);
  e21smio(ju)=e_cartesian(2,1);
  e22smio(ju)=e_cartesian(2,2);
  e23smio(ju)=e_cartesian(2,3);
  e31smio(ju)=e_cartesian(3,1);
  e32smio(ju)=e_cartesian(3,2);
  e33smio(ju)=e_cartesian(3,3); 
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  
  e_polar=[err1(ju),ert1(ju),0;ert1(ju),ett1(ju),0;0,0,eff1(ju)];
  e_cartesianA=matrixR*e_polar*matrixR';
  e_cartesian=matrixR2*e_cartesianA*matrixR2';
  
  e11smio1(ju)=e_cartesian(1,1);
  e12smio1(ju)=e_cartesian(1,2);
  e13smio1(ju)=e_cartesian(1,3);
  e21smio1(ju)=e_cartesian(2,1);
  e22smio1(ju)=e_cartesian(2,2);
  e23smio1(ju)=e_cartesian(2,3);
  e31smio1(ju)=e_cartesian(3,1);
  e32smio1(ju)=e_cartesian(3,2);
  e33smio1(ju)=e_cartesian(3,3); 
  
  
  e_polar=[err2(ju),ert2(ju),0;ert2(ju),ett2(ju),0;0,0,eff2(ju)];
  e_cartesianA=matrixR*e_polar*matrixR';
  e_cartesian=matrixR2*e_cartesianA*matrixR2';
  
  e11smio2(ju)=e_cartesian(1,1);
  e12smio2(ju)=e_cartesian(1,2);
  e13smio2(ju)=e_cartesian(1,3);
  e21smio2(ju)=e_cartesian(2,1);
  e22smio2(ju)=e_cartesian(2,2);
  e23smio2(ju)=e_cartesian(2,3);
  e31smio2(ju)=e_cartesian(3,1);
  e32smio2(ju)=e_cartesian(3,2);
  e33smio2(ju)=e_cartesian(3,3); 
  
  u_polar=[urE(ju),utE(ju),0];
  
  u_polar1=[urt1(ju),utt1(ju),0];
  u_polar2=[urt2(ju),utt2(ju),0];
  
  u_cartesianA=matrixR*u_polar';
  u_cartesian=matrixR2*u_cartesianA;
  

  ux(ju)=u_cartesian(1);
  uy(ju)=u_cartesian(2);
  uz(ju)=u_cartesian(3);
  
   u_cartesianA=matrixR*u_polar1';
  u_cartesian=matrixR2*u_cartesianA;
  
  ux1(ju)=u_cartesian(1);
  uy1(ju)=u_cartesian(2);
  uz1(ju)=u_cartesian(3);
  
  u_cartesianA=matrixR*u_polar2';
  u_cartesian=matrixR2*u_cartesianA;
  
   ux2(ju)=u_cartesian(1);
   uy2(ju)=u_cartesian(2);
   uz2(ju)=u_cartesian(3);
%   
%   tau_polar=[taurr(ju),0,taurt(ju);0,tauff(ju),0;taurt(ju),0,tautt(ju)];
%   tau_cartesian=matrixR'*tau_polar*matrixR;  
  
  tau_polar=[taurr(ju),taurt(ju),0;taurt(ju),tautt(ju),0;0,0,tauff(ju)];
  tau_cartesianA=matrixR*tau_polar*matrixR';
  tau_cartesian=matrixR2*tau_cartesianA*matrixR2';
  
  
  tau11smio(ju)=tau_cartesian(1,1);
  tau12smio(ju)=tau_cartesian(1,2);
  tau13smio(ju)=tau_cartesian(1,3);
  tau21smio(ju)=tau_cartesian(2,1);
  tau22smio(ju)=tau_cartesian(2,2);
  tau23smio(ju)=tau_cartesian(2,3);
  tau31smio(ju)=tau_cartesian(3,1);
  tau32smio(ju)=tau_cartesian(3,2);
  tau33smio(ju)=tau_cartesian(3,3); 
  
  
  tau_polar=[taurr1(ju),taurt1(ju),0;taurt1(ju),tautt1(ju),0;0,0,tauff1(ju)];
  tau_cartesianA=matrixR*tau_polar*matrixR';
  tau_cartesian=matrixR2*tau_cartesianA*matrixR2';
  
  
  tau11smio1(ju)=tau_cartesian(1,1);
  tau12smio1(ju)=tau_cartesian(1,2);
  tau13smio1(ju)=tau_cartesian(1,3);
  tau21smio1(ju)=tau_cartesian(2,1);
  tau22smio1(ju)=tau_cartesian(2,2);
  tau23smio1(ju)=tau_cartesian(2,3);
  tau31smio1(ju)=tau_cartesian(3,1);
  tau32smio1(ju)=tau_cartesian(3,2);
  tau33smio1(ju)=tau_cartesian(3,3); 
  
  tau_polar=[taurr2(ju),taurt2(ju),0;taurt2(ju),tautt2(ju),0;0,0,tauff2(ju)];
  tau_cartesianA=matrixR*tau_polar*matrixR';
  tau_cartesian=matrixR2*tau_cartesianA*matrixR2';
  
 
  tau11smio2(ju)=tau_cartesian(1,1);
  tau12smio2(ju)=tau_cartesian(1,2);
  tau13smio2(ju)=tau_cartesian(1,3);
  tau21smio2(ju)=tau_cartesian(2,1);
  tau22smio2(ju)=tau_cartesian(2,2);
  tau23smio2(ju)=tau_cartesian(2,3);
  tau31smio2(ju)=tau_cartesian(3,1);
  tau32smio2(ju)=tau_cartesian(3,2);
  tau33smio2(ju)=tau_cartesian(3,3); 
  
end


ux1ve=ux1+uxns1;
uz1ve=uz1+uzns1;

 ux2ve=ux2+uxns2;
 uz2ve=uz2+uzns2;


 u1=(ux+u1ns);
 %u2=(uy+u2ns);
 u3=(uz+u3ns);


  if abs(z)<=high/2
       
 tau11smio(x<a)=tau11smio(x<a)+gfAlambda(0)*2*a/high;
 tau22smio(x<a)=tau22smio(x<a)+gfAlambda(0)*2*a/high;
 tau33smio(x<a)=tau33smio(x<a)+2*gfAmu(0)*2*a/high+gfAlambda(0)*2*a/high;  
 
 tau11smio1(x<a)=tau11smio1(x<a)+gfAlambda(t)*2*a/high;
 tau22smio1(x<a)=tau22smio1(x<a)+gfAlambda(t)*2*a/high;
 tau33smio1(x<a)=tau33smio1(x<a)+2*gfAmu(t)*2*a/high+gfAlambda(t)*2*a/high;  
 
 tau11smio2(x<a)=tau11smio2(x<a)+gfAlambda(t2)*2*a/high;
 tau22smio2(x<a)=tau22smio2(x<a)+gfAlambda(t2)*2*a/high;
 tau33smio2(x<a)=tau33smio2(x<a)+2*gfAmu(t2)*2*a/high+gfAlambda(t2)*2*a/high;    
       
 e33smio(x<a)=e33smio(x<a)+2*A*a/high;
 e33smio1(x<a)=e33smio1(x<a)+2*gf(t)*a/high;
 e33smio2(x<a)=e33smio2(x<a)+2*gf(t2)*a/high;
 
 end

   e11=e11smio+e11ns;
   e22=e22smio+e22ns;
   e33=e33smio+e33ns;
   e13=e13smio+e13ns;
  
   e111ve=e11smio1+e11ns1;
   e221ve=e22smio1+e22ns1;
   e331ve=e33smio1+e33ns1;
   e131ve=e13smio1+e13ns1;
%    
   e112ve=e11smio2+e11ns2;
   e222ve=e22smio2+e22ns2;
   e332ve=e33smio2+e33ns2;
   e132ve=e13smio2+e13ns2;
   
   tau11ve=tau11smio+tau11ns;
   tau22ve=tau22smio+tau22ns;
   tau33ve=tau33smio+tau33ns;
   tau13ve=tau13smio+tau13ns;
 
   
   tau111ve=tau11smio1+tau11ns1;
   tau221ve=tau22smio1+tau22ns1;
   tau331ve=tau33smio1+tau33ns1;
   tau131ve=tau13smio1+tau13ns1;
   
   tau112ve=tau11smio2+tau11ns2;
   tau222ve=tau22smio2+tau22ns2;
   tau332ve=tau33smio2+tau33ns2;
   tau132ve=tau13smio2+tau13ns2;
   
   
 for i=1:length(x)
  if abs(z)<=high/2 && x(i)<=a
 tau11ve(i)=tau11ve(i)-3*K*e0;
 tau22ve(i)=tau22ve(i)-3*K*e0;
 tau33ve(i)=tau33ve(i)-3*K*e0;
 
 tau111ve(i)=tau111ve(i)-3*K*e0;
 tau221ve(i)=tau221ve(i)-3*K*e0;
 tau331ve(i)=tau331ve(i)-3*K*e0;

 tau112ve(i)=tau112ve(i)-3*K*e0;
 tau222ve(i)=tau222ve(i)-3*K*e0;
 tau332ve(i)=tau332ve(i)-3*K*e0;

 end
end

   

for i=1:length(x)
  if abs(z)<=high/2 
 if x(i)<=a
 tau11(i)=2*mu.*(e11(i))+lambda*(e11(i)+e22(i)+e33(i))-3*K*e0;
 tau22(i)=2*mu.*(e22(i))+lambda*(e11(i)+e22(i)+e33(i))-3*K*e0;
 tau33(i)=2*mu.*(e33(i))+lambda*(e11(i)+e22(i)+e33(i))-3*K*e0;
 tau13(i)=2*mu.*(e13smio(i));
 
 tau11s(i)=2*mu.*(e11smio(i))+lambda*(e11smio(i)+e22smio(i)+e33smio(i))-3*K*e0;
 tau22s(i)=2*mu.*(e22smio(i))+lambda*(e11smio(i)+e22smio(i)+e33smio(i))-3*K*e0;
 tau33s(i)=2*mu.*(e33smio(i))+lambda*(e11smio(i)+e22smio(i)+e33smio(i))-3*K*e0;


elseif x(i)>a
 tau11(i)=2*mu.*(e11(i))+lambda*(e11(i)+e22(i)+e33(i));
 tau22(i)=2*mu.*(e22(i))+lambda*(e11(i)+e22(i)+e33(i));
 tau33(i)=2*mu.*(e33(i))+lambda*(e11(i)+e22(i)+e33(i));
 tau13(i)=2*mu.*(e13smio(i));

tau11s(i)=2*mu.*(e11smio(i))+lambda*(e11smio(i)+e22smio(i)+e33smio(i));
 tau22s(i)=2*mu.*(e22smio(i))+lambda*(e11smio(i)+e22smio(i)+e33smio(i));
 tau33s(i)=2*mu.*(e33smio(i))+lambda*(e11smio(i)+e22smio(i)+e33smio(i));
 end
  else
  %TAU OUT    
 tau11(i)=2*mu.*(e11(i))+lambda*(e11(i)+e22(i)+e33(i));
 tau22(i)=2*mu.*(e22(i))+lambda*(e11(i)+e22(i)+e33(i));
 tau33(i)=2*mu.*(e33(i))+lambda*(e11(i)+e22(i)+e33(i));
 tau13(i)=2*mu.*(e13(i));

 tau11s(i)=2*mu.*(e11smio(i))+lambda*(e11smio(i)+e22smio(i)+e33smio(i));
 tau22s(i)=2*mu.*(e22smio(i))+lambda*(e11smio(i)+e22smio(i)+e33smio(i));
 tau33s(i)=2*mu.*(e33smio(i))+lambda*(e11smio(i)+e22smio(i)+e33smio(i));
 
 end
end

end

