# TPVE_solutions
Matlab scripts to compute Thermo-Poro-Visco-Elastic displacement, stress and strain due to a disk-shaped inclusion

&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
&&&&&&&                                                                                                               &&&&&&&
&&&&&&&    The two script generate visco-elastic displacements, strain and stresses for a disk-shaped TPE inclusion   &&&&&&&
&&&&&&&    Input parameters are in S.I. The scripts are written in MATLAB (2022b).                                    &&&&&&&
&&&&&&&                                                                                                               &&&&&&&
&&&&&&&                    Massimo Nespoli 05/05/2023                                                                 &&&&&&&
&&&&&&&                    E-mail: massimo.nespoli2@unibo.it                                                          &&&&&&&
&&&&&&&                                                                                                               &&&&&&&
&&&&&&&    If you use these scripts please cite the paper:                                                            &&&&&&&
&&&&&&&                                                                                                               &&&&&&&
&&&&&&&    Nespoli Massimo, Belardinelli Maria Elina, Bonafede Maurizio (2023). Thermo-poro-visco-elastic             &&&&&&&
&&&&&&&    response of a disk shaped inclusion. Geophysical Journal International.                                    &&&&&&&
&&&&&&&                                                                                                               &&&&&&&
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
&&&&&&&                                                                                                               &&&&&&&
&&&&&&&     The first script "Compute_solutions.m" must be executed to call the function "TPVE_SOLUTIONS_DISK.m"      &&&&&&&
&&&&&&&                                                                                                               &&&&&&&
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

%%%Compute_solutions.m%%%

The script computes displacements, strain and stresses fields for a disk-shaped TPE inclusion with an assigned geometry.
The fields are computed along different "x" profiles from x=0 to x=limiteplot with a step of k(m), for different "z",
 from z=Zmin to z=Zmax with a step of Zstep (m). If Zmin=Zmax only one x-profile will be computed 
(e.g. if Zmin=Zmax=0 the script computes the fields at the surface). The input parameters are reported at the beginning 
of the script.
Their description is given below:

H:                         Constant of Biot (Pa)
alfa:                      Thermal expansion (1/K)
dp:                        Pore pressure change inside the disk (Pa)
dT:                        Temperature change inside the disk (°C or K)
a:                         Disk radius (m)
db:                        Disk height, thickness (m)
mu:                        Shear modulus (Pa)
lambda:                    Lamè constant (Pa)
MedianPlane:               Depth of the median plane of the TPE inclusion (m)
limiteplot:                Limit to compute the fields along a x profile (m)
k:                         Step interval for plot in x (m)
Zmin:                      minimum z for the computation of fields (m)
Zmax:                      maximum z for the computation of fields (m)
Zstep:                     Step interval for plot in z (m)
eta:                       Viscosity of the Maxwell body (Pa*s)
TV:                        t1/tau (a first, normalized, instant of time to compute fields as a function of time)
TV2:                       t2/tau (a second, normalized, instant of time to compute fields as a function of time)
plotshow:                  plotshow=n, where Zetav(n) is the depth to plot in the three following figures

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EXAMPLE: Below there is the choice of parameters which was used in Nespoli et al. (2023). 
Figures of displacement, strain and stress obtained using such paramters are located in the directory "Example_figures".

H=10*10^9;               
alfa=3*10^(-5);         
dp=10e6;                 
dT=100;                  
a=2500;                  
db=500;                 
mu=6*10^9;               
lambda=4*10^9;           
MedianPlane=2000;        
limiteplot=8000;        
k=25;                    
Zmin=0;                  
Zmax=0;                  
Zstep=25;                
eta=10^16;               
TV=2;                    
TV2=0.5;                 
plotshow=1;              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
&&&&&&&                                                                                                               &&&&&&&
&&&&&&&       The second script "TPVE_SOLUTIONS_DISK.m" is the function, it does not have to be executed directly !   &&&&&&&
&&&&&&&                                                                                                               &&&&&&&
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

%%%"TPVE_SOLUTIONS_DISK.m"%%%

The function solves the equations presented in Nespoli et al. (2023) which includes some equations previously 
derived by Belardinelli et al. (2022) and Mantiloni et al (2020). See the main paper references.

The function takes as input the parameters specified in the "compute_solutions.m" file and gives as output the different
 components of the computed displacement, strain and stress fields (Elastic solution, Visco-elastic solution computed
 at TV=t1/tau, Visco-elastic solution computed at TV2=t2/tau, where tau=eta/mu).

[x,e11,e22,e33,u1,u3,e111ve,e221ve,e331ve,ux1ve,uz1ve,e112ve,e222ve,e332ve,ux2ve,uz2ve,tau11ve,tau22ve,tau33ve,
tau111ve,tau221ve,tau331ve,tau112ve,tau222ve,tau332ve,e13,e131ve,e132ve,tau13ve,tau131ve,tau132ve]=TPVE_SOLUTIONS_DISK(H,alfa,dp,dT,a,db,ni,mu,lambda,MedianPlane,limiteplot,zlm,k,eta,TV,TV2)

The Symbolic Math Toolbox is required to compute the Laplace transforms: https://it.mathworks.com/products/symbolic.html

Description of OUTPUT variables:

x,                                  x coordinates (m)
e11,e22,e33,                        Elastic strain components (11,22,33 i.e. xx, yy, zz)
u1,u3,                              Elastic displacements components (1,3 i.e. x, z) (m)
e111ve,e221ve,e331ve,               Visco-Elastic strain components (11,22,33 i.e. xx, yy, zz) computed at TV
ux1ve,uz1ve,                        Visco-Elastic displacement components (1,3 i.e. x, z) (m) computed at TV
e112ve,e222ve,e332ve,               Visco-Elastic strain components (11,22,33 i.e. xx, yy, zz) computed at TV2
ux2ve,uz2ve,                        Visco-Elastic displacement components (1,3 i.e. x, z) (m) computed at TV2
tau11ve,tau22ve,tau33ve,            Elastic stress components (11,22,33 i.e. xx, yy, zz) (Pa) computed at TV
tau111ve,tau221ve,tau331ve,         Visco-Elastic stress components (11,22,33 i.e. xx, yy, zz) (Pa) computed at TV
tau112ve,tau222ve,tau332ve,         Visco-Elastic stress components (11,22,33 i.e. xx, yy, zz) (Pa) computed at TV2
e13,                                Elastic strain component (13 i.e. xz)
e131ve,                             Visco-Elastic strain component 13 (13 i.e. xz) computed at TV
e132ve,                             Visco-Elastic strain component 13 (13 i.e. xz) computed at TV2
tau13ve,                            Elastic stress component 13 (13 i.e. xz) (Pa)
tau131ve,                           Visco-Elastic stress component 13 (13 i.e. xz) (Pa) computed at TV
tau132ve                            Visco-Elastic stress component 13 (13 i.e. xz) (Pa) computed at TV2

&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
