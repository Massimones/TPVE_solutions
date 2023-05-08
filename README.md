# Description of the scripts (TPVE solutions)



    The two scripts compute visco-elastic displacement, strain and stress due to a disk-shaped Thermo-Poro-Visco-Elastic
    (TPVE) inclusion   
    
    Input parameters are in S.I. The scripts are written in MATLAB (2022b).                                    
                                                                                                               
                   Massimo Nespoli 05/05/2023                                                                 
                   E-mail: massimo.nespoli2@unibo.it                                                          
                                                                                                               
    If you use these scripts please cite the paper:                                                            
                                                                                                              
    Nespoli Massimo, Belardinelli Maria Elina, Bonafede Maurizio (2023). Thermo-poro-visco-elastic             
    response of a disk shaped inclusion. Geophysical Journal International.                                    

                                                                                                        
## The first script "Compute_solutions.m" must be executed to call the function "TPVE_SOLUTIONS_DISK.m"                                                                                                             

### "Compute_solutions.m"

The script computes displacements, strain and stresses fields for a disk-shaped TPE inclusion with an assigned geometry.
The fields are computed along different "x" profiles from x=0 to x=limiteplot with a step of k(m), for different "z",
 from z=Zmin to z=Zmax with a step of Zstep (m). If Zmin=Zmax only one x-profile will be computed 
(e.g. if Zmin=Zmax=0 the script computes the fields at the surface). The input parameters are reported at the beginning 
of the script.
Their description is given below:

**H**:     &nbsp; &nbsp; &nbsp;          Constant of Biot (Pa)  <br />
**alfa**:    &nbsp; &nbsp; &nbsp;        	 Thermal expansion (1/K) <br />
**dp**:      &nbsp; &nbsp; &nbsp;          Pore pressure change inside the disk (Pa) <br />
**dT**:     &nbsp; &nbsp; &nbsp;           Temperature change inside the disk (°C or K) <br />
**a**:      &nbsp; &nbsp; &nbsp;           Disk radius (m) <br />
**db**:     &nbsp; &nbsp; &nbsp;           Disk height, thickness (m) <br />
**mu**:     &nbsp; &nbsp; &nbsp;           Shear modulus (Pa) <br />
**lambda**:   &nbsp; &nbsp; &nbsp;         Lamè constant (Pa) <br />
**MedianPlane**:  &nbsp; &nbsp; &nbsp;     Depth of the median plane of the TPE inclusion (m) <br />
**limiteplot**:   &nbsp; &nbsp; &nbsp;     Limit to compute the fields along a x profile (m) <br />
**k**:      &nbsp; &nbsp; &nbsp;           Step interval for plot in x (m) <br />
**Zmin**:    &nbsp; &nbsp; &nbsp;          minimum z for the computation of fields (m) <br />
**Zmax**:     &nbsp; &nbsp; &nbsp;         maximum z for the computation of fields (m) <br />
**Zstep**:    &nbsp; &nbsp; &nbsp;         Step interval for plot in z (m) <br />
**eta**:  &nbsp; &nbsp; &nbsp;             Viscosity of the Maxwell body (Pa x s) <br />
**TV**:   &nbsp; &nbsp; &nbsp;             t1/tau (a first, normalized, instant of time to compute fields) <br />
**TV2**:    &nbsp; &nbsp; &nbsp;           t2/tau (a second, normalized, instant of time to compute fields) <br />
**plotshow**:   &nbsp; &nbsp; &nbsp;       plotshow=n, where Zetav(n) is the depth to plot in the three following figures <br />

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

###  EXAMPLE: Below there is the choice of parameters which was used in Nespoli et al. (2023). 
Figures of displacement, strain and stress obtained using such paramters are located in the directory "Example_figures".

**H**=10x10^9;               
**alfa**=3x10^(-5);         
**dp**=10e6;                 
**dT**=100;                  
**a**=2500;                  
**db**=500;                 
**mu**=6x10^9;               
**lambda**=4x10^9;           
**MedianPlane**=2000;        
**limiteplot**=8000;        
**k**=25;                    
**Zmin**=0;                  
**Zmax**=0;                  
**Zstep**=25;                
**eta**=10^16;               
**TV**=2;                    
**TV2**=0.5;                 
**plotshow**=1;       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                                        
## The second script "TPVE_SOLUTIONS_DISK.m" is the function, it does not have to be executed directly !                                                                                                               

### "TPVE_SOLUTIONS_DISK.m"

The function solves the equations presented in Nespoli et al. (2023) which includes some equations previously 
derived by Belardinelli et al. (2022) and Mantiloni et al (2020). See the main paper references.

The function takes as input the parameters specified in the "compute_solutions.m" file and gives as output the different
 components of the computed displacement, strain and stress fields (Elastic solution, Visco-elastic solution computed
 at TV=t1/tau, Visco-elastic solution computed at TV2=t2/tau, where tau=eta/mu).

***[x,e11,e22,e33,u1,u3,e111ve,e221ve,e331ve,ux1ve,uz1ve,e112ve,e222ve,e332ve,ux2ve,uz2ve,tau11ve,tau22ve,tau33ve,
tau111ve,tau221ve,tau331ve,tau112ve,tau222ve,tau332ve,e13,e131ve,e132ve,tau13ve,tau131ve,tau132ve]=TPVE_SOLUTIONS_DISK(H,alfa,dp,dT,a,db,ni,mu,lambda,MedianPlane,limiteplot,zlm,k,eta,TV,TV2)***

The Symbolic Math Toolbox is required to compute the Laplace transforms: https://it.mathworks.com/products/symbolic.html

### Description of the OUTPUT variables:

**x:**   &nbsp; &nbsp; &nbsp;                                 x coordinates (m) <br />
**e11,e22,e33:**        &nbsp; &nbsp; &nbsp;                  Elastic strain components (11,22,33 i.e. xx, yy, zz) <br />
**u1,u3:**       &nbsp; &nbsp; &nbsp;                         Elastic displacements components (1,3 i.e. x, z) (m) <br />
**e111ve,e221ve,e331ve:**   &nbsp; &nbsp; &nbsp;              Visco-Elastic strain components (11,22,33 i.e. xx, yy, zz) computed at TV <br />
**ux1ve,uz1ve:**        &nbsp; &nbsp; &nbsp;                  Visco-Elastic displacement components (1,3 i.e. x, z) (m) computed at TV <br />
**e112ve,e222ve,e332ve:**   &nbsp; &nbsp; &nbsp;              Visco-Elastic strain components (11,22,33 i.e. xx, yy, zz) computed at TV2 <br />
**ux2ve,uz2ve:**         &nbsp; &nbsp; &nbsp;                 Visco-Elastic displacement components (1,3 i.e. x, z) (m) computed at TV2 <br />
**tau11ve,tau22ve,tau33ve:**   &nbsp; &nbsp; &nbsp;           Elastic stress components (11,22,33 i.e. xx, yy, zz) (Pa) computed at TV <br />
**tau111ve,tau221ve,tau331ve:**   &nbsp; &nbsp; &nbsp;        Visco-Elastic stress components (11,22,33 i.e. xx, yy, zz) (Pa) computed at TV <br />
**tau112ve,tau222ve,tau332ve:**  &nbsp; &nbsp; &nbsp;         Visco-Elastic stress components (11,22,33 i.e. xx, yy, zz) (Pa) computed at TV2 <br />
**e13:**         &nbsp; &nbsp; &nbsp;                         Elastic strain component (13 i.e. xz) <br />
**e131ve:**   &nbsp; &nbsp; &nbsp;                            Visco-Elastic strain component (13 i.e. xz) computed at TV <br />
**e132ve:**    &nbsp; &nbsp; &nbsp;                           Visco-Elastic strain component (13 i.e. xz) computed at TV2 <br />
**tau13ve:**    &nbsp; &nbsp; &nbsp;                          Elastic stress component (13 i.e. xz) (Pa) <br />
**tau131ve:**     &nbsp; &nbsp; &nbsp;                        Visco-Elastic stress component (13 i.e. xz) (Pa) computed at TV <br />
**tau132ve:**     &nbsp; &nbsp; &nbsp;                        Visco-Elastic stress component (13 i.e. xz) (Pa) computed at TV2 <br />
