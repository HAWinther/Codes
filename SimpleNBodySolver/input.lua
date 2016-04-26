-----------------------------------------
--            Inputfile                --
-----------------------------------------
                
-- Particlefile
particlefile = "ic/part.ascii" 

-- Starting redshift
aexp_ini = 0.1

-- Static test of solvers
static   = false

-- Boxsize in Mpc/h
BoxSize  = 200.0   

-- Omega_m
omega_m  = 0.267       

-- Size of grid used to solve the Poisson and scalar field eq. on
ngrid    = 128

-- Courant-fac for time-step control
courant_fac = 0.05

-- Modified gravity run?
modified_gravity = false

-- f(R) parameters
nfofr    = 1.0      
fofr0    = 1.0e-5 

-- Symmetron parameters
assb_symm   = 0.5  
beta_symm   = 1.0 
lambda_symm = 5.0 

-- Multigrid parameters
NGS_FINE     = 5
NGS_COARSE   = 5   
EPS_CONVERGE = 1e-5  
verboseMG    = false  

