! This is the namelist file for ice part

&ice_stress
Pstar = 15000.0			!20000. 27500. 30000.
c_pressure = 20.0		!20.0
/

&ice_fric
Cd_oce_ice = 5.0e-3    		!5.5e-3  3.0e-3
Kh_ice=0.0 
/

&ice_rheology
EVP_rheology=.true.
evp_rheol_steps=40	 
evp_Tdamp_ratio=3		!ratio dt/T_damp
vp_rheol_steps=500   
/

&ice_scheme
ice_gamma_fct=0.2
/

&ice_therm
Sice=5.0
h0=1.0
emiss_ice=0.97
emiss_wat=0.97
albsn=0.85  
albsnm=0.75   
albi=0.75
albim=0.66    
albw=0.1
/
