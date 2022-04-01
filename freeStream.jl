

function calcSutherlandViscosityOF(T::Float64)::Float64

	As = 1.4792e-6;
	Ts = 116;
	return As*sqrt(T)/(1.0 + Ts/T);

end


## U.S. Standard Atmosphere Air Properties - SI Units 
## h = 5000  (m) Geo potential Altitude above Sea Level 
## T = -17.47 degC 
## P = 5.405e+4 N/m2 
## rho = 0.7364 kg/m3 
## mu = 1.628e-5 Pa/s

D = 0.1
L = 1.0
gamma  = 1.4
PInf = 2.5605e+4
R = 286.7
TInf = 273.0-50.47
UInf = 35.0

muInf = calcSutherlandViscosityOF(TInf); # inlet dynamic viscosity [Pa*s]
##muInf = 2e-5;
rhoInf = PInf/R/TInf

nuInf = muInf/rhoInf

gamma = 1.4;
aInf = sqrt(gamma*PInf/rhoInf);
MInf = UInf/aInf;
ReInf = rhoInf*UInf*D/muInf ;

#Aref =  D*L; # [m2]  
A = rhoInf*UInf*UInf;



println("The turbulent flow over a circular cylinder");
println("Free-stream conditions:");
@printf("===================================\n");
@printf("PInf:\t%.2f\n", PInf);
@printf("TInf:\t%.2f\n", TInf);
@printf("rhoInf:\t%.4f\n", rhoInf);
@printf("muInf:\t%.7f\n", muInf);
@printf("UInf:\t%.3f\n", UInf);
@printf("nutInf:\t%.10f\n", 5.0*nuInf);
@printf("Re:\t%.0f\n", ReInf);
@printf("Ma:\t%.2f\n", MInf);
@printf("D[m]:\t%.2f\n", D);

