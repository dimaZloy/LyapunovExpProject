
using PyPlot;
using Dates;
using Printf;
using DelimitedFiles;
using Statistics; 
using HDF5;
using FFTW;
using DSP;

include("utils.jl")
include("freeStream.jl")
include("computeFFT.jl")

###################################################################################

debug1 = true;
if (debug1)

	HC00_T, HC00_CD, HC00_CL =   readForcesDataOF3("hcyl2d_00.dat",374530);
	HC01_T, HC01_CD, HC01_CL =   readForcesDataOF3("hcyl2d_01.dat",425509);
	HC02_T, HC02_CD, HC02_CL =   readForcesDataOF3("hcyl2d_02.dat",357704);
	HC01_BL_T, HC01_BL_CD, HC01_BL_CL =   readForcesDataOF3("hcyl2d_01_BL.dat",240589);
	HC01_WN1_T, HC01_WN1_CD, HC01_WN1_CL =   readForcesDataOF3("hcyl2d_01_WN1.dat",234779);
	
end

N00 = 100000;
L00 = 0.107703;
Area00 = D*L00;
A00 = A*Area00;
A01 = A*Area00;
A02 = A*Area00;

HC00_T1 = deepcopy(HC00_T[N00:end]); 
HC00_CD1 = deepcopy(HC00_CD[N00:end])./A00; 
HC00_CL1 = deepcopy(HC00_CL[N00:end])./A00;


HC01_T1 = deepcopy(HC01_T[N00:end]); 
HC01_CD1 = deepcopy(HC01_CD[N00:end])./A01; 
HC01_CL1 = deepcopy(HC01_CL[N00:end])./A01;


HC01_BL_T1 = deepcopy(HC01_BL_T[56609:end]); 
HC01_BL_CD1 = deepcopy(HC01_BL_CD[56609:end])./A01; 
HC01_BL_CL1 = deepcopy(HC01_BL_CL[56609:end])./A01;

HC01_WN1_T1 = deepcopy(HC01_WN1_T[56609:end]); 
HC01_WN1_CD1 = deepcopy(HC01_WN1_CD[56609:end])./A01; 
HC01_WN1_CL1 = deepcopy(HC01_WN1_CL[56609:end])./A01;


HC02_T1 = deepcopy(HC02_T[N00:end]); 
HC02_CD1 = deepcopy(HC02_CD[N00:end])./A02; 
HC02_CL1 = deepcopy(HC02_CL[N00:end])./A02;


freqHC00, pHC00, maxFreqHC00, ST_HC00 = computeFFT(HC00_T1,HC00_CL1,D,UInf);
freqHC01, pHC01, maxFreqHC01, ST_HC01 = computeFFT(HC01_T1,HC01_CL1,D,UInf);
freqHC02, pHC02, maxFreqHC02, ST_HC02 = computeFFT(HC02_T1,HC02_CL1,D,UInf);
freqHC01_BL, pHC01_BL, maxFreqHC01_BL, ST_HC01_BL = computeFFT(HC01_BL_T1,HC01_BL_CL1,D,UInf);
freqHC01_WN1, pHC01_WN1, maxFreqHC01_WN1, ST_HC01_WN1 = computeFFT(HC01_WN1_T1,HC01_WN1_CL1,D,UInf);

meanCD1_H00 = mean(HC00_CD1);
meanCL1_H00 = mean(HC00_CL1);

meanCD1_H01 = mean(HC01_CD1);
meanCL1_H01 = mean(HC01_CL1);

meanCD1_H02 = mean(HC02_CD1);
meanCL1_H02 = mean(HC02_CL1);

meanCD1_H01_BL = mean(HC01_BL_CD1);
meanCL1_H01_BL = mean(HC01_BL_CL1);

meanCD1_H01_WN1 = mean(HC01_WN1_CD1);
meanCL1_H01_WN1 = mean(HC01_WN1_CL1);


debug2 = true;





if (debug2)


	yminCL =  -1.0;
	ymaxCL  =  2.0;

	yminCD =   1.0;
	ymaxCD  =  2.0;

	tmin = 0.1;
	tmax = 0.374516;

	

    lw2 = 0.5;
	lw1 = 1.0;
	

	figure(10);
	clf();

	subplot(2,1,1)

	plot(HC00_T1, HC00_CL1, "-.r", label="HC00", linewidth = lw2);
	plot(HC01_T1, HC01_CL1, "-.b", label="HC01", linewidth = lw2);
	plot(HC02_T1, HC02_CL1, "-.g", label="HC02", linewidth = lw2);
	plot(HC01_BL_T1, HC01_BL_CL1, "-.c", label="HC01_BL", linewidth = lw2);
	plot(HC01_WN1_T1, HC01_WN1_CL1, "-.m", label="HC01_WN1", linewidth = lw2);

	xlabel("t[s]");
	ylabel("Cl");
	#xlim(tmin, tmax);
	#ylim(yminCL, ymaxCL);
	legend(fontsize=8);
	grid();
	
	subplot(2,1,2)
	plot(HC00_T1, HC00_CD1, "-.r", label="HC00", linewidth = lw2);
	plot(HC01_T1, HC01_CD1, "-.b", label="HC01", linewidth = lw2);
	plot(HC02_T1, HC02_CD1, "-.g", label="HC02", linewidth = lw2);
	plot(HC01_BL_T1, HC01_BL_CD1, "-.c", label="HC01_BL", linewidth = lw2);
	plot(HC01_WN1_T1, HC01_WN1_CD1, "-.m", label="HC01_WN1", linewidth = lw2);

	xlabel("t[s]");
	ylabel("Cd");
	#xlim(tmin, tmax);
	#ylim(yminCD, ymaxCD);
	legend(fontsize=8);
	grid();


end


println("Mean aerodynamic characteristics");
println("EXP (Re=50000)", "Cd:\t", 0.500);
println("EXP           ", "Cd:\t", 0.370);
println("EXP           ", "Cd:\t", 0.467);
println("EXP           ", "Cd:\t", 0.476);
println("HC00          ", "Cd:\t", round(meanCD1_H00,digits=3));
println("HC01          ", "Cd:\t", round(meanCD1_H01,digits=3));
println("HC02          ", "Cd:\t", round(meanCD1_H02,digits=3));
println("HC01_BL       ", "Cd:\t", round(meanCD1_H01_BL,digits=3));
println("HC01_WN1      ", "Cd:\t", round(meanCD1_H01_WN1,digits=3));

println("=====================================");
println("EXP           ", "Cl:\t", -1.100);
println("EXP           ", "Cl:\t", -0.540);
println("EXP           ", "Cl:\t", -0.634);
println("EXP           ", "Cl:\t", -0.816);
println("HC00          ", "Cl:\t", round(meanCL1_H00,digits=3));
println("HC01          ", "Cl:\t", round(meanCL1_H01,digits=3));
println("HC02          ", "Cl:\t", round(meanCL1_H02,digits=3));
println("HC01_BL       ", "Cl:\t", round(meanCL1_H01_BL,digits=3));
println("HC01_WN1      ", "Cl:\t", round(meanCL1_H01_WN1,digits=3));


println("=====================================");
println("EXP           ", "St :\t", 0.37);
println("EXP           ", "St :\t", 0.41);
println("EXP           ", "St :\t", 0.44);
println("HC00          ", "St :\t", round(ST_HC00,digits=2));
println("HC01          ", "St :\t", round(ST_HC01,digits=2));
println("HC02          ", "St :\t", round(ST_HC02,digits=2));
println("HC01_BL       ", "St :\t", round(ST_HC01_BL,digits=2));
println("HC01_WN1      ", "St :\t", round(ST_HC01_WN1,digits=2));
