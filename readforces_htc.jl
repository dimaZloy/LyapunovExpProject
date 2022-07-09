
using PyPlot;
using Dates;
using Printf;
using DelimitedFiles;
using Statistics; 
using HDF5;
using FFTW;
using DSP;
using CSV;
using DataFrames;

include("utils.jl")
include("freeStream.jl")
include("computeFFT.jl")


###################################################################################

# DL> 9-July-2022
# 2D URANS of the flow past semicircular cylinder at Re = 100,000 and M = 0.12
# simple SA (w/o strain-worticity correction) is used!!!
# Goal is to investigate heat transfer effects 
# for cases are considered, when the difference between cylinder surface temperature 
# and free-stream temperature>
# a) delta T = 0 [K]
# b) delta T = 100 [K]
# c) delta T = 200 [K]
# d) delta T = 300 [K]
# Basic resuls and conclusions with dT -> from 0 to 300 K: 
# a) lift coeff goes down
# b) drag coeff goes up 
# c) Lr recirculation zone is increasing
# d) wake dynamics is decreasing> St goes down 

#youtube> https://youtu.be/kfIDH4F_ulo 

###################################################################################

debug1 = true;
if (debug1)

	# HC00_T, HC00_CD, HC00_CL =   readForcesDataOF3("hcyl2d_00.dat",374530);
	HC01_T, HC01_CD, HC01_CL =   readForcesDataOF3("hcyl2d_01.dat",425509);
	# HC02_T, HC02_CD, HC02_CL =   readForcesDataOF3("hcyl2d_02.dat",357704);
	# HC01_BL_T, HC01_BL_CD, HC01_BL_CL =   readForcesDataOF3("hcyl2d_01_BL.dat",240589);
	# HC01_CY_T, HC01_CY_CD, HC01_CY_CL =   readForcesDataOF3("hcyl2d_01_clarkY.dat",256235);
	# HC01_WN1_T, HC01_WN1_CD, HC01_WN1_CL =   readForcesDataOF3("hcyl2d_01_WN1.dat",300014);
	# HC01_WN5_T, HC01_WN5_CD, HC01_WN5_CL =   readForcesDataOF3("hcyl2d_01_WN5.dat",300014);
	# HC01_WN10_T, HC01_WN10_CD, HC01_WN10_CL =   readForcesDataOF3("hcyl2d_01_WN10.dat",600014);
	
	HC01_isoT_T, HC01_isoT_CD, HC01_isoT_CL =   readForcesDataOF3("hcyl2d_01_isoT.dat",541072);
	HC01_isoTd100K_T, HC01_isoTd100K_CD, HC01_isoTd100K_CL =   readForcesDataOF3("hcyl2d_01_isoTd100K.dat",400014);
	HC01_isoTd200K_T, HC01_isoTd200K_CD, HC01_isoTd200K_CL =   readForcesDataOF3("hcyl2d_01_isoTd200K.dat",400014);
	HC01_isoTd300K_T, HC01_isoTd300K_CD, HC01_isoTd300K_CL =   readForcesDataOF3("hcyl2d_01_isoTd300K.dat",400114);
	
end

debug2 = true;
if (debug2)

	HC01_isoTd000K_U = CSV.read("dt000k_uaxial.csv", DataFrame);
	HC01_isoTd100K_U = CSV.read("dt100k_uaxial.csv", DataFrame);
	HC01_isoTd200K_U = CSV.read("dt200k_uaxial.csv", DataFrame);
	HC01_isoTd300K_U = CSV.read("dt300k_uaxial.csv", DataFrame);

	HC01_isoTd000K_Uwake = CSV.read("dt000k_uwake.csv", DataFrame);
	HC01_isoTd100K_Uwake = CSV.read("dt100k_uwake.csv", DataFrame);
	HC01_isoTd200K_Uwake = CSV.read("dt200k_uwake.csv", DataFrame);
	HC01_isoTd300K_Uwake = CSV.read("dt300k_uwake.csv", DataFrame);

	
end

HC01_isoTd000K_t = HC01_isoTd000K_Uwake[:,1];
HC01_isoTd000K_Uy = HC01_isoTd000K_Uwake[:,3];
HC01_isoTd100K_t = HC01_isoTd100K_Uwake[:,1];
HC01_isoTd100K_Uy = HC01_isoTd100K_Uwake[:,3];
HC01_isoTd200K_t = HC01_isoTd200K_Uwake[:,1];
HC01_isoTd200K_Uy = HC01_isoTd200K_Uwake[:,3];
HC01_isoTd300K_t = HC01_isoTd300K_Uwake[:,1];
HC01_isoTd300K_Uy = HC01_isoTd300K_Uwake[:,3];





HC01_isoT_X = HC01_isoTd000K_U[:,29]./D;

HC01_isoTd000K_Ux = HC01_isoTd000K_U[:,5]./UInf;
HC01_isoTd000K_UxRMS = sqrt.(HC01_isoTd000K_U[:,8])./UInf;

HC01_isoTd100K_Ux = HC01_isoTd100K_U[:,5]./UInf;
HC01_isoTd100K_UxRMS = sqrt.(HC01_isoTd100K_U[:,8])./UInf;

HC01_isoTd200K_Ux = HC01_isoTd200K_U[:,5]./UInf;
HC01_isoTd200K_UxRMS = sqrt.(HC01_isoTd200K_U[:,8])./UInf;

HC01_isoTd300K_Ux = HC01_isoTd300K_U[:,5]./UInf;
HC01_isoTd300K_UxRMS = sqrt.(HC01_isoTd300K_U[:,8])./UInf;



N00 = 100000;
L00 = 0.107703;
Area00 = D*L00;
A00 = A*Area00;
A01 = A*Area00;
A02 = A*Area00;

# HC00_T1 = deepcopy(HC00_T[N00:end]); 
# HC00_CD1 = deepcopy(HC00_CD[N00:end])./A00; 
# HC00_CL1 = deepcopy(HC00_CL[N00:end])./A00;


HC01_T1 = deepcopy(HC01_T[N00:end]); 
HC01_CD1 = deepcopy(HC01_CD[N00:end])./A01; 
HC01_CL1 = deepcopy(HC01_CL[N00:end])./A01;


# HC01_BL_T1 = deepcopy(HC01_BL_T[56609:end]); 
# HC01_BL_CD1 = deepcopy(HC01_BL_CD[56609:end])./A01; 
# HC01_BL_CL1 = deepcopy(HC01_BL_CL[56609:end])./A01;

HC01_isoT_T1 = deepcopy(HC01_isoT_T[N00:end]); 
HC01_isoT_CD1 = deepcopy(HC01_isoT_CD[N00:end])./A01; 
HC01_isoT_CL1 = deepcopy(HC01_isoT_CL[N00:end])./A01;


HC01_isoTd100K_T1 = deepcopy(HC01_isoTd100K_T[N00:end]); 
HC01_isoTd100K_CD1 = deepcopy(HC01_isoTd100K_CD[N00:end])./A01; 
HC01_isoTd100K_CL1 = deepcopy(HC01_isoTd100K_CL[N00:end])./A01;


HC01_isoTd200K_T1 = deepcopy(HC01_isoTd200K_T[N00:end]); 
HC01_isoTd200K_CD1 = deepcopy(HC01_isoTd200K_CD[N00:end])./A01; 
HC01_isoTd200K_CL1 = deepcopy(HC01_isoTd200K_CL[N00:end])./A01;

HC01_isoTd300K_T1 = deepcopy(HC01_isoTd300K_T[N00:end]); 
HC01_isoTd300K_CD1 = deepcopy(HC01_isoTd300K_CD[N00:end])./A01; 
HC01_isoTd300K_CL1 = deepcopy(HC01_isoTd300K_CL[N00:end])./A01;


# HC01_CY_T1 = deepcopy(HC01_CY_T[N00:end]); 
# HC01_CY_CD1 = deepcopy(HC01_CY_CD[N00:end])./A01; 
# HC01_CY_CL1 = deepcopy(HC01_CY_CL[N00:end])./A01;


# HC01_WN1_T1 = deepcopy(HC01_WN1_T[1:end]); 
# HC01_WN1_CD1 = deepcopy(HC01_WN1_CD[1:end])./A01; 
# HC01_WN1_CL1 = deepcopy(HC01_WN1_CL[1:end])./A01;

# HC01_WN5_T1 = deepcopy(HC01_WN5_T[1:end]); 
# HC01_WN5_CD1 = deepcopy(HC01_WN5_CD[1:end])./A01; 
# HC01_WN5_CL1 = deepcopy(HC01_WN5_CL[1:end])./A01;

# HC01_WN10_T1 = deepcopy(HC01_WN10_T[1:end]); 
# HC01_WN10_CD1 = deepcopy(HC01_WN10_CD[1:end])./A01; 
# HC01_WN10_CL1 = deepcopy(HC01_WN10_CL[1:end])./A01;


# HC02_T1 = deepcopy(HC02_T[N00:end]); 
# HC02_CD1 = deepcopy(HC02_CD[N00:end])./A02; 
# HC02_CL1 = deepcopy(HC02_CL[N00:end])./A02;


# freqHC00, pHC00, maxFreqHC00, ST_HC00 = computeFFT(HC00_T1,HC00_CL1,D,UInf);
freqHC01, pHC01, maxFreqHC01, ST_HC01 = computeFFT(HC01_T1,HC01_CL1,D,UInf);
# freqHC02, pHC02, maxFreqHC02, ST_HC02 = computeFFT(HC02_T1,HC02_CL1,D,UInf);
# freqHC01_BL, pHC01_BL, maxFreqHC01_BL, ST_HC01_BL = computeFFT(HC01_BL_T1,HC01_BL_CL1,D,UInf);

#freqHC01_isoT, pHC01_isoT, maxFreqHC01_isoT, ST_HC01_isoT = computeFFT(HC01_isoT_T1,HC01_isoT_CL1,D,UInf);
freqHC01_isoT, pHC01_isoT, maxFreqHC01_isoT, ST_HC01_isoT = computeFFT(HC01_isoTd000K_t, HC01_isoTd000K_Uy,D,UInf);
freqHC01_isoTd100K, pHC01_isoTd100K, maxFreqHC01_isoTd100K, ST_HC01_isoTd100K = computeFFT(HC01_isoTd100K_t, HC01_isoTd100K_Uy,D,UInf);
freqHC01_isoTd200K, pHC01_isoTd200K, maxFreqHC01_isoTd200K, ST_HC01_isoTd200K = computeFFT(HC01_isoTd200K_t, HC01_isoTd200K_Uy,D,UInf);
freqHC01_isoTd300K, pHC01_isoTd300K, maxFreqHC01_isoTd300K, ST_HC01_isoTd300K = computeFFT(HC01_isoTd300K_t, HC01_isoTd300K_Uy,D,UInf);


T_isoT = 0.407193- 0.400001;
freq_isoT = 1.0/T_isoT;
ST_HC01_isoT = freq_isoT*D/UInf;




# freqHC01_CY, pHC01_CY, maxFreqHC01_CY, ST_HC01_CY = computeFFT(HC01_CY_T1,HC01_CY_CL1,D,UInf);

# freqHC01_WN1, pHC01_WN1, maxFreqHC01_WN1, ST_HC01_WN1 = computeFFT(HC01_WN1_T1,HC01_WN1_CL1,D,UInf);
# freqHC01_WN5, pHC01_WN5, maxFreqHC01_WN5, ST_HC01_WN5 = computeFFT(HC01_WN5_T1,HC01_WN5_CL1,D,UInf);
# freqHC01_WN10, pHC01_WN10, maxFreqHC01_WN10, ST_HC01_WN10 = computeFFT(HC01_WN10_T1,HC01_WN10_CL1,D,UInf);

meanCD1_H00 = mean(HC00_CD1);
meanCL1_H00 = mean(HC00_CL1);

meanCD1_H01 = mean(HC01_CD1);
meanCL1_H01 = mean(HC01_CL1);

# meanCD1_H02 = mean(HC02_CD1);
# meanCL1_H02 = mean(HC02_CL1);

# meanCD1_H01_BL = mean(HC01_BL_CD1);
# meanCL1_H01_BL = mean(HC01_BL_CL1);

meanCD1_H01_isoT = mean(HC01_isoT_CD1);
meanCL1_H01_isoT = mean(HC01_isoT_CL1);

meanCD1_H01_isoTd100K = mean(HC01_isoTd100K_CD1);
meanCL1_H01_isoTd100K = mean(HC01_isoTd100K_CL1);

meanCD1_H01_isoTd200K = mean(HC01_isoTd200K_CD1);
meanCL1_H01_isoTd200K = mean(HC01_isoTd200K_CL1);

meanCD1_H01_isoTd300K = mean(HC01_isoTd300K_CD1);
meanCL1_H01_isoTd300K = mean(HC01_isoTd300K_CL1);


# meanCD1_H01_CY = mean(HC01_CY_CD1);
# meanCL1_H01_CY = mean(HC01_CY_CL1);


# meanCD1_H01_WN1 = mean(HC01_WN1_CD1);
# meanCL1_H01_WN1 = mean(HC01_WN1_CL1);

# meanCD1_H01_WN5 = mean(HC01_WN5_CD1);
# meanCL1_H01_WN5 = mean(HC01_WN5_CL1);

# meanCD1_H01_WN10 = mean(HC01_WN10_CD1);
# meanCL1_H01_WN10 = mean(HC01_WN10_CL1);


debug2 = true;


deltaT = [
0
100
200
300
];

cLiftAll = [
	meanCL1_H01_isoT
	meanCL1_H01_isoTd100K
	meanCL1_H01_isoTd200K
	meanCL1_H01_isoTd300K

];

cDragAll = [
	meanCD1_H01_isoT
	meanCD1_H01_isoTd100K
	meanCD1_H01_isoTd200K
	meanCD1_H01_isoTd300K
];

LrAll =[
	1.145
	1.1
	1.2
	1.464
];

StAll = [
	ST_HC01_isoT
	ST_HC01_isoTd100K
	ST_HC01_isoTd200K
	ST_HC01_isoTd300K
];

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

	plot(HC01_T1, HC01_CL1, "-.m", label="HC01", linewidth = lw2);
	
	plot(HC01_isoT_T1, HC01_isoT_CL1, "--k", label="HC01_isoT", linewidth = lw2);
	plot(HC01_isoTd100K_T1, HC01_isoTd100K_CL1, "-.r", label="HC01_isoTd100K", linewidth = lw2);
	plot(HC01_isoTd200K_T1, HC01_isoTd200K_CL1, "-.g", label="HC01_isoTd200K", linewidth = lw2);
	plot(HC01_isoTd300K_T1, HC01_isoTd300K_CL1, "-.b", label="HC01_isoTd300K", linewidth = lw2);
	

	xlabel("t[s]");
	ylabel("Cl");
	#xlim(tmin, tmax);
	#ylim(yminCL, ymaxCL);
	legend(fontsize=8);
	grid();
	
	subplot(2,1,2)

	plot(HC01_T1, HC01_CD1, "-.m", label="HC01", linewidth = lw2);
	
	plot(HC01_isoT_T1, HC01_isoT_CD1, "--k", label="HC01_isoT", linewidth = lw2);
	plot(HC01_isoTd100K_T1, HC01_isoTd100K_CD1, "-.r", label="HC01_isoTd100K", linewidth = lw2);
	plot(HC01_isoTd200K_T1, HC01_isoTd200K_CD1, "-.g", label="HC01_isoTd200K", linewidth = lw2);
	plot(HC01_isoTd300K_T1, HC01_isoTd300K_CD1, "-.b", label="HC01_isoTd300K", linewidth = lw2);	
	

	xlabel("t[s]");
	ylabel("Cd");
	#xlim(tmin, tmax);
	#ylim(yminCD, ymaxCD);
	legend(fontsize=8);
	grid();


	# figure(11);
	# clf();
	# loglog(freqHC01_WN1, pHC01_WN1,"-.m",label="HC01_WN1");
	# loglog(freqHC01_WN5, pHC01_WN5,"-.y",label="HC01_WN5");
	# loglog(freqHC01_WN10, pHC01_WN10,"-.k",label="HC01_WN10");
	# xlabel("freq [1/s]");
	# ylabel("PWD");
	# #xlim(tmin, tmax);
	# #ylim(yminCD, ymaxCD);
	# legend(fontsize=8);
	# grid();

	figure(11);
	clf();
	subplot(2,2,1)
	loglog(freqHC01_isoT, pHC01_isoT,"-.m",label="HC01_isoT");
	loglog(freqHC01_isoTd100K, pHC01_isoTd100K,"-.r",label="HC01_isoTd100K");
	loglog(freqHC01_isoTd200K, pHC01_isoTd200K,"-.g",label="HC01_isoTd200K");
	loglog(freqHC01_isoTd300K, pHC01_isoTd300K,"-.b",label="HC01_isoTd300K");
	xlabel("freq [1/s]");
	ylabel("PWD");
	#xlim(tmin, tmax);
	#ylim(yminCD, ymaxCD);
	legend(fontsize=8);
	grid();

	subplot(2,2,2)
	plot(HC01_isoTd000K_t, HC01_isoTd000K_Uy,"-.m",label="HC01_isoT");
	plot(HC01_isoTd100K_t, HC01_isoTd100K_Uy,"-.r",label="HC01_isoT100K");
	plot(HC01_isoTd200K_t, HC01_isoTd200K_Uy,"-.g",label="HC01_isoT200K");
	plot(HC01_isoTd300K_t, HC01_isoTd300K_Uy,"-.b",label="HC01_isoT300K");
	xlabel("t [s]");
	ylabel("Uy [m/s]");
	#xlim(tmin, tmax);
	#ylim(yminCD, ymaxCD);
	legend(fontsize=8);
	grid();


	# figure(12);
	# clf();
	subplot(2,2,3)
	plot(HC01_isoT_X, HC01_isoTd000K_Ux,"-.m",label="HC01_isoT");
	plot(HC01_isoT_X, HC01_isoTd100K_Ux,"-.r",label="HC01_isoTd100K");
	plot(HC01_isoT_X, HC01_isoTd200K_Ux,"-.g",label="HC01_isoTd200K");
	plot(HC01_isoT_X, HC01_isoTd300K_Ux,"-.b",label="HC01_isoTd300K");
	xlabel("x/D");
	ylabel(L"U/U_\infty");
	#xlim(tmin, tmax);
	#ylim(yminCD, ymaxCD);
	legend(fontsize=8);
	grid();

	subplot(2,2,4)
	plot(HC01_isoT_X, HC01_isoTd000K_UxRMS,"-.m",label="HC01_isoT");
	plot(HC01_isoT_X, HC01_isoTd100K_UxRMS,"-.r",label="HC01_isoTd100K");
	plot(HC01_isoT_X, HC01_isoTd200K_UxRMS,"-.g",label="HC01_isoTd200K");
	plot(HC01_isoT_X, HC01_isoTd300K_UxRMS,"-.b",label="HC01_isoTd300K");
	xlabel("x/D");
	ylabel(L"U'/U_\infty");
	#xlim(tmin, tmax);
	#ylim(yminCD, ymaxCD);
	legend(fontsize=8);
	grid();


	figure(13);
	clf();
	subplot(2,2,1)
	plot(deltaT, cLiftAll, "ok", label="Cl", markersize = ms);
	plot(deltaT, cLiftAll, "-r", label="Cl", markersize = ms, linewidth=lw1);
	xlabel(L"\Delta T");
	ylabel(L"C_l");
	#legend(fontsize=8);
	grid();
	subplot(2,2,2)
	plot(deltaT, cDragAll, "ok", label="Cd", markersize = ms);
	plot(deltaT, cDragAll, "-r", label="Cd", markersize = ms, linewidth=lw1);
	xlabel(L"\Delta T");
	ylabel(L"C_d");
	#legend(fontsize=8);
	grid();
	subplot(2,2,3)
	plot(deltaT, LrAll, "ok", label="Lr", markersize = ms);
	plot(deltaT, LrAll, "-r", label="Lr", markersize = ms, linewidth = lw1);
	xlabel(L"\Delta T");
	ylabel(L"L_r");
	#legend(fontsize=8);
	grid();
	subplot(2,2,4)
	plot(deltaT, StAll, "ok", label="St", markersize = ms);
	plot(deltaT, StAll, "-r", label="St", markersize = ms, linewidth = lw1);
	xlabel(L"\Delta T");
	ylabel(L"St");
	#legend(fontsize=8);
	grid();
	


end


println("Mean aerodynamic characteristics");
println("EXP (Re=50000)", "Cd:\t", 0.500);
println("EXP           ", "Cd:\t", 0.370);
println("EXP           ", "Cd:\t", 0.467);
println("EXP           ", "Cd:\t", 0.476);
#println("HC00          ", "Cd:\t", round(meanCD1_H00,digits=3));
println("HC01          ", "Cd:\t", round(meanCD1_H01,digits=3));
#println("HC02          ", "Cd:\t", round(meanCD1_H02,digits=3));
#println("HC01_BL       ", "Cd:\t", round(meanCD1_H01_BL,digits=3));
println("HC01_isoT     ", "Cd:\t", round(meanCD1_H01_isoT,digits=3));
println("HC01_isoTd100K", "Cd:\t", round(meanCD1_H01_isoTd100K,digits=3));
println("HC01_isoTd200K", "Cd:\t", round(meanCD1_H01_isoTd200K,digits=3));
println("HC01_isoTd300K", "Cd:\t", round(meanCD1_H01_isoTd300K,digits=3));
#println("HC01_CL       ", "Cd:\t", round(meanCD1_H01_CY,digits=3));

# println("HC01_WN1      ", "Cd:\t", round(meanCD1_H01_WN1,digits=3));
# println("HC01_WN5      ", "Cd:\t", round(meanCD1_H01_WN5,digits=3));
# println("HC01_WN10     ", "Cd:\t", round(meanCD1_H01_WN10,digits=3));

println("=====================================");
println("EXP           ", "Cl:\t", -1.100);
println("EXP           ", "Cl:\t", -0.540);
println("EXP           ", "Cl:\t", -0.634);
println("EXP           ", "Cl:\t", -0.816);
#println("HC00          ", "Cl:\t", round(meanCL1_H00,digits=3));
println("HC01          ", "Cl:\t", round(meanCL1_H01,digits=3));
# println("HC02          ", "Cl:\t", round(meanCL1_H02,digits=3));
# println("HC01_BL       ", "Cl:\t", round(meanCL1_H01_BL,digits=3));
# println("HC01_CY       ", "Cl:\t", round(meanCL1_H01_CY,digits=3));
println("HC01_isoT     ", "Cl:\t", round(meanCL1_H01_isoT,digits=3));
println("HC01_isoTd100K", "Cl:\t", round(meanCL1_H01_isoTd100K,digits=3));
println("HC01_isoTd200K", "Cl:\t", round(meanCL1_H01_isoTd200K,digits=3));
println("HC01_isoTd300K", "Cl:\t", round(meanCL1_H01_isoTd300K,digits=3));
# println("HC01_WN1      ", "Cl:\t", round(meanCL1_H01_WN1,digits=3));
# println("HC01_WN5      ", "Cl:\t", round(meanCL1_H01_WN5,digits=3));
# println("HC01_WN10     ", "Cl:\t", round(meanCL1_H01_WN10,digits=3));

println("=====================================");
println("EXP           ", "St :\t", 0.37);
println("EXP           ", "St :\t", 0.41);
println("EXP           ", "St :\t", 0.44);
#println("HC00          ", "St :\t", round(ST_HC00,digits=2));
println("HC01          ", "St :\t", round(ST_HC01,digits=2));
#println("HC02          ", "St :\t", round(ST_HC02,digits=2));
#println("HC01_BL       ", "St :\t", round(ST_HC01_BL,digits=2));
println("HC01_isoT     ", "St :\t", round(ST_HC01_isoT,digits=2));
println("HC01_isoTd100K", "St :\t", round(ST_HC01_isoTd100K,digits=2));
println("HC01_isoTd200K", "St :\t", round(ST_HC01_isoTd200K,digits=2));
println("HC01_isoTd300K", "St :\t", round(ST_HC01_isoTd300K,digits=2));
# println("HC01_CY       ", "St :\t", round(ST_HC01_CY,digits=2));
# println("HC01_WN1      ", "St :\t", round(ST_HC01_WN1,digits=2));
# println("HC01_WN5      ", "St :\t", round(ST_HC01_WN1,digits=2)); ## freq WN1 == WN5 
# println("HC01_WN10     ", "St :\t", round(ST_HC01_WN1,digits=2)); ## freq WN1 == WN10 
