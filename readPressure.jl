
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
include("atan2.jl")
include("hc_exp_cp.jl")


function computeHCnormalizedCoords(radialCoords, D::Float64)

	cp = pi*D/2.0/D;
	cs = 1.0;
	ct = cp + cs;
	cpn = cp/ct;

	cn = cpn/0.5;

	L = size(radialCoords,1);
	
	normCoords = zeros(Float64,L,1);

	for i =1:L
		if (radialCoords[i] <= 0.5)
			normCoords[i,1] = radialCoords[i]*cn;
		else
			normCoords[i,1] = radialCoords[i]*(1.0/cn) + 0.2;
		end
	end

	return normCoords;
end

function computeRadialCoords(x,y, D::Float64)

	L = length(x);

	s = zeros(Float64,L);
	
	epsilon_::Float64 = 1e-4; 

	for i =1:L
		if ( y[i] <= epsilon_) ## (HC7_x[i] <= eps && HC7_y[i] <=eps)
			s[i] =  0.75 - x[i]/D/2.0;	
		else
			theta = atan2(y[i],-x[i]);
			s[i] =  theta/360.0; 
		end
	
	end
	
	return s; 
	
end


###################################################################################

debug1 = true;
if (debug1)

	HC01_psurface = CSV.read("hcyl2d_01_psurf.csv", DataFrame);
	
end


HC01_x = HC01_psurface[:,1];
HC01_y = HC01_psurface[:,2];

L = length(HC01_x);
HC01_CpSpan = zeros(Float64,L);
HC01_s = zeros(Float64,L);

for i=1:L
	HC01_CpSpan[i] = 2.0*(HC01_psurface[i,5] - PInf)/rhoInf/UInf/UInf;
end

HC01_sSpan = zeros(Float64,L);
HC01_ssSpan = zeros(Float64,L);

HC01_sSpan = computeRadialCoords(HC01_x,HC01_y, D::Float64);
HC01_ssSpan = computeHCnormalizedCoords(HC01_sSpan,D);


	figure(11);
	clf();

	ms = 6;
	lw = 1.5;
	fs = 10;
	
	TL = pi*D/2.0 + D;
	Stedge = pi*D/2.0 / TL;


	plot(Isaev2018_exp_4D[:,1],Isaev2018_exp_4D[:,2],"ow",  markersize=ms, markeredgecolor="gray", label = "HC, Exp (4xD), Isaev et al. 2018 ");
	plot(Isaev2018_exp_16D[:,1],Isaev2018_exp_16D[:,2],"sw",  markersize=ms, markeredgecolor="gray", label = "HC, Exp (16xD), Isaev et al. 2018 ");
	plot(Sluchanovskaya1973[:,1]*0.99,Sluchanovskaya1973[:,2],"^c",  markersize=ms+1.0, markeredgecolor="gray", label = "HC, Exp (16xD), Sluchanovskaya 1973");

	plot(HC01_ssSpan, HC01_CpSpan, "or", label="HC01", linewidth = lw2,markersize=ms-1.0);

	xlabel("s");
	ylabel("Cp");
	xlim([0.0 , 1.0]);
	ylim([ -2.0 , 1.1]);
	xticks([ 0.0, Stedge, 1.0]);
	yticks([ -2.0, -1.0, 0.0, 1.0 ]);

	legend(fontsize=8);
	grid();

