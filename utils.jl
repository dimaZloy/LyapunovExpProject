
#rev 4 July 2020
#rev 26 July 2020

function trapz(x::T1, y::T2) where {fT,T1<:AbstractVector{fT},T2<:AbstractVector{fT}}
    n = length(x)
    r = zero(fT)
    
    if n == 1; return r; end
    for i in 2:n-1
       @inbounds r += (x[i+1]-x[i-1]) * y[i]
    end
    @inbounds r += (x[end]-x[end-1]) * y[end] + (x[2] - x[1])*y[1]
    trapz_int = r/2
    return trapz_int
end

function calcBulkVelocity(x::T1, y::T2) where {fT,T1<:AbstractVector{fT},T2<:AbstractVector{fT}}

    area = trapz(x,y);
    mass = trapz(x,y.*y);
    return mass/area;  		

end


function addItem(x,dx)
	return [x; dx];
end


function calcSutherlandViscosityOF(T::Float64)

	As = 1.4792e-6;
	Ts = 116;
	return As*sqrt(T)/(1.0 + Ts/T);

end

function readForcesAF(fname, Nsize::Int64)
	
	zzzT = zeros(Float64,Nsize);
	zzzF = zeros(Float64,Nsize);
	

	io = open(fname,"r");

	for i=1:2
		head = readline(io);
	end

	counter = 0;

	while(true)

		line = readline(io);
		counter = counter+1;
	
		if line == "" #eof(io)
			break;
		else
			#println(line);
			z = split(line);

			timeStep = parse(Float64,z[1]);
			Coeff = parse(Float64,z[2]);
			
			zzzT[counter,1] = timeStep;
			zzzF[counter,1] = Coeff;
			
		end
		
		end # end while
	

	close(io);

	return zzzT, zzzF
	
end

function readForcesDataOF(fname, Nsize::Int64)
	
	zzzT = zeros(Float64,Nsize);
	zzzCD = zeros(Float64,Nsize);
	zzzCL = zeros(Float64,Nsize);


	io = open(fname,"r");

	for i=1:10
		head = readline(io);
	end

	counter = 0;

	while(true)

		line = readline(io);
		counter = counter+1;
	
		if line == "" #eof(io)
			break;
		elseif 	(counter == Nsize)
			break;			
		else
			#println(line);
			z = split(line);

			timeStep = parse(Float64,z[1]);
			dragCoeff = parse(Float64,z[3]);
			liftCoeff = parse(Float64,z[4]);
			#global counter += 1;
		
			#zzzT = [zzzT; timeStep];
			#zzzCD = [zzzCD; dragCoeff];
			#zzzCL = [zzzCL; liftCoeff];

			zzzT[counter,1] = timeStep;
			zzzCD[counter,1] = dragCoeff;
			zzzCL[counter,1] =  liftCoeff;
			
		end
		
		end # end while
	

	close(io);

	return zzzT, zzzCD, zzzCL
	
end




function readDataFile(fname)
	
	zzzT = zeros(Float64,0);
	zzzCD = zeros(Float64,0);
	zzzCL = zeros(Float64,0);


	io = open(fname,"r");

	for i=1:10
		head = readline(io);
	end

	counter = 0;

	while(true)

		line = readline(io);
	
				
		if line == "" #eof(io)
			break;
		else
			#println(line);
			z = split(line);

			timeStep = parse(Float64,z[1]);
			dragCoeff = parse(Float64,z[3]);
			liftCoeff = parse(Float64,z[4]);
			#global counter += 1;
		
			zzzT = [zzzT; timeStep];
			zzzCD = [zzzCD; dragCoeff];
			zzzCL = [zzzCL; liftCoeff];
		end
		
		end # end while
	

	close(io);

	return zzzT, zzzCD, zzzCL
	
end

function readDataFileUprobe(fname)
	
	
	UDATA = zeros(Float64,190910,Int(54*3 + 1)) 
	

	io = open(fname,"r");
	

	for i=1:56
		head = readline(io);
	end

	counter = 0;
	
	println("starting parse file: ", fname);

	while(true)

		line = readline(io);
		#display(line)
		#pause(10)
		
		#global UDATA;
		#global counter;
		
		counter = counter+1;
	
		if line == "" #eof(io)
			break;
		else
			#println(line);
			z = split(line);
			N = length(z);
			#display(N)
			#display(z)
			#pause(10)
			
			for i=1:N
				tmp =z[i];
				n = length(z[i]);
				if ( occursin("(",z[i]))
					tmp = SubString(z[i],2,n);
				elseif ( occursin(")",z[i]))
					tmp = SubString(z[i],1,n-1);
				elseif (n != 0)
					tmp = z[i];
				end
				UDATA[counter,i] = parse(Float64,tmp);
			end
			
			if ( mod(counter,10000) == 0 )
				display(counter);
				pause(1e-4);
			end
		
		end
		
		end # end while
	

	close(io);
	
	println("done");
	println("saving U data to X_Udata.h5");
	

	h5open("X_Udata.h5","w") do file
		write(file, "UDATA", UDATA);
	end
	
	
end



function readDataFileUprobe(fnameRead, fnameWrite, offset::Int, L::Int, N::Int)
	
	
	UDATA = zeros(Float64,L,N) 
	

	io = open(fnameRead,"r");
	

	for i=1:offset
		head = readline(io);
	end

	counter = 0;
	
	println("starting parse file: ", fnameRead);

	while(true)

		line = readline(io);
		#display(line)
		#pause(10)
		
		#global UDATA;
		#global counter;
		
		counter = counter+1;
	
		if line == "" #eof(io)
			break;
		elseif 	(counter == L)
			break;		
		else
			#println(line);
			z = split(line);
			N = length(z);
			#display(N)
			#display(z)
			#pause(10)
			
			for i=1:N
				tmp =z[i];
				n = length(z[i]);
				if ( occursin("(",z[i]))
					tmp = SubString(z[i],2,n);
				elseif ( occursin(")",z[i]))
					tmp = SubString(z[i],1,n-1);
				elseif (n != 0)
					tmp = z[i];
				end
				UDATA[counter,i] = parse(Float64,tmp);
			end
			
			if ( mod(counter,10000) == 0 )
				display(counter);
				pause(1e-4);
			end
		
		end
		
		end # end while
	

	close(io);
	
	println("done");
	println("Total lines read: ", counter);
	println("saving U data to h5 format: ", fnameWrite );
	

	h5open(fnameWrite,"w") do file
		write(file, "UDATA", UDATA);
	end
	
	
end

function readForcesDataOF2(fname, Nsize::Int64)
	


	io = open(fname,"r");

	for i=1:9
		head = readline(io);
		println(head);
	end

	counter = 0;
	Nsize = Nsize-9-1;

	zzzT = zeros(Float64,Nsize);
	zzzCD = zeros(Float64,Nsize);
	zzzCL = zeros(Float64,Nsize);


	while(true)

		line = readline(io);
		counter = counter+1;
		
	
		if line == "" #eof(io)
			println(counter);			
			println(Nsize);			

			break;
		#elseif 	(counter == Nsize)
		#	println(counter);			
		#	println(Nsize);			
		#	break;			
		else
			#println(line);
			z = split(line);

			timeStep = parse(Float64,z[1]);
			dragCoeff = parse(Float64,z[3]);
			liftCoeff = parse(Float64,z[4]);
			#global counter += 1;
		
			#zzzT = [zzzT; timeStep];
			#zzzCD = [zzzCD; dragCoeff];
			#zzzCL = [zzzCL; liftCoeff];

			zzzT[counter,1] = timeStep;
			zzzCD[counter,1] = dragCoeff;
			zzzCL[counter,1] =  liftCoeff;
			
		end

		
	end # end while
		
	

	close(io);

	return zzzT, zzzCD, zzzCL
	
end

## blyd, for OF2020 or later 
function readForcesDataOF3(fname, Nsize::Int64)
	


	io = open(fname,"r");

	for i=1:13
		head = readline(io);
		println(head);
	end

	counter = 0;
	Nsize = Nsize-13-1;

	zzzT = zeros(Float64,Nsize);
	zzzCD = zeros(Float64,Nsize);
	zzzCL = zeros(Float64,Nsize);


	while(true)

		line = readline(io);
		counter = counter+1;
		
	
		if line == "" #eof(io)
			println(counter);			
			println(Nsize);			

			break;
		#elseif 	(counter == Nsize)
		#	println(counter);			
		#	println(Nsize);			
		#	break;			
		else
			#println(line);
			z = split(line);

			timeStep = parse(Float64,z[1]);
			dragCoeff = parse(Float64,z[2]);
			liftCoeff = parse(Float64,z[4]);
			#global counter += 1;
		
			#zzzT = [zzzT; timeStep];
			#zzzCD = [zzzCD; dragCoeff];
			#zzzCL = [zzzCL; liftCoeff];

			zzzT[counter,1] = timeStep;
			zzzCD[counter,1] = dragCoeff;
			zzzCL[counter,1] =  liftCoeff;
			
		end

		
	end # end while
		
	

	close(io);

	return zzzT, zzzCD, zzzCL
	
end

