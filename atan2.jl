
function atan2(y::Float64, x::Float64)::Float64
	
	atan2::Float64 = 0.0;	
	if (x == 0.0 && y == 0.0)
	elseif (x == 0.0 && y < 0.0)
		atan2 = -pi/2.0;
	elseif (x == 0.0 && y > 0.0)
		atan2 =  pi/2.0;
	elseif (x < 0.0 && y < 0.0)
		atan2 =  atan(y/x)-pi;
	elseif (x < 0.0 && y >= 0.0)
		atan2 =  atan(y/x) + pi;
	elseif (x > 0.0)
		atan2 =  atan(y/x) ;
	end
	
	if ((atan2*180.0/pi) < 1e-5) 
		return 0.0;
	else
		return  (atan2*180.0/pi) ;
	end
	
	
end