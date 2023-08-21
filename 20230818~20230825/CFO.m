function [hat_delta_f,t] = CFO(frame_data,y,Fs,CFO_ignore)
	index  			= 1;
	CFO_sum	 		= 0;
	CP_num			= 0;
	t 				= 0:(10e-3/1228800):10e-3;
	t 				= t(2:1228801);
    thread          = parfeval(backgroundPool,@reshape,1,t(frame_data.CPdataPos),2048,560);
	for symbol = 1:560
		if (mod(symbol,28) - 1)%size 144
			CFO_head = y( : , index+CFO_ignore     : index+143      );
			CFO_tail = y( : , index+CFO_ignore+2048: index+143+2048 );
			CFO_sum	 = CFO_sum+ sum(CFO_head .* (CFO_tail').','all');
			CP_num	 = CP_num + 144 - CFO_ignore;
			index    = index + 144 + 2048;
		else				 %size 208
			CFO_head = y( : , index+CFO_ignore     : index+207      );
			CFO_tail = y( : , index+CFO_ignore+2048: index+207+2048 );
			CFO_sum	 = CFO_sum+ sum(CFO_head .* (CFO_tail').','all');
			CP_num	 = CP_num + 208 - CFO_ignore;
			index    = index  + 208 + 2048;
		end
    end
    t = fetchOutputs(thread);
    hat_delta_f = -(Fs/2048) * angle( CFO_sum/CP_num )/(2*pi) ;
    
end

