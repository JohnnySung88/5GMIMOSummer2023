function [y_out] = syc(frame_data,y)
    thread1 =parfeval(backgroundPool,@xcorr,2,y(1,1:1228800/5),frame_data.PSS_t);
	thread2 =parfeval(backgroundPool,@xcorr,2,y(2,1:1228800/5),frame_data.PSS_t);
    [corr_1,corr_pos_1] = fetchOutputs(thread1);
    [corr_2,corr_pos_2] = fetchOutputs(thread2);
	[~,max_pos_1]=max(corr_1);
	[~,max_pos_2]=max(corr_2);
	y_pos(1) = corr_pos_1(max_pos_1) - frame_data.syc_pos+1;
	y_pos(2) = corr_pos_2(max_pos_2) - frame_data.syc_pos+1;
    if( (y_pos(1)+1228800-1 <= length(y)) & (y_pos(2)+1228800-1 <= length(y)) )
        y_out    = zeros(2,1228800);
        y_out(1,:) 	  = y(1,y_pos(1):y_pos(1)+1228800-1);
	    y_out(2,:)	  = y(2,y_pos(2):y_pos(2)+1228800-1);
    else
        y_out    = zeros(2,1228800);
    end
    
end