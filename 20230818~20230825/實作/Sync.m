function [y_out] = Sync(frame_data,y)
    [corr_1,corr_pos_1] = xcorr(y(1,1:floor(length(y(1,:))/2)),frame_data.PSS_t);
    [corr_2,corr_pos_2] = xcorr(y(2,1:floor(length(y(2,:))/2)),frame_data.PSS_t);
    
    [max1,max_pos_1]=max(abs(corr_1));
    [max2,max_pos_2]=max(abs(corr_2));
    y_pos(1) = corr_pos_1(max_pos_1) - frame_data.syc_pos+1 - 3;  % 從找到的位置往前抓三個位元
    y_pos(2) = corr_pos_2(max_pos_2) - frame_data.syc_pos+1 - 3;  % 從找到的位置往前抓三個位元

    if( (y_pos(1)+1228800-1 <= length(y)) & (y_pos(2)+1228800-1 <= length(y)) & (y_pos(1) > 0) & (y_pos(2) > 0) & max1>0.05 & max2>0.05)
        y_out = zeros(2,1228800);
        y_out(1,:) = y(1,y_pos(1):y_pos(1)+1228800-1);
        y_out(2,:) = y(2,y_pos(2):y_pos(2)+1228800-1);
    else
        y_out = zeros(2,1228800);
    end
end
