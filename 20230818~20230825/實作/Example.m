%Sync
%{
function [y_out] = Sync(frame_data,y)
    thread1 =parfeval(backgroundPool,@xcorr,2,y(1,1:1228800/5),frame_data.PSS_t);
	thread2 =parfeval(backgroundPool,@xcorr,2,y(2,1:1228800/5),frame_data.PSS_t);
    [corr_1,corr_pos_1] = fetchOutputs(thread1);
    [corr_2,corr_pos_2] = fetchOutputs(thread2);
	[~,max_pos_1]=max(abs(corr_1));
	[~,max_pos_2]=max(abs(corr_2));
	y_pos(1) = corr_pos_1(max_pos_1) - frame_data.syc_pos+1 - 3;
	y_pos(2) = corr_pos_2(max_pos_2) - frame_data.syc_pos+1 - 3;
    if( (y_pos(1)+1228800-1 <= length(y)) & (y_pos(2)+1228800-1 <= length(y)) )
        y_out    = zeros(2,1228800);
        y_out(1,:) 	  = y(1,y_pos(1):y_pos(1)+1228800-1);
	    y_out(2,:)	  = y(2,y_pos(2):y_pos(2)+1228800-1);
    else
        y_out    = zeros(2,1228800);
    end    
end
%}


%Rx_fun
%{
function [Time,Biterror,Capacity_sum,Rx1_SNR,Rx2_SNR,JPG_RGB] = Rx_fun(frame_data,Rx_signal,DTinfo,CFO_ignore)
	Fs 	= frame_data.Fs;
	Tx	= frame_data.Tx;
	Rx  = frame_data.Rx;
	[X1,Y1]  = meshgrid(3:14:549,1:1644);
	[Xq,Yq]	 = meshgrid(3: 1:549,1:1644);
	y_rmCP		= zeros(2048 , 560, Tx);
	H_INTER		= zeros(1644,560,2,2);	
	DMRS_DATA 	= +0.7071 + 0.7071*1i ;
	H_LMMSE		= zeros(1644,40,2,2);
	y_in 		= Rx_signal;%2*1228800
	tic;
	%同步
	y = Sync(frame_data,y_in);

	
	%CFO估測(需多重路徑猜測) 背景運算
	f(1) = parfeval(backgroundPool,@CFO,2,frame_data,y,Fs,CFO_ignore);
	%移除CP 0.072sec
	y_rmCP(:,:,1) = reshape(y(1,frame_data.CPdataPos),2048,560);
	y_rmCP(:,:,2) = reshape(y(2,frame_data.CPdataPos),2048,560);
	%CFO補償
	[hat_delta_f,t] = fetchOutputs(f(1));
	y_rmCP = y_rmCP ./ exp( 1i * 2 * pi * hat_delta_f * t);
	%FFT----以下為頻域 0.02sec
	Y_fft 	= fftshift( fft( y_rmCP/sqrt(2048) ) ,1);
	%rm Guard Band 0.03sec
	Y	  		= [ Y_fft( 203:1024,:,:) ; Y_fft( 1026:1847,:,:) ];
	%LMMSE估測 0.025sec
	Y_DMRS 	= Y(2:2:1644 , 3:14:560,:);
	H_LMMSE	= pagemtimes(frame_data.LMMSE , DMRS_DATA' .*Y_DMRS);
	%雜訊估測(背景執行) 0.10sec(雜訊估測 + 線性內插)
	f(5) = parfeval(backgroundPool,@SNR_EST,1,frame_data,H_LMMSE,Y_DMRS);
	%線性內差(LMMSE only) 
	%連續
	f(1) = parfeval(backgroundPool,@interp2,1,X1,Y1,H_LMMSE(:,:,1,1),Xq,Yq);
	f(2) = parfeval(backgroundPool,@interp2,1,X1,Y1,H_LMMSE(:,:,1,2),Xq,Yq);
	f(3) = parfeval(backgroundPool,@interp2,1,X1,Y1,H_LMMSE(:,:,2,1),Xq,Yq);
	f(4) = parfeval(backgroundPool,@interp2,1,X1,Y1,H_LMMSE(:,:,2,2),Xq,Yq);
	%邊界
	for symbol  = 549:560 
		head_dist	= symbol - 549;
		back_dist	= 563    - symbol;
		H_INTER(:,symbol,:,:) = ( back_dist * H_LMMSE(:,40,:,:) + head_dist * H_LMMSE(:,1,:,:)  )  /14;
	end
	%邊界
	for symbol  = 1:2     
		head_dist	= 11 + symbol;
		back_dist	= 3  - symbol;
		H_INTER(:,symbol,:,:) = ( back_dist * H_LMMSE(:,40,:,:) + head_dist * H_LMMSE(:,1,:,:)  )  /14;
	end
	H_INTER(:,3:549,1,1) = fetchOutputs(f(1));
	H_INTER(:,3:549,1,2) = fetchOutputs(f(2));
	H_INTER(:,3:549,2,1) = fetchOutputs(f(3));
	H_INTER(:,3:549,2,2) = fetchOutputs(f(4));
	%雜訊估測 (回傳)
	Rx_No = fetchOutputs(f(5));	
	%detector
	norm_Y = Y 		 ./ Rx_No;
	norm_H = H_INTER ./ Rx_No;
	switch DTinfo
		case 'LMMSE'
			X_hat = LMMSE(norm_Y,norm_H,Tx);
		case 'ZF'
			X_hat = ZFD(norm_Y,norm_H,Tx);
	end
	X_hat = X_hat/frame_data.NF;
	%反解資料
	LDPC_mod_L_hat	= X_hat(frame_data.DATA_Pos(1:1770336));
	%解碼 0.17sec
	LDPC_dec_L_hat = qamdemod(LDPC_mod_L_hat,frame_data.QAM  ,'gray');		
	LDPC_bin_L_hat = reshape(dec2bin (LDPC_dec_L_hat,frame_data.q_bit).' - '0',[],1) ;
	LDPC_bin_part  = permute(reshape(LDPC_bin_L_hat,[],1296) ,[2,1]);
	JPG_bin_hat   = reshape(LDPC_bin_part(1:648,:).',[],8);
	%Time
	Time = toc;
	%noLDPC decode(image)
	bin_table	= 2 .^ (7:-1:0);
	JPG_row = frame_data.JPG_row;
	JPG_col = frame_data.JPG_col;
	JPG_size= frame_data.JPG_size;
	JPG_bin_hat = JPG_bin_hat(1:JPG_size,:);
	JPG_dec_hat = sum(JPG_bin_hat.*bin_table,2);%bin2dec
	JPG_Csize	 = JPG_row*JPG_col;
	JPG_RGB	 = zeros(JPG_row,JPG_col,3);
	JPG_RGB(:,:,1) = reshape( JPG_dec_hat(             1:JPG_Csize  ) ,JPG_row,JPG_col);
	JPG_RGB(:,:,2) = reshape( JPG_dec_hat(JPG_Csize  +1:JPG_Csize*2) ,JPG_row,JPG_col);
	JPG_RGB(:,:,3) = reshape( JPG_dec_hat(JPG_Csize*2+1:JPG_Csize*3) ,JPG_row,JPG_col);
	JPG_RGB		= uint8(JPG_RGB);
	%BER
	Biterror = sum(JPG_bin_hat ~= frame_data.JPG_bin,'all');
	%SNR
	Rx1_SNR  = -10*log10(Rx_No(1) / mean(abs(Y(:,:,1).^2),'all'));
	Rx2_SNR  = -10*log10(Rx_No(2) / mean(abs(Y(:,:,2).^2),'all'));
	%Capacity
	No = (Rx_No(1)+Rx_No(2))/2;
	SNR = 1/No;
	Capacity_sum = 0;
	parfor SC = 1:1644
		for slot = 1:560
			unit_H 	 = reshape(H_INTER(SC,slot,:,:),Rx,Tx);
			Capacity_sum = Capacity_sum + abs( log2( det( eye(Tx) + (1/Tx) .* SNR .* (unit_H*unit_H'))));
		end
	end
end
%}

%Rx_fun_LDPC
%{
function [Time,Biterror,Capacity_sum,Rx1_SNR,Rx2_SNR,JPG_RGB] = Rx_fun_LDPC(frame_data,Rx_signal,DTinfo,CFO_ignore,iteration)
	Fs 	= frame_data.Fs;
	Tx	= frame_data.Tx;
	Rx  = frame_data.Rx;
	[X1,Y1]  = meshgrid(3:14:549,1:1644);
	[Xq,Yq]	 = meshgrid(3: 1:549,1:1644);
	y_rmCP		= zeros(2048 , 560, Tx);
	H_INTER		= zeros(1644,560,2,2);	
	DMRS_DATA 	= +0.7071 + 0.7071*1i ;
	H_LMMSE		= zeros(1644,40,2,2);
	y 	= Rx_signal;%2*1228800
	tic;
	%同步
	y = Sync(frame_data,y);
	%CFO估測(需多重路徑猜測) 背景運算
	f(1) = parfeval(backgroundPool,@CFO,2,frame_data,y,Fs,CFO_ignore);
	%移除CP 0.072sec
	y_rmCP(:,:,1) = reshape(y(1,frame_data.CPdataPos),2048,560);
	y_rmCP(:,:,2) = reshape(y(2,frame_data.CPdataPos),2048,560);
	%CFO補償
	[hat_delta_f,t] = fetchOutputs(f(1));
	y_rmCP = y_rmCP ./ exp( 1i * 2 * pi * hat_delta_f * t);
	%FFT----以下為頻域 0.02sec
	Y_fft 	= fftshift( fft( y_rmCP/sqrt(2048) ) ,1);
	%rm Guard Band 0.03sec
	Y	  		= [ Y_fft( 203:1024,:,:) ; Y_fft( 1026:1847,:,:) ];
	%LMMSE估測 0.025sec
	Y_DMRS 	= Y(2:2:1644 , 3:14:560,:);
	H_LMMSE	= pagemtimes(frame_data.LMMSE , DMRS_DATA' .*Y_DMRS);
	%雜訊估測(背景執行) 0.10sec(雜訊估測 + 線性內插)
	f(5) = parfeval(backgroundPool,@SNR_EST,1,frame_data,H_LMMSE,Y_DMRS);
	%線性內差(LMMSE only) 
	%連續
	f(1) = parfeval(backgroundPool,@interp2,1,X1,Y1,H_LMMSE(:,:,1,1),Xq,Yq);
	f(2) = parfeval(backgroundPool,@interp2,1,X1,Y1,H_LMMSE(:,:,1,2),Xq,Yq);
	f(3) = parfeval(backgroundPool,@interp2,1,X1,Y1,H_LMMSE(:,:,2,1),Xq,Yq);
	f(4) = parfeval(backgroundPool,@interp2,1,X1,Y1,H_LMMSE(:,:,2,2),Xq,Yq);
	%邊界
	for symbol  = 549:560 
		head_dist	= symbol - 549;
		back_dist	= 563    - symbol;
		H_INTER(:,symbol,:,:) = ( back_dist * H_LMMSE(:,40,:,:) + head_dist * H_LMMSE(:,1,:,:)  )  /14;
	end
	%邊界
	for symbol  = 1:2     
		head_dist	= 11 + symbol;
		back_dist	= 3  - symbol;
		H_INTER(:,symbol,:,:) = ( back_dist * H_LMMSE(:,40,:,:) + head_dist * H_LMMSE(:,1,:,:)  )  /14;
	end
	H_INTER(:,3:549,1,1) = fetchOutputs(f(1));
	H_INTER(:,3:549,1,2) = fetchOutputs(f(2));
	H_INTER(:,3:549,2,1) = fetchOutputs(f(3));
	H_INTER(:,3:549,2,2) = fetchOutputs(f(4));
	%雜訊估測 (回傳)
	Rx_No = fetchOutputs(f(5));	
	%detector
	norm_Y = Y 		 ./ Rx_No;
	norm_H = H_INTER ./ Rx_No;
	switch DTinfo
		case 'LMMSE'
			X_hat = LMMSE(norm_Y,norm_H,Tx);
		case 'ZF'
			X_hat = ZFD(norm_Y,norm_H,Tx);
	end
	X_hat = X_hat/frame_data.NF;
	%反解資料
	LDPC_mod_L_hat	= X_hat(frame_data.DATA_Pos(1:1770336));
	%LLR
	yI = real(LDPC_mod_L_hat).';
	yQ = imag(LDPC_mod_L_hat).';
	No = (Rx_No(1) + Rx_No(2))/2;
	LLR= zeros( 4*length(LDPC_mod_L_hat) ,1);
	f(1) = parfeval(backgroundPool,@LLRcul,1,yI,No, 1, 3,-1,-3);
	f(2) = parfeval(backgroundPool,@LLRcul,1,yI,No,-1, 1,-3, 3);
	f(3) = parfeval(backgroundPool,@LLRcul,1,yQ,No,-1,-3, 1, 3);
	LLR(4:4:length(LLR)) = (1/(2*No)) .* ( min( (yQ-[-1; 1]).^2 ) - min( (yQ-[-3; 3]).^2 ));;
	LLR(1:4:length(LLR)) = fetchOutputs(f(1));
	LLR(2:4:length(LLR)) = fetchOutputs(f(2));
	LLR(3:4:length(LLR)) = fetchOutputs(f(3));
	LLR_part			 = reshape(LLR,[],1296);
	LLR_part			 = permute(LLR_part,[2,1]);
	
	LLR_MSA			  	 = LDPC_MSA(LLR_part,iteration,frame_data.H_row_master,frame_data.H_row_master_size,5464,0.5);
	JPG_bin_hat 		 = reshape(LLR_MSA(1:648,:).'<0,[],8);
	%Time
	Time = toc;
	%LDPC decode(image)
	bin_table	= 2 .^ (7:-1:0);
	JPG_row = frame_data.JPG_row;
	JPG_col = frame_data.JPG_col;
	JPG_size= frame_data.JPG_size;
	JPG_bin_hat = JPG_bin_hat(1:JPG_size,:);
	JPG_dec_hat = sum(JPG_bin_hat.*bin_table,2);%bin2dec
	JPG_Csize	 = JPG_row*JPG_col;
	JPG_RGB	 = zeros(JPG_row,JPG_col,3);
	JPG_RGB(:,:,1) = reshape( JPG_dec_hat(             1:JPG_Csize  ) ,JPG_row,JPG_col);
	JPG_RGB(:,:,2) = reshape( JPG_dec_hat(JPG_Csize  +1:JPG_Csize*2) ,JPG_row,JPG_col);
	JPG_RGB(:,:,3) = reshape( JPG_dec_hat(JPG_Csize*2+1:JPG_Csize*3) ,JPG_row,JPG_col);
	JPG_RGB		= uint8(JPG_RGB);
	%BER
	Biterror = sum(JPG_bin_hat ~= frame_data.JPG_bin,'all');
	%SNR
	Rx1_SNR  = -10*log10(Rx_No(1) / mean(abs(Y(:,:,1).^2),'all'));
	Rx2_SNR  = -10*log10(Rx_No(2) / mean(abs(Y(:,:,2).^2),'all'));
	%Capacity
	SNR = 1/No;
	Capacity_sum = 0;
    parfor SC = 1:1644
		for slot = 1:560
			unit_H 	 = reshape(H_INTER(SC,slot,:,:),Rx,Tx);
			Capacity_sum = Capacity_sum + abs( log2( det( eye(Tx) + (1/Tx) .* SNR .* (unit_H*unit_H'))));
		end
	end
end
%}

%CFO
%{
function [hat_delta_f,t] = CFO(frame_data,y,Fs,CFO_ignore)
	index  			= 1;
	CFO_sum	 		= 0;
	CP_num			= 0;
	t 				= 0:(20e-3/1228800):20e-3;
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
%}