function [Time,Biterror,Capacity_sum,Rx1_SNR,Rx2_SNR,Lena_RGB] = Rx_fun_LDPC(frame_data,Rx_signal,DTinfo,CFO_ignore,iteration)
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
	y = syc(frame_data,y);
	%CFO估測(需多重路徑猜測) 背景運算
	f(1) = parfeval(backgroundPool,@CFO_EST,2,frame_data,y,Fs,CFO_ignore);
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
			X_hat = ZF_detector(norm_Y,norm_H,Tx);
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
	
	LLR_MSA			  	 = LDPCMSAcustom(LLR_part,iteration,frame_data.H_row_master,frame_data.H_row_master_size,5464,0.5);
	Lena_bin_hat 		 = reshape(LLR_MSA(1:648,:).'<0,[],8);
	%Time
	Time = toc;
	%LDPC decode(image)
	bin_table	= 2 .^ (7:-1:0);
	Lena_row = frame_data.Lena_row;
	Lena_col = frame_data.Lena_col;
	Lena_size= frame_data.Lena_size;
	Lena_bin_hat = Lena_bin_hat(1:Lena_size,:);
	Lena_dec_hat = sum(Lena_bin_hat.*bin_table,2);%bin2dec
	Lena_Csize	 = Lena_row*Lena_col;
	Lena_RGB	 = zeros(Lena_row,Lena_col,3);
	Lena_RGB(:,:,1) = reshape( Lena_dec_hat(             1:Lena_Csize  ) ,Lena_row,Lena_col);
	Lena_RGB(:,:,2) = reshape( Lena_dec_hat(Lena_Csize  +1:Lena_Csize*2) ,Lena_row,Lena_col);
	Lena_RGB(:,:,3) = reshape( Lena_dec_hat(Lena_Csize*2+1:Lena_Csize*3) ,Lena_row,Lena_col);
	Lena_RGB		= uint8(Lena_RGB);
	%BER
	Biterror = sum(Lena_bin_hat ~= frame_data.Lena_bin,'all');
	%SNR
	Rx1_SNR  = -10*log10(Rx_No(1) );
	Rx2_SNR  = -10*log10(Rx_No(2) );
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

