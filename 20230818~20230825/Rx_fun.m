function [Time,Biterror,Capacity_sum,Rx1_SNR,Rx2_SNR,Pict_RGB] = Rx_fun(frame_data,Rx_signal,DTinfo,CFO_ignore)
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
	Pict_bin_hat   = reshape(LDPC_bin_part(1:648,:).',[],8);
	%Time
	Time = toc;
	%noLDPC decode(image)
	bin_table	= 2 .^ (7:-1:0);
	Pict_row = frame_data.Pict_row;
	Pict_col = frame_data.Pict_col;
	Pict_size= frame_data.Pict_size;
	Pict_bin_hat = Pict_bin_hat(1:Pict_size,:);
	Pict_dec_hat = sum(Pict_bin_hat.*bin_table,2);%bin2dec
	Pict_Csize	 = Pict_row*Pict_col;
	Pict_RGB	 = zeros(Pict_row,Pict_col,3);
	Pict_RGB(:,:,1) = reshape( Pict_dec_hat(             1:Pict_Csize  ) ,Pict_row,Pict_col);
	Pict_RGB(:,:,2) = reshape( Pict_dec_hat(Pict_Csize  +1:Pict_Csize*2) ,Pict_row,Pict_col);
	Pict_RGB(:,:,3) = reshape( Pict_dec_hat(Pict_Csize*2+1:Pict_Csize*3) ,Pict_row,Pict_col);
	Pict_RGB		= uint8(Pict_RGB);
	%BER
	Biterror = sum(Pict_bin_hat ~= frame_data.Pict_bin,'all');
	%SNR
	Rx1_SNR  = -10*log10(Rx_No(1) );
	Rx2_SNR  = -10*log10(Rx_No(2) );
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

