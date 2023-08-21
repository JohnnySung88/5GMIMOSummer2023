y = Rx_signal;%2*1228800
tic;
%CFO估測(需多重路徑猜測)
CFO_ignore		= 5;%尚待自訂益
index  			= 1;
hat_delta_f_sym = zeros(560,2);
CFO_sum	 		= 0;
CP_num			= 0;
t 				= 0:(10e-3/1228800):10e-3;
t 				= t(2:1228801);				%移除第一個位置資料
for symbol = 1:560
	if (mod(symbol,28)-1)%size 144
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
		index    = index + 208 + 2048;
	end
end
hat_delta_f = -(Fs/2048) * angle( CFO_sum/CP_num )/(2*pi) ;
%CFO補償
y = y ./ exp( 1i * 2 * pi * hat_delta_f * t);
%移除CP
y_rmCP		= zeros(2048 , 560, Tx);
CP_tpos		= [];
index  = 1;
for symbol = 1:560
	if (mod(symbol,28)-1)
		for ram = 1:Rx
			y_rmCP(:,symbol,ram) = y(ram,index+144:index+144+2048-1);
		end
		CP_tpos = [CP_tpos,index+144:index+144+2048-1];
		index  = index+144+2048;
	else
		for ram = 1:Rx
			y_rmCP(:,symbol,ram) = y(ram,index+208:index+208+2048-1);
		end
		CP_tpos = [CP_tpos,index+208:index+208+2048-1];
		index  = index +208+2048;
	end
end
%FFT----以下為頻域
Y_fft 	= fftshift( fft( y_rmCP/sqrt(2048) ) ,1);
%rm Guard Band
Y	  	= [ Y_fft( 203:1024,:,:) ; Y_fft( 1026:1847,:,:) ];
Y_DMRS 	= Y(2:2:1644 , 3:14:560,:);
H_LMMSE			 = zeros(1644,40,2,2);
H_LMMSE(:,:,1,1) = LMMSE_1_1 * (DMRS_DATA' .*Y_DMRS(:,:,1));
H_LMMSE(:,:,2,1) = LMMSE_2_1 * (DMRS_DATA' .*Y_DMRS(:,:,2));
H_LMMSE(:,:,1,2) = LMMSE_1_2 * (DMRS_DATA' .*Y_DMRS(:,:,1));
H_LMMSE(:,:,2,2) = LMMSE_2_2 * (DMRS_DATA' .*Y_DMRS(:,:,2));
%線性內差(LMMSE only)
DMRS_Spos	= 3:14:560;
H_INTER		= zeros(1644,560,2,2);	
for symbol  = 549:560  %邊界
	head_dist	= symbol - 549;
	back_dist	= 563    - symbol;
	H_INTER(:,symbol,:,:) = ( back_dist * H_LMMSE(:,40,:,:) + head_dist * H_LMMSE(:,1,:,:)  )  /14;
end
for symbol  = 1:2     %邊界
	head_dist	= 11 + symbol;
	back_dist	= 3  - symbol;
	H_INTER(:,symbol,:,:) = ( back_dist * H_LMMSE(:,40,:,:) + head_dist * H_LMMSE(:,1,:,:)  )  /14;
end
for symbol 	= 3:548  %連續
	pos  		= floor((symbol-3)/14) + 1;
	head_dist	= symbol 			- DMRS_Spos(pos);
	back_dist	= DMRS_Spos(pos+1) 	- symbol;
	H_INTER(:,symbol,:,:) = ( back_dist * H_LMMSE(:,pos,:,:) + head_dist * H_LMMSE(:,pos+1,:,:)  )  /14;
end
%雜訊估測
DMRS  		= LDPC_mod(2:2:1644 , 3:14:560,:);%%%%%%%
DMRS_H		= H_LMMSE (2:2:1644,:,:,:);
DMRS_hat_Y  = zeros(822,40,2);
for SC = 1:822
	for symbol = 1:40
		DMRS_hat_Y(SC,symbol,:) = reshape(DMRS_H(SC,symbol,:,:),2,2) * reshape(DMRS(SC,symbol,:),2,1);
	end
end
Rx1_No = sum( abs( Y_DMRS(:,:,1)  - DMRS_hat_Y(:,:,1) ).^2 ,'all' )/(40*822);
Rx2_No = sum( abs( Y_DMRS(:,:,2)  - DMRS_hat_Y(:,:,2) ).^2 ,'all' )/(40*822);
%detector
switch DTinfo
	case 'LMMSE'
		norm_Y = zeros(1644,560,2);
		norm_H = zeros(1644,560,2,2);
		norm_Y(:,:,1)   = Y(:,:,1)          ./ Rx1_No;
		norm_Y(:,:,2)   = Y(:,:,2)          ./ Rx2_No;
		norm_H(:,:,1,:) = H_INTER(:,:,1,:)  ./ Rx1_No;
		norm_H(:,:,2,:) = H_INTER(:,:,2,:)  ./ Rx2_No;
		X_hat = LMMSE(norm_Y,norm_H,Tx);
		X_hat = X_hat/NF;
	case 'ZF'
		X_hat = ZF_detector(Y,H_frame,Tx);
        X_hat = X_hat/NF;
end
%反解資料
LDPC_mod_L_hat	= X_hat(DATA_Pos);
LDPC_mod_L_hat  = LDPC_mod_L_hat(1:1770336);%magic num
%LLR
yI = real(LDPC_mod_L_hat).';
yQ = imag(LDPC_mod_L_hat).';
LLR= zeros( 4*length(LDPC_mod_L_hat) ,1);
LLR(1:4:length(LLR)) = (1/(2*No)) .* ( min( (yI-[ 1; 3]).^2 ) - min( (yI-[-1;-3]).^2 ));
LLR(2:4:length(LLR)) = (1/(2*No)) .* ( min( (yI-[-1; 1]).^2 ) - min( (yI-[-3; 3]).^2 ));
LLR(3:4:length(LLR)) = (1/(2*No)) .* ( min( (yQ-[-1;-3]).^2 ) - min( (yQ-[ 1; 3]).^2 ));
LLR(4:4:length(LLR)) = (1/(2*No)) .* ( min( (yQ-[-1; 1]).^2 ) - min( (yQ-[-3; 3]).^2 ));
LLR_part			 = reshape(LLR,5464,1296);
LLR_part			 = permute(LLR_part,[2,1]);

%LDPC decode(image)
bin_table	= 2 .^ (7:-1:0);
LLR_MSA			  = LDPCMSAcustom(LLR_part,iteration,H_row_master,H_row_master_size,5464,0.5);
Lena_bin_hat_LDPC = reshape(LLR_MSA(1:648,:).'<0,[],8);
Lena_bin_hat_LDPC = Lena_bin_hat_LDPC(1:Lena_size,:);
Lena_bin_hat_LDPC = sum(Lena_bin_hat_LDPC.*bin_table,2);%bin2dec
Lena_Csize		  = Lena_row*Lena_col;
Lena_RGB_LDPC	  = zeros(Lena_row,Lena_col,3);
Lena_RGB_LDPC(:,:,1) = reshape( Lena_bin_hat_LDPC(             1:Lena_Csize  ) ,Lena_row,Lena_col);
Lena_RGB_LDPC(:,:,2) = reshape( Lena_bin_hat_LDPC(Lena_Csize  +1:Lena_Csize*2) ,Lena_row,Lena_col);
Lena_RGB_LDPC(:,:,3) = reshape( Lena_bin_hat_LDPC(Lena_Csize*2+1:Lena_Csize*3) ,Lena_row,Lena_col);
%imshow( uint8(Lena_RGB_LDPC) );

%C
for SC = 1:1644
	for slot = 1:560
		unit_H 	= reshape(H_INTER(SC,slot,:,:),Rx,Tx);
		C		= C + abs( log2( det( eye(Tx) + (1/Tx) .* SNR .* (unit_H*unit_H'))));
	end
end