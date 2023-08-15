%MIMO system
clear
clc
rng(123);
%假設規範
frame_num 	= 500;
SNR_in_dB 	= 0:5:40;        % 自己設訂雜訊大小
SNR_weight 	= 45;
window 		= 10;
DMRS_DATA 	= +0.7071 + 0.7071*1i ;
CFO_ignore	= 5;
Fs			= 122.88e6;
delta_f		= 1000;

%Data pos
DMRS_SC = [2:2:1644];
Data_SC = [1:1:1644];
DMRS_sym= [3:14:560];

%LMMSE init
N_fft 	=2048;
R_HD_HD	=zeros( 822,822);
R_H_HD	=zeros(1644,822);
DMRS_pos=[204:2:1024,1027:2:1847];
Real_pos=[203:1:1024,1026:1:1847];
for p=1:822
	for k=1:822%R_HD_HD 822*822
		if DMRS_pos(k)==DMRS_pos(p)
			R_HD_HD  (k,p)=1;
		else
			R_HD_HD  (k,p)=( 1 - exp( -1i*2*pi*window*(DMRS_pos(k)-DMRS_pos(p))/N_fft ) )/( 1i*2*pi*window*(DMRS_pos(k)-DMRS_pos(p))/N_fft );
		end
	end
	for k=1:1644
		if Real_pos(k)==DMRS_pos(p)
			R_H_HD(k,p)=1;
		else
			R_H_HD(k,p)=( 1 - exp( -1i*2*pi*window*(Real_pos(k)-DMRS_pos(p))/N_fft ) )/( 1i*2*pi*window*(Real_pos(k)-DMRS_pos(p))/N_fft );
		end
	end
end
SNR_W 	= 10^( SNR_weight /10);


% choose QAM= 4/16/64;
QAM 	= 16;
Eavg 	= (qammod([0:QAM-1],QAM) * qammod([0:QAM-1],QAM)') / QAM;
NF 		= 1 / sqrt(Eavg);
q_bit 	= log2(QAM);        % 一個symbol可以傳幾個bit
Tx		  	= 2;
Rx		  	= 2;

%output
BER_SNR=zeros(1,length(SNR_in_dB));

for a=1:length(SNR_in_dB)
	SNR = 10^( SNR_in_dB(a)/10);
	No  = 10^(-SNR_in_dB(a)/10);
	BER = 0;                         % 算error rate要作平均
	fprintf("SNR : %d/%d \n\n",a,length(SNR_in_dB));
    parfor_progress(frame_num);
	parfor frame=1:frame_num	
		%輸入資料(含DMRS)
		data_dec_L	= randi  ([0,QAM-1], Tx*(1644*560-822*40),1);
		data_mod_L	= qammod (data_dec_L,  QAM,'gray')*NF;
		data_bin_L 	= dec2bin(data_dec_L,q_bit);
		data_mod	= zeros(1644,560,2);
		index		= 1;
		for symbol=1:560
			if mod( (symbol - 3) , 14) == 0
				data_range 		     = reshape(data_mod_L(index:index+ 822*Tx-1),1, 822,Tx );
				data_mod(:,symbol,:) = reshape([ data_range  ; DMRS_DATA * ones(1,822,Tx)],1644,2);
				index				 = index + 822*Tx;
			else
				data_mod(:,symbol,:) = reshape(data_mod_L(index:index+1644*Tx-1),1644,Tx );
				index				 = index + 1644*Tx;
			end
		end
		%CDM
		data_mod(2:4:1644,3:14:560,2:2:Tx) = -data_mod(2:4:1644,3:14:560,2:2:Tx);
		%Guard Band
		DC =   zeros(   1,560,Tx);
		X  = [ zeros( 202,560,Tx) ;data_mod(1:822,:,:) ;DC ;data_mod(823:end,:,:) ;zeros(201,560,Tx) ];	
		%IFFT
		x  = ifft(ifftshift(X,1))*sqrt(2048);	%	2048 x 14*4*10						
		%CP
		x_CP = zeros(1,1228800,Tx);
		index=1;
		for symbol=1:14*4*10
			if	mod(symbol,28)-1
				x_CP(1, index:index+2048+144-1,:)=[ x(2048-144+1:2048,symbol,:) ; x(:,symbol,:)];
				index = index+2048+144;
			else
				x_CP(1, index:index+2048+208-1,:)=[ x(2048-208+1:2048,symbol,:) ; x(:,symbol,:)];
				index = index+2048+208;
			end
		end
		%通道
		PowerdB 	= [ -2 -8 -10 -12 -15 -18].';
		PowerdB_MIMO= repmat(PowerdB,1,Rx,Tx);
		H_Channel 	= sqrt(10.^(PowerdB_MIMO./10)).* sqrt( 1/ Tx );
		Ntap		= length(PowerdB);
		H_Channel   = H_Channel .* ( sqrt( 1/2 ) .* ( randn(Ntap,Rx,Tx) + 1i*randn(Ntap,Rx,Tx) ) );
		%捲積
		pure_y		= zeros(Rx,1228800+5);
		for Tx_n = 1:Tx
			for Rx_n = 1:Rx
				pure_y(Rx_n,:) = pure_y(Rx_n,:) + conv( x_CP(1,:,Tx_n) , H_Channel(:,Rx_n,Tx_n) );
			end
		end
		pure_y(:,1228801:end) = [];	%刪除最後5筆資料
		%產生訊號
		n 			= sqrt(No/2) *( randn(Tx,1228800) + randn(Tx,1228800)*1i );% randn產生noise variance=No
		y			= pure_y + n;
		%頻偏
		t 			= 0:(10e-3/1228800):10e-3;
		t 			= t(2:1228801);				%移除第一個位置資料
		y			= exp( 1i * 2 * pi * delta_f * t) .* y;
		%CFO估測(需多重路徑猜測)
		index  			= 1;
		hat_delta_f_sym = zeros(560,2);
		CFO_sum	 		= 0;
		CP_num			= 0;
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
		%
		y = y ./ exp( 1i * 2 * pi * hat_delta_f * t);
		%移除CP
		y_rmCP		= zeros(2048 , 560, Tx);
		index  = 1;
		for symbol = 1:560
			if (mod(symbol,28)-1)
				for ram = 1:Rx
					y_rmCP(:,symbol,ram) = y(ram,index+144:index+144+2048-1);
				end
				index  = index+144+2048;
			else
				for ram = 1:Rx
					y_rmCP(:,symbol,ram) = y(ram,index+208:index+208+2048-1);
				end
				index  = index +208+2048;
			end
		end
		%FFT----以下為頻域
		Y_fft 	= fftshift( fft( y_rmCP/sqrt(2048) ) ,1);
		%rm Guard Band
		Y	  	= [ Y_fft( 203:1024,:,:) ; Y_fft( 1026:1847,:,:) ];
		%取得LMMSE估測結果
		Y_DMRS 			 = Y(2:2:1644 , 3:14:560,:);
		P0 				 = eye(822);
		P1 				 = ones(822,1);P1(1:2:822) = -1;
		P1 				 = diag(P1);
		H_LMMSE			 = zeros(1644,40,2,2);
		H_LMMSE(:,:,1,1) = R_H_HD* P0' * inv( P0*R_HD_HD*P0' + P1*R_HD_HD*P1' + (1/SNR_W)*eye(822) ) * (DMRS_DATA' .*Y_DMRS(:,:,1));
		H_LMMSE(:,:,2,1) = R_H_HD* P0' * inv( P0*R_HD_HD*P0' + P1*R_HD_HD*P1' + (1/SNR_W)*eye(822) ) * (DMRS_DATA' .*Y_DMRS(:,:,2));
		H_LMMSE(:,:,1,2) = R_H_HD* P1' * inv( P0*R_HD_HD*P0' + P1*R_HD_HD*P1' + (1/SNR_W)*eye(822) ) * (DMRS_DATA' .*Y_DMRS(:,:,1));
		H_LMMSE(:,:,2,2) = R_H_HD* P1' * inv( P0*R_HD_HD*P0' + P1*R_HD_HD*P1' + (1/SNR_W)*eye(822) ) * (DMRS_DATA' .*Y_DMRS(:,:,2));
		%線性內差
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
		DMRS  		= data_mod(2:2:1644 , 3:14:560,:);
		DMRS_H		= H_LMMSE (2:2:1644,:,:,:);
		DMRS_hat_Y  = zeros(822,40,2);
		for SC = 1:822
			for symbol = 1:40
				DMRS_hat_Y(SC,symbol,:) = reshape(DMRS_H(SC,symbol,:,:),2,2) * reshape(DMRS(SC,symbol,:),2,1);
			end
		end
		Rx1_No = sum( abs( Y_DMRS(:,:,1)  - DMRS_hat_Y(:,:,1) ).^2 ,'all' )/(40*822);
		Rx2_No = sum( abs( Y_DMRS(:,:,2)  - DMRS_hat_Y(:,:,2) ).^2 ,'all' )/(40*822);
		%LMMSE
        norm_Y = zeros(1644,560,2);
		norm_H = zeros(1644,560,2,2);
		norm_Y(:,:,1)   = Y(:,:,1)          ./ Rx1_No;
        norm_Y(:,:,2)   = Y(:,:,2)          ./ Rx2_No;
        norm_H(:,:,1,:) = H_INTER(:,:,1,:)  ./ Rx1_No;
        norm_H(:,:,2,:) = H_INTER(:,:,2,:)  ./ Rx2_No;
        X_hat_L = LMMSE(norm_Y,norm_H,Tx);
        X_hat_L = X_hat_L/NF;
		%反解資料
		data_mod_L_hat 	= zeros(Tx*(1644*560-822*40),1);
		index			= 1;
		for symbol=1:560
			if mod( (symbol - 3) , 14) == 0
				data_mod_L_hat(index:index+ 822*Tx-1)	 = reshape(X_hat_L(1:2:1644,symbol,:), 822*Tx,1 );
				index				 = index + 822*Tx;
			else
				data_mod_L_hat(index:index+1644*Tx-1)	 = reshape(X_hat_L(1:1:1644,symbol,:),1644*Tx,1 );
				index				 = index + 1644*Tx;
			end
		end
		
		%demod
		data_dec_L_hat = qamdemod(data_mod_L_hat,  QAM,'gray');
		data_bin_L_hat = dec2bin (data_dec_L_hat,q_bit);     	%0~3 to '00'~'11'
		%BER計算
		BER 		 =  BER + sum(sum(data_bin_L_hat ~= data_bin_L ));
		parfor_progress;
	end
	parfor_progress(0);
	BER_SNR(1,a) = 	BER/(1644 * 560 * frame_num * q_bit * Tx );
end

%輸入SNR_in_dB和BER
figure(1)
semilogy(SNR_in_dB,BER_SNR(1,:),'k--','LineWidth',2)
hold on
grid on
title('BER of 5G NR MIMO OFDM(CFO)')
xlabel('SNR (dB)')
ylabel('BER')
legend('16QAM CFO')