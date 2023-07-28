clear
clc

%假設規範
SNR_in_dB 	= 0:5:40;	% 自己設訂雜訊大小
SNR_weight 	= 40;
window 		= 6;
frame_num 	= 50;		%10ms 10subframe 
DMRS_DATA 	= -0.7071 - 0.7071*1i ;

% choose QAM= 4/16/64;
QAM 	= 16;
Eavg 	= (qammod([0:QAM-1],QAM) * qammod([0:QAM-1],QAM)') / QAM;
NF 		= 1 / sqrt(Eavg);
q_bit 	= log2(QAM);        % 一個symbol可以傳幾個bit
Tx 		= 1;

%output
MSE_dB_LS		=zeros(1,length(SNR_in_dB));
MSE_dB_LMMSE	=zeros(1,length(SNR_in_dB));

%LMMSE 權重矩陣
W = make_W_S(window,SNR_weight);

for a=1:length(SNR_in_dB)
	MSE_LS			= 0;
	MSE_LMMSE		= 0;
	
	parfor frame=1:frame_num
        fprintf("SNR : %d/%d \t frame : %d/%d\n",a,length(SNR_in_dB),frame,frame_num);
		SNR = 10^( SNR_in_dB(a)/10);
		No  = 10^(-SNR_in_dB(a)/10);
		%生產數據
		data_dec_RB	= randi([0,QAM-1],12,14); 		% 隨機產生QAM
		data_mod_RB	= qammod(data_dec_RB,QAM)*NF;   % 0~3 to complex (Modulation); remember to normalize
		%安置DMRS
		data_mod_RB(2:2:12,3) = DMRS_DATA;
		%擴展資料
		data_mod = repmat(data_mod_RB,137,40);
		%Guard Band
		DC =     zeros(  1,14*4*10);
		X  =   [ zeros(202,14*4*10) ;data_mod(1:822,:) ;DC ;data_mod(823:end,:) ;zeros(201,14*4*10) ];	
		%IFFT
		x  = ifft(ifftshift(X,1))*sqrt(2048);	%	2048 x 14*4*10						
		%CP
		x_CP = zeros(1,1228800);
		index=1;
		for symbol=1:14*4*10
			if	mod(symbol,28)-1
				x_CP(1, index:index+2048+144-1)=[ x(2048-144+1:2048,symbol) ; x(:,symbol)];
				index = index+2048+144;
			else
				x_CP(1, index:index+2048+208-1)=[ x(2048-208+1:2048,symbol) ; x(:,symbol)];
				index = index+2048+208;
			end
		end
		%通道與雜訊
		PowerdB 		= [ -2 -8 -10 -12 -15 -18];
		Total_H_Power 	= sum(10.^(PowerdB/10)); 	%總通道能量 = 1
		Ntap = 6;									%通道數量
		H_Channel 	=  sqrt(10.^(PowerdB/10));
		H_Channel   =  H_Channel .* ( sqrt( 1/(2*Tx) ) * ( randn(1,Ntap) + 1i*randn(1,Ntap) ) );
		%訊號通過通道
		H_y					= conv( x_CP, H_Channel );
		H_y(:,1228801:end)  = [];			%刪除最後5筆資料		
		%產生訊號
		n 			= sqrt(No/2) *( randn(1,1228800) + randn(1,1228800)*1i );% randn產生noise variance=No
		y			= H_y+ n;
		%移除CP
		y_rmCP		= zeros(2048 , 560);
		index  = 1;
		for symbol = 1:560
			if mod(symbol,28)-1;
				y_rmCP(:,symbol) = y(1,index+144:index+144+2048-1);
				index  = index+144+2048;
			else
				y_rmCP(:,symbol) = y(1,index+208:index+208+2048-1);
				index  = index +208+2048;
			end
		end
		%FFT----以下為頻域 (test OK
		Y_fft 	= fftshift( fft( y_rmCP/sqrt(2048) ) ,1);
		%rm Guard Band(test OK
		Y	  	= [ Y_fft( 203:1024,:) ; Y_fft( 1026:1847,:) ];
		%實際通道  (test OK
		h		= [H_Channel,zeros(1,2042)];
		H		= fftshift(fft(h));
		H_Data 	= [H(1,203:1024),H(1,1026:1847)].';	%rmGB
		H_frame = repmat(H_Data,1,560);				%allframe
		%取得LS估測結果，20個slot 每個都要估
		X_LS 	= DMRS_DATA * eye(1644/2);
		Y_LS 	= Y(2:2:1644 , 3:14:560);
		H_LS 	= inv(X_LS)*Y_LS;
		%取得LMMSE估測結果
		H_LMMSE 	= W * H_LS;
		%產生真實通道
		H_LS_R 		= H_frame(2:2:1644 , 3:14:560);
		H_LMMSE_R 	= H_frame(1:1:1644 , 3:14:560);
		%累加
		MSE_LS		= MSE_LS    + sum( abs( H_LS_R    - H_LS    ).^2,'all');
		MSE_LMMSE	= MSE_LMMSE + sum( abs( H_LMMSE_R - H_LMMSE ).^2,'all');
	end
	%計算MSE數值
	MSE_dB_LS   (a)	= 10*log10( MSE_LS    / ( 822*40*frame_num) );
	MSE_dB_LMMSE(a)	= 10*log10( MSE_LMMSE / (1644*40*frame_num) );
end

figure(1);
plot(SNR_in_dB,MSE_dB_LS   ,'r-','LineWidth',2);
hold on;
plot(SNR_in_dB,MSE_dB_LMMSE,'b-','LineWidth',2);
hold on;
title('LS 與 LMMSE 的MSE比較');
grid on;
xlabel('SNR (dB)');
ylabel('MSE (dB)');
legend('LS MSE','LMMSE MSE');