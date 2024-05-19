% SISO system
clear
clc
warning('off')

SNR_in_dB = 0:1:16;        % 自己設訂雜訊大小
BER_SNR=zeros(1,length(SNR_in_dB));

QAM 	= 4;
Eavg 	= (qammod([0:QAM-1],QAM) * qammod([0:QAM-1],QAM)') / QAM;
NF 		= 1 / sqrt(Eavg);
q_bit 	= 2;        % 一個symbol可以傳幾個bit

for a=1:length(SNR_in_dB)
	SNR = 10^( SNR_in_dB(a)/10);
	No  = 10^(-SNR_in_dB(a)/10);
	BER = 0;                         % 算error rate要作平均
	%mod
	data_dec	= randi([0,QAM-1],1644,560); 		% 隨機產生0~3 for 4QAM
    data_bin = de2bi(data_dec, q_bit, 'left-msb');	% 將 0~3 轉為 '00'~'11'
	data_mod = NF*qammod(data_dec, QAM, 'gray');    % 0~3 to complex (Modulation); remember to normalize
	%Guard Band
	DC =   zeros(1,560);
	X  = [ zeros(202,560) ;data_mod(1:822,:) ;DC ;data_mod(823:end,:) ;zeros(201,560) ];	
	%IFFT
	x  = ifft(ifftshift(X))*sqrt(2048);	%	2048 x 560								
	%CP
	x_CP = zeros(1,1228800);
	index=1;
	for symbol=1:560
		if	mod(symbol,28)-1
			x_CP(1, index:index+2048+144-1)=[ x(2048-144+1:2048,symbol) ; x(:,symbol)];
			index = index+2048+144;
		else
			x_CP(1, index:index+2048+208-1)=[ x(2048-208+1:2048,symbol) ; x(:,symbol)];
			index = index+2048+208;
		end
	end
	%通道與雜訊
	PowerdB 	= [ -2 -8 -10 -12 -15 -18];
	H_Channel 	= sqrt(10.^(PowerdB/10));
	Total_H_Power = sum(10.^(PowerdB/10)); %總通道能量 = 1
	pure_y		= conv( x_CP, H_Channel );
	pure_y(:,1228801:end) = [];			%刪除最後5筆資料
	n 			= sqrt(No/2) *( randn(1,1228800) + randn(1,1228800)*1i );% randn產生noise variance=No
	%產生訊號
	y			= pure_y + n;
	%移除CP
	y_rmCP		= zeros(2048 , 560);
	index  = 1;
	for symbol = 1:560
		if mod(symbol,28)-1;
			y_rmCP(:,symbol) = y(1,index+144:index+144+2048-1);
			index  = index+144+2048;
			symbol = symbol+1;
		else
			y_rmCP(:,symbol) = y(1,index+208:index+208+2048-1);
			index  = index+208+2048;
			symbol = symbol+1;
		end
	end
	%FFT
	Y_fft = fftshift( fft( y_rmCP/sqrt(2048) ) );%%%%%%%%
	%rm Guard Band
	Y	  = [ Y_fft( 203:1024,:) ; Y_fft( 1026:1847,:) ];
	%ZF detector
	h=[H_Channel,zeros(1,2042)];
	H = fftshift(fft(h));
	H_Data = [H(1,203:1024),H(1,1026:1847)].';
	H_frame = repmat(H_Data,1,560);
	X_hat = (Y ./ H_frame)/NF;
	%demod
	data_dec_hat = qamdemod(X_hat,QAM,'gray');
	data_bin_hat = de2bi(data_dec_hat, q_bit, 'left-msb');     	%0~3 to '00'~'11'
	%BER計算
    BER 		 =  sum(sum(data_bin ~= data_bin_hat ));
	BER_SNR(1,a) = 	BER/(920640 * 2);
end

%輸入SNR_in_dB和BER
figure(1)
semilogy(SNR_in_dB,BER_SNR(1,:),'r-x', 'LineWidth',2)
hold on
grid on
axis square
axis tight
title('BER of SISO')
xlabel('SNR (dB)')
ylabel('BER')
legend('4QAM ZF MATLAB')