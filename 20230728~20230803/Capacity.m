% MIMO system
clear
clc

MIMO_range= 1:6;
frame_num = 200;
SNR_in_dB = 0:5:40;        % 自己設訂雜訊大小
C_in_SNR  = zeros(length(MIMO_range),length(SNR_in_dB));

%% choose QAM= 4/16/64;
QAM 	= 16;
Eavg 	= (qammod([0:QAM-1],QAM) * qammod([0:QAM-1],QAM)') / QAM;
NF 		= 1 / sqrt(Eavg);
q_bit 	= log2(QAM);        % 一個symbol可以傳幾個bit

for smui=1:length(MIMO_range)
	Tx	= MIMO_range(smui);
	Rx	= Tx;
	for a=1:length(SNR_in_dB)
		SNR = 10^( SNR_in_dB(a)/10);
		No  = 10^(-SNR_in_dB(a)/10);
		C 	= 0;
		parfor frame = 1:frame_num	
			fprintf("MIMO : %d/%d\t SNR : %d/%d \t frame : %d/%d\n",smui,length(MIMO_range),a,length(SNR_in_dB),frame,frame_num);
			%mod
			data_dec	= randi([0,QAM-1],1644,560,Tx); 		% 隨機產生0~3 for 4QAM
			data_bin 	= dec2bin(data_dec,q_bit);   	% 將 0~3 轉為 '00'~'11'
			data_mod	= qammod(data_dec,QAM    	 ,'gray'    	)*NF;       % 0~3 to complex (Modulation); remember to normalize
			%Guard Band
			DC =   zeros(   1,560,Tx);
			X  = [ zeros( 202,560,Tx) ;data_mod(1:822,:,:) ;DC ;data_mod(823:end,:,:) ;zeros(201,560,Tx) ];	
			%mod
			data_dec	= randi([0,QAM-1],1644,560,Tx); 		% 隨機產生0~3 for 4QAM
			data_bin 	= dec2bin(data_dec,q_bit);   	% 將 0~3 轉為 '00'~'11'
			data_mod	= qammod(data_dec,QAM    	 ,'gray'    	)*NF;       % 0~3 to complex (Modulation); remember to normalize
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
			%通道與雜訊
			PowerdB 	= [ -2 -8 -10 -12 -15 -18].';
			PowerdB_MIMO= repmat(PowerdB,1,Rx,Tx);
			H_Channel 	= sqrt(10.^(PowerdB_MIMO./10));
			Ntap		= length(PowerdB);
			H_Channel   = H_Channel .* ( sqrt( 1/(2*Tx) ) .* ( randn(Ntap,Rx,Tx) + 1i*randn(Ntap,Rx,Tx) ) );
			%捲積
			pure_y		= zeros(Rx,1228800+5);
			for Tx_n = 1:Tx
				for Rx_n = 1:Rx
					pure_y(Rx_n,:) = pure_y(Rx_n,:) + conv( x_CP(1,:,Tx_n) , H_Channel(:,Rx_n,Tx_n) );
				end
			end
			pure_y(:,1228801:end) = [];	%刪除最後5筆資料
			n 			= sqrt(No/(2/Tx)) *( randn(Tx,1228800) + randn(Tx,1228800)*1i );% randn產生noise variance=No
			%產生訊號
			y			= pure_y + n;
			%移除CP
			y_rmCP		= zeros(2048 , 560, Tx);
			index  = 1;
			for symbol = 1:560
				if mod(symbol,28)-1;
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
			%ZF detector
			h = [H_Channel ; zeros(2042,Rx,Tx)];
			H = fftshift(fft(h,[],1),1);
			H_Data 	= [H(203:1024,:,:);H(1026:1847,:,:)];	%rmGB
			H_frame = permute(repmat( H_Data(:,:,:),1,1,1,560),[1 4 2 3]  );
			for SC = 1:1644
				for slot = 1:560
					unit_H 	= reshape(H_frame(SC,slot,:,:),Rx,Tx);
					C		= C + abs( log2( det( eye(Tx) + (1/Tx) .* SNR .* (unit_H*unit_H'))));
				end
			end
		end
		C_in_SNR(smui,a) = C / (1644 * 560 * frame_num);
	end

end

% Plotting
figure(1)
hold on
for i=1:6
    semilogy(SNR_in_dB, C_in_SNR(i,:), 'LineWidth',2)
end
grid on
axis tight
axis square
title('Capacity of MIMO OFDM')
xlabel('SNR (dB)')
ylabel('Capacity(bits/secHz)')
legend('1x1','2X2','3x3','4x4','5x5','6x6')
