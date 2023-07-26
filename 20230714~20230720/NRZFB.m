clear
clc
warning('off')

SNR_in_dB = 0:1:16;
BER_SNR_MATLAB=zeros(1,length(SNR_in_dB));
BER_SNR_C=zeros(1,length(SNR_in_dB));

QAM 	= 4;
Eavg 	= (qammod([0:QAM-1],QAM) * qammod([0:QAM-1],QAM)') / QAM;
NF 		= 1 / sqrt(Eavg);
q_bit 	= 2;

for a=1:length(SNR_in_dB)
	SNR = 10^( SNR_in_dB(a)/10);
	No  = 10^(-SNR_in_dB(a)/10);
	BER_MATLAB = 0;
    BER_C = 0;

	data_dec	= randi([0,QAM-1],1644,560); 
    data_bin = de2bi(data_dec, q_bit, 'left-msb');
	data_mod = NF*qammod(data_dec, QAM, 'gray');
	DC =   zeros(1,560);
	X  = [ zeros(202,560) ;data_mod(1:822,:) ;DC ;data_mod(823:end,:) ;zeros(201,560) ];
	x  = ifft(ifftshift(X))*sqrt(2048);										
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

	PowerdB 	= [ -2 -8 -10 -12 -15 -18];
	H_Channel 	= sqrt(10.^(PowerdB/10));
	Total_H_Power = sum(10.^(PowerdB/10)); 
	n 			= sqrt(No/2) *( randn(1,1228800) + randn(1,1228800)*1i );
    
    % Using MATLAB conv
	pure_y_MATLAB	= conv( x_CP, H_Channel );
	pure_y_MATLAB(:,1228801:end) = [];
	y_MATLAB			= pure_y_MATLAB + n;
    
    % Using myConv
    pure_y_C	= myConv( x_CP, H_Channel );
	pure_y_C(:,1228801:end) = [];
	y_C			= pure_y_C + n;

	y_rmCP_MATLAB		= zeros(2048 , 560);
    y_rmCP_C		= zeros(2048 , 560);
	index  = 1;
	for symbol = 1:560
		if mod(symbol,28)-1;
			y_rmCP_MATLAB(:,symbol) = y_MATLAB(1,index+144:index+144+2048-1);
            y_rmCP_C(:,symbol) = y_C(1,index+144:index+144+2048-1);
			index  = index+144+2048;
			symbol = symbol+1;
		else
			y_rmCP_MATLAB(:,symbol) = y_MATLAB(1,index+208:index+208+2048-1);
            y_rmCP_C(:,symbol) = y_C(1,index+208:index+208+2048-1);
			index  = index+208+2048;
			symbol = symbol+1;
		end
	end

	Y_fft_MATLAB = fftshift( fft( y_rmCP_MATLAB/sqrt(2048) ) );
    Y_fft_C = fftshift( fft( y_rmCP_C/sqrt(2048) ) );
	Y_MATLAB	  = [ Y_fft_MATLAB( 203:1024,:) ; Y_fft_MATLAB( 1026:1847,:) ];
    Y_C	  = [ Y_fft_C( 203:1024,:) ; Y_fft_C( 1026:1847,:) ];
    
	h=[H_Channel,zeros(1,2042)];
	H = fftshift(fft(h));
	H_Data = [H(1,202:1023),H(1,1025:1846)].';
	H_frame = repmat(H_Data,1,560);

	X_hat_MATLAB = (Y_MATLAB ./ H_frame)/NF;
    X_hat_C = (Y_C ./ H_frame)/NF;
	data_dec_hat_MATLAB = qamdemod(X_hat_MATLAB,QAM,'gray');
    data_dec_hat_C = qamdemod(X_hat_C,QAM,'gray');
	data_bin_hat_MATLAB = de2bi(data_dec_hat_MATLAB, q_bit, 'left-msb');
    data_bin_hat_C = de2bi(data_dec_hat_C, q_bit, 'left-msb');
    
    BER_MATLAB 		 =  sum(sum(data_bin ~= data_bin_hat_MATLAB ));
	BER_SNR_MATLAB(1,a) = 	BER_MATLAB/(920640 * 2);
    
    BER_C 		 =  sum(sum(data_bin ~= data_bin_hat_C ));
	BER_SNR_C(1,a) = 	BER_C/(920640 * 2);
end

figure(1)
semilogy(SNR_in_dB,BER_SNR_MATLAB(1,:),'r-X', 'LineWidth',2)
hold on
semilogy(SNR_in_dB,BER_SNR_C(1,:),'b--O', 'LineWidth',2)
grid on
axis square
axis tight
title('BER of SISO')
xlabel('SNR (dB)')
ylabel('BER')
legend('4QAM ZF MATLAB', '4QAM ZF C')
