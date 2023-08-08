%MIMO system
clear;
clc;

Tx = 2; % number of transmit antennas
Rx = 2; % number of receive antennas
Data_num = 200;
SNR_in_dB = 0:5:40;
% BER_SNR=zeros(1,length(SNR_in_dB));
SNR_weight  = 45;
window      = 10;
DMRS_DATA   = -0.7071 - 0.7071*1i ;

QAM     = 16;
Eavg    = (qammod([0:QAM-1],QAM) * qammod([0:QAM-1],QAM)') / QAM;
NF      = 1 / sqrt(Eavg);
q_bit   = log2(QAM);        % 一個symbol可以傳幾個bit

%output
MSE_dB_LS               =zeros(1,length(SNR_in_dB));
MSE_dB_LMMSE            =zeros(1,length(SNR_in_dB));
MSE_dB_INTERPOLATION    =zeros(1,length(SNR_in_dB));

%LMMSE 權重矩陣
W = make_W_S(window,SNR_weight);

for a=1:length(SNR_in_dB)
    MSE_LS			= 0;
	MSE_LMMSE		= 0;
    MSE_dB_INTERPOLATION = 0;
    SNR = 10^(SNR_in_dB(a)/10);
    No  = 1/SNR;
    % BER = 0;
    for Data=1:Data_num
    fprintf("SNR : %d/%d \t Data : %d/%d\n",a,length(SNR_in_dB),Data,Data_num);

        % TX生產數據
        data_dec    = randi([0,QAM-1],1644,560,Tx);
        data_bin = de2bi(data_dec, q_bit, 'left-msb');  
        data_mod = NF*qammod(data_dec, QAM, 'gray');
        
        %安置DMRS
		data_mod_RB(2:2:12,3) = DMRS_DATA;

        %擴展資料
		data_mod = repmat(data_mod_RB,137,40);

        %Guard Band
        DC =   zeros(1,560,Tx);
        X  = [ zeros(202,560,Tx) ;data_mod(1:822,:,:) ;DC ;data_mod(823:end,:,:) ;zeros(201,560,Tx) ];    

        %IFFT
        x  = ifft(ifftshift(X))*sqrt(2048); %   2048 x 560

        %CP
        x_CP = zeros(1,1228800,Tx);
        index=1;
        for symbol=1:560
            if  mod(symbol,28)-1
                x_CP(1, index:index+2048+144-1,:)=[ x(2048-144+1:2048,symbol,:) ; x(:,symbol,:)];
                index = index+2048+144;
            else
                x_CP(1, index:index+2048+208-1,:)=[ x(2048-208+1:2048,symbol,:) ; x(:,symbol,:)];
                index = index+2048+208;
            end
        end

        % Channel & Noise
        PowerdB     = [ -2 -8 -10 -12 -15 -18].';
        PowerdB_MIMO= repmat(PowerdB,1,Rx,Tx);
        Ntap        = length(PowerdB);
        H_Channel   = sqrt(10.^(PowerdB_MIMO./10));
        H_Channel   = H_Channel .* ( sqrt( 1/(2*Tx) ) .* ( randn(Ntap,Rx,Tx) + 1i*randn(Ntap,Rx,Tx) ) );
        
        %捲積
        pure_y      = zeros(Rx,1228800+5);
        for Tx_n = 1:Tx
            for Rx_n = 1:Rx
                pure_y(Rx_n,:) = pure_y(Rx_n,:) + conv( x_CP(1,:,Tx_n) , H_Channel(:,Rx_n,Tx_n) );
            end
        end
        pure_y(:,1228801:end) = []; %刪除最後5筆資料
        n           = sqrt(No/(2)) *( randn(Tx,1228800) + randn(Tx,1228800)*1i );% randn產生noise variance=No
        
        %產生訊號
        y           = pure_y + n;

        % RX
        %移除CP
        y_rmCP      = zeros(2048 , 560, Tx);
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
        Y_fft   = fftshift( fft( y_rmCP/sqrt(2048) ) ,1);
        
        %rm Guard Band
        Y       = [ Y_fft( 203:1024,:,:) ; Y_fft( 1026:1847,:,:) ];
        
        % ZF Detector
        % h = [H_Channel ; zeros(2042,Rx,Tx)];
        % H = fftshift(fft(h,[],1),1);
        % H_Data  = [H(203:1024,:,:);H(1026:1847,:,:)];   %rmGB
        % H_frame = permute(repmat( H_Data(:,:,:),1,1,1,560),[1 4 2 3]  );
        % X_hat = zeros(1644,560,Tx);
        % X_hat = ZFDC(Y,H_frame,Tx);
        % X_hat = X_hat/NF;

        %實際通道
        h		= [H_Channel ; zeros(2042,Rx,Tx)];
		H = fftshift(fft(h,[],1),1);
        H_Data  = [H(203:1024,:,:);H(1026:1847,:,:)];   %rmGB
        H_frame = permute(repmat( H_Data(:,:,:),1,1,1,560),[1 4 2 3]  );

        % LS Detector
        X_LS 	= DMRS_DATA * eye(1644/2);
		Y_LS 	= Y(2:2:1644 , 3:14:560);
		H_LS 	= inv(X_LS)*Y_LS;

        % LMMSE Detector
        X_LMMSE = zeros(1644,560,Tx);
        X_LMMSE = LMMSE(Y,H_frame,Tx);
        X_LMMSE = X_LMMSE/NF;

        % Interpolation Detector
        X_Interpolation = zeros(1644,560,Tx);
        X_Interpolation = InterpolationDetector(Y, H_frame, Tx, W);
        X_Interpolation = X_Interpolation / NF;

        % De-Mod
        % data_dec_hat = qamdemod(X_hat,QAM,'gray');
        % data_bin_hat = de2bi(data_dec_hat, q_bit, 'left-msb');
        
        %產生真實通道
		H_LS_R 		= H_frame(2:2:1644 , 3:14:560);
		H_LMMSE_R 	= H_frame(1:1:1644 , 3:14:560);
        H_Interpolation_R = H_frame(1:1:1644, 3:14:560);
		%累加
		MSE_LS		= MSE_LS    + sum( abs( H_LS_R    - H_LS    ).^2,'all');
		MSE_LMMSE	= MSE_LMMSE + sum( abs( H_LMMSE_R - H_LMMSE ).^2,'all');
        MSE_INTERPOLATION = sum(abs(H_Interpolation_R - X_Interpolation).^2, 'all');

        % BER
        % BER = BER + sum(sum(data_bin ~= data_bin_hat ));
    end
    %計算MSE數值
	MSE_dB_LS   (a)	= 10*log10( MSE_LS    / ( 822*560*Data_num) );
	MSE_dB_LMMSE(a)	= 10*log10( MSE_LMMSE / (1644*560*Data_num) );
    MSE_dB_INTERPOLATION(a) = 10*log10( MSE_INTERPOLATION / (1644*40*Data_num) );
    %BER_SNR(1,a) = BER/(1644 * 560 * Data_num * q_bit * Tx );
end

figure(1)
semilogy(SNR_in_dB,MSE_dB_LS(1,:),'r-x', 'LineWidth',2)
semilogy(SNR_in_dB,MSE_dB_LMMSE(1,:),'b-o', 'LineWidth',2)
semilogy(SNR_in_dB,MSE_dB_INTERPOLATION(1,:),'g-o', 'LineWidth',2)
hold on
grid on
axis square;
axis tight;
title('5G NR MIMO CHANNEL ESTIMATION')
xlabel('SNR (dB)')
ylabel('MSE (dB)')
legend('LS','LMMSE','Interpolation')
