%MIMO system
clear;
clc;

Tx = 4; % number of transmit antennas
Rx = 4; % number of receive antennas
Data_num = 200;
SNR_in_dB = 0:5:40;
BER_SNR=zeros(1,length(SNR_in_dB));

QAM 	= 16;
Eavg 	= (qammod([0:QAM-1],QAM) * qammod([0:QAM-1],QAM)') / QAM;
NF 		= 1 / sqrt(Eavg);
q_bit 	= log2(QAM);        % 一個symbol可以傳幾個bit

for a=1:length(SNR_in_dB)
	SNR = 10^(SNR_in_dB(a)/10);
	No  = 1/SNR;
	BER = 0;
    for Data=1:Data_num
    fprintf("SNR : %d/%d \t Data number : %d/%d\n",a,length(SNR_in_dB),Data,Data_num);

        % TX
        data_dec	= randi([0,QAM-1],1644,560,Tx);
        data_bin = de2bi(data_dec, q_bit, 'left-msb');	
        data_mod = NF*qammod(data_dec, QAM, 'gray');
        
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
        H_Channel   = sqrt(10.^(PowerdB_MIMO./10));
        Ntap        = length(PowerdB);
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
        h = [H_Channel ; zeros(2042,Rx,Tx)];
        H = fftshift(fft(h,[],1),1);
        H_Data  = [H(203:1024,:,:);H(1026:1847,:,:)];   %rmGB
        H_frame = permute(repmat( H_Data(:,:,:),1,1,1,560),[1 4 2 3]  );
        X_hat = zeros(1644,560,Tx);
        for i=1:560
            for k=1:1644
                H_eq = reshape(H_frame(k,i,:,:),[Rx Tx]);
                X_hat(k,i,:) = reshape(Y(k,i,:),[Rx 1])'/(H_eq'*H_eq)*H_eq';  %ZF detection
            end
        end

        % De-Mod
        data_dec_hat = qamdemod(X_hat,QAM,'gray');
        data_bin_hat = de2bi(data_dec_hat, q_bit, 'left-msb');
        
        % BER
        BER = BER + sum(sum(data_bin ~= data_bin_hat ));
    end
    
	BER_SNR(1,a) = BER/(1644 * 560 * frame_num * q_bit * Tx );
end

figure(1)
semilogy(SNR_in_dB,BER_SNR(1,:),'r-x', 'LineWidth',2)
hold on
grid on
axis square;
axis tight;
title('BER of 4x4 MIMO with 16-QAM and ZF Detector')
xlabel('SNR (dB)')
ylabel('BER')
legend('16QAM ZF MATLAB')
