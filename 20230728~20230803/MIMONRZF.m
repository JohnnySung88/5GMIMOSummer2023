% 5G 4x4 MIMO system
clear
clc

Tx = 4;     % 傳送端個數
Rx = 4;     % 接收端個數

one_frame  = 1644*560;
frame_num  = 1000;
SNR_in_dB  = 0:5:40;        % 自己設訂雜訊大小
BER_SNR_ZF = zeros(1,length(SNR_in_dB));

QAM   = 16;
Eavg  = (qammod((0:QAM-1),QAM)*qammod((0:QAM-1),QAM)')/QAM;
NF    = 1/sqrt(Eavg);
q_bit = log2(QAM);      % 一個symbol可以傳幾個bit

for a      = 1:length(SNR_in_dB)
    Es     = 1;
    SNR    = 10^(SNR_in_dB(a)/10);
    No     = Es/SNR;
    BER_ZF = 0; 
	for frame = 1:frame_num	
		fprintf("SNR : %d/%d \t frame : %d/%d\n",a,length(SNR_in_dB),frame,frame_num);
        % Mod
        dec_data = randi([0 QAM-1], 1644, 560, Tx);
        bin_data = de2bi(dec_data, q_bit, 'left-msb');
        mod_data = NF*qammod(dec_data, QAM, 'gray');
        % Guard Band
        DC = zeros(1, 560, Tx);
        X  = [zeros(202, 560, Tx); mod_data(1:822, :, :); DC; mod_data(823:end, :, :); zeros(201, 560, Tx)];
        % IFFT
        x = ifft(ifftshift(X))*sqrt(2048);
        % ADD Cyclic Prefix
        x_CP  = zeros(1, 1228800, Tx);
        index = 1;
        for symbol=1:560
            if mod(symbol, 28) - 1        % 為0的話會不成立走else
                x_CP(1, index:index+2048+144-1, :) = [x(2048-144+1:2048, symbol, :); x(:, symbol, :)];
                index = index + 2048 + 144;
            else
                x_CP(1, index:index+2048+208-1, :) = [x(2048-208+1:2048, symbol, :); x(:, symbol, :)];
                index = index + 2048 + 208;
            end
        end
        % Channel & Noise
        PowerdB       = [-2 -8 -10 -12 -15 -18].';
        PowerdB_MIMO  = repmat(PowerdB,1,Rx,Tx);
        H_Channel     = sqrt(10.^(PowerdB_MIMO./10));
        Ntap	      = length(PowerdB);
        H_Channel     = H_Channel .* (sqrt(1/2*Tx) * ( randn(1,Ntap,Rx,Tx) + 1i*randn(1,Ntap,Rx,Tx) ) );
        Total_H_Power = sum(10.^(PowerdB/10));      % 總通道能量 = 1
        % Convolution
        pure_y = zeros(Rx,1228800+5);
        for Tx_n = 1:Tx
            for Rx_n = 1:Rx
                pure_y(Rx_n,:) = pure_y(Rx_n,:) + conv( x_CP(1,:,Tx_n) , H_Channel(:,Rx_n,Tx_n) );
            end
        end
        pure_y(:,1228801:end) = []; %刪除最後5筆資料
        n = sqrt(No/(2)) *( randn(Tx,1228800) + randn(Tx,1228800)*1i );% randn產生noise variance=No
        % 產生接收訊號
        y = pure_y + n;
        % Remove Cyclic Prefix
        y_rmCP = zeros(2048, 560, Tx);
        index = 1;
        for symbol = 1:560
            if mod(symbol, 28) - 1;
                for ram = 1:Rx
                    y_rmCP(:,symbol,ram) = y(ram,index+144:index+144+2048-1);
                end                
                index = index + 144 + 2048;
            else
                for ram = 1:Rx
                    y_rmCP(:,symbol,ram) = y(ram,index+208:index+208+2048-1);
                end                
                index = index + 208 + 2048;
            end
        end
        % FFT
        Y_fft = fftshift(fft(y_rmCP/sqrt(2048)),1);
        % Remove Guard Band
        Y = [Y_fft(203:1024, :, :); Y_fft(1026:1847, :, :)];
        % ZF Detector
        h = [H_Channel, zeros(2042,Rx,Tx)];
        H = fftshift(fft(h, [], 1),1);
        H_Data = [H(203:1024,:,:);H(1026:1847,:,:)]; ;
        H_frame = permute(repmat( H_Data(:,:,:),1,1,1,560),[1 4 2 3]  );
        % X_hat_ZF = (Y ./ H_frame)/NF;
        X_hat_ZF = zeros(1644, 560, Tx);
        X_hat_ZF = ZFD_4x4(Y, H_frame, Tx);
        X_hat_ZF = X_hat_ZF/NF;
        % Demod
        dec_data_hat_ZF = qamdemod(X_hat_ZF, QAM, 'gray');
        bin_data_hat_ZF = de2bi(dec_data_hat_ZF, q_bit, 'left-msb');
        % BER Calculation
        BER_ZF = BER_ZF + sum(bin_data_hat_ZF ~= bin_data, 'all');
    end
	BER_SNR_ZF(1,a) = BER_ZF/(one_frame * frame_num * q_bit * Tx );
end

figure(1)
semilogy(SNR_in_dB,BER_SNR_ZF(1,:), 'r-X', 'LineWidth',2)
hold on
grid on
axis square
axis tight
title('BER of 4x4 MIMO')
xlabel('SNR (dB)')
ylabel('BER')
legend('16QAM ZF C')
