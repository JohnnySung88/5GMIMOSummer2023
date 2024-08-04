% 5G MIMO Capacity
clear
clc

MIMO_range = 1:6;
one_frame  = 1644*560;
frame_num  = 200;
SNR_in_dB  = 0:5:40;        % 自己設訂雜訊大小
Capacity_in_SNR  = zeros(length(MIMO_range),length(SNR_in_dB));

QAM   = 16;
Eavg  = (qammod((0:QAM-1),QAM)*qammod((0:QAM-1),QAM)')/QAM;
NF    = 1/sqrt(Eavg);
q_bit = log2(QAM);      % 一個symbol可以傳幾個bit

for i = 1:length(MIMO_range)
    Tx = MIMO_range(i);
    Rx = Tx; 
    for a      = 1:length(SNR_in_dB)
        Es     = 1;
        SNR    = 10^(SNR_in_dB(a)/10);
        No     = Es/SNR;
        Capacity = 0; 
        for frame = 1:frame_num	
			fprintf("MIMO : %d/%d\t SNR : %d/%d \t frame : %d/%d\n",i,length(MIMO_range),a,length(SNR_in_dB),frame,frame_num);
            % Mod
            dec_data = randi([0 QAM-1], 1644, 560, Tx);
            bin_data = de2bi(dec_data, q_bit, 'left-msb');
            mod_data = NF*qammod(dec_data, QAM, 'gray');
            % Guard Band
            DC = zeros(1, 560, Tx);
            X  = [zeros(202, 560, Tx); mod_data(1:822, :, :); DC; mod_data(823:end, :, :); zeros(201, 560, Tx)];
            % IFFT
            x = ifft(ifftshift(X, 1))*sqrt(2048);
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
            H_Channel     = H_Channel .* ( sqrt( 1/(2*Tx) ) .* ( randn(Ntap,Rx,Tx) + 1i*randn(Ntap,Rx,Tx) ) );
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
            %% Capacity Calculation
            h = [H_Channel; zeros(2042,Rx,Tx)];
            H = fftshift(fft(h, [], 1),1);
            H_Data = [H(203:1024,:,:);H(1026:1847,:,:)]; ;
            H_frame = permute(repmat( H_Data(:,:,:),1,1,1,560),[1 4 2 3]);
            % Capacity
            for SC = 1:1644
				for slot = 1:560
					unit_H 	 = reshape(H_frame(SC,slot,:,:),Rx,Tx);
					Capacity = Capacity + abs( log2( det( eye(Tx) + (1/Tx) .* SNR .* (unit_H*unit_H'))));
				end
			end
        end
        Capacity_in_SNR(i,a) = Capacity/(one_frame * frame_num);
    end
end

figure(1)
hold on
for i=1:6
    semilogy(SNR_in_dB, Capacity_in_SNR(i,:), 'LineWidth',2)
end
grid on
axis tight
axis square
title('Capacity of MIMO OFDM')
xlabel('SNR (dB)')
ylabel('Capacity(bits/sec/Hz)')
legend('1x1','2X2','3x3','4x4','5x5','6x6')
