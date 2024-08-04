% LMMSE Window 數量不同比較
clear
clc

% Tx = 1;     % 傳送端個數
% Rx = 1;     % 接收端個數

% one_frame  = 1644*560;
frame_num  = 50;
SNR_in_dB  = 0:5:30;        % 自己設訂雜訊大小
SNR_weight = 40;
Window     = [4 6 8 14];
DMRS_DATA  = -0.7071 - 0.7071*1i;

QAM   = 16;
Eavg  = (qammod((0:QAM-1),QAM)*qammod((0:QAM-1),QAM)')/QAM;
NF    = 1/sqrt(Eavg);
q_bit = log2(QAM);      % 一個symbol可以傳幾個bit

MSE_dB_LS    = zeros(length(Window), length(SNR_in_dB));
MSE_dB_LMMSE = zeros(length(Window), length(SNR_in_dB));

for w = 1:length(Window)
    %LMMSE 權重矩陣
    W = make_W_S(Window(w),SNR_weight);

    for a      = 1:length(SNR_in_dB)
        Es     = 1;
        SNR    = 10^(SNR_in_dB(a)/10);
        No     = Es/SNR;
        MSE_LS    = 0;
        MSE_LMMSE = 0;
        for frame = 1:frame_num	
            fprintf("Window: %d/%d \t SNR : %d/%d \t frame : %d/%d\n", w, length(Window), a,length(SNR_in_dB),frame,frame_num);
            % Mod
            dec_data = randi([0 QAM-1], 12, 14);           % One Resource Block
            % bin_data = de2bi(dec_data, q_bit, 'left-msb');
            mod_data = NF*qammod(dec_data, QAM, 'gray');
            % Add DMRS
            mod_data(2:2:12, 3) = DMRS_DATA;
            % Extend Data
            mod_data = repmat(mod_data, 137, 40);           % One Frame, 137 Subcarriers, 40 OFDM Symbols
            % Guard Band
            DC = zeros(1, 560);
            X  = [zeros(202, 560); mod_data(1:822, :); DC; mod_data(823:end, :); zeros(201, 560)];
            % IFFT
            x = ifft(ifftshift(X,1))*sqrt(2048);
            % ADD Cyclic Prefix
            x_CP  = zeros(1, 1228800);
            index = 1;
            for symbol=1:560
                if mod(symbol, 28) - 1        % 為0的話會不成立走else
                    x_CP(1, index:index+2048+144-1) = [x(2048-144+1:2048, symbol); x(:, symbol)];
                    index = index + 2048 + 144;
                else
                    x_CP(1, index:index+2048+208-1) = [x(2048-208+1:2048, symbol); x(:, symbol)];
                    index = index + 2048 + 208;
                end
            end
            % Channel & Noise
            PowerdB   = [-2 -8 -10 -12 -15 -18];
            Ntap	  = length(PowerdB);
            H_Channel = sqrt(10.^(PowerdB/10));
            H_Channel = H_Channel .* (sqrt(1/2) * ( randn(1,Ntap) + 1i*randn(1,Ntap) ) );
            Total_H_Power = sum(10.^(PowerdB/10));      % 總通道能量 = 1
            pure_y = conv(x_CP, H_Channel);
            pure_y(:, 1228801:end) = [];                % 刪除最後5筆資料
            n = sqrt(No/2)*(randn(1, 1228800) + randn(1, 1228800)*1i);      % randn產生noise variance = No
            % 產生接收訊號
            y = pure_y + n;
            % Remove Cyclic Prefix
            y_rmCP = zeros(2048, 560);
            index = 1;
            for symbol = 1:560
                if mod(symbol, 28) - 1;
                    y_rmCP(:, symbol) = y(1, index+144:index+144+2048-1);
                    index = index + 144 + 2048;
                else
                    y_rmCP(:, symbol) = y(1, index+208:index+208+2048-1);
                    index = index + 208 + 2048;
                end
            end
            % FFT
            Y_fft = fftshift(fft(y_rmCP),1)/sqrt(2048);
            % Remove Guard Band
            Y = [Y_fft(203:1024, :); Y_fft(1026:1847, :)];
            %% 實際通道估測
            h = [H_Channel, zeros(1, 2042)];
            H = fftshift(fft(h, 2048));
            H_Data = [H(1, 203:1024),H(1, 1026:1847)].';        % Remove Guard Band
            H_frame = repmat(H_Data, 1, 560);                   % One Frame, 560 OFDM Symbols
            % 取得LS 估測結果，20 個Slot 都要估測
            X_LS 	= DMRS_DATA * eye(1644/2);
            Y_LS 	= Y(2:2:1644 , 3:14:560);
            H_LS 	= inv(X_LS)*Y_LS;
            % 取得LMMSE 估測結果
            H_LMMSE 	= W * H_LS;
            % 產生真實通道
            H_LS_R 		= H_frame(2:2:1644 , 3:14:560);
            H_LMMSE_R 	= H_frame(1:1:1644 , 3:14:560);
            % 累加
            MSE_LS		= MSE_LS    + sum( abs( H_LS_R    - H_LS    ).^2,'all');
            MSE_LMMSE	= MSE_LMMSE + sum( abs( H_LMMSE_R - H_LMMSE ).^2,'all');
        end
        MSE_dB_LS(w,a)    = 10*log10( MSE_LS    / (822  * 40 * frame_num) );
        MSE_dB_LMMSE(w,a) = 10*log10( MSE_LMMSE / (1644 * 40 * frame_num) );
    end
end

figure(1);
plot(SNR_in_dB, MSE_dB_LS(1,:), 'r-', 'LineWidth', 2);
hold on;
colors = ['g', 'b', 'm', 'k'];
for w=1:length(Window)
    plot(SNR_in_dB, MSE_dB_LMMSE(w,:), [colors(w) '-'], 'LineWidth', 2);
    hold on;
end
title('不同Window 數量大小比較');
grid on;
xlabel('SNR (dB)');
ylabel('MSE (dB)');
legend('LS MSE', 'LMMSE MSE (W=4)', 'LMMSE MSE (W=6)', 'LMMSE MSE (W=8)', 'LMMSE MSE (W=14)');
