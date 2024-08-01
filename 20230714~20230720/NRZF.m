% 5G SISO system
clear
clc

% Tx = 1;     % 傳送端個數
% Rx = 1;     % 接收端個數

one_frame  = 1644*560;
SNR_in_dB  = 0:1:16;        % 自己設訂雜訊大小
BER_SNR_ZF = zeros(1,length(SNR_in_dB));

QAM   = 4;
Eavg  = (qammod((0:QAM-1),QAM)*qammod((0:QAM-1),QAM)')/QAM;
NF    = 1/sqrt(Eavg);
q_bit = log2(QAM);      % 一個symbol可以傳幾個bit

for a      = 1:length(SNR_in_dB)
    Es     = 1;
    SNR    = 10^(SNR_in_dB(a)/10);
    No     = Es/SNR;
    BER_ZF = 0; 
    % Mod
    dec_data = randi([0 QAM-1], 1644, 560);
    bin_data = de2bi(dec_data, q_bit, 'left-msb');
    mod_data = NF*qammod(dec_data, QAM, 'gray');
    % Guard Band
    DC = zeros(1, 560);
    X  = [zeros(202, 560); mod_data(1:822, :); DC; mod_data(823:end, :); zeros(201, 560)];
    % IFFT
    x = ifft(ifftshift(X))*sqrt(2048);
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
    PowerdB = [-2 -8 -10 -12 -15 -18];
    H_Channel = sqrt(10.^(PowerdB/10));
    Total_H_Power = sum(10.^(PowerdB/10));      % 總通道能量 = 1
    pure_y = conv(x_CP, H_Channel);
    pure_y(:, 1228801:end) = [];                % 刪除最後5筆資料
    n = sqrt(No/2)*(randn(1, 1228800) + randn(1, 1228800)*1i);      % randn產生noise variance = No
    % 產生接收訊號
    y = pure_y + n;
    % Remove Cyclic Prefix
    y_rmCP = zeros(2048, 560);
    index = 1;
    for symbol=1:560
        if mod(symbol, 28) - 1;
            y_rmCP(:, symbol) = y(1, index+144:index+144+2048-1);
            index = index + 144 + 2048;
            symbol = symbol + 1;
        else
            y_rmCP(:, symbol) = y(1, index+208:index+208+2048-1);
            index = index + 208 + 2048;
            symbol = symbol + 1;
        end
    end
    % FFT
    Y_fft = fftshift(fft(y_rmCP))/sqrt(2048);
    % Remove Guard Band
    Y = [Y_fft(203:1024, :); Y_fft(1026:1847, :)];
    % ZF Detector
    h = [H_Channel, zeros(1, 2042)];
    H = fftshift(fft(h, 2048));
    H_Data = [H(1, 203:1024),H(1, 1026:1847)].';
	H_frame = repmat(H_Data, 1, 560);
	X_hat_ZF = (Y ./ H_frame)/NF;
    % Demod
    dec_data_hat_ZF = qamdemod(X_hat_ZF, QAM, 'gray');
    bin_data_hat_ZF = de2bi(dec_data_hat_ZF, q_bit, 'left-msb');
    % BER Calculation
    BER_ZF = sum(bin_data_hat_ZF ~= bin_data, 'all');
    BER_SNR_ZF(1, a) = BER_ZF/(one_frame*q_bit);
end

%% 輸入SNR_in_dB
figure(1)
semilogy(SNR_in_dB,BER_SNR_ZF(1,:), 'r-X', 'LineWidth',2)
hold on
grid on
axis square
axis tight
title('BER of SISO')
xlabel('SNR (dB)')
ylabel('BER')
legend('4QAM ZF MATLAB')
