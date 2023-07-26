% SISO system with synchronization
clear
clc
warning('off')

Data_number = 10000;  % 資料筆數
SNR_in_dB = 20;       % SNR設定為20dB
Channel_dB = [-2, -8, -10, -12, -15, -18];
Ntap = length(Channel_dB);
N_ID_Cell = 86;       % 同步訊號的長度
sync_offset_count = zeros(1, 21);  % 偏移量範圍：-10~10

for a = 1:Data_number
    SNR = 10^(SNR_in_dB/10);
    No = 10^(-SNR_in_dB/10);
    
    %mod
    data_dec = randi([0, 3], 1644, 560);  % 隨機產生0~3 for 4QAM
    data_bin = de2bi(data_dec, 2, 'left-msb');  % 將0~3轉換為'00'~'11'
    data_mod = qammod(data_dec, 4, 'gray') * sqrt(1/2);  % 0~3 to complex (Modulation); 請記得進行正規化

    %Guard Band
    DC = zeros(1, 560);
    X = [zeros(202, 560); data_mod(1:822, :); DC; data_mod(823:end, :); zeros(201, 560)];

    %IFFT
    x = ifft(ifftshift(X)) * sqrt(2048);

    %CP
    x_CP = zeros(1, 1228800);
    index = 1;
    for symbol = 1:560
        if mod(symbol, 28) - 1
            x_CP(1, index:index+2048+144-1) = [x(2048-144+1:2048, symbol); x(:, symbol)];
            index = index + 2048 + 144;
        else
            x_CP(1, index:index+2048+208-1) = [x(2048-208+1:2048, symbol); x(:, symbol)];
            index = index + 2048 + 208;
        end
    end

    %通道與雜訊
    PowerdB = Channel_dB;
    H_Channel = sqrt(10.^(PowerdB/10));
    h = H_Channel .* sqrt((1/(2*Ntap))*(randn(1, Ntap) + 1i*randn(1, Ntap)));
    
    pure_y = conv(x_CP, h);
    pure_y(:, 1228801:end) = [];  % 刪除最後5筆資料
    
    n = sqrt(No/2) * (randn(1, 1228800) + 1i*randn(1, 1228800));  % randn產生雜訊，其變異數為No
    y = pure_y + n;  % 產生最終傳輸信號，考慮通道和雜訊
    
    % 加入同步訊號和空白資料
    data_sync = [zeros(1, N_ID_Cell), data_mod(1, :), zeros(1, 100)];  % 假設同步訊號長度為 N_ID_Cell，前後加上 100 個空白資料
    rx_signal = [zeros(1, N_ID_Cell), y, zeros(1, 100)];  % 前後加上 100 個空白資料
    
    % 將 data_sync 長度補齊至 rx_signal 長度
    data_sync_padded = [data_sync, zeros(1, length(rx_signal) - length(data_sync))];
    
    % 找到同步偏移量
    correlation_result = xcorr(rx_signal, fliplr(data_sync_padded));
    [~, max_index] = max(abs(correlation_result));  % Find the index of the peak
    sync_offset = max_index - (length(rx_signal) + N_ID_Cell - 1);  % 轉換為同步偏移量
    
    % 檢查同步偏移量是否在範圍內
    if sync_offset >= -10 && sync_offset <= 10
        sync_offset_count(sync_offset + 11) = sync_offset_count(sync_offset + 11) + 1;
    else
        % Debugging: print the sync_offset if it's outside the range
        fprintf('Trial %d: sync_offset = %d\n', a, sync_offset);
    end
end

% 輸入同步偏移量次數並繪製柱狀圖
figure(1)
bar(-10:10, sync_offset_count)
grid on
axis square
axis tight
title('Sync Offset Count of SISO with Synchronization')
xlabel('Sync Offset')
ylabel('Sync Offset Count')
