% 5G Synchronization
clear
clc

% Default Settings
Tx = 1;
Rx = 1;
SNR_in_dB = 20;
frame_num = 1e4;
Blank_num = 150;
SSB_5G_NR;      %生成PSS SSS PBCH資料

% 16QAM init
QAM = 16;
Eavg  = (qammod((0:QAM-1),QAM)*qammod((0:QAM-1),QAM)')/QAM;
NF    = 1/sqrt(Eavg);
q_bit = log2(QAM);      % 一個symbol可以傳幾個bit

% Format Sync Signal
PBCH = [zeros(702, 1); qammod(randi([0, QAM-1], 240, 1), QAM, 'gray')*NF; zeros(702, 1)];
PSS  = [zeros(758, 1); PSS; zeros(759, 1)];
SSS  = [PBCH(1:750); zeros(8, 1); SSS; zeros(9, 1); PBCH(895:end)];

% Add Guard Band
g_PBCH = [zeros(202, 1); PBCH(1:822); 0; PBCH(823:end); zeros(201, 1)];
g_PSS  = [zeros(202, 1);  PSS(1:822); 0;  PSS(823:end); zeros(201, 1)];
g_SSS  = [zeros(202, 1);  SSS(1:822); 0;  SSS(823:end); zeros(201, 1)];

% Sync Signal IFFT
syn_t  = ifft(ifftshift(g_PSS))*sqrt(2048);

% Init Frame Shift
frame_shift_ans = zeros(1, frame_num);

for frame  = 1:frame_num
    Es     = 1;
    SNR    = 10^(SNR_in_dB/10);
    No     = Es/SNR;

    % Mod
    data     = randi([0 QAM-1], 1644, 560);
    % bin_data = de2bi(data, q_bit, 'left-msb');
    mod_data = NF*qammod(data, QAM, 'gray');

    % Sync Signal
    mod_data(:, 5) = PSS;
    mod_data(:, 6) = PBCH;
    mod_data(:, 7) = SSS;
    mod_data(:, 8) = PBCH;

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
    PowerdB 		= [ -2 -8 -10 -12 -15 -18];
	Ntap 			= 6;                            %通道數量
	H_Channel 		= sqrt(10.^(PowerdB/10));
	H_Channel   	= H_Channel .* (sqrt( 1/(2*Tx) ) * ( randn(1,Ntap) + 1i*randn(1,Ntap) ) );
    Total_H_Power 	= sum(10.^(PowerdB/10));        %總通道能量 = 1

	% 訊號通過通道
	pure_y				  = conv( x_CP, H_Channel );
	pure_y(:,1228801:end) = [];     %刪除最後5筆資料

    % 產生接收訊號
    % n = sqrt(No/2)*(randn(1, 1228800) + randn(1, 1228800)*1i);      % randn產生noise variance = No
    % y = pure_y + n;
	n = sqrt(No/2)*(randn(1, 1228800 + 2*Blank_num) + randn(1, 1228800 + 2*Blank_num)*1i);      % randn產生noise variance = No
	y = [zeros(1, Blank_num), pure_y, zeros(1, Blank_num)] + n;     %noise off

	% 同步偵測(Synchronization)
	syn_corr            = xcorr(syn_t,y);
	[syn_max, syn_pos]	= max(syn_corr);
	[ram_a, ram_b]		= size(syn_corr);

	% 紀錄目前偏移
	element_pos = 623376;
	frame_shift_ans(frame) = ceil(ram_a/2) - element_pos - syn_pos - Blank_num - 1;
end

% Bar Chart Parameters
bar_x = -12:12;
bar_y = zeros(1, 25);
for i   = 1:25
	num = i-13;
	bar_y(i) = length(find(frame_shift_ans==num));
end

% Plot Figure
bar(bar_x,bar_y);
title('Synchronization 不同的偏移量次數')
xlabel('偏移量')
ylabel('偏移次數')
