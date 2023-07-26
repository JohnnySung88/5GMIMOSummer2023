%生成PSS SSS PBCH資料
SSB_5G_NR;

%基礎設置
Tx = 1;
Rx = 1;
SNR_in_dB = 20;
frame_num = 1e4;
Blank_num = 150;

%4QAM init
QAM 	= 16;
Eavg 	= (qammod([0:QAM-1],QAM) * qammod([0:QAM-1],QAM)') / QAM;
NF 		= 1 / sqrt(Eavg);
q_bit 	= 4;        % 一個symbol可以傳幾個bit

%格式化 同步訊號
PBCH = [ zeros(702,1) ; qammod(randi([0,QAM-1],240,1),QAM,'gray')*NF ; zeros(702,1)];
PSS  = [ zeros(758,1) ; PSS 		; zeros(759,1) ];
SSS  = [ PBCH(1:750);zeros(8,1);SSS;zeros(9,1);PBCH(895:end) ];

%加上guard band
g_PBCH = [ zeros(202,1) ;PBCH(1:822) ;0 ;PBCH(823:end) ;zeros(201,1) ];
g_PSS  = [ zeros(202,1) ; PSS(1:822) ;0 ; PSS(823:end) ;zeros(201,1) ];
g_SSS  = [ zeros(202,1) ; SSS(1:822) ;0 ; SSS(823:end) ;zeros(201,1) ];


%同步訊號IFFT
syn_t  = ifft(ifftshift( g_PSS ))*sqrt(2048);

%長條圖參數
bar_x = -12:12;
bar_y = zeros(1,25);

%ans zone
frame_shift_ans = zeros(1,frame_num);

for frame=1:frame_num
	SNR = 10^( SNR_in_dB/10);
	No  = 10^(-SNR_in_dB/10);
	%生產數據
	data_dec	= randi([0,QAM-1],1644,14*4*10); 		% 隨機產生0~3 for 4QAM
	%data_bin 	= dec2bin(bin2gray(data_dec	 ,'qam'	,QAM 	),q_bit);   	% 將 0~3 轉為 '00'~'11'
	data_mod	= qammod(data_dec,QAM)*NF;       % 0~3 to complex (Modulation); remember to normalize
	%Synchronization Signal
	data_mod(:,5) = PSS;
	data_mod(:,6) = PBCH;
	data_mod(:,7) = SSS;
	data_mod(:,8) = PBCH;
	
	%Guard Band
	DC =   zeros(1,14*4*10);
	X  =   [ zeros(202,14*4*10) ;data_mod(1:822,:) ;DC ;data_mod(823:end,:) ;zeros(201,14*4*10) ];	
	%IFFT
	x  = ifft(ifftshift(X))*sqrt(2048);	%	2048 x 14*4*10	%%%%%%%										
	%CP
	x_CP = zeros(1,1228800);
	index=1;
	for symbol=1:14*4*10
		if	mod(symbol,28)-1
			x_CP(1, index:index+2048+144-1)=[ x(2048-144+1:2048,symbol) ; x(:,symbol)];
			index = index+2048+144;
		else
			x_CP(1, index:index+2048+208-1)=[ x(2048-208+1:2048,symbol) ; x(:,symbol)];
			index = index+2048+208;
		end
	end
	%通道與雜訊
	PowerdB 		= [ -2 -8 -10 -12 -15 -18];
	Total_H_Power 	= sum(10.^(PowerdB/10)); %總通道能量 = 1
	Ntap 			= 6;%通道數量
	H_Channel 		= sqrt(10.^(PowerdB/10));
	H_Channel   	= H_Channel .* (sqrt( 1/(2*Tx) ) * ( randn(1,Ntap) + 1i*randn(1,Ntap) ) );
	%訊號通過通道
	H_y					= conv( x_CP, H_Channel );
	H_y(:,1228801:end)  = [];			%刪除最後5筆資料
	
	%產生訊號
	n = sqrt(No/2) *( randn(1,1228800+2*Blank_num) + randn(1,1228800+2*Blank_num)*1i );% randn產生noise variance=No
	y = [zeros(1,Blank_num) ,H_y ,zeros(1,Blank_num)]+n;%noise off
	
	%接收後處理開始
	%同步偵測(Synchronization)
	syn_corr = xcorr(syn_t,y);
	[syn_max , syn_pos]	= max(syn_corr);
	[ram_a,ram_b]		= size(syn_corr);
	%紀錄目前偏移
	element_pos = 623378;
	frame_shift_ans(frame) = ceil(ram_a/2)- element_pos -syn_pos - Blank_num;
end

%畫圖
for i=1:25
	num = i-13;
	bar_y(i) = length(find(frame_shift_ans==num));
end
bar(bar_x,bar_y);
title('Synchronization 不同的偏移量次數')
xlabel('偏移量')
ylabel('偏移次數')
