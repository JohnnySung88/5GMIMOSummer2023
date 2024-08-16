% 5G 2x2 MIMO Channel Estimation
clear
clc

Tx = 2;     % 傳送端個數
Rx = 2;     % 接收端個數

one_frame  = 1644*560;
frame_num  = 200;
SNR_in_dB  = 0:5:40;        % 自己設訂雜訊大小
SNR_weight = 45;
window     = 10;
DMRS_DATA  = +0.7071 + 0.7071*1i;

QAM   = 16;
Eavg  = (qammod((0:QAM-1),QAM)*qammod((0:QAM-1),QAM)')/QAM;
NF    = 1/sqrt(Eavg);
q_bit = log2(QAM);      % 一個symbol可以傳幾個bit

MSE_dB_LS    = zeros(1,length(SNR_in_dB));
MSE_dB_LMMSE = zeros(1,length(SNR_in_dB));
MSE_dB_INTER = zeros(1,length(SNR_in_dB));

% LMMSE Init
N_fft 	=2048;
R_HD_HD	=zeros( 822,822);
R_H_HD	=zeros(1644,822);
DMRS_pos=[204:2:1024,1027:2:1847];
Real_pos=[203:1:1024,1026:1:1847];
for p=1:822
	for k=1:822%R_HD_HD 822*822
		if DMRS_pos(k)==DMRS_pos(p)
			R_HD_HD  (k,p)=1;
		else
			R_HD_HD  (k,p)=( 1 - exp( -1i*2*pi*window*(DMRS_pos(k)-DMRS_pos(p))/N_fft ) )/( 1i*2*pi*window*(DMRS_pos(k)-DMRS_pos(p))/N_fft );
		end
	end
	for k=1:1644
		if Real_pos(k)==DMRS_pos(p)
			R_H_HD(k,p)=1;
		else
			R_H_HD(k,p)=( 1 - exp( -1i*2*pi*window*(Real_pos(k)-DMRS_pos(p))/N_fft ) )/( 1i*2*pi*window*(Real_pos(k)-DMRS_pos(p))/N_fft );
		end
	end
end
SNR_W 	= 10^( SNR_weight /10);

for a      = 1:length(SNR_in_dB)
    Es     = 1;
    SNR    = 10^(SNR_in_dB(a)/10);
    No     = Es/SNR;
    MSE_LS    = 0;
    MSE_LMMSE = 0;
    MSE_INTER = 0; 
	for frame = 1:frame_num	
		fprintf("SNR : %d/%d \t frame : %d/%d\n",a,length(SNR_in_dB),frame,frame_num);
        % Mod
        dec_data = randi([0 QAM-1], 12, 14, Tx);
        % bin_data = de2bi(dec_data, q_bit, 'left-msb');
        mod_data = NF*qammod(dec_data, QAM, 'gray');
        % 安置DMRS
		mod_data(2:2:12,3,:) = DMRS_DATA;
		mod_data = repmat(mod_data,137,40,1);
		% CDM
		mod_data(2:4:1644,3:14:560,2:2:Tx) = -mod_data(2:4:1644,3:14:560,2:2:Tx);
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
        %% 實際通道估測
        h = [H_Channel; zeros(2042,Rx,Tx)];
        H = fftshift(fft(h, [], 1),1);
        H_Data  = [H(203:1024,:,:);H(1026:1847,:,:)];
        H_frame = permute(repmat( H_Data(:,:,:),1,1,1,560),[1 4 2 3]  );
        % 取得LS 估測結果，20 個Slot 都要估測 (For 2x2 MIMO Only)
		Y_DMRS 	= Y(2:2:1644,3:14:560,:);
		H_LS	= zeros(822,40,2,2);
		H_LS(1:2:822,:,1,1) 	= DMRS_DATA' * (  Y_DMRS(1:2:822,:,1) + Y_DMRS(2:2:822,:,1) )/2;
		H_LS(1:2:822,:,1,2) 	= DMRS_DATA' * ( -Y_DMRS(1:2:822,:,1) + Y_DMRS(2:2:822,:,1) )/2;
		H_LS(1:2:822,:,2,1) 	= DMRS_DATA' * (  Y_DMRS(1:2:822,:,2) + Y_DMRS(2:2:822,:,2) )/2;
		H_LS(1:2:822,:,2,2) 	= DMRS_DATA' * ( -Y_DMRS(1:2:822,:,2) + Y_DMRS(2:2:822,:,2) )/2;
		H_LS(2:2:822,:,:,:) 	= H_LS(1:2:822,:,:,:);
		% 取得LMMSE 估測結果
		P0 				 = eye(822);
		P1 				 = ones(822,1);P1(1:2:822) = -1;
		P1 				 = diag(P1);
		H_LMMSE			 = zeros(1644,40,2,2);
		H_LMMSE(:,:,1,1) = R_H_HD* P0' * inv( P0*R_HD_HD*P0' + P1*R_HD_HD*P1' + (1/SNR_W)*eye(822) ) * (DMRS_DATA' .*Y_DMRS(:,:,1));
		H_LMMSE(:,:,2,1) = R_H_HD* P0' * inv( P0*R_HD_HD*P0' + P1*R_HD_HD*P1' + (1/SNR_W)*eye(822) ) * (DMRS_DATA' .*Y_DMRS(:,:,2));
		H_LMMSE(:,:,1,2) = R_H_HD* P1' * inv( P0*R_HD_HD*P0' + P1*R_HD_HD*P1' + (1/SNR_W)*eye(822) ) * (DMRS_DATA' .*Y_DMRS(:,:,1));
		H_LMMSE(:,:,2,2) = R_H_HD* P1' * inv( P0*R_HD_HD*P0' + P1*R_HD_HD*P1' + (1/SNR_W)*eye(822) ) * (DMRS_DATA' .*Y_DMRS(:,:,2));
		% Interpolation 線性內差
		DMRS_Spos	= 3:14:560;
		H_INTER		= zeros(1644,560,2,2);	
		for symbol  = 549:560   %邊界
			head_dist	= symbol - 549;
			back_dist	= 563    - symbol;
			H_INTER(:,symbol,:,:) = ( back_dist * H_LMMSE(:,40,:,:) + head_dist * H_LMMSE(:,1,:,:)  )  /14;
		end
		for symbol  = 1:2       %邊界
			head_dist	= 11 + symbol;
			back_dist	= 3  - symbol;
			H_INTER(:,symbol,:,:) = ( back_dist * H_LMMSE(:,40,:,:) + head_dist * H_LMMSE(:,1,:,:)  )  /14;
		end
		for symbol 	= 3:548     %連續
			pos  		= floor((symbol-3)/14) + 1;
			head_dist	= symbol 			- DMRS_Spos(pos);
			back_dist	= DMRS_Spos(pos+1) 	- symbol;
			H_INTER(:,symbol,:,:) = ( back_dist * H_LMMSE(:,pos,:,:) + head_dist * H_LMMSE(:,pos+1,:,:)  )  /14;
		end
		% 產生真實通道
		H_LS_R 		= H_frame(2:2:1644 , 3:14:560,:,:);
		H_LMMSE_R 	= H_frame(1:1:1644 , 3:14:560,:,:);
		% 累加
		MSE_LS		= MSE_LS    + sum( abs( H_LS_R    - H_LS    ).^2,'all');
		MSE_LMMSE	= MSE_LMMSE + sum( abs( H_LMMSE_R - H_LMMSE ).^2,'all');
		MSE_INTER	= MSE_INTER + sum( abs( H_frame   - H_INTER ).^2,'all');
	end
	MSE_dB_LS   (a)	= 10*log10( MSE_LS    / ( 822* 40*Tx*Rx*frame_num) );
	MSE_dB_LMMSE(a)	= 10*log10( MSE_LMMSE / (1644* 40*Tx*Rx*frame_num) );
	MSE_dB_INTER(a)	= 10*log10( MSE_INTER / (1644*560*Tx*Rx*frame_num) );
end

figure(1)
plot(SNR_in_dB,MSE_dB_LS   ,'r-','LineWidth',2);
hold on;
plot(SNR_in_dB,MSE_dB_LMMSE,'b-','LineWidth',2);
hold on;
plot(SNR_in_dB,MSE_dB_INTER,'g-','LineWidth',2);
hold on;
grid on;
axis square;
axis tight;
title('LS 與 LMMSE 與Interpolation 的MSE 比較');
xlabel('SNR (dB)');
ylabel('SNR (dB) from MSE');
legend('LS MSE','LMMSE MSE','Interpolation MSE');
