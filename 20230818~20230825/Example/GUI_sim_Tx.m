SNR = 10^( SNR_in_dB/10);
No  = 10^(-SNR_in_dB/10);

%LDPC數據生成  totle 1772232(扣除DMRS、PSS) 5469*1296+1104
%load('LDPC_11nD2_1296b_R12.mat');

Lena_origin  = imread('ngc6543a.jpg');
Lena_origin  = imresize(Lena_origin,[384,384]);
Lena_row	 = size(Lena_origin,1);
Lena_col	 = size(Lena_origin,2);
Lena_bin	 = [dec2bin(Lena_origin(:,:,1));
				dec2bin(Lena_origin(:,:,2));
				dec2bin(Lena_origin(:,:,3))];
Lena_size	 = size(Lena_bin,1);
Lena_bin_RE	 = [Lena_bin;zeros(648-mod(Lena_row*Lena_col*3,648),8)];
Lena_bin_RE	 = reshape(Lena_bin_RE,[],648);
encode 		 = mod(Lena_bin_RE * double(LDPC.G.x),2);
LDPC_num	 = size(encode,1);

LDPC_bin 	= reshape(encode,q_bit ,[]);%資料順序不論(準備進調變)
LDPC_dec_L	= sum(LDPC_bin.'.*binTable,2);%bin2dec
LDPC_mod_L	= QAMMOD (LDPC_dec_L,  QAM)*NF;
LDPC_mod	= zeros(1644,560,Tx);
LDPC_mod(DATA_Pos) = [ LDPC_mod_L ; zeros( length(DATA_Pos)-length(LDPC_mod_L) ,1)];
LDPC_mod(DMRS1_Pos) = DMRS_DATA;
LDPC_mod(DMRS2_Pos) =-DMRS_DATA;
LDPC_mod(  PSS_pos) = [PSS;PSS];

%Guard Band
DC =   zeros(   1,560,Tx);
X  = [ zeros( 202,560,Tx) ;LDPC_mod(1:822,:,:) ;DC ;LDPC_mod(823:end,:,:) ;zeros(201,560,Tx) ];	
%IFFT
x  = ifft(ifftshift(X,1))*sqrt(2048);	%	2048 x 14*4*10						
%CP
x_CP = zeros(1,1228800,Tx);
index=1;
for symbol=1:14*4*10
	if	mod(symbol,28)-1
		x_CP(1, index:index+2048+144-1,:)=[ x(2048-144+1:2048,symbol,:) ; x(:,symbol,:)];
		index = index+2048+144;
	else
		x_CP(1, index:index+2048+208-1,:)=[ x(2048-208+1:2048,symbol,:) ; x(:,symbol,:)];
		index = index+2048+208;
	end
end
%通道
PowerdB 	= [ -2 -8 -10 -12 -15 -18].';
PowerdB_MIMO= repmat(PowerdB,1,Rx,Tx);
H_Channel 	= sqrt(10.^(PowerdB_MIMO./10)).* sqrt( 1/ Tx );
Ntap		= length(PowerdB);
H_Channel   = H_Channel .* ( sqrt( 1/2 ) .* ( randn(Ntap,Rx,Tx) + 1i*randn(Ntap,Rx,Tx) ) );
%捲積
pure_y		= zeros(Rx,1228800+5);
for Tx_n = 1:Tx
	for Rx_n = 1:Rx
		pure_y(Rx_n,:) = pure_y(Rx_n,:) + conv( x_CP(1,:,Tx_n) , H_Channel(:,Rx_n,Tx_n) );
	end
end
pure_y(:,1228801:end) = [];	%刪除最後5筆資料
%產生訊號
n 			= sqrt(No/2) *( randn(Tx,1228800) + randn(Tx,1228800)*1i );% randn產生noise variance=No
Rx_signal	= pure_y + n;