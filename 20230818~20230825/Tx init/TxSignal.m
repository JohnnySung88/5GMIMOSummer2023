load('LDPC_11nD2_1296b_R12.mat');
Pict_origin  = imread('Cat.jpg');
Pict_origin  = imresize(Pict_origin,[384,384]);
Pict_row	 = size(Pict_origin,1);
Pict_col	 = size(Pict_origin,2);
Pict_bin	 = [dec2bin(Pict_origin(:,:,1),8);
				dec2bin(Pict_origin(:,:,2),8);
				dec2bin(Pict_origin(:,:,3),8)];
JPG_bin      = Pict_bin-'0';
Pict_size	 = size(Pict_bin,1);
Pict_bin_RE	 = [Pict_bin;zeros(648-mod(Pict_row*Pict_col*3,648),8)];   %補0 for LDPC
Pict_bin_RE	 = reshape(Pict_bin_RE,[],648);
%%% random order %%%
[m, n] = size(Pict_bin_RE);
NumElements = m*n;
RandOrder = randperm(NumElements);
Pict_bin_RE_flatten = reshape(Pict_bin_RE, [1 NumElements]);
Pict_bin_RE_flatten_rand = Pict_bin_RE_flatten(RandOrder);
Pict_bin_RE_rand = reshape(Pict_bin_RE_flatten_rand,[m, n]);
%%%%%%%%%%%%%%%%%%%%
encode 		 = mod(Pict_bin_RE_rand * double(LDPC.G.x),2);

SNR = 10^( SNR_in_dB/10);
No  = 10^(-SNR_in_dB/10);
q_bit        = 4;
LDPC_num	 = size(encode,1);
LDPC_bin 	= reshape(encode,q_bit ,[]);%資料順序不論(準備進調變)
LDPC_dec_L	= sum(LDPC_bin.'.*binTable,2);%bin2dec
LDPC_mod_L	= qammod (LDPC_dec_L,  QAM)*NF;
LDPC_mod	= zeros(1644,560,Tx);
LDPC_mod(DATA_Pos) = [ LDPC_mod_L ; zeros( length(DATA_Pos)-length(LDPC_mod_L) ,1)];   %補0 for FrameData
LDPC_mod(DMRS1_Pos) = DMRS_DATA;
LDPC_mod(DMRS2_Pos) =-DMRS_DATA;
LDPC_mod(  PSS_pos) = [PSS;PSS]*40;

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

figure
plot(abs(x_CP(1,:,1)));
title("Data signal in Time domain",FontSize=14)
xlabel("Sample",FontSize=12)
ylabel("Voltage",FontSize=12)
