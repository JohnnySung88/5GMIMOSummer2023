SNR = 10^( SNR_in_dB/10);
No  = 10^(-SNR_in_dB/10);

%LDPC數據生成  totle 1772232(扣除DMRS、PSS) 5469*1296+1104
%load('LDPC_11nD2_1296b_R12.mat');

figure_origin = imread("Cat.jpg");
figure_origin = imresize(figure_origin, [384,384]);
figure_row	  = size(figure_origin,1);
figure_col	  = size(figure_origin,2);
figure_bin1	  = dec2bin(figure_origin(:,:,1));
figure_bin2	  = dec2bin(figure_origin(:,:,2));
figure_bin3	  = dec2bin(figure_origin(:,:,3));
figure_bin	  = [figure_bin1,figure_bin2,figure_bin3];
figure_size	 = size(figure_bin,1);
figure_bin_RE	 = [figure_bin;zeros(648-mod(figure_row*figure_col*3,648),8)];
figure_bin_RE	 = reshape(figure_bin_RE,[],648);

encode 		 = mod(figure_bin_RE * double(LDPC.G.x),2);
q_bit        = 4;
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

