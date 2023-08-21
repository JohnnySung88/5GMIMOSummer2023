load('LDPC_11nD2_1296b_R12.mat');
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

bin_table	= 2 .^ (7:-1:0);
Lena_bin_hat_LDPC = reshape(encode(:,1:648),[],8);
Lena_bin_hat_LDPC = Lena_bin_hat_LDPC(1:Lena_size,:);
Lena_bin_hat_LDPC = sum(Lena_bin_hat_LDPC.*bin_table,2);%bin2dec
Lena_Csize		  = Lena_row*Lena_col;
Lena_RGB_LDPC	  = zeros(Lena_row,Lena_col,3);
Lena_RGB_LDPC(:,:,1) = reshape( Lena_bin_hat_LDPC(             1:Lena_Csize  ) ,Lena_row,Lena_col);
Lena_RGB_LDPC(:,:,2) = reshape( Lena_bin_hat_LDPC(Lena_Csize  +1:Lena_Csize*2) ,Lena_row,Lena_col);
Lena_RGB_LDPC(:,:,3) = reshape( Lena_bin_hat_LDPC(Lena_Csize*2+1:Lena_Csize*3) ,Lena_row,Lena_col);
imshow( uint8(Lena_RGB_LDPC) );