load('LDPC_11nD2_1296b_R12.mat');
Pict_origin  = imread('Cat.jpg');
Pict_origin  = imresize(Pict_origin,[384,384]);
Pict_row	 = size(Pict_origin,1);
Pict_col	 = size(Pict_origin,2);
Pict_bin	 = [dec2bin(Pict_origin(:,:,1));
				dec2bin(Pict_origin(:,:,2));
				dec2bin(Pict_origin(:,:,3))];
Pict_size	 = size(Pict_bin,1);
Pict_bin_RE	 = [Pict_bin;zeros(648-mod(Pict_row*Pict_col*3,648),8)];
Pict_bin_RE	 = reshape(Pict_bin_RE,[],648);
encode 		 = mod(Pict_bin_RE * double(LDPC.G.x),2);

bin_table	= 2 .^ (7:-1:0);
Pict_bin_hat_LDPC = reshape(encode(:,1:648),[],8);
Pict_bin_hat_LDPC = Pict_bin_hat_LDPC(1:Pict_size,:);
Pict_bin_hat_LDPC = sum(Pict_bin_hat_LDPC.*bin_table,2);%bin2dec
Pict_Csize		  = Pict_row*Pict_col;
Pict_RGB_LDPC	  = zeros(Pict_row,Pict_col,3);
Pict_RGB_LDPC(:,:,1) = reshape( Pict_bin_hat_LDPC(             1:Pict_Csize  ) ,Pict_row,Pict_col);
Pict_RGB_LDPC(:,:,2) = reshape( Pict_bin_hat_LDPC(Pict_Csize  +1:Pict_Csize*2) ,Pict_row,Pict_col);
Pict_RGB_LDPC(:,:,3) = reshape( Pict_bin_hat_LDPC(Pict_Csize*2+1:Pict_Csize*3) ,Pict_row,Pict_col);
imshow( uint8(Pict_RGB_LDPC) );