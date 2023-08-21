figure_origin = imread("Figure.png");
figure_reshape = imresize(figure_origin, [270, 271]);
imshow(figure_reshape,'Parent',app.org_pit)

figure_pirolbin1 = dec2bin(figure_reshape(:,:,1));
figure_pirolbin2 = dec2bin(figure_reshape(:,:,2));
figure_pirolbin3 = dec2bin(figure_reshape(:,:,3));
figure_pirolbin  = [figure_pirolbin1,figure_pirolbin2,figure_pirolbin3];
figure_pirolbit  = [figure_pirolbin;zeros(648-mod(270*271*3,648),8)];
figure_pirolbit  = reshape(figure_pirolbit',[],648,);

figure_decode_bit_ldpc = reshape(figure_decode_bit_ldpc',[],8);
figure_decode_bit_ldpc = double(figure_decode_bit_ldpc(1:figure_size,:));

figure_decode_bin_ldpc = bin2dec(num2str(figure_decode_bit_ldpc));
figure_red_ldpc        = figure_decode_bin_ldpc(1:73170);
pika_green_Idpc        = figure_decode_bin_ldpc(73170+1:73170*2); 
pika_blue_Idpc         = figure_decode_bin_Idpc(73170*2+1:73170*3);

figure_encoded_ALL(:,:,1) = reshape (figure_red_1dpc, 270,271); 
figure_encoded_ALL(:,:,2) = reshape (figure_green_Idpc, 270,271); 
figure_encoded_ALL(:,:,3) = reshape (figure_blue_Idpc, 270,271);
figure_encoded_bin = uint8(figure_encoded_ALL);
imshow (figure_encoded_bin, 'Parent' ,app.LDPC_Pit);