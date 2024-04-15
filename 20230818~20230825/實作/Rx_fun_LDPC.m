function [Time,Biterror,Capacity_sum,Rx1_SNR,Rx2_SNR,JPG_RGB] = Rx_fun_LDPC(frame_data,Rx_signal,DTinfo,CFO_ignore,iteration)
    Fs  = frame_data.Fs;
    Tx  = frame_data.Tx;
    Rx  = frame_data.Rx;
    [X1,Y1]  = meshgrid(3:14:549,1:1644);
    [Xq,Yq]  = meshgrid(3: 1:549,1:1644);
    y_rmCP   = zeros(2048 , 560, Tx);
    H_INTER  = zeros(1644,560,2,2);    
    DMRS_DATA = +0.7071 + 0.7071*1i ;
    H_LMMSE  = zeros(1644,40,2,2);
    y  = Rx_signal;%2*1228800

    tic;
    % 同步
    y = Sync(frame_data,y);
    if(sum(y,'all')~=0)
        % CFO估測(需多重路徑猜測)
        [hat_delta_f,t] = CFO(frame_data,y,Fs,CFO_ignore);
        
        % 移除CP
        y_rmCP(:,:,1) = reshape(y(1,frame_data.CPdataPos),2048,560);
        y_rmCP(:,:,2) = reshape(y(2,frame_data.CPdataPos),2048,560);
        
        % CFO補償
        y_rmCP = y_rmCP ./ exp( 1i * 2 * pi * hat_delta_f * t);

        % FFT
        Y_fft  = fftshift( fft( y_rmCP/sqrt(2048) ) ,1);

        % rm Guard Band
        Y      = [ Y_fft( 203:1024,:,:) ; Y_fft( 1026:1847,:,:) ];

        % LMMSE估測
        Y_DMRS = Y(2:2:1644 , 3:14:560,:);
        H_LMMSE = pagemtimes(frame_data.LMMSE , DMRS_DATA' .*Y_DMRS);
        
        % 雜訊估測
        Rx_No = SNR_EST(frame_data,H_LMMSE,Y_DMRS);

        % 線性內差(LMMSE only) 
        H_INTER(:,3:549,1,1) = interp2(X1,Y1,H_LMMSE(:,:,1,1),Xq,Yq);
        H_INTER(:,3:549,1,2) = interp2(X1,Y1,H_LMMSE(:,:,1,2),Xq,Yq);
        H_INTER(:,3:549,2,1) = interp2(X1,Y1,H_LMMSE(:,:,2,1),Xq,Yq);
        H_INTER(:,3:549,2,2) = interp2(X1,Y1,H_LMMSE(:,:,2,2),Xq,Yq);

        % 邊界
        for symbol  = 549:560 
            head_dist  = symbol - 549;
            back_dist  = 563 - symbol;
            H_INTER(:,symbol,:,:) = ( back_dist * H_LMMSE(:,40,:,:) + head_dist * H_LMMSE(:,1,:,:)  )  /14;
        end
        for symbol  = 1:2     
            head_dist  = 11 + symbol;
            back_dist  = 3 - symbol;
            H_INTER(:,symbol,:,:) = ( back_dist * H_LMMSE(:,40,:,:) + head_dist * H_LMMSE(:,1,:,:)  )  /14;
        end

        % detector
        norm_Y = Y ./ sqrt(Rx_No);
        norm_H = H_INTER ./ sqrt(Rx_No);

        switch DTinfo
            case 'LMMSE'
                X_hat = LMMSE(norm_Y,norm_H,Tx);
            case 'ZF'
                X_hat = ZFD(norm_Y,norm_H,Tx);
        end
        X_hat = X_hat/frame_data.NF;

        % 反解資料
        LDPC_mod_L_hat = X_hat(frame_data.DATA_Pos(1:1770336));

        % LLR
        yI = real(LDPC_mod_L_hat).';
        yQ = imag(LDPC_mod_L_hat).';
        No = (Rx_No(1) + Rx_No(2))/2;
        LLR = zeros( 4*length(LDPC_mod_L_hat) ,1);
        LLR(1:4:length(LLR)) = LLRcul(yI,No, 1, 3,-1,-3);
        LLR(2:4:length(LLR)) = LLRcul(yI,No,-1, 1,-3, 3);
        LLR(3:4:length(LLR)) = LLRcul(yQ,No,-1,-3, 1, 3);
        LLR(4:4:length(LLR)) = (1/(2*No)) .* ( min( (yQ-[-1; 1]).^2 ) - min( (yQ-[-3; 3]).^2 ));

        LLR_part = reshape(LLR,[],1296);
        LLR_part = permute(LLR_part,[2,1]);

        LLR_MSA = LDPC_MSA(LLR_part,iteration,frame_data.H_row_master,frame_data.H_row_master_size,5464,0.5);
        JPG_bin_hat = reshape(LLR_MSA(1:648,:).'<0,[],8);

        % Time
        Time = toc;

        % LDPC decode(image)
        bin_table = 2 .^ (7:-1:0);
        JPG_row = frame_data.JPG_row;
        JPG_col = frame_data.JPG_col;
        JPG_size= frame_data.JPG_size;
        JPG_bin_hat = JPG_bin_hat(1:JPG_size,:);
        JPG_dec_hat = sum(JPG_bin_hat.*bin_table,2);%bin2dec

        JPG_Csize = JPG_row*JPG_col;
        JPG_RGB = zeros(JPG_row,JPG_col,3);
        JPG_RGB(:,:,1) = reshape( JPG_dec_hat(             1:JPG_Csize  ) ,JPG_row,JPG_col);
        JPG_RGB(:,:,2) = reshape( JPG_dec_hat(JPG_Csize  +1:JPG_Csize*2) ,JPG_row,JPG_col);
        JPG_RGB(:,:,3) = reshape( JPG_dec_hat(JPG_Csize*2+1:JPG_Csize*3) ,JPG_row,JPG_col);
        JPG_RGB = uint8(JPG_RGB);

        % BER
        Biterror = sum(JPG_bin_hat ~= frame_data.JPG_bin,'all');

        % SNR
        Rx1_SNR = -10*log10(Rx_No(1)/ mean(abs(Y(:,:,1).^2),'all'));
        Rx2_SNR = -10*log10(Rx_No(2)/ mean(abs(Y(:,:,2).^2),'all'));

        % Capacity
        SNR = 1/No;
        Capacity_sum = 0;
        for SC = 1:1644
            for slot = 1:560
                unit_H = reshape(H_INTER(SC,slot,:,:),Rx,Tx);
                Capacity_sum = Capacity_sum + abs( log2( det( eye(Tx) + (1/Tx) .* SNR .* (unit_H*unit_H'))));
            end
        end
    else
        Time = -1;
        JPG_RGB = 0;
        Rx1_SNR = 0;
        Rx2_SNR = 0;
        Biterror = 0;
        Capacity_sum = 0;
    end
end
