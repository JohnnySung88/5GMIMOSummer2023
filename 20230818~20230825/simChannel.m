function [Rx_signal] = simChannel(Tx_signal,SNR_in_dB)
    Tx = 2;
    Rx = 2;
    PowerdB 	= [ -2 -8 -10 -12 -15 -18].';
    PowerdB_MIMO= repmat(PowerdB,1,Rx,Tx);
    H_Channel 	= sqrt(10.^(PowerdB_MIMO./10)).* sqrt( 1/ Tx );
    Ntap		= length(PowerdB);
    H_Channel   = H_Channel .* ( sqrt( 1/2 ) .* ( randn(Ntap,Rx,Tx) + 1i*randn(Ntap,Rx,Tx) ) );
    %捲積
    pure_y		= zeros(Rx,1228800+5);
    for Tx_n = 1:Tx
	    for Rx_n = 1:Rx
		    pure_y(Rx_n,:) = pure_y(Rx_n,:) + conv( Tx_signal(1,:,Tx_n) , H_Channel(:,Rx_n,Tx_n) );
	    end
    end
    %pure_y(:,1228801:end) = [];	%刪除最後5筆資料
    %產生訊號
    No          = 10^(-SNR_in_dB/10);
    Blank_num   = 50;
    n = sqrt(No/2) *( randn(2,1228800+2*Blank_num) + randn(2,1228800+2*Blank_num)*1i );
    Rx_signal	= [zeros(2,Blank_num) ,pure_y ,zeros(2,Blank_num-5)]+n;
end

