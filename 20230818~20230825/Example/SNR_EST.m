function [Rx_No] = SNR_EST(frame_data,H_LMMSE,Y_DMRS)
    DMRS  		= frame_data.DMRS;
	DMRS_H		= H_LMMSE (2:2:1644,:,:,:);
	DMRS        = DMRS(:,:,:,ones(1,2));
    DMRS        = permute(DMRS  ,[3,4,1,2]);
	DMRS_H      = permute(DMRS_H,[3,4,1,2]);
	DMRS_ram    = pagemtimes(DMRS_H,DMRS);
	DMRS_ram    = permute(DMRS_ram,[3,4,1,2]);
	DMRS_hat_Y  = DMRS_ram(:,:,:,1);
    Rx_No = mean( abs(Y_DMRS - DMRS_hat_Y).^2,[1,2] );
end

