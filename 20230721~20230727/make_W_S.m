%LMMSE 權重矩陣
function W = make_W(window_loc,SNR_weight_loc)
	N_fft 	=2048;
	R_H_H	=zeros( 822,822);
	R_H_HLS	=zeros(1644,822);
	DMRS_pos=[204:2:1024,1027:2:1847];
	Real_pos=[203:1:1024,1026:1:1847];
	for p=1:822
		for k=1:822%R_H_H 822*822
			if DMRS_pos(k)==DMRS_pos(p)
				R_H_H  (k,p)=1;
			else
				R_H_H  (k,p)=( 1 - exp( -1i*2*pi*window_loc*(DMRS_pos(k)-DMRS_pos(p))/N_fft ) )/( 1i*2*pi*window_loc*(DMRS_pos(k)-DMRS_pos(p))/N_fft );
			end
		end
		for k=1:1644
			if Real_pos(k)==DMRS_pos(p)
				R_H_HLS(k,p)=1;
			else
				R_H_HLS(k,p)=( 1 - exp( -1i*2*pi*window_loc*(Real_pos(k)-DMRS_pos(p))/N_fft ) )/( 1i*2*pi*window_loc*(Real_pos(k)-DMRS_pos(p))/N_fft );
			end
		end
	end
	SNR_W 	= 10^( SNR_weight_loc /10);
	W		= R_H_HLS * inv(R_H_H + (1/SNR_W) * eye(822) );%權重矩陣
end