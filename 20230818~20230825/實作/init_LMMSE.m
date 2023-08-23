function [LMMSE] = init_LMMSE(window,SNR_weight)
 %初始化通道估測
 N_fft  =2048;
 R_HD_HD =zeros( 822,822);
 R_H_HD =zeros(1644,822);
 DMRS_pos=[204:2:1024,1027:2:1847];
 Real_pos=[203:1:1024,1026:1:1847];
 %parfor p=1:822
 for p=1:822
  for k=1:822%R_HD_HD 822*822
   if DMRS_pos(k)==DMRS_pos(p)
    R_HD_HD  (k,p)=1;
   else
    R_HD_HD  (k,p)=( 1 - exp( -1i*2*pi*window*(DMRS_pos(k)-DMRS_pos(p))/N_fft ) )/( 1i*2*pi*window*(DMRS_pos(k)-DMRS_pos(p))/N_fft );
   end
  end
  for k=1:1644
   if Real_pos(k)==DMRS_pos(p)
    R_H_HD(k,p)=1;
   else
    R_H_HD(k,p)=( 1 - exp( -1i*2*pi*window*(Real_pos(k)-DMRS_pos(p))/N_fft ) )/( 1i*2*pi*window*(Real_pos(k)-DMRS_pos(p))/N_fft );
   end
  end
 end
 SNR_W  = 10^( SNR_weight /10);

 P0    = eye(822);
 P1    = ones(822,1);P1(1:2:822) = -1;
 P1    = diag(P1);
 LMMSE  = zeros(1644,822,2,2);
 LMMSE(:,:,1,1) = R_H_HD* P0' * inv( P0*R_HD_HD*P0' + P1*R_HD_HD*P1' + (1/SNR_W)*eye(822) );
 LMMSE(:,:,2,1) = R_H_HD* P0' * inv( P0*R_HD_HD*P0' + P1*R_HD_HD*P1' + (1/SNR_W)*eye(822) );
 LMMSE(:,:,1,2) = R_H_HD* P1' * inv( P0*R_HD_HD*P0' + P1*R_HD_HD*P1' + (1/SNR_W)*eye(822) );
 LMMSE(:,:,2,2) = R_H_HD* P1' * inv( P0*R_HD_HD*P0' + P1*R_HD_HD*P1' + (1/SNR_W)*eye(822) );
end