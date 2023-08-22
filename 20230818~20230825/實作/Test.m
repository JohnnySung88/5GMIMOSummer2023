load('Tx_signal.mat')
load('frame_data.mat')
Rx_signal = simChannel(Tx_signal,20);
DTinfo='LMMSE'
CFO_ignore = 5;
[Time,Biterror,Capacity_sum,Rx1_SNR,Rx2_SNR,JPG_RGB] = Rx_fun_LDPC(frame_data,Rx_signal,DTinfo,CFO_ignore,5)
imshow(JPG_RGB);