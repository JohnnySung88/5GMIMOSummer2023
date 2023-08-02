clear;
clc;

SNR_in_dB = 0:5:40;
BER_SNR=zeros(1,length(SNR_in_dB));

QAM 	= 16;
Eavg 	= (qammod([0:QAM-1],QAM) * qammod([0:QAM-1],QAM)') / QAM;
NF 		= 1 / sqrt(Eavg);
q_bit 	= log2(QAM);        % 一個symbol可以傳幾個bit

Nt = 4; % number of transmit antennas
Nr = 4; % number of receive antennas

for a=1:length(SNR_in_dB)
	SNR = 10^(SNR_in_dB(a)/10);
	No  = 1/SNR;
	BER = 0;

    for iterations=1:1e5
        % TX
        data_dec	= randi([0,QAM-1],Nt,1);
        data_bin = de2bi(data_dec, q_bit, 'left-msb');	
        data_mod = NF*qammod(data_dec, QAM, 'gray');

        % Channel
        H = 1/sqrt(2)*(randn(Nr,Nt) + 1i*randn(Nr,Nt));
        
        % Noise
        n = sqrt(No/2) *( randn(Nr,1) + randn(Nr,1)*1i );
        
        % RX
        y = H*data_mod + n;
        
        % ZF Detector
        X_hat = H\y;
        
        % De-Mod
        data_dec_hat = qamdemod(X_hat/NF,QAM,'gray');
        data_bin_hat = de2bi(data_dec_hat, q_bit, 'left-msb');
        
        % BER
        BER = BER + sum(sum(data_bin ~= data_bin_hat ));
    end
    
	BER_SNR(1,a) = 	BER/(1e5 * Nt * q_bit);
end

figure(1)
semilogy(SNR_in_dB,BER_SNR(1,:),'r-x', 'LineWidth',2)
hold on
grid on
axis square;
axis tight;
title('BER of 4x4 MIMO with 16-QAM and ZF Detector')
xlabel('SNR (dB)')
ylabel('BER')
legend('16QAM ZF MATLAB')
