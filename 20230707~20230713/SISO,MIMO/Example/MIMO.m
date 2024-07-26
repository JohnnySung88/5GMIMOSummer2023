% MIMO system
clear
clc

Tx =4;         %傳送端個數
Rx =4;         %接收端個數

%可以先自己設
data_num = 100000;           % data量(注意BER要跑到10^(-3)!!)
SNR_in_dB = 0:5:40;          % 自己設訂雜訊大小

SER_SNR_ZF=zeros(3,length(SNR_in_dB));
BER_SNR_ZF=zeros(3,length(SNR_in_dB));
SER_SNR_LMMSE=zeros(3,length(SNR_in_dB));
BER_SNR_LMMSE=zeros(3,length(SNR_in_dB));

%% choose QAM= 4/16/64;
for v=1:3
    QAM = 4^v;
    Eavg = ;
    NF = ;
    q_bit = ;        % 一個symbol可以傳幾個bit
    N = Tx;

    for  a=1:length(SNR_in_dB)
        SNR = 10^(SNR_in_dB(a)/10);
        No = 10^(-SNR_in_dB(a)/10);
        Es = 1;
        
        SER_ZF = 0;     
        BER_ZF = 0; 
        SER_LMMSE = 0;  
        BER_LMMSE = 0;
        
        for seperate = 1:data_num
            data = ;           % 隨機產生0~3 for 4QAM
            bin_data = ;       % 將 0~3 轉為 '00'~'11'
            X = ;              % 0~3 to complex (Modulation); remember to normalize
            H = ;              % randn產生channel(注意正規化的問題)
            n = ;              % randn產生noise variance=No
            
            Y = ;
            
            %% type 1 = ZF
            X_hat_ZF = ;
            data_hat_ZF = ;          %complex to 0~3  ; remember to inverse-normalize
            bin_data_hat =;          %0~3 to '00'~'11'
            % 算SER/BER
            SER_ZF = ;
            BER_ZF = ;
            
            %% type 2 = LMMSE
            X_hat_LMMSE = ;
            data_hat_LMMSE = ;          % complex to 0~3 ; remember to inverse-normalize
            bin_data_hat_LMMSE = ;      %0~3 to '00'~'11'
            % SER/BER
            SER_LMMSE = ;
            BER_LMMSE = ;
           
        end
%         
%         按照SNR把算好的SER/BER存在矩陣裡
        SER_SNR_ZF(v,a) = SER_ZF/(data_num*Tx);
        BER_SNR_ZF(v,a) = BER_ZF/(data_num*q_bit*Tx);
        
        SER_SNR_LMMSE(v,a) = SER_LMMSE/(data_num*Tx);
        BER_SNR_LMMSE(v,a) = BER_LMMSE/(data_num*q_bit*Tx);
        
    end
end
% 
% 輸入SNR_in_dB和SER
figure(1)
semilogy(SNR_in_dB,SER_SNR_ZF(1,:),'r-X')
hold on
semilogy(SNR_in_dB,SER_SNR_ZF(2,:),'r-diamond')
hold on
semilogy(SNR_in_dB,SER_SNR_ZF(3,:),'r-O')
hold on
semilogy(SNR_in_dB,SER_SNR_LMMSE(1,:),'b-X')
hold on
semilogy(SNR_in_dB,SER_SNR_LMMSE(2,:),'b-diamond')
hold on
semilogy(SNR_in_dB,SER_SNR_LMMSE(3,:),'b-O')
hold on
grid on
title('SER of MIMO')
xlabel('SNR (dB)')
ylabel('SER')
legend('4QAM ZF','16QAM ZF','64QAM ZF','4QAM LMMSE','16QAM LMMSE','64QAM LMMSE')

%輸入SNR_in_dB和BER
figure(2)
semilogy(SNR_in_dB,BER_SNR_ZF(1,:),'r-X')
hold on
semilogy(SNR_in_dB,BER_SNR_ZF(2,:),'r-diamond')
hold on
semilogy(SNR_in_dB,BER_SNR_ZF(3,:),'r-O')
hold on
semilogy(SNR_in_dB,BER_SNR_LMMSE(1,:),'b-X')
hold on
semilogy(SNR_in_dB,BER_SNR_LMMSE(2,:),'b-diamond')
hold on
semilogy(SNR_in_dB,BER_SNR_LMMSE(3,:),'b-O')
hold on
grid on
title('BER of MIMO')
xlabel('SNR (dB)')
ylabel('BER')
legend('4QAM ZF','16QAM ZF','64QAM ZF','4QAM LMMSE','16QAM LMMSE','64QAM LMMSE')

