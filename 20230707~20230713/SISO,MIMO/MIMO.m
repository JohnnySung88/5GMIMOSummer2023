% MIMO system
clear
clc

Tx = 4;     % 傳送端個數
Rx = 4;     % 接收端個數

% 可以先自己設
data_num  = 1000000;        % data量(注意BER要跑到10^(-3)!!)
SNR_in_dB = 0:5:40;         % 自己設訂雜訊大小

SER_SNR_ZF    = zeros(3,length(SNR_in_dB));
BER_SNR_ZF    = zeros(3,length(SNR_in_dB));
SER_SNR_LMMSE = zeros(3,length(SNR_in_dB));
BER_SNR_LMMSE = zeros(3,length(SNR_in_dB));

%% choose QAM= 4/16/64;
for v     = 1:3
    QAM   = 4^v;
    Eavg  = (qammod((0:QAM-1),QAM)*qammod((0:QAM-1),QAM)')/QAM;
    NF    = 1/sqrt(Eavg);
    q_bit = log2(QAM);      % 一個symbol可以傳幾個bit
    N     = Tx;

    parfor a   = 1:length(SNR_in_dB)
        Es  = 1;
        SNR = 10^(SNR_in_dB(a)/10);
        No  = Es/SNR;
        
        SER_ZF    = 0;     
        BER_ZF    = 0;
        SER_LMMSE = 0;  
        BER_LMMSE = 0;
        
        for seperate = 1:data_num
            data     = randi([0 QAM-1], Tx, 1);                             % 隨機產生0~3 for 4QAM
            bin_data = de2bi(data, q_bit, 'left-msb');                      % 將 0~3 轉為 '00'~'11'
            X        = NF*qammod(data, QAM, 'gray');                        % 0~3 to complex (Modulation); remember to normalize
            H        = (randn(Rx, Tx) + 1j*randn(Rx, Tx)) /sqrt(2*Tx);      % randn產生channel(注意正規化的問題)
            n        = (randn(Rx, 1)  + 1j*randn(Rx, 1))  *sqrt(No/2);      % randn產生noise variance=No
            Y        = H*X + n;
            
            %% type 1 = ZF
            X_hat_ZF        = (inv(H)*Y)/NF;
            data_hat_ZF     = qamdemod(X_hat_ZF, QAM, 'gray');              % complex to 0~3; remember to inverse-normalize
            bin_data_hat_ZF = de2bi(data_hat_ZF, q_bit, 'left-msb');        % 0~3 to '00'~'11'
            % SER/BER
            SER_ZF = SER_ZF + sum(data_hat_ZF     ~= data);
            BER_ZF = BER_ZF + sum(bin_data_hat_ZF ~= bin_data, 'all');
            %% type 2 = LMMSE
            X_hat_LMMSE        = (inv(H'*H + No*eye(N))*H'*Y)/NF;
            data_hat_LMMSE     = qamdemod(X_hat_LMMSE, QAM, 'gray');            % complex to 0~3; remember to inverse-normalize
            bin_data_hat_LMMSE = de2bi(data_hat_LMMSE, q_bit, 'left-msb');      % 0~3 to '00'~'11'
            % SER/BER
            SER_LMMSE = SER_LMMSE + sum(data_hat_LMMSE     ~= data);
            BER_LMMSE = BER_LMMSE + sum(bin_data_hat_LMMSE ~= bin_data, 'all');
           
        end
        
        % 按照SNR把算好的SER/BER存在矩陣裡
        SER_SNR_ZF(v,a)    = SER_ZF/(data_num*Tx);
        BER_SNR_ZF(v,a)    = BER_ZF/(data_num*q_bit*Tx);
        SER_SNR_LMMSE(v,a) = SER_LMMSE/(data_num*Tx);
        BER_SNR_LMMSE(v,a) = BER_LMMSE/(data_num*q_bit*Tx);
        
    end
end

%% 輸入SNR_in_dB和SER
figure(1)
semilogy(SNR_in_dB,SER_SNR_ZF(1,:),    'r-X',       'LineWidth',2)
hold on
semilogy(SNR_in_dB,SER_SNR_ZF(2,:),    'r-diamond', 'LineWidth',2)
hold on
semilogy(SNR_in_dB,SER_SNR_ZF(3,:),    'r-O',       'LineWidth',2)
hold on
semilogy(SNR_in_dB,SER_SNR_LMMSE(1,:), 'b-X',       'LineWidth',2)
hold on
semilogy(SNR_in_dB,SER_SNR_LMMSE(2,:), 'b-diamond', 'LineWidth',2)
hold on
semilogy(SNR_in_dB,SER_SNR_LMMSE(3,:), 'b-O',       'LineWidth',2)
hold on
grid on
axis square, title('SER of MIMO')
xlabel('SNR (dB)')
ylabel('SER')
legend('4QAM ZF','16QAM ZF','64QAM ZF','4QAM LMMSE','16QAM LMMSE','64QAM LMMSE')

% 輸入SNR_in_dB和BER
figure(2)
semilogy(SNR_in_dB,BER_SNR_ZF(1,:),    'r-X',       'LineWidth',2)
hold on
semilogy(SNR_in_dB,BER_SNR_ZF(2,:),    'r-diamond', 'LineWidth',2)
hold on
semilogy(SNR_in_dB,BER_SNR_ZF(3,:),    'r-O',       'LineWidth',2)
hold on
semilogy(SNR_in_dB,BER_SNR_LMMSE(1,:), 'b-X',       'LineWidth',2)
hold on
semilogy(SNR_in_dB,BER_SNR_LMMSE(2,:), 'b-diamond', 'LineWidth',2)
hold on
semilogy(SNR_in_dB,BER_SNR_LMMSE(3,:), 'b-O',       'LineWidth',2)
hold on
grid on
axis square, title('BER of MIMO')
xlabel('SNR (dB)')
ylabel('BER')
legend('4QAM ZF','16QAM ZF','64QAM ZF','4QAM LMMSE','16QAM LMMSE','64QAM LMMSE')
