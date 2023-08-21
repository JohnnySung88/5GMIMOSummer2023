SSB_5G_NR;


load('LDPC_11nD2_1296b_R12.mat');
code_rate	= 1/2;
iteration	= 5;
SNR_in_dB 	= 20;
SNR_weight 	= 45;
window 		= 10;
DMRS_DATA 	= +0.7071 + 0.7071*1i ;
CFO_ignore	= 5;
Fs			= 122.88e6;
delta_f		= 0;
DTinfo 		= 'LMMSE';

% choose QAM= 4/16/64;
QAM 	= 16;
Eavg 	= (qammod([0:QAM-1],QAM) * qammod([0:QAM-1],QAM)') / QAM;
NF 		= 1 / sqrt(Eavg);
q_bit 	= log2(QAM);        % 一個symbol可以傳幾個bit
Tx		  	= 2;
Rx		  	= 2;
binTable  	= 2.^(q_bit-1:-1:0);%[8 4 2 1]
%生成PSS SSS PBCH資料
PSS  = [ zeros(758,1);PSS;zeros(759,1) ];
%初始化通道估測
N_fft 	=2048;
R_HD_HD	=zeros( 822,822);
R_H_HD	=zeros(1644,822);
DMRS_pos=[204:2:1024,1027:2:1847];
Real_pos=[203:1:1024,1026:1:1847];
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
SNR_W 	= 10^( SNR_weight /10);

P0 			= eye(822);
P1 			= ones(822,1);P1(1:2:822) = -1;
P1 			= diag(P1);
LMMSE_1_1	= R_H_HD* P0' * inv( P0*R_HD_HD*P0' + P1*R_HD_HD*P1' + (1/SNR_W)*eye(822) );
LMMSE_2_1	= R_H_HD* P0' * inv( P0*R_HD_HD*P0' + P1*R_HD_HD*P1' + (1/SNR_W)*eye(822) );
LMMSE_1_2	= R_H_HD* P1' * inv( P0*R_HD_HD*P0' + P1*R_HD_HD*P1' + (1/SNR_W)*eye(822) );
LMMSE_2_2	= R_H_HD* P1' * inv( P0*R_HD_HD*P0' + P1*R_HD_HD*P1' + (1/SNR_W)*eye(822) );


%LDPC init
[row,col]	= find(double(LDPC.H.x)==1);
[A,I] = sort(row);
H_row_master 	  = zeros(648,8); 
H_row_master_size = ones(648,1);
for i=1:4644
	H_row_master( A(i) , H_row_master_size( A(i)) ) = col(I(i));
	H_row_master_size(A(i)) = H_row_master_size(A(i)) + 1;
end
H_row_master_size = H_row_master_size - 1;

%frame架構 0=data 1=DMRS+ -1=DMRS- PSS=2
frameTable 		= zeros(1644,560,Tx);
frameTable(2:2:1644,3:14:560, :    ) =  1;%DMRS frame pos
frameTable(2:4:1644,3:14:560,2:2:Tx) = -1;%DMRS frame pos(CDM)
frameTable(   :    ,   5    ,   :  ) =  2;%PSS  frame pos

DATA_Pos = find(frameTable==  0);
DMRS1_Pos= find(frameTable==  1);
DMRS2_Pos= find(frameTable== -1);
PSS_pos  = find(frameTable==  2);
