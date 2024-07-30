clear
clc
%% NR frame structure
SCS=60000;         % subcarrier spacing (Hz): Type A
Nfft=2048;         % FFT size
Trf=1e-2;          % 10 ms radio frame (sec)
Ts=1/(SCS*Nfft);   % Sample time (sec)
FS=1/Ts;           % Sample Frequency (Hz)
Nts=Trf/Ts;        % Number of samples per 10 ms radio frame

NSymbolPerSlot=14; % Number of OFDM symbols per slot
NSlotPerFrame=40;  % Number of slots per 10 ms radio frame
NSlotPerHalfFrame=5*SCS/15000; % Number of slots per 5 ms half radio frame
Ncp1=208;          % Number of samples for the first CP in a slot
Ncp2=144;          % Number of samples for the 2nd to the 14th CP in a slot
NSamplePerSlot=(Ncp1+Nfft)+(NSymbolPerSlot-1)*(Ncp2+Nfft); % Number of samples per slot

Nrb=137;           % Number of resource blocks
Nrb_sc=12;         % Number of subcarriers per resource block
%% NR physical cell ID (0,1,2,...1007)
NRPCI=86;          % NR physical cell ID (0,1,2,...1007)
N_id_2=mod(NRPCI,3);
N_id_1=floor((NRPCI-N_id_2)/3);
%% NR-SS sequence generation
PSSLength=127;     % Length of NR-PSS sequence
SSSLength=127;     % Length of NR-SSS sequence
PSS=PSS_Seq(N_id_2,PSSLength).';
SSS=SSS_Seq(N_id_1,N_id_2,SSSLength).';
% %% SSB Structure 
% SSB=zeros(240,4);
% SSB(58:58+PSSLength-1,1)=PSS;
% SSB(58:58+PSSLength-1,3)=SSS;









