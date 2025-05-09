% 复现：MPI 多径串扰的模拟仿真
clear;clc;close all;
% addpath('D:\PhD\Codebase\')
addpath("D:\BIT_PhD\Base_Code\Codebase_using")
addpath('Fncs\')
addpatn('Decode\')
% paramer
sps = 6;
Rs  = 40e9;
Ts  = 1/Rs ;
fs  = sps*Rs;
Ta  = 1/fs;

% mzm
Vpi = 10;
Vb = -Vpi/2;


% fiber
param=struct();
param.Ltotal = 20; %km
param.Lspan =10;
param.hz= 0.1;
param.alpha=0.2;
param.D = 16;
param.gamma = 1.3;
param.Fc = 193.1e12;
param.NF = 4.5;
param.amp='ideal';
param.Fs=fs;


paramPD=struct();
paramPD.B =Rs;
paramPD.R =1;
paramPD.type = 'ideal';
paramPD.Fs=fs;




rng(123);
%PAM
M=4;
data_2bit=randi([0,1],log2(M),80000);
% 相当于四个电平
symbols = 2.^(0:log2(M)-1)*data_2bit;
% Mapeia bits para pulsos
symbTx = pammod(symbols,M,0,'gray');
% 传输信号归一化
symbTx = pnorm(symbTx);
% 参考向量,转换为行向量
label=symbTx.';


% 装载DSP-CDR 和解码模块
SignalRec  =  DspSyncDecoding(...
    fs,...    % 接收信号的采样率
    Rs,...    % 接收信号的波特率
    4,...     % 接收信号的调制格式
    sps,...   % 上采样率
    2*sps,...   % 时钟信号的上采样率
    1e5,...   % 误码计算起始位置
    label);      % 参考信号



% Upsampling
symbolsUp = upsample(symbTx, sps);

% Puls
hsqrt = rcosdesign(0.1,256,sps,'sqrt');
% pulse shaping
sigTx=conv(symbolsUp,hsqrt,'same');

% EA
amp_factor=0.5;
Sig_Tx=sigTx*amp_factor;

% mod index
m=Modulation_index(Sig_Tx,Vpi,'pam');
fprintf(' the module index =%.3f \n', m);


% Laser
% 输入光功率
type='laser';
Pi_dBm = 10;
Pi = 10^(Pi_dBm/10)*1e-3; %W
Ai= sqrt(Pi);
lw      = 1e3;    % laser linewidth
phi_pn_lo = phaseNoise(lw, length(sigTx), Ta);
sigLO = exp(1i * phi_pn_lo);
if strcmp(type,'laser_phase')
    Pin=Ai*sigLO;
else
    Pin=Ai;
end

%mzm
sigTxo = mzm(Pin, Sig_Tx, Vb,Vpi);
power=signalpower(sigTxo);
fprintf(' after module signal power: %.2f dBm\n', 10 * log10(power / 1e-3));

% 1:9 copulter
aerfa=0.1;
[sigTxo_train,sigTxo_reflect]=Coupler(aerfa,sigTxo);

% [sig1,sigTotal1]=Coupler(1-aerfa,sigTxo_train,sigTxo_reflect); % 两路信号通过Coupler进行合并
% (耦合比：后面耦合比加上前面矩阵耦合比为1)

power3=signalpower(sigTxo_train);
fprintf(' after Coupler signal power: %.2f dBm\n', 10 * log10(power3 / 1e-3));




% PAM_tran
sigRxo=ssfm(sigTxo_train,param);
power2=signalpower(sigRxo);
fprintf(' after ssfm signal power: %.2f dBm\n', 10 * log10(power2 / 1e-3));

%pd
ipdTrain = pd(sigRxo, paramPD);
%
% match and Dc remove
sigI=ipdTrain.';
sigI=sigI-mean(sigI);
% match
sigMatch=conv(sigI,hsqrt,'same');
% norm
sigMatch=pnorm(sigMatch);

% 降采样
sigRx_E=downsample(sigMatch,sps);
% 解码
[decodedData,ber]=SignalRec.PAM_ExecuteDecoding(sigRx_E);


% PAM_refection
% fiber-reflect
param_reflect=struct();
param_reflect.Ltotal = 1; %km
param_reflect.Lspan =0.5;
param_reflect.hz= 0.1;
param_reflect.alpha=0.2;
param_reflect.D = 16;
param_reflect.gamma = 1.3;
param_reflect.Fc = 193.1e12;
param_reflect.NF = 4.5;
param_reflect.amp='none';
param_reflect.Fs=fs;

sigRxo_reflect=ssfm(sigTxo_reflect,param_reflect);

% delay tao
delay_ps=5.6;
% sigRxo_delay1= delay_signal(sigRxo_reflect, delay_ps*1e-12);
sigRxo_delay = iqdelay(sigRxo_reflect, fs, delay_ps*1e-12);

% Attenuator
Att_dB=2;
sigRxo_ref=Attenuator(sigRxo_delay,Att_dB);


% Combine
[~,sigTotal]=Coupler(1-aerfa,sigRxo,sigRxo_ref); % 两路信号通过Coupler进行合并

power5=signalpower(sigTotal);
fprintf(' after MPI signal power: %.2f dBm\n', 10 * log10(power5 / 1e-3));





%pd
ipd = pd(sigTotal, paramPD);

% match and Dc Remove
sigMPI=ipd.';
sigMPI=sigMPI-mean(sigMPI);
% match
sigMatchMPI=conv(sigMPI,hsqrt,'same');
% norm
sigMatchMPI=pnorm(sigMatchMPI);


% downsample
sigRxMPI=downsample(sigMatchMPI,sps);
% decode
[decodedData1,ber1]=SignalRec.PAM_ExecuteDecoding(sigRxMPI);




