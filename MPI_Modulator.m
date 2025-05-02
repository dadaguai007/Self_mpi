clc;clear;close all;
addpath("D:\BIT_PhD\Base_Code\Codebase_using")
% PAM transimision

sps = 2;
Rs  = 40e9;
Ts  = 1/Rs ;
fs  = sps*Rs;
Ta  = 1/fs;

Vpi = 2;
Vb = -Vpi/2;

% ssfm
param=struct();
param.Ltotal = 10; %km
param.Lspan =5;
param.hz= 0.1;
param.alpha=0.2;
param.D = 16;
param.gamma = 1.3;
param.Fc = 193.1e12;
param.NF = 4.5;
param.amp='ideal';
param.Fs=fs;

% phase noise
data_num=1e6;
lw      = 1e3;    % laser linewidth
phi_pn_lo = phaseNoise(lw, data_num, Ta);
sigLO = exp(1j * phi_pn_lo);

%LO
Pi_dBm = 10;
Pi = 10^(Pi_dBm/10)*1e-3; %W
Ai= sqrt(Pi);
Pin=Ai;

%PAM
M=4;
data_2bit=randi([0,1],log2(M),80000);
% 相当于四个电平
symbols = 2.^(0:log2(M)-1)*data_2bit;

% Mapeia bits para pulsos elétricos
symbTx = pammod(symbols,M,0,'gray');
symbTx = pnorm(symbTx);

% MDRC 编码
res1=[];
N=20;
for j=1:length(symbTx)/N
    x=symbTx((j-1)*N+1:j*N);
    % 计算RDS
    for i=1:N
        if i==1
            z(1)=x(1);
        else
            z(i)=z(i-1)+x(i);
        end
    end
    % 确定反转点
    z_k(j)=(z(N)/2);
    [~,k_index(j)]=min(abs(z-z_k(j)));

    % 进行反转
    res=[x(1:k_index(j)),-x(k_index(j)+1:end)];

    % 排列传输信号
    res1=[res1,res];
end


% Upsampling
symbolsUp_reverse = upsample(res1, sps);

% %Pulso
hsqrt = rcosdesign(0.2,256,sps,'sqrt');
hsqrt=hsqrt./max(abs(hsqrt));
% pulse shaping
sigTx=conv(symbolsUp_reverse,hsqrt,'same');




% MDRC的反转编码
for i=1:length(sigTx)
    sigTx_reverse(i)=(-1).^i*sigTx(i);
end

% 对PAM信号脉冲成型
symbolsUp = upsample(symbTx, sps);
sigTx_alamouti_pluse=conv(symbolsUp,hsqrt,'same');


% Alamouti Code
codedData = alamouti_code_pam(sigTx_alamouti_pluse);


% 符号相反的Alamouti Code
codedData2 = alamouti_code_pam_reverse(sigTx_alamouti_pluse);

% MDRC码进行Alamouti编码
codedData_MDRC=alamouti_code_pam(sigTx);

% Alamouti——相加编码
codedData3 = alamouti_code_pam_linear(sigTx_alamouti_pluse);

% EA放大
amp_factor=0.1;


% 测试 Alamouti 编码 在 DD-mzm 和 mzm  中的性能

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       mzm     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sig_Tx=amp_factor*codedData3(1,:);
sigTxo = mzm(Pin, Sig_Tx, Vb,Vpi);

% 通过光纤
E_x=ssfm(sigTxo,param);
I_x=(E_x).*conj(E_x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     DD-mzm     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vbias1=Vpi/4;
Vbias2=-Vpi/4;

VRF1=amp_factor*codedData3(1,:);
VRF2=amp_factor*codedData3(2,:);
[Eout] = MZM_DualDrive(Pin,VRF1,VRF2,Vbias1,Vbias2,Vpi,1);

% 通过光纤
Eo=ssfm(Eout,param);
I_pd= (Eo).*conj(Eo);

figure;
mon_ESA(I_pd,fs);
title("DD-MZM")
figure;
mon_ESA(I_x,fs);
title("MZM")