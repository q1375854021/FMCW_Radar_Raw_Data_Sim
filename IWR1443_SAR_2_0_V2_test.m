%%%This script is used to test the signal sampled by the IWR1443 ,if it's
%%%expressed by complex number.
%%%�ο�����Kdil
%%%���������420��
%��ʹ�þ�����ѹ���ͷ�λ��ѹ��
clear all;close all;clc;
%% The parameter of IWR1443.
c=3.0e8;             %Light Velocity.
B=1798.92e6;             %Bandwidth    change   Hz
K=29.982e12;             %Frequency Slop  change  Hz/S
T=B/K;                   %Ramp End Time  change S
Tc=160e-6;               %Chirp Cycle Time = Ramp End Time + Idle Time  =��100+60��e-6  change S
fs=10e6;                 %Sample rate change sps��How many times ADC sampling per second�� ��λsps  1sps=1Hz
f0=77e9;                 %Start Freq   change Hz
lambda=c/f0;             %Radar signal wavelength
d=lambda/2;              %Antenna array spacing   
nChirps = 12;
nSamples = 512;      %ADC samples per chirp
allFrame = 3750;       %�ܹ��ɼ�����֡  ������汾���õ�allFrame��ʱ����mFFT
%%%
V = 0.09166667 ;           %��λ m/s �����ٶ� 5500mm/min  
Lsar = 1400e-3;
Nr = 512;%��ʱ������� 
Na = 3750;%��ʱ�������
%%%

%%%
PRF = 250;%��λ�����Ƶ��
PRT=1/PRF;
% Ta=Na*PRT;
% Tr=Nr/fs;
% % fr=linspace(-fs/2,fs/2,Nr); 
% fr=linspace(0,fs,Nr); %%%ȥ��������FFT��Ƶ��
% ta=linspace(Ta,0,Na);  %��ʱ��   %�о������������ط�    %��ʱ��
% fa=linspace(-PRF/2,PRF/2,Na);
% tr=linspace(-Tr/2,Tr/2,Nr);%��ʱ�� 
fr = (0:Nr-1)./Nr.*fs;                  %������Ƶ��
tr = fr/K;                          %������ʱ��
fa = (-Na/2:Na/2-1)./Na.*PRF;       %��λ�������Ƶ��
ta = (-Na/2:Na/2-1).*PRT;           %��λ��ʱ��
range_f_axis = fr;
range_r_axis = (c*range_f_axis)./(2*K);      %fft ��ʱ������������ 
range_fang_axis = V.*ta;   %��λ������߶Ȼ��� ���֣��������е�ʱ�䣻���߳��ȣ��������        ����ʱ���Ϊ���ĶԳƣ�0~N-1��Ϊ -N/2~N/2-1

d_res = c/(2*B);  %����ֱ��ʣ��ο�PPT������ȷ��x������

load sif;

%һ��512�㣬��Ӧ��ʱ�䣬һ��Chirp�Ĳ�����
%һ��3750�㣬��Ӧ��ʱ�䣬��Ӧ�����λ��
figure(1)
subplot(2,1,1)
plot(0:511,real(sif((fix(Na/2)),:)));  %NaΪ��ʱ�����������Ҳ��3750  fixΪ��0�������ȡ����Na/2��ǡ��Ϊ���쵽���ĵ���ʱ��Ŀ�ʱ�������
xlabel('��ʱ�������');
ylabel('�ź�')
title('��Ƶ�ź�ʱ����ʵ������');
subplot(2,1,2)
plot(0:511,imag(sif((fix(Na/2)),:)));  %NaΪ��ʱ�����������Ҳ��3750  fixΪ��0�������ȡ����Na/2��ǡ��Ϊ���쵽���ĵ���ʱ��Ŀ�ʱ�������
xlabel('��ʱ�������');
ylabel('�ź�')
title('��Ƶ�ź�ʱ������������');


%% ���ȶԾ��������FFT
nr = 0:Nr-1;        %range �����������
Rvp = exp(-1i*pi*K*(nr./fs).^2);
range_win = hamming(Nr);   %�Ӹ�������   Ч��Ҳ������
for i = 1:Na
    sif(i,:) = (sif(i,:)-mean(sif(i,:)));   %����ֱ����Ӱ�첻��
%     sif(i,:) = fftshift(fft(fftshift(sif(i,:))));    %��ʱ��fft
    sif(i,:) = fft(sif(i,:).*range_win').*Rvp;    %��ʱ��fft
end

%ȷ������������
% range_f_axis = ((0:Nr-1)-Nr/2)/Nr*fs;    %fft ��ʱ���Ƶ��������     fftshift��



figure(2)
imagesc(range_r_axis,range_fang_axis,abs(sif));
xlabel('���������');
ylabel('��λ��ʱ���Ӧ�ľ���');
title('������FFT����');

%% ��λ��FFT
doppler_win = hamming(Na);        %���ٽض������Ƶ��й¶  
for ii=1:Nr    %��ʱ��ѭ�����Է�λ����д���  
             sif(:,ii) = (sif(:,ii)-mean(sif(:,ii)));       %�ź�ά��3750*512   //һ�����ݴ���һ֡��chirp(1 Frame�е�3��)�����㣬������3750��
                                                            %һ��512�㣬��Ӧ��ʱ�䣬һ��Chirp�Ĳ�����
                                                            %һ��3750�㣬��Ӧ��ʱ�䣬��Ӧ�����λ��
                                                            %����ȥ�����������𣬻���ȥ��ֱ������
            sif(:,ii)=fftshift(fft((fftshift(sif(:,ii).*doppler_win))));    %����ʱ��ά�ȼӴ�����������FFT
                                                            %һ��512�㣬��Ӧ��ʱ�䣬һ��Chirp�Ĳ����㣬
                                                            %һ��3750�㣬��Ӧ��ʱ���FFT����Ӧ����(��λ��FFT
end

%% �����㶯����
Rd=c*fr/2/K;
delta_R = c.*fs./(2*K*Nr);
for i = 1:Na
    for j = 1:Nr
        delta_R_eta_feta = (lambda.^2.*Rd(j).*fa(i).^2)./(8*V.^2);
        delta_n = round(delta_R_eta_feta./delta_R);
        if(j+delta_n<=Nr)
            sif(i,j) = sif(i,j+delta_n);
        else
            sif(i,j) = 0;
        end
    end
end
figure(3)
imagesc(range_r_axis,range_fang_axis,abs(ifft(sif,Na)));
xlabel('���������');
ylabel('��λ��ʱ���Ӧ�ľ���');
title('�����㶯��������');


%% ���sif
figure(4)
m = size(sif);
x = 1:m(2);
y = 1:m(1);
meshc(x,y,abs(sif));
xlabel('x  ��ʱ��  ����');
ylabel('y  ��ʱ��  ��λ');

%% ƥ���˲�����з�λ��IFFT
for i=1:Nr                   
%     Rd=-c*fr(ii)/2/K+Rc;                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Rd=-c*fr(i)/2/K;                
%     Kd=-2*V^2/lambda/Rd;
    Rd=c*fr(i)/2/K;                
    Kd=2*V^2/lambda/Rd;
    %sif(:,ii)=ifftshift(ifft(ifftshift(sif(:,ii).*exp(1i*pi*fa.^2/Kd).'.*exp(1i*4*pi/lamda*(Rd)))));
    sif(:,i)=ifftshift(ifft(ifftshift(sif(:,i).*exp(1i*pi*fa.^2/Kd).')));
end


% sif(abs(sif)>1700 & abs(sif)<6000)=sif(abs(sif)>1700 & abs(sif)<6000)*10;     %������Ŀ�������ԣ����ǻ��ڸǵ������Χ���Ŀ��
% sif(abs(sif)<1000 )=sif(abs(sif)<1000)*0;
% sif(:,220:511) = 0;
% sif(1:150,:) = 0;%%C

figure(5)
%% ���sif
m = size(sif);
x = 1:m(2);
y = 1:m(1);
meshc(x,y,abs(sif));
xlabel('x  ��ʱ��  ����');
ylabel('y  ��ʱ��  ��λ');


range_r_axis = (c*range_f_axis)./(2*K);      %fft ��ʱ������������ 
range_a_axis = linspace(0,Lsar,Na) - Lsar/2;     %������ʱ�������з���

figure(6)
% imagesc(range_r_axis,range_a_axis,255-abs(sif));
imagesc(range_r_axis,range_a_axis,abs(sif));
axis equal
xlabel('������ m');
ylabel('��λ�� m');
title('������');
image_bad = sif;
save image_bad image_bad;
