% ����һ�������ź�
Fs = 1000;
T = 1/Fs;
time = 2;
L = Fs* time;
t = (0: L-1)* T;
x1 = 4* sin(2* pi* 10* t);
x2 = 2* sin(2* pi* 20* t);
x3 = 1* sin(2* pi* 30* t);
x =  x1+ x2+ x3;
figure(1);
subplot(4,1,1); plot(t, x1);
subplot(4,1,2); plot(t, x2);
subplot(4,1,3); plot(t, x3);
subplot(4,1,4); plot(t, x);

% ���ź���ʱ�������
TE = sum(x.^2);

% �������ź�x�����ٸ���Ҷ�任
% �󵥱��ף���������
fx = fft(x);
ft = Fs* (0: L/2-1)/L;
fxP2 = (fx./ L).* (conj(fx)./ L);
fxP1 = fxP2(1: L/2);
fxP1(2:end) = 2* fxP1(2:end);
figure(2);
plot(ft(1:50/Fs*L), fxP1(1:50/Fs*L));

% ��Ƶ��������
ftE = Fs* (0: L-1)/L;
fE = fx.*conj(fx)/L;
figure(3);
plot(ftE, fE);

% ���ź���Ƶ�������
FE = sum(fE);

% ʹ��Welch�����ƹ������ܶ�
nfft = 2^nextpow2(L);
[fxPw, ftW] = pwelch(x, 1000, 500, nfft, Fs);
figure(4);
plot(ftW(1:50/Fs*L), fxPw(1:50/Fs*L));





