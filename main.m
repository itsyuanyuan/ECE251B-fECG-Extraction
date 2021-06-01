clear all; close all;
%%
D_name = ["r01.csv"; "r04.csv"; "r07.csv"; "r08.csv"; "r10.csv"]';
for i = 1:5
    file = "./dset/" + D_name(i);
    if i == 1
        D = readmatrix(file);
    else
        tmp = readmatrix(file);
        D = [D; tmp];
    end
end
D = D * 1000;
%%
qrsfile = "dset/r07.edf";
[header, S1]=edfread(qrsfile);
%%
Fs = 1000;
channel = 2;
u = D(300000*3+1:300000*4,5);
d = D(300000*3+1:300000*4,1);
% d = TWRemoveBL(d,Fs,5);
% u = TWRemoveBL(u,Fs,5);\
%%
up = D(300000*2+1:300000*3,5);
dp = D(300000*2+1:300000*3,1);
u_f = medfilt1(up, 100);
u_clean_p = u-u_f;
d_f = medfilt1(dp, 100);
d_clean_p = dp - d_f;
aaa = conv(u_clean_p, w,'same');
%%
aaaa = aaa*10000;
figure;
hold on
plot(up(10000:16000))
plot(aaaa(10000:16000))
%%
figure;
plot(d(10000:20000))
xlabel('time(ms)');
ylabel('Voltage(mV)');
ylim([-0.2,0.2]);
title('Signal observed from V0')
%%
figure;
plot(u(10000:20000))
xlabel('time(ms)');
ylabel('Voltage(mV)');
ylim([-0.2,0.2]);
title('Signal observed from A4')
u_f = medfilt1(u, 100);
u_clean = u-u_f;
figure;
plot(u_clean(10000:20000))
%%
d_f = medfilt1(d, 100);
d_clean = d - d_f;
%%
figure;
hold on 
plot(d_clean(10000:12000));
plot(u_clean(10000:12000))

xlabel('time(ms)');
ylabel('Voltage(mV)');
legend('V0','AB4');
title('Maternal v.s. Abnomed');
d_clean = d_clean / 4.0;
%%
figure;
hold on 
plot(u_clean(10000:20000));
plot(d_clean(10000:20000));
%%
d_low = lowpass(d_clean, 5 , 1000);
u_low = lowpass(u_clean, 3 , 1000);

figure;
hold on 
plot(d_low(10000:20000));
plot(u_low(10000:20000));
%% LMS
lms = dsp.LMSFilter(20, 'StepSize', 0.00001, 'Method', 'LMS', 'WeightsOutputPort', true);
[mu_max, mu_min] = maxstep(lms, u_clean);

[y,e,w] = step(lms, u_clean, d_clean);
figure;
hold on
plot(e(10000:18000),'-');
plot(u_clean(10000:18000));
legend('e', 'u')
%% NLMS
lms = dsp.LMSFilter(20, 'StepSize', 0.00001, 'Method', 'Normalized LMS', 'WeightsOutputPort', true);
[y,e,w] = step(lms, u_clean, d_clean);
figure;
hold on
plot(e(10000:18000),'-');
plot(u_clean(10000:18000));
legend('e', 'u')
%% SNR
SNR =10 * log(sum(e.^2)/sum(d_clean.^2))
%%
k = y*100000;
figure;
hold on
plot(k(10000:12000));
plot(d_clean(10000:12000),'--');
legend('y', 'dclean');
%%
figure;
hold on
plot(e(10000:12000));
plot(d_clean(10000:12000),'--');
xlabel('time(ms)');
ylabel('Voltage(mV)');
legend('Fetal(LMS)', 'Maternal(V0)');
title('LMS output compares to Maternal');

%%
l = (u_clean - e);
figure;
hold on 
plot(l(10000:12000));
plot(u_clean(10000:12000),'--');
xlabel('time(ms)');
ylabel('Voltage(mV)');

legend('Fetal(LMS)', 'Abdomen4(AB4)');
title('LMS filtered signal');
%%
figure
plot(e(10000:12000));
figure;
hold on
plot(e(10000:12000));
plot(d_clean(10000:12000));
legend('e', 'd');
%%
Fech = e-u_clean;
figure;
plot(Fech(10000:12000))
%%
figure;
freqz(w);
%% RLS
rls = dsp.RLSFilter('Length', 10, 'ForgettingFactor', 1, 'Method', 'Conventional RLS');
[y,e] = rls(u_clean, d_clean);
%%
figure;
hold on
plot(e(10000:12000));
plot(d_clean(10000:12000),'--');

xlabel('time(ms)');
ylabel('Voltage(mV)');
legend('RLS output', 'Maternal(V0)')
title('RLS filtered signal')
%%
l = u_clean - y;
figure;
hold on
plot(l(10000:12000))
plot(y(10000:12000),'--');
legend('Fetal(RLS)', 'Abdomen4(AB4)');
xlabel('time(ms)');
ylabel('Voltage(mV)');
title('RLS filtered signal')

%%
figure
%%
[qrs_amp_raw , qrs_i_raw , delay] = pan_tompkin(d_clean,Fs,1);
[qrs_amp_raw2 , qrs_i_raw2 , delay2] = pan_tompkin(e,Fs,1);
[qrs_amp_raw3 , qrs_i_raw3 , delay3] = pan_tompkin(l,Fs,1);

total_loc = 0;
for i = 2:length(qrs_i_raw)
    range = abs(qrs_i_raw(1,i) - qrs_i_raw(1,i-1));
    total_loc = total_loc + range;
end

total_loc2 = 0;
for i = 2:length(qrs_i_raw2)
    range2 = abs(qrs_i_raw2(1,i) - qrs_i_raw2(1,i-1));
    total_loc2 = total_loc2 + range2;
end

total_loc3 = 0;
for i = 2:length(qrs_i_raw3)
    range3 = abs(qrs_i_raw3(1,i) - qrs_i_raw3(1,i-1));
    total_loc3 = total_loc3 + range3;
end

mean_loc = total_loc/(length(qrs_i_raw) - 1); 
mean_loc2 = total_loc2/(length(qrs_i_raw2) - 1); 
mean_loc3 = total_loc3/(length(qrs_i_raw3) - 1);

bpm_abdomen = (60*1000)/mean_loc           
bpm_direct = (60*1000)/mean_loc2
bpm_fetus = (60*1000)/mean_loc3

%%
function [ECG_Clean] = TWRemoveBL(ecgy, Fs, Fc)
%  ecgy:        the contamined signal
%  Fs:          sample frequiency
%  Fc:          cut-off frequency
%  ECG_Clean :  processed signal without BLW

    % Wavelet type
    w='sym10';
    % Threshold type
    thr_met='s';

    % Perform decomposition at level 5
    %[c,l] = wavedec(s,5,w);

    % The decomposition leve is estimated based in the cut-off frequency
    lev = ceil(log2(Fs/Fc));

    BL = wden(ecgy,'heursure',thr_met,'one',lev, w);

    ECG_Clean = ecgy - BL;

end

