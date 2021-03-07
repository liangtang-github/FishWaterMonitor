%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  药物名称：铜离子（硫酸铜）,中华青鱂鱼
%  药物浓度: 10TU,  1TU=7.9mg/L
%  开始时间：2017.8.3.9:33
%  加药时间：2017.8.3.13:34
%  结束时间：2017.8.3.16:33
%  按照一分钟一次：241分钟时下毒
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;
close all;

dosingInv = 241;

h = figure;
set(h,'position',[300 150 900 580]);

% read 10TU TCP data
load('CuSO410TU20170803.mat');
IntegSignal = VarName3 + VarName4 + VarName5 + VarName6 + VarName7 + ...
              VarName8 + VarName9 + VarName10;
IntegSignal = IntegSignal ./ 8;
signal = GetSignalByMinute(IntegSignal, 1200);

% orignal signal
subplot(3,2,[1, 2]);
set(gca,'position',[0.055,0.72,0.925,0.24]);
plot(IntegSignal,'k','LineWidth',1.5);
hold on;
line([dosingInv*1200,dosingInv*1200],[0,1],'linestyle',':','LineWidth',2,'color','k');
grid on;
annotation(h,'textarrow',[0.560 0.586], [0.928 0.928],'String',{'13:34'});
xlabel('Sequentially');
ylabel('Signal intensity');
axis([0 size(IntegSignal,1) 0 1])
title('(a). 10TU CuSO_4 exposure experiment original signal');

% EMD method
subplot(3,2,3);
set(gca,'position',[0.055,0.40,0.42,0.24]);
signal1 = emd(signal);
signal1 = (signal1 - min(signal1)) ./ (max(signal1) - min(signal1));
plot(signal1(:,4),'k','LineWidth',1.5);
hold on;
line([dosingInv,dosingInv],[0,1],'linestyle',':','LineWidth',2,'color','k');
grid on;
annotation(h,'textarrow',[0.271 0.297], [0.614 0.614],'String',{'241'});
xlabel('Time span/(minute)');
ylabel('Signal intensity');
xticks([0 50 100 150 200 250 300 350 400]);
axis([0 size(signal1,1) 0 1])
title('(b). EMD method');

% SG smooth method
subplot(3,2,4);
set(gca,'position',[0.53,0.40,0.448,0.24]);
signal2 = sgolayfilt(signal,3,61);
signal2 = (signal2 - min(signal2)) ./ (max(signal2) - min(signal2));
plot(signal2,'k','LineWidth',1.5);
hold on;
line([dosingInv,dosingInv],[0,1],'linestyle',':','LineWidth',2,'color','k');
grid on;
annotation(h,'textarrow',[0.761 0.788], [0.519 0.519],'String',{'241'});
annotation(h,'textarrow',[0.898 0.874], [0.429 0.407],'String',{'332'});
xlabel('Time span/(minute)');
ylabel('Signal intensity');
xticks([0 50 100 150 200 250 300 350 400]);
axis([0 size(signal2,1) 0 1])
title('(c). Savitzky Golay method');

% wavelet method
subplot(3,2,5);
set(gca,'position',[0.055,0.0748,0.42,0.24]);
[c,l] = wavedec(signal,3,'db3');
approx = appcoef(c,l,'db3');
approx = (approx - min(approx)) ./ (max(approx) - min(approx));
signal3 = interp(approx, 8);
plot(signal3,'k','LineWidth',1.5);
hold on;
line([dosingInv,dosingInv],[0,1],'linestyle',':','LineWidth',2,'color','k');
grid on;
annotation(h,'textarrow',[0.254 0.281], [0.195 0.195],'String',{'241'});
annotation(h,'textarrow',[0.392 0.368], [0.102 0.073],'String',{'337'});
xlabel('Time span/(minute)');
ylabel('Signal intensity');
xticks([0 50 100 150 200 250 300 350 400]);
axis([0 size(signal3,1) 0 1])
title('(d). Wavelet method');

% unwinding method
subplot(3,2,6);
set(gca,'position',[0.53,0.0748,0.448,0.24]);
t=linspace(0,2*pi,length(signal));
L=15;
[F, err, S1, a, amp, ~, blaschke_z] = Unwinding_Blaschke(signal.',L,t);
[~, amplitude] = disp_uw(a,F);
wavesignal = amplitude(10,:);
wavesignal = (wavesignal - min(wavesignal)) ./ (max(wavesignal) - min(wavesignal));
signal4 = wavesignal;
plot(signal4,'k','LineWidth',1.5);
hold on;
line([dosingInv,dosingInv],[0,1],'linestyle',':','LineWidth',2,'color','k');
grid on;
annotation(h,'textarrow',[0.761 0.788], [0.283 0.283],'String',{'241'});
annotation(h,'textarrow',[0.907 0.882], [0.294 0.314],'String',{'330'});
xlabel('Time span/(minute)');
ylabel('Signal intensity');
xticks([0 50 100 150 200 250 300 350 400]);
axis([0 size(signal4,2) 0 1])
title('(e). Unwinding method');