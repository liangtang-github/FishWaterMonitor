clear;
clc;
close all;

h = figure;
set(h,'position',[500 300 450 200]);

% Behavioral signals in a safe state.
load('CuSO410TU20170803.mat');
IntegSignal = VarName3 + VarName4 + VarName5 + VarName6 + VarName7 + ...
              VarName8 + VarName9 + VarName10;
IntegSignal = IntegSignal ./ 8;
saftySignal = GetSignalByMinute(IntegSignal, 1200);
saftySignal = saftySignal(1:200);
saftySignal = (saftySignal - min(saftySignal)) ./ (max(saftySignal) - min(saftySignal));
plot(saftySignal,'k','LineWidth',1.5);
grid on;
yticks([0 0.25 0.50 0.75 1.0]);
%xticks([0 100 200 300 400 500 600]);
axis([0 size(saftySignal,1) 0 1])
xlabel('Time span/(minute)');
ylabel('Behavior intensity');
title('Medaka Behavior Signal');
set(gca,'LooseInset',get(gca,'TightInset'))
set(gca,'looseInset',[0 0 0 0])