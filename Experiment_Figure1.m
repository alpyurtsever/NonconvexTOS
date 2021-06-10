%% Test setup for QAP experiment 1
% Replicates the experiment and generates Figure 1 in [YMS21].
%
% [YMS21] A. Yurtsever, V. Mangalick, and S. Sra,
% "Three Operator Splitting with a Nonconvex Loss Function"
% International Conference on Machine Learning, 2021
% 
% contact information: https://github.com/alpyurtsever

%% Preamble
clearvars
close all
rng(0,'twister');
addpath(genpath('utils'));
addpath methods;

%% Load data
% NOTE: You need to download data from QAPLIB and locate them to under the
% "./data/" folder, together with the solution files
% Links for QAPLIB:
%   http://anjos.mgi.polymtl.ca/qaplib/inst.html
%   https://www.opt.math.tugraz.at/qaplib/inst.html
%	https://coral.ise.lehigh.edu/data-sets/qaplib/qaplib-problem-instances-and-solutions/

dataname = 'chr12a';
% dataname = 'esc128';

[A,B,OPT,P] = qapread(['./data/',dataname]);
n = size(A,1);

%% Choose initial point
% Generate a random direction with Gaussian iid entries, then project
% (approximately) onto the Birkhoff polytope
X0 = ApproxProjBirkhoff(randn(n,n)./n,1e3);

%% Frank-Wolfe
[X, infoFW] = FW(A,B,'maxit',1e5,'x0',X0);

%% TOS (Split1)
[XTOSM1, infoTOSM1] = TOSM1(A,B,'maxit',1e5,'x0',X0);

%% TOS (Split2)
[XTOSM2, infoTOSM2] = TOSM2(A,B,'maxit',1e5,'x0',X0);

%% Plot the graphs
close all

cTOS1 = [247,179,34]/255;
cTOS2 = [162,30,51]/255;
cFW = [79,79,79]/255;

hfig = figure('Position',[100,100,1300,240]);
set(hfig,'name',dataname(9:end),'numbertitle','off');
set(hfig,'name',dataname,'numbertitle','off');

subplot(141)
loglog(infoTOSM1.iter,abs(infoTOSM1.gap)./max(infoTOSM1.obj,1),'-','color',cTOS1);
hold on;
loglog(infoTOSM2.iter,abs(infoTOSM2.gap)./max(infoTOSM2.obj,1),'--','color',cTOS2);
loglog(infoFW.iter , abs(infoFW.gap)./max(infoFW.obj,1),':','color',cFW);
xlabel('iteration','Interpreter','latex','fontsize',14)
ylabel('nonstationarity err.','Interpreter','latex','fontsize',14)

subplot(142)
loglog(infoTOSM1.iter,infoTOSM1.feas./sqrt(n),'-','color',cTOS1);
hold on;
loglog(infoTOSM2.iter,infoTOSM2.feas./sqrt(n),'--','color',cTOS2);
xlabel('iteration','Interpreter','latex','fontsize',14)
ylabel('infeasibility err.','Interpreter','latex','fontsize',14)

subplot(143)
loglog(infoTOSM1.time,abs(infoTOSM1.gap)./max(infoTOSM1.obj,1),'-','color',cTOS1);
hold on;
loglog(infoTOSM2.time,abs(infoTOSM2.gap)./max(infoTOSM2.obj,1),'--','color',cTOS2);
loglog(infoFW.time , abs(infoFW.gap)./max(infoFW.obj,1),':','color',cFW);
xlabel('time (sec)','Interpreter','latex','fontsize',14)
ylabel('nonstationarity err.','Interpreter','latex','fontsize',14)

subplot(144)
loglog(infoTOSM1.time,infoTOSM1.feas./sqrt(n),'-','color',cTOS1);
hold on;
loglog(infoTOSM2.time,infoTOSM2.feas./sqrt(n),'--','color',cTOS2);
xlabel('time (sec)','Interpreter','latex','fontsize',14)
ylabel('infeasibility err.','Interpreter','latex','fontsize',14)

for t = 1:4
    
    subplot(1,4,t)
    axis tight;
    ax = gca;
    set(findall(ax, 'Type', 'line'),'LineWidth',2.5);
    ax.FontSize = 14;
    ax.TickLabelInterpreter = 'latex';
    ax.TickDir = 'out';
    grid on; grid minor; grid minor;
    set(gca,'TickDir','out')
    set(gca,'LineWidth',0.75,'TickLength',[0.02 0.02]);
    ax.XTick = 10.^(-10:10);
    ax.YTick = 10.^(-100:2:100);
    ax.Box = 'on';
        
end

subplot(1,4,1)
xlim([1,1e5])
subplot(1,4,2)
xlim([1,1e5])

subplot(1,4,3)
XLIM = xlim();
subplot(1,4,4)
xlim(XLIM);
