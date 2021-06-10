%% Test setup for QAP experiment 2
% Replicates the experiment and generates Figure 2 in [YMS21].
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

% Change the directory to your data folder
dataDir = './data/';
QAPDATA = {'bur26a','bur26b','bur26c','bur26d','bur26e','bur26f','bur26g','bur26h',...
    'chr12a','chr12b','chr12c','chr15a','chr15b','chr15c','chr18a','chr18b','chr20a','chr20b','chr20c','chr22a','chr22b','chr25a',...
    'els19',...
    'esc16a','esc16b','esc16c','esc16d','esc16e','esc16f','esc16g','esc16h','esc16i','esc16j','esc32a','esc32b','esc32c','esc32d','esc32e','esc32g','esc32h','esc64a','esc128',...
    'had12','had14','had16','had18','had20',...
    'kra30a','kra30b','kra32',...
    'lipa20a','lipa20b','lipa30a','lipa30b','lipa40a','lipa40b','lipa50a','lipa50b','lipa60a','lipa60b','lipa70a','lipa70b','lipa80a','lipa80b','lipa90a','lipa90b',...
    'nug12','nug14','nug15','nug16a','nug16b','nug17','nug18','nug20','nug21','nug22','nug24','nug25','nug27','nug28','nug30',...
    'rou12','rou15','rou20',...
    'scr12','scr15','scr20',...
    'sko42','sko49','sko56','sko64','sko72','sko81','sko90','sko100a','sko100b','sko100c','sko100d','sko100e','sko100f',...
    'ste36a','ste36b','ste36c',...
    'tai12a','tai12b','tai15a','tai15b','tai17a','tai20a','tai20b','tai25a','tai25b','tai30a','tai30b','tai35a','tai35b','tai40a','tai40b','tai50a','tai50b','tai60a','tai60b','tai64c','tai80a','tai80b','tai100a','tai100b','tai150b','tai256c',...
    'tho30','tho40','tho150',...
    'wil50','wil100'};

T = length(QAPDATA);

errFW.assignmentErr = nan(T,1); 
errTOS2.assignmentErr = nan(T,1); 

for t = 1:T

fprintf(['\n\n****',QAPDATA{t},'****\n']);

%% Load data

% Change the directory to your data folder
[A,B,OPT,P] = qapread(['./data/',QAPDATA{t}]);

n = size(A,1);

%% Choose initial point
% Generate a random direction with Gaussian iid entries, then project
% (approximately) onto the Birkhoff polytope
X0 = ApproxProjBirkhoff(randn(n,n)./n,1e4);

%% Frank-Wolfe
fprintf('>>FW\n');
[~, infoFW] = FW(A,B,'maxit',1e7,'x0',X0,'tol',1e-5);
errFW.assignmentErr(t) = (infoFW.rddobj(end) - OPT)/max(OPT,1); 

%% TOS (Split2)
fprintf('\n--TOS2\n');
[~, infoTOS2] = TOSM2(A,B,'maxit',1e7,'x0',X0,'tol',1e-5);
errTOS2.assignmentErr(t) = (infoTOS2.rddobj(end) - OPT)/max(OPT,1); 

end

%% save 
if ~exist('./results','dir'), mkdir('results'); end
save('./results/Test2Results.mat','errFW','errTOS2','QAPDATA');

%% Plot Results

hfig = {};

for t = 1:3
    
    tt = ((t-1)*45+1):(t*45);
    tt(tt>T) = [];
    
    Xpart = categorical(QAPDATA(tt));
    Ypart = [errFW.assignmentErr(tt), errTOS2.assignmentErr(tt)];
    
    hfig{t} = figure;
    hfig{t}.Position = [0,0,1200,250];
    hfig{t}.Name = ['qapcost',num2str(t)];
    hfig{t}.NumberTitle = 'off';
    
    hBar = bar(Xpart,Ypart,'grouped');
    
    ax = gca;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.FontSize = 11;
    ax.TickDir = 'out';
    ax.YScale = 'log';
    ax.TickLabelInterpreter = 'latex';
    ylabel('assignment err.','Interpreter','latex','fontsize',15)
    
    hBar(1).FaceColor = [0.6,0.3,0.3];
    hBar(2).FaceColor = [1,0.9,0.9];

    if t == 1
        hl = legend('FW','TOS (Split 2)');
        hl.Interpreter = 'latex';
        hl.FontSize = 11;
        hl.Location = 'NorthWest';
    end
    
    ylim([1e-3,2])
    
end
