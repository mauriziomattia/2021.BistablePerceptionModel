%% Set option and parameters for the analysis of the local set of data...
%
setOptionsAndDataset

% This is the best minima found for the logaritmic contrast mapping.
%
% Best min. 1 fit error = 0.1214 +- 0.0011 (*16/13 = 14.9% average error)
Params = [1.949  0.01767  15.21  33.38  152.2  2.34  32.1  -4.938  0.07088  0.08202  0.6555];

TIME_WINDOW = [-1 1];
Options.FitError.Repetitions = 50; % To improve the signal to noise ratio...


%% Plots an example inner dynamics for the model...
%
X1 = 1; %Population indexing
X2 = 2;
Y1 = 3;
Y2 = 4;

% Contrasts = DataSet.ContrastValues(1)+[0 0]; % C = 6.25 %
Contrasts = DataSet.ContrastValues(3)+[0 0];   % C = 25 %
% Contrasts = DataSet.ContrastValues(5)+[0 0]; % C = 100 %
   
Net = createEmptyNetwork(Options.Net.NumOfPools,Options.Net.ModulesPerPool);
Net = SetSimulationParams(Net,Params);
   
Net.SPParam.Theta(X1) = Params(11)*log(Contrasts(1) + Params(9)) + Params(10);
Net.SPParam.Theta(Y1) = Params(11)*log(Contrasts(2) + Params(9)) + Params(10);
   
% Integration
[~, History] = integrateNetworkDynamics(Net, Options.Net.Life);

% Transitions
TS = rebuildTimeSerie(Options.Net.Life, Options.Dt, History); % TS means TimeSeries
Pconsi = FindConsistentRegions(TS, Options.StateDetectionThreshold);
Triggers = TS.t(Pconsi(Pconsi(:,4)<-0.5,end)); % from Y to X...
Triggers = [Triggers TS.t(Pconsi(Pconsi(:,4)>0.5,end))]; % from X to Y...

% Plotting network history
figure

subplot(2,5,1:3)
hold('on');
Options.PopsToPlot = [X2 Y2];
Options.PopsColors = [1 0 0; 0 0 1];
Options.PopsLabels = {'r(t)','r''(t)'};
plotNetworkHistory(Net, History, Options);
set(gca,'XLim',[0 Options.Net.Life],'YLim',[0.0 1.0]);
set(gca,'Layer','top','Box','on','TickDir','out');
xlabel('Time, t [s]')
ylabel('Activity level')
title('Decision pools')

subplot(2,5,6:8)
hold('on');
Options.PopsToPlot = [X1 Y1];
Options.PopsColors = [1 0 0; 0 0 1];
Options.PopsLabels = {'e(t)','e''(t)'};
plot([Triggers; Triggers],repmat([0 1]',1,numel(Triggers)),'k--')
plotNetworkHistory(Net, History, Options);
set(gca,'XLim',[0 Options.Net.Life],'YLim',[0.0 1.0]);
set(gca,'Layer','top','Box','on','TickDir','out');
xlabel('Time, t [s]')
ylabel('Activity level')
title('Evidence pools')

subplot(2,5,4:5)
hold('on');
YRange = [-0.35 0.35];
X = TS.t(1:10:end); % Downsampling...
Y = TS.Y1dt(1:10:end)-TS.X1dt(1:10:end);
% patch([X(1) X X(end)],[0 Y 0], ([1 0 1]+3)/4,'EdgeColor','none')
plot([0 Options.Net.Life],[0 0],'k-')
plot(X, Y,'m-','LineWidth',0.75)
plot([Triggers; Triggers],repmat(YRange',1,numel(Triggers)),'k--')
set(gca,'XLim',[0 Options.Net.Life],'YLim',YRange);
set(gca,'Layer','top','Box','on','TickDir','out');
xlabel('Time (s)')
ylabel('e''(t) - e(t)')

subplot(2,5,9:10)
hold('on');
XRange = [-0.35 0.35];
YRange = [-1.2 1.2];
X = TS.Y1dt(1:10:end)-TS.X1dt(1:10:end); % Downsampling...
Y = TS.Y2dt(1:10:end)-TS.X2dt(1:10:end);
X = X + randn(size(X))*1/Net.SPParam.N(X1)*0.5;
Y = Y + randn(size(Y))*1/Net.SPParam.N(X2)*0.5;
plot(XRange,[0 0],'k-')
plot([0 0],YRange,'k-')
plot(X, Y, 'k.')
set(gca,'XLim',XRange,'YLim',YRange);
set(gca,'Layer','top','Box','on','TickDir','out');
xlabel('e''(t) - e(t)')
ylabel('r''(t) - r(t)')

FigSize = [10 6];
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 FigSize],'PaperSize',FigSize);
print('-dpdf','NetworkDynamicsPlot')