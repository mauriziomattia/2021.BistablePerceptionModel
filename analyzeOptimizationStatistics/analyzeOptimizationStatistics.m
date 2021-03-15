%% Set option and parameters for the analysis of the local set of data...
%
setOptionsAndDataset


%% Plot statistics of optimization process...
%
PlotBorder = 0.1;
FCORangeLog = 0;
FCOThreshold = 0.4;

load('OptimalModel.mat');
NumOfAllMinima = numel(FCostOpt);
if ~isempty(FCOThreshold)
   ndx = find(FCostOpt<=FCOThreshold);
   InitParams = InitParams(ndx,:);
   FParamsOpt = FParamsOpt(ndx,:);
   FCostOpt = FCostOpt(ndx);
   FEvalNum = FEvalNum(ndx);
end
[~,ndxMin] = min(FCostOpt);
Params = FParamsOpt(ndxMin,:);

% Remap parameter \alpha into u_0^e and \beta into w_{vis}
Gamma = FParamsOpt(:,9);
Beta = FParamsOpt(:,10);
Alpha = FParamsOpt(:,11);
FParamsOpt(:,10) = Alpha.*log((1+Gamma)./Gamma); % w_{vis}
FParamsOpt(:,11) = Alpha.*log(Gamma) + Beta;     % u_0^e
Options.Param.Labels{10} = 'w_{vis}';
Options.Param.Labels{11} = 'u_0^e';

% Select range of the explored degrees of freedom.
FCORange = [min(FCostOpt) max(FCostOpt)];
if FCORangeLog
   FCORange = log(FCORange);
   FCORange = FCORange + diff(FCORange)*PlotBorder*[-1 +1];
   FCORange = exp(FCORange);
else
   FCORange = FCORange + diff(FCORange)*PlotBorder*[-1 +1];
end

figure

for k = 1:size(FParamsOpt,2)
   subplot(4,3,k); 
   plot(FCostOpt,FParamsOpt(:,k),'.',FCostOpt(ndxMin),FParamsOpt(ndxMin,k),'or'); 
   xlabel('Fit Error'); 
   ylabel(Options.Param.Labels{k});
   title(sprintf('%s^{(opt)} = %.3g',Options.Param.Labels{k},FParamsOpt(ndxMin,k)));

   YRange = [min(FParamsOpt(:,k)) max(FParamsOpt(:,k))];
   if diff(YRange) > 0
      YRange = YRange + diff(YRange)*PlotBorder*[-1 +1];
   else
      YRange = YRange + [-0.5 +0.5];
   end
   set(gca,'XLim',FCORange,'YLim',YRange,'TickDir','out');
   if FCORangeLog
      set(gca,'XScale','log');
   end
end

k = k + 1;
subplot(4,3,k)
FCORange = [min(FCostOpt) max(FCostOpt)];
FCORange = FCORange + diff(FCORange)*PlotBorder*[-1 +1];
hist(FCostOpt,20);
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[1 1 1]*0.85,'EdgeColor','k')
set(gca,'XLim',FCORange,'TickDir','out');
xlabel('Fit Error'); 
ylabel('Samples');
title(['n = ' num2str(numel(FCostOpt)) ' of ' num2str(NumOfAllMinima)]);

FigSize = [6 8];
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 FigSize],'PaperSize',FigSize);
print('-dpdf','OptimizationStatistics')


%% PLOTS the correlation between fitted parameters and between cost
%  function...
%
if ~isempty(FCOThreshold)
   FE_THRESHOLD_FOR_PCA = FCOThreshold;
end
FEThreshold = FE_THRESHOLD_FOR_PCA; % to catch only the first two local minima...
ndxGlobMin = find(FCostOpt < FEThreshold);
[R,P]=corrcoef([FParamsOpt(ndxGlobMin,:) log(FCostOpt(ndxGlobMin)')]);
figure
subplot(1,2,1)
imagesc(R,[-1 1])
set(gca,'YDir','norm');
title('Correlation coefficient \rho')
for r = 1:size(R,1)
   for c = 1:size(R,2)
      text(c,r,num2str(round(100*R(r,c))/100),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',7);
   end
end
colormap(gca,gradedColormap([0 0 1],[1 0 0]))
colorbar
xlabel('Params')
ylabel('Params')
subplot(1,2,2)
imagesc(log(P)/log(10),[-3 0])
for r = 1:size(P,1)
   for c = 1:size(P,2)
      text(c,r,num2str(log(P(r,c))/log(10),1),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',7);
   end
end
set(gca,'YDir','norm');
title('log_{10} P')
colormap(gca,1-gray())
colorbar
xlabel('Params')
ylabel('Params')

FigSize = [10 4];
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 FigSize],'PaperSize',FigSize);
print('-dpdf','FittedParamsCorrelation')


%% PCA to evaluate the number of degrees of freedom.
%
ndxLowCF = find(FCostOpt<FEThreshold);
% ndxLowCF = find(FCostOpt<prctile(FCostOpt,50));
LogFCostOpt = log(FCostOpt(ndxLowCF))/log(10);
FCORange = [min(LogFCostOpt) max(LogFCostOpt)];

[pc,score,latent] = pca(FParamsOpt(ndxLowCF,:));
VarExplained = cumsum(latent)./sum(latent);
Coords = FParamsOpt(ndxLowCF,:) * pc;

figure

subplot(1,3,1)
plot(1:size(FParamsOpt,2),VarExplained*100,'k+-','LineWidth',0.75)
set(gca,'XLim',[0.5 size(FParamsOpt,2)+0.5],'YLim',[0 100])
set(gca,'Layer','top','Box','on','TickDir','out');
grid('on')
xlabel('PCs')
ylabel('Explained variance [%]')

subplot(1,3,2)
hold on
clear hlgn
hlgn(1) = stem((1:size(FParamsOpt,2))-0.2,pc(:,1)./std(FParamsOpt(ndxLowCF,:))','ro');
hlgn(2) = stem((1:size(FParamsOpt,2))-0.0,pc(:,2)./std(FParamsOpt(ndxLowCF,:))','gs');
hlgn(3) = stem((1:size(FParamsOpt,2))+0.2,pc(:,3)./std(FParamsOpt(ndxLowCF,:))','bd');
set(gca,'XLim',[0.5 size(FParamsOpt,2)+0.5])
set(gca,'Layer','top','Box','on','TickDir','out');
legend(hlgn,{'PC_1','PC_2','PC_3'})
legend('boxoff')
xlabel('Params')
ylabel('PC_x/St.Dev.Params')

subplot(1,3,3)
hold on
clr = jet();
for j = 1:numel(ndxLowCF)
   ndxClr = round((LogFCostOpt(j)-FCORange(1))/diff(FCORange)*(size(clr,1)-1))+1;
   plot3(Coords(j,1),Coords(j,2),Coords(j,3),'ko',...
        'MarkerFaceColor',clr(ndxClr,:),'MarkerSize',3.); 
end
[val,ndxMin]=min(LogFCostOpt);
plot(Coords(ndxMin,1),Coords(ndxMin,2),'ro'); 
set(gca,'Layer','top','Box','on','TickDir','out');
xlabel('PC_1')
ylabel('PC_2')
zlabel('PC_3')
grid('on')
view([0 90])
title(['Variance explained = ' num2str(VarExplained(3)*100,3) '%'])

FigSize = [9 4];
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 FigSize],'PaperSize',FigSize);
print('-dpdf','PCAofOptimizedParamSpace')


%% Probe the fit error surface around a parameter point. If it is a minimum
%  provides the standard error of the parameter found... 
%
FERepetitions = 30;
Options.Param.Labels{10} = '\beta';
Options.Param.Labels{11} = '\alpha';

tic
estimateParameterSEofFitErrorMinima(Params, FERepetitions, Options)
toc
