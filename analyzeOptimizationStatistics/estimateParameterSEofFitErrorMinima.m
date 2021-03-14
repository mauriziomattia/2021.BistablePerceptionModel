function estimateParameterSEofFitErrorMinima(Params, FERepetitions, Options)
%
%  estimateParameterSEofFitErrorMinima(Params, FERepetitions)
%
%     Params: Parameters of the fit error minimum to test.
%     FERepetitions: Number of Fit Error computation to carry out averages 
%                    and gradients.
%     Options: Information about data set and parameters needed to compute
%              fit error.

% global Options DataSet
% global Options

PlotResults = 1;


%% COMPUTES FIT ERROR and its GRADIENT for the best optima...
%
if isfield(Options.Param,'TestRange')
   DParam = Options.Param.TestRange;
else
   if isfield(Options.Param,'RelDerStep')
      % RelDerStep = [0.2 0.4 0.4 0.1 0.05 0.2 0.1 0.1 0.4 0.1];
      RelDerStep = Options.Param.RelDerStep;
   else
      RelDerStep = ones(size(Params))*0.1;
   end
   DParam = abs(Params.*RelDerStep);
end

ParamsSE = zeros(1,numel(Params));
FESurf = struct('Pvals',{},'FEs',{},'p',{},'psem',{},'xOpt',{},'yOpt',{},'Y',{},'DY',{});
FESurf(numel(Params)).psem = 0; % this is to set the array of structures...

TestedParams = setdiff(1:numel(Params),Options.ParamToFix.ndx);
% for p = TestedParams
%    [ParamsSE(p), FESurf(p)] = CoreLoop(p, Params, DParam, FERepetitions, Options);
%
% To use 'parfor', the following command has to be launched:
%   matlabpool local 4
%
parfor p = 1:numel(Params)
   if ismember(p,TestedParams)
      [ParamsSE(p), FESurf(p)] = CoreLoop(p, Params, DParam, FERepetitions, Options);
   end
end
fprintf('\n');

FEMean = 0;
FESEM = 0;
for p = TestedParams
   FEMean = FEMean + FESurf(p).p(3);
   FESEM = FESEM + FESurf(p).psem(3)^2;
end   
FEMean = FEMean/numel(TestedParams);
FESEM = sqrt(FESEM)/numel(TestedParams);


% Save computed data...
save('FEMinimumHighRes.mat', ...
   'Params','ParamsSE','FESurf','FEMean','FESEM',...
   'FERepetitions','DParam','TestedParams');


% Plot the FE surface estation...
if PlotResults
   figure

   for p = TestedParams
      subplot(4,3,p);
      plot(FESurf(p).Pvals,FESurf(p).FEs,'.',FESurf(p).Pvals,FESurf(p).Y,'r-',...
           FESurf(p).Pvals,FESurf(p).Y+FESurf(p).DY,'r:',FESurf(p).Pvals,FESurf(p).Y-FESurf(p).DY,'r:',...
           Params(p),polyval(FESurf(p).p,0),'rx',FESurf(p).xOpt,FESurf(p).yOpt,'go');
      [valMin,ndxMin] = min(FESurf(p).Y);
      text(Params(p),polyval(FESurf(p).p,0)*0.85,sprintf('(%.4g,%.4g)',FESurf(p).Pvals(ndxMin),valMin),...
           'Color',[0 0.5 0],'HorizontalAlignment','Center')
      xlabel(Options.Param.Labels{p});
      ylabel('Fit error');
      title(sprintf('%s = %.4g \\pm %.2g',Options.Param.Labels{p},Params(p),ParamsSE(p)))
      set(gca,'XLim',FESurf(p).Pvals([1 end]));
   end

   text(max(xlim())+diff(xlim()),mean(ylim())+diff(ylim())/6,...
      sprintf('Fit error = %.4g \\pm %.2g',FEMean,FESEM),...
      'HorizontalAlignment','center');
   text(max(xlim())+diff(xlim()),mean(ylim())-diff(ylim())/6,...
      sprintf('Repetitions = %d',FERepetitions),...
      'HorizontalAlignment','center');
   
   set(gcf, 'PaperUnit', 'inch', 'PaperPosition', [0 0 7 8])
   print('-deps2c','FEMinimumHighRes.eps');
end



%% This function is defined to allow distributed computing...
%
function [ParamsSE, FESurf] = CoreLoop(p, Params, DParam, FERepetitions, Options)

   fprintf('\n%d ',p);
   
   ParamsToTest = Params;
   Prange = [Params(p)-DParam(p) Params(p)+DParam(p)];
   Prange = [max([Prange(1) Options.Param.LowerBound(p)]) ...
             min([Prange(2) Options.Param.UpperBound(p)])];
   FESurf.Pvals = linspace(Prange(1),Prange(2),FERepetitions);
   FESurf.FEs = zeros(1,FERepetitions);
   for r = 1:FERepetitions
      if mod(r-1,10) == 0
         fprintf('.');
      end
      ParamsToTest(p) = FESurf.Pvals(r);
      FESurf.FEs(r) = FitErrorWithOptions(ParamsToTest,Options);
      pause(0.1); % To allow Ctrl-C...
   end
   [FESurf.p,S] = polyfit(FESurf.Pvals-Params(p),FESurf.FEs,2);
   FESurf.psem = sqrt(diag(inv(S.R)*inv(S.R')).*S.normr.^2./S.df)';
   FESurf.xOpt = Params(p)-FESurf.p(2)/(2*FESurf.p(1));
   FESurf.yOpt = FESurf.p(3)-FESurf.p(2)^2/(4*FESurf.p(1));
   
   % This is the range of values which do not increase Fit Error more than
   % the standard error of the minimum...
   ParamsSE = sqrt(FESurf.psem(3)/FESurf.p(1));

   % Fitted polinomial and its confidence interval for plot.
   [FESurf.Y,FESurf.DY]=polyval(FESurf.p,FESurf.Pvals-Params(p),S);

return
