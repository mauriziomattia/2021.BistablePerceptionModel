function FE = FitErrorWithOptions(Params,Options)
%
%  FE = FitErrorWithOptions(Params,Options)
%
% Params(1) = 2; %Time-constant, integration layer
% Params(2) = 0.01; %Time-constant, decision layer
% Params(3) = 15; %WE_self
% Params(4) = 38.5; % WI_top
% Params(5) = 140; %WE_ff
% Params(6) = 1.4; %WI_fb
% Params(7) = 32; %WI_AS
% Params(8) = Params(6)/2; %Offset, decision layer
% Params(9) = ...;  % is \gamma
% Params(10) = ...; % is \beta
% Params(11) = ...; % is \alpha of the logarithmic mapping:
%                   % \alpha log(contrast + \gamma) + \beta

% This is the case of a square root mapping of the contrast...
% Params(9) = ...;  % is \beta
% Params(10) = ...; % is \alpha of the square root mapping:
%                   % \alpha \sqrt(contrast) + \beta

% Options.Net.NumOfPools = 4;
% Options.Net.ModulesPerPool = [30, 50] ;
% Options.Net.Life = 10000;
% Options.Dt = 0.001;
% Options.StateDetectionThreshold = 0.4;

% DataSet.Tmean
% DataSet.Tcv
% DataSet.Tgamma1
% DataSet.ContrastValues

% global Options
DataSet = Options.DataSet;

PrintDebug = 0;

X1 = 1; %Population indexing
X2 = 2;
Y1 = 3;
Y2 = 4;


Tmean = zeros([Options.FitError.Repetitions 2 size(DataSet.Tmean)]);
Tcv = zeros([Options.FitError.Repetitions 2 size(DataSet.Tcv)]);
Tratio = zeros([Options.FitError.Repetitions 2 size(DataSet.Tratio)]);
Tcc = zeros([Options.FitError.Repetitions 2 size(DataSet.Tcc)]);

for nr = 1:Options.FitError.Repetitions
   
   Net = createEmptyNetwork(Options.Net.NumOfPools,Options.Net.ModulesPerPool);
   Net = SetSimulationParams(Net,Params);
   
   for ni = 1:numel(DataSet.ContrastValues)
      for nc = 1:numel(DataSet.ContrastValues)
         if PrintDebug
            fprintf('(%d,%d) %g %g\n',ni, nc, DataSet.ContrastValues(ni), DataSet.ContrastValues(nc))
         end
         % Logarithmic mapping of the contrast
         Net.SPParam.Theta(X1) = Params(11)*log(DataSet.ContrastValues(ni) + Params(9)) + Params(10);
         Net.SPParam.Theta(Y1) = Params(11)*log(DataSet.ContrastValues(nc) + Params(9)) + Params(10);

%          % Square root mapping of the contrast
%          Net.SPParam.Theta(X1) = Params(10)*sqrt(DataSet.ContrastValues(ni)) + Params(9);
%          Net.SPParam.Theta(Y1) = Params(10)*sqrt(DataSet.ContrastValues(nc)) + Params(9);

         [~, History] = integrateNetworkDynamics(Net, Options.Net.Life);

         % Look for pure dominance periods and mistakes...
         TS = rebuildTimeSerie(Options.Net.Life, Options.Dt,History); % TS means TimeSeries
         Pconsi = FindConsistentRegions(TS, Options.StateDetectionThreshold);
         % Pconsi structured as follows [ start, end, duration, size of difference, Threshold Val, Z value at threshold ]
         
         % Dominance periods and mixed states statistics
         % Fields are: *.Mix, *.Percept1, *.Percept2
         EnoughDomTime = 0;
         if ~isempty(Pconsi) % At least one transition...
            [DD_T1, DD_CV, DD_Gamma1OverCV, ~, ~] = extractDominancePeriods(Pconsi);
            
            Ti1 = Pconsi(Pconsi(:,4)>0.5,3);       % Dominance periods of Percept1
            Ti1Start = Pconsi(Pconsi(:,4)>0.5,1);  % Start time of Percept1 dominance
            Ti2 = Pconsi(Pconsi(:,4)<-0.5,3);      % Dominance periods of Percept2
            Ti2Start = Pconsi(Pconsi(:,4)<-0.5,1); % Start time of Percept2 dominance
            
            if numel(Ti1)>=3 && numel(Ti2)>=3 % Three dominance period per percept...
               EnoughDomTime = 1;

               ndx1 = 1:numel(Ti1);
               ndx2 = 1:numel(Ti2);
               if Ti1Start(1) > Ti2Start(1) % Percept 1 has to be the first for 1 -> 2 CC
                  ndx2 = ndx2(2:end);
               end
               if numel(ndx1) < numel(ndx2)
                  ndx2 = ndx2(1:numel(ndx1));
               else
                  ndx1 = ndx1(1:numel(ndx2));
               end
               CC = corrcoef(Ti1(ndx1),Ti2(ndx2));
               Tcc(nr,2,ni,nc) = CC(1,2);
               
               ndx1 = 1:numel(Ti1);
               ndx2 = 1:numel(Ti2);
               if Ti2Start(1) > Ti1Start(1) % Percept 2 has to be the first for 2 -> 1 CC
                  ndx1 = ndx1(2:end);
               end
               if numel(ndx1) < numel(ndx2)
                  ndx2 = ndx2(1:numel(ndx1));
               else
                  ndx1 = ndx1(1:numel(ndx2));
               end
               CC = corrcoef(Ti2(ndx2),Ti1(ndx1));
               Tcc(nr,1,ni,nc) = CC(1,2);
            end
         end
         if ~EnoughDomTime
            if ~isempty(Pconsi)
               Tmean(nr,1,ni,nc) = DD_T1.Percept2;
               Tmean(nr,2,ni,nc) = DD_T1.Percept1;
            else
               Tmean(nr,1,ni,nc) = 0;
               Tmean(nr,2,ni,nc) = 0;
            end
            Tcv(nr,1,ni,nc) = 1; % Assumes exponential distribution
            Tcv(nr,2,ni,nc) = 1;
            Tratio(nr,1,ni,nc) = 2; % Assumes exponential distribution
            Tratio(nr,2,ni,nc) = 2;
            Tcc(nr,1,ni,nc) = 0; % Assumes no correlation
            Tcc(nr,2,ni,nc) = 0;
         else
            Tmean(nr,1,ni,nc) = DD_T1.Percept2;
            Tmean(nr,2,ni,nc) = DD_T1.Percept1;
            Tcv(nr,1,ni,nc) = DD_CV.Percept2;
            Tcv(nr,2,ni,nc) = DD_CV.Percept1;
            Tratio(nr,1,ni,nc) = DD_Gamma1OverCV.Percept2;
            Tratio(nr,2,ni,nc) = DD_Gamma1OverCV.Percept1;
         end
         
      end % for nc = ...
   end % for ni = ...
end % for nr = ...

% Mean over repetitions and symmetrization...
if Options.FitError.Repetitions > 1
   Tmean = squeeze(nanmean(Tmean));
   Tcv = squeeze(nanmean(Tcv));
   Tratio = squeeze(nanmean(Tratio));
   Tcc = squeeze(nanmean(Tcc));
else
   Tmean = squeeze(Tmean);
   Tcv = squeeze(Tcv);
   Tratio = squeeze(Tratio);
   Tcc = squeeze(Tcc);
end
Tmean = (squeeze(Tmean(1,:,:)) + squeeze(Tmean(2,:,:))')/2;
Tcv = (squeeze(Tcv(1,:,:)) + squeeze(Tcv(2,:,:))')/2;
Tratio = (squeeze(Tratio(1,:,:)) + squeeze(Tratio(2,:,:))')/2;
Tcc = (squeeze(Tcc(1,:,:)) + squeeze(Tcc(2,:,:))')/2;

Residuals = zeros(1,4);
Residuals(1) = mean(abs(Tmean(:) - DataSet.Tmean(:))/mean(DataSet.Tmean(:)));
Residuals(2) = mean(abs(Tcv(:) - DataSet.Tcv(:))/mean(DataSet.Tcv(:)));
Residuals(3) = mean(abs(Tratio(:) - DataSet.Tratio(:))/mean(DataSet.Tratio(:)));
% Residuals(4) = mean(abs(Tcc(:) - DataSet.Tcc(:))/mean(DataSet.Tcc(:)));
% This is to reduce the large variability in the estimate of the CC.
% in this way are assuming a flat CC surface.
Residuals(4) = abs(mean(Tcc(:)) - mean(DataSet.Tcc(:)))/mean(DataSet.Tcc(:)); 

Weights = ones(1,4);
if isfield(Options.FitError,'Weights')
   Weights = Options.FitError.Weights;
end
Residuals = Residuals.*Weights;

FE = mean(Residuals);
