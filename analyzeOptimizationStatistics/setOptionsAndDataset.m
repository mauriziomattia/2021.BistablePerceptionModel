% setOptionsAndDataset

global Options DataSet

Options.Net.NumOfPools = 4;
Options.Net.ModulesPerPool = [25, 25] ;
% Options.Net.Life = 200;
Options.Net.Life = 120;
Options.Dt = 0.001;
Options.StateDetectionThreshold = 0.4;

Options.ParamToFix.ndx = [];
Options.ParamToFix.val = [];
% Options.Param.Labels = {'\tau_{int}','\tau_{dec}','W_{E,dec}','W_{I,dec}',...
%                         'W_{E,ff}','W_{I,fb}','W_{I,ff}','\theta_{dec}',...
%                         '\gamma','\beta','\alpha'};
Options.Param.Labels = {'1/\nu_e','1/\nu_r','w_{coop}','w_{comp}',...
                        'w_{exc}','w_{supp}','w_{inh}','u_0^r',...
                        '\gamma','\beta','\alpha'};
% Options.Param.RelDerStep = [0.2   0.6   0.6   0.1   0.05  0.2   0.2   0.2   0.4   0.4   0.4];
Options.Param.TestRange  = [0.4   2.5   2.0   9.0   20.0   1.0  6.0   2.0   0.06   0.2   0.075]; % If this is defined .RelDerStep is not used.
Options.Param.UpperBound = [Inf   Inf   Inf   Inf   Inf   Inf   Inf   Inf   Inf   Inf   Inf];
Options.Param.LowerBound = [0.010 0.010 0.0   0.0   0.0   0.0   0.0   -Inf  0.001 -Inf  0];

Options.SDParam.UpperBound = [10   10   30  80  205 5  80  -5   0.1   +1 1.01];
Options.SDParam.LowerBound = [0.01 0.01 5   5   105 0  5   -30  0.001 -1 0.01];

Options.FitError.Repetitions = 10;
Options.FitError.Weights = [1 1 1 1/4];

% Load the dataset, after symmetrization and smoothing...
% 'ContrastValues', 'nMU', 'nSIGMA', 'nCV', 'nRATIO', 'fCCL2R', 'nCCL2R', 'nCCR2L'
DS = load('ObservablesToFit.mat');
DataSet.ContrastValues = DS.ContrastValues;
DataSet.Tmean = DS.nMU;
DataSet.Tcv = DS.nCV;
DataSet.Tratio = DS.nRATIO;
DataSet.Tcc = DS.fCCL2R;
clear('DS');

Options.DataSet = DataSet;
