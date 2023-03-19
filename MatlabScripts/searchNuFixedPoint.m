function NuFP = searchNuFixedPoint(Net, NuGuess, PreMFDynamics, PlotMFDynamics)
%
%  NuFP = searchNuFixedPoint(Net, NuGuess[, PreMFDynamics[, PlotMFDynamics]])
%
%  Searches for a fixed point of the module dynamics. Two modes are allowed:
%  with or without a pre-search of the fixed point using the relaxation 
%  properties of the system, then following the mean-field (MF) dynamics for a
%  small time period starting from <NuGuess> initial condition. The last (and 
%  the only, for the first mode) step in search is searching a minimum around 
%  <NuGuess> or the resulting activity from the dynamics.
%  Pays attention that if you use the second mode, possibly only the stable
%  fixed points are found. In addition if the dynamics does not relaxes (perhaps 
%  showing a limit cycle) the guess provided to the minimization algorithm
%  can be worst than the one suggested as external parameter.
%
%  <Net>: the network structure for which the fixed point is searched.
%  <NuGuess>: an array of emission rates of the populations of the module,
%     used as initial condition for the MF dynamics or for the minimization
%     algorithm, respectively in the second and first mode of search.
%  <PreMFDynamics> (OPTIONAL, 0): a flag to switch from the first mode (0)
%     and the second one (1).
%  <PlotMFDynamics> (OPTIONAL, 0): if in the second mode of search, it is
%     a flag to show a plot of the MF dynamics (1) or not (0). Of course if
%     <PreMFDynamics> == 0, the value of this parameter is neglected.
%
%  <NuFP>: is the activity of the network in the found fixed point (if it is
%     not a fixed point several warnings appear).
%  

COST_FUNCTION_TOLERANCE = 5e-4;
MAX_FMINSEARCH_REPETITION = 3;
MAX_DYNAMICS_REPETITION = 3;
MIN_DNUDT = 1e-2;

%
% Parameters analysis...
%
if exist('PreMFDynamics') == 0
   PreMFDynamics = 0;
else
   if exist('PlotMFDynamics') == 0
      PlotMFDynamics = 0;
   end
end

%
% If required, searches for the closest stable activity to <NuGuess>...
%
if PreMFDynamics
   Life = 2*max(Net.SNParam.TauD);
   Dt = min(Net.SNParam.TauD)/10;
   [t, Nu, NuN] = MFTDynamicsNoNoise(NuGuess, 0, Life, Dt, Net); 
   if PlotMFDynamics
      figure;
      plot(t,Nu);
      xlabel('time (seconds)');
      ylabel('\nu (Hz)');
   end
   Nu0 = Nu(:,end);
else
   Nu0 = NuGuess;
end

for n = 1:MAX_FMINSEARCH_REPETITION
   [NuFP,fval] = fminsearch(@NuPhiCostFunction, Nu0, [], Net);
   if fval > COST_FUNCTION_TOLERANCE
      Nu0 = NuFP;
      disp(['[searchNuFixedPoint] Warning: the found minimum could not be a fixed point, fval = ' num2str(fval)]);
   else
      break
   end
end

if fval > COST_FUNCTION_TOLERANCE
   Life = 2*max(Net.SNParam.TauD);
   Dt = min(Net.SNParam.TauD)/10;
   n = 1;
   while n <= MAX_DYNAMICS_REPETITION
      [t, Nu, NuN] = MFTDynamicsNoNoise(Nu0, 0, Life, Dt, Net);
      Nu0 = Nu(:,end);

      % Evaluate the stability of the end-point...
      MaxDNuDt = max(abs((Nu(:,end) - Nu(:,end-1)) / Dt));
      if MaxDNuDt < MIN_DNUDT
         break
      end
   end
   
   for n = 1:MAX_FMINSEARCH_REPETITION
      [NuFP,fval] = fminsearch(@NuPhiCostFunction, Nu0, [], Net);
      if fval > COST_FUNCTION_TOLERANCE
         Nu0 = NuFP;
         disp(['[searchNuFixedPoint] Warning: the found minimum could not be a fixed point, fval = ' num2str(fval)]);
      else
         break
      end
   end
end