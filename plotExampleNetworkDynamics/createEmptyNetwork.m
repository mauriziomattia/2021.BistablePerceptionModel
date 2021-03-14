function Net = createEmptyNetwork(NumOfPools,ModulePerPools,ModuleType)
%
% Net = createEmptyNetwork(NumOfPools,ModulePerPools[,ModuleType])
%
% - A cortical module, or simply a 'module', is a stochastic switch in which
%   only two preferred states are allowed (0 and 1, Down and Up, low and
%   high firing rates, respectively). Transition rates of state flips
%   completely determine the Poissonian statistics of the inter-event
%   intervals. The NuPlus and NuMinus rates determining the average
%   frequencies of  Down-to-Up and Up-to-Down transitions, respectively,
%   are state dependent, because function of the number of active switches
%   (those in Up state) contributing to the input received.
% - A 'pool' of modules is a homogeneous network of interacting modules
%   with same properties and input statistics.
% - A network of pools, is a set of interacting pools of modules with
%   possibly heterogeneous features.


if ~exist('ModuleType','var')
   ModuleType = 'EXP';
end

%% Set the module type...
%
Net.Constants.MT_EXP = 0;
Net.Constants.MT_EXP_KICK = 1;
switch upper(ModuleType)
   case 'EXP'
      mtype = Net.Constants.MT_EXP;
      tpup = @(field,tau0)exp(field/2)./(2*tau0);
      tpdown = @(field,tau0)exp(-field/2)./(2*tau0);
   case 'EXP_KICK'
      mtype = Net.Constants.MT_EXP_KICK;
      tpup = @(field,tau0)exp(-field/2)./(2*tau0);
      tpdown = @(field,tau0)exp(+field/2)./(2*tau0);
   otherwise
      disp('[setNetworkParams] Error: Unknown module type. Accepted values: EXP.');
%       disp('[setNetworkParams] Error: Unknown module type. Accepted values: EXP, EXP_SFA.');
      Net = [];
      return;
end


%% Set parameters for single pools of modules...
%
Net.P = NumOfPools; % Number of module pools in the network.
Net.ModType = mtype;     % Type of modules in the pools...
Net.TPup = tpup;      % Down-to-Up transition rate...
Net.TPdown = tpdown;  % Up-to-Down transition rate...

Net.SPParam.N(1) = ModulePerPools(1);
Net.SPParam.N(3) = ModulePerPools(1);
Net.SPParam.N(2) = ModulePerPools(2);
Net.SPParam.N(4) = ModulePerPools(2);
Net.SPParam.N = Net.SPParam.N';
%Net.SPParam.N = repmat(ModulePerPools,Net.P,1);   % Number of stochastic modules in the pools...
Net.SPParam.Tau0 = ones(Net.P,1); % Time scale of the pool in absence of input...
Net.SPParam.Theta = zeros(Net.P,1); % External input...

Net.SPParam.NUP = zeros(Net.P,1); % NuPlus0
Net.SPParam.NUM = zeros(Net.P,1); % NuMinus0

% if mtype == Net.Constants.MT_EXP_SFA
% % TO DO...
% end


%% Set parameters for pool coupling (postsyn,presyn)...
%
Net.CParam.W = zeros(Net.P,Net.P);
% Net.CParam.Delay(postsyn,presyn) = num(7);


%% Set initial condition for the network state...
%
Net.n = zeros(Net.P,1); % Number ot active switches.
Net.t = 0;              % Simulation time.
Net.LastTrans = zeros(Net.P,1); % This is the time when last transition 
                                % occurred in each pool. It is needed for 
                                % the event-driven numerical integration of 
                                % the stochastic dynamics.
