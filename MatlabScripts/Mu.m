function MuOut = Mu(NuIn, Net)
%
% MuOut = Mu(NuIn, Net)
%
%   Copyright 2008-2019 Maurizio Mattia 
%   Version: 1.1 - Jul. 10, 2019
%   Version: 1.0 - Oct. 18, 2008
%

% Verifies input params 
MuOut = [];
if size(NuIn,1) ~= Net.P
   disp('Error [Mu]: NuIn rows different from Net.P');
   return
end
if ~isfield(Net.SNParam, 'IExt')
   Net.SNParam.IExt = zeros(Net.P,1);
end

% Works out Mu...
CJ = Net.CParam.c .* Net.CParam.J;

if Net.SNParam.Type == Net.Constants.NT_LIFCA || ...
   Net.SNParam.Type == Net.Constants.NT_VIFCA
   MuOut = CJ * (Net.SNParam.N .* NuIn) + Net.SNParam.IExt + ...
           Net.SNParam.NExt .* Net.SNParam.JExt .* Net.SNParam.NuExt - ...
           Net.SNParam.AlphaC .* Net.SNParam.TauC .* Net.SNParam.GC .*  NuIn;
else
   MuOut = CJ * (Net.SNParam.N .* NuIn) + Net.SNParam.IExt + ...
           Net.SNParam.NExt .* Net.SNParam.JExt .* Net.SNParam.NuExt;
end