function NuOut = EffectivePhi(NuIn, Ca, Net)
%
% NuOut = EffectivePhi(NuIn, Ca, Net)
%
%   Copyright 2013-2019 Maurizio Mattia 
%   Version: 1.2 - Oct. 22, 2019
%

global CFNet
CFNet = Net;

if numel(Net.ndxEFg) > 1
   error('[EffectivePhi] More than 1 foreground population');
end

if isfield(Net.SNParam,'GC')
   GC = Net.SNParam.GC(Net.ndxEFg);
   Net.SNParam.GC(Net.ndxEFg) = 0;
   IExt = Net.SNParam.IExt;
   Net.SNParam.IExt(Net.ndxEFg) = IExt(Net.ndxEFg) - GC*Ca;
else
   GC = 0;
end

ndxNotEFg = setdiff((1:Net.P)',Net.ndxEFg);

NuOut = zeros(Net.P,numel(NuIn));
Nu = Net.SNParam.Nu;
Options = optimset('MaxFunEvals',400*Net.P);
for nt = numel(NuIn):-1:1
   fprintf('%g -> ',NuIn(nt));
   
   CFNet.SNParam.Nu(Net.ndxEFg) = NuIn(nt);
   Nu(ndxNotEFg) = fminsearch(@CostFunc,Nu(ndxNotEFg),Options);
   Nu(Net.ndxEFg) = NuIn(nt);
   Nu = Phi(Nu,Net);
   
   fprintf('%g ',Nu)
   fprintf('\n')
   
   NuOut(:,nt) = Nu;
end % for nt = ...

if isfield(Net.SNParam,'GC')
   Net.SNParam.GC = GC;
   Net.SNParam.IExt = IExt;
end

end

%%
function out = CostFunc(NuNotEFg)
global CFNet

ndxNotEFg = setdiff((1:CFNet.P)',CFNet.ndxEFg);
cfNu = zeros(CFNet.P,1);
cfNu(ndxNotEFg) = NuNotEFg;
cfNu(CFNet.ndxEFg) = CFNet.SNParam.Nu(CFNet.ndxEFg);
NuOut = Phi(cfNu,CFNet);
out = norm(NuOut(ndxNotEFg)-NuNotEFg);

end
