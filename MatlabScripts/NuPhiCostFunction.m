function f = NuPhiCostFunction(Nu, Net)
%
%  f = NuPhiCostFunction(Nu, Net)
%

f = norm(Nu - Phi(Nu, Net));

OFFSET = 0.01;
f = log(f + OFFSET) - log(OFFSET);