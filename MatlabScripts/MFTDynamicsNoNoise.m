function [t, Nu, NuN] = MFTDynamicsNoNoise(Nu0, t0, tMax, Dt, Net, DNuExt)
%
%  [t, Nu, NuN] = MFTDynamicsNoNoise(Nu0, t0, tMax, Dt, Net[, DNuExt])
%

ONLINE_PLOT = 0;
PRINT_ELAPSED_TIME = 10; % Seconds between two consecutive prints of the
                         % numerical integration state.

%
% Parameter processing.
%
if size(Nu0,1)==1
   Nu0 = Nu0';
end
if exist('DNuExt') == 0
   DNuExt.Freq = zeros(Net.P, floor((tMax-t0)/Dt)+1);
end

%
% Support vector processing.
%
t = t0:Dt:tMax;
Nu = zeros(Net.P, floor((tMax-t0)/Dt)+1);
NuN = zeros(Net.P, floor((tMax-t0)/Dt)+1);
NuNMean = zeros(Net.P, floor((tMax-t0)/Dt)+1);
Nu(:,1) = Nu0;
NuN(:,1) = Nu0;
NuNMean(:,1) = Nu0;

OriginalNuExt = Net.SNParam.NuExt;
NuExt = repmat(Net.SNParam.NuExt, 1, length(t)) + DNuExt.Freq;


if ONLINE_PLOT
   %
   % Color plot initialization.
   %
   clrbase = 'rbmcygkk';
   clrndx = (mean(Net.CParam.J,1)<0)+1;
   inhndx = find(clrndx==2);
   for j = 1:(length(inhndx)-1)
      clrndx((inhndx(j)+1):length(clrndx)) = clrndx((inhndx(j)+1):length(clrndx)) + 2;
   end
   clr = clrbase(clrndx);
   
   %
   % Lines initialization.
   %
   for p = 1:Net.P
      h(p) = line('Color', clr(p), 'Marker','.', 'markersize',5, 'erase','none', 'xdata',[],'ydata',[]);
   end
   
end

%
% Dynamic integration.
%
StartPrintStatus = cputime;
ElapsedTime = PRINT_ELAPSED_TIME;
for n = 2:length(t)
   NuNMean(:,n) = (1 - Dt ./ Net.SNParam.TauD) .* NuNMean(:,n-1) + Dt * NuN(:,n-1) ./ Net.SNParam.TauD;
   
   Net.SNParam.NuExt = NuExt(:,n);
   Nu(:,n) = Phi(NuNMean(:,n), Net);
   NuN(:,n) = Nu(:,n);
   
   if ONLINE_PLOT
      for p = 1:Net.P
         set(h(p),'xdata',t(n),'ydata',Nu(p,n));
      end
      drawnow;
   end
   
   MuVsTime(:,n) = Mu(NuNMean(:,n), Net);
   Sigma2VsTime(:,n) = Sigma2(NuNMean(:,n), Net);
   
   %
   % Print status
   %
   if cputime > ElapsedTime + StartPrintStatus
      if ElapsedTime > 60
         if ElapsedTime > 3600
            str = sprintf('%g h', ElapsedTime/3600);
         else
            str = sprintf('%g min', ElapsedTime/60);
         end
      else
         str = sprintf('%g s', ElapsedTime);
      end
      ElapsedTime = ElapsedTime + PRINT_ELAPSED_TIME;
      disp(sprintf('Integrated time %g s, Processing status %3g %% (Real time elapsed: %s)', t(n)-t(1), 100*(t(n)-t(1))/(t(end)-t(1)), str));
   end

end

%
% Restore original values...
%
Net.SNParam.NuExt = OriginalNuExt;