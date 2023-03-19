%% Loading and setting of the network params of the bistable module.
%  Before running this Matlab script be sure to have in the path the folder
%  MatlabScripts

NuExtEFg = 3.06;  % With this external input the cortical module is bistable.
% NuExtEFg = 3.00;
% NuExtEFg = 3.0;
% NuExtEFg = 2.94;

NuE = 3.0; % Guess for the firing rate of excitatory and inhibitory neurons
NuI = 6.0;

Net = loadPerseusParams('modules_CM.ini', 'connectivity_CM.ini','LIFCA');

NuGuess(Net.ndxE) = NuE;
NuGuess(Net.ndxI) = NuI;
Net.SNParam.Nu = searchNuFixedPoint(Net, NuGuess, 1);

Net.ndxEFg = Net.ndxE(1);
Net.SNParam.NuExt(Net.ndxEFg) = NuExtEFg;

%% Compute the effective gain function for the foreground excitatory population.
NuInRange = [0 60];
NuInSamples = 60;
NuIn = exp(linspace(log(NuInRange(1)+1),log(NuInRange(2)+1),NuInSamples))-1;
NuOut = EffectivePhi(NuIn, 0, Net);

fnForce = csapi(NuIn,NuOut(1,:)-NuIn);
fnEnergy = fnint(fncmb(fnForce,-1));
fnEffPhi = csapi(NuIn,NuOut(Net.ndxEFg,:));
fnNuEBg = csapi(NuIn,mean(NuOut(setdiff(Net.ndxE,Net.ndxEFg),:)));
fnNuI = csapi(NuIn,NuOut(Net.ndxI,:));
FPNu = fnzeros(fnForce);
FPNu = FPNu(1,:);
FPStability = fnval(fnder(fnForce),FPNu) < 0; % 1 = stable; 0 = unstable.

%%
%
EnergyRange = [-10 40];

figure

subplot(4,1,1)
% plot(NuIn,NuIn*0,'k--',NuIn,NuOut(1,:)-NuIn,'r','LineWidth',1.)
plot(NuIn,NuIn*0,'k--')
hold('on')
fnplt(fnForce,'r-',1.)
for k = 1:numel(FPNu)
   if FPStability(k)
      plot(FPNu(k),0,'ko','LineWidth',1.)
      text(FPNu(k),diff(ylim())*0.1, sprintf(' %.3g Hz',FPNu(k)),'Color','k')
   else
      plot(FPNu(k),0,'kx','LineWidth',1.)
      text(FPNu(k),-diff(ylim())*0.1, sprintf(' %.3g Hz',FPNu(k)),'Color','k')
   end
end
grid('on')
set(gca,'XLim',NuIn([1 end]))
set(gca,'Layer','top','Box','on','TickDir','out')
ylabel('\Phi_{eff}(\nu_E^{(in)}) - \nu_E^{(in)} [1/s]')
% title(sprintf('J_{EE} = %g mV, Recurrency = %g%%',JEE,RecurFrac*100))

subplot(4,1,2)
hold('on')
fnplt(fnEnergy,'r-',1.)
for k = 1:numel(FPNu)
   if FPStability(k)
      plot(FPNu(k),fnval(fnEnergy,FPNu(k)),'ko','LineWidth',1.)
   else
      plot(FPNu(k),fnval(fnEnergy,FPNu(k)),'kx','LineWidth',1.)
   end
end
set(gca,'XLim',NuIn([1 end]),'YLim',EnergyRange)
set(gca,'Layer','top','Box','on','TickDir','out')
ylabel('Energy [1/s^2]')

subplot(4,1,3:4)
% plot(NuIn,NuIn,'k--',NuIn,NuOut(1,:),'r','LineWidth',1.)
plot(NuIn,NuIn,'k--')
hold('on')
fnplt(fnEffPhi,'r-',1.);
for k = 1:numel(FPNu)
   if FPStability(k)
      plot(FPNu(k),FPNu(k),'ko','LineWidth',1.)
   else
      plot(FPNu(k),FPNu(k),'kx','LineWidth',1.)
   end
end
fnplt(fnNuEBg,'m-',1.);
fnplt(fnNuI,'b-',1.);

hlgn(1) = findobj(gca,'Type','line','-and','color','r');
hlgn(2) = findobj(gca,'Type','line','-and','color','m');
hlgn(3) = findobj(gca,'Type','line','-and','color','b');
legend(hlgn,{'E','E_{bg}','I'},'Location','NorthWest')
legend('boxoff')

set(gca,'XLim',NuIn([1 end]))
set(gca,'Layer','top','Box','on','TickDir','out')
xlabel('\nu_E^{(in)} [1/s]')
ylabel('\Phi_{eff}(\nu_E^{(in)}) [1/s]')

fp = get(gcf,'Position');
DH = round(fp(4)*0.5);
set(gcf,'Position', [fp(1) fp(2)-DH fp(3) fp(4)+DH]);

FigSize = [4 6];
set(gcf, 'PaperUnits', 'inch', 'PaperSize',FigSize, 'PaperPosition', [0 0 FigSize]);
print('-dpdf', '-painters', 'EffGainFuncAndEnerLand');
