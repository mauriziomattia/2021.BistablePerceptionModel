function plotNetworkHistory(Net, History, Options)
%
%  plotNetworkHistory(Net, History[, Options])
%

ndxPops = 1:Net.P; 
clrPops = jet(numel(ndxPops));
labPops = [];

% figure
% hold on

if exist('Options', 'var')
   if isfield(Options,'PopsToPlot')
      ndxPops = Options.PopsToPlot;
   end
   if isfield(Options,'PopsColors')
      clrPops = Options.PopsColors;
   end
   if isfield(Options,'PopsLabels')
      labPops = Options.PopsLabels;
   end
end

nc = 0;
for np = ndxPops
   ndx = find(History.P == np);
   t = History.t(ndx);
   X = History.X(ndx);
   plt.X = [t(2:end-1); t(2:end-1)];
   plt.Y = [X(1:end-2); X(2:end-1)];
   plt.X = [t(1); plt.X(:); t(end)];
   plt.Y = [X(1); plt.Y(:); X(end)];
   nc = nc + 1;
   hlgn(nc) = plot(plt.X,plt.Y,'-','Color',clrPops(nc,:),'LineWidth',0.75);
end
xlabel('Time (s)')
ylabel('Activity X')
if isempty(labPops) == 0
   legend(hlgn,labPops,'Location','NorthEastOutside')
end