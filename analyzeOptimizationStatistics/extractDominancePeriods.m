function [DD_T1 ,  DD_CV , DD_Gamma1 , Ti , Pdomi] = extractDominancePeriods(Pconsi)

%% Statistics of mixed dominance periods 
%
Mi  = Pconsi(Pconsi(:,4) == 0 ,3);            % list of mixed periods
DD_T1.Mix = mean( Mi );
DD_CV.Mix = std( Mi ) / mean( Mi );
DD_Gamma1.Mix = skewness( Mi ) / ( std( Mi ) / mean( Mi ) );


%% Statistics of dominance periods.
%
Ti  = Pconsi(Pconsi(:,4) ~= 0,3);            % list of dominance durations
Pdomi = Pconsi( Pconsi(:,4) ~= 0 , :);    %remove mixed states from statistic

Ti1 = Pdomi(Pdomi(:,4) == +1 ,3);            % list of dominance periods for 1
if ~isempty(Ti1)
   DD_T1.Percept1 = mean( Ti1 );
   DD_CV.Percept1 = std( Ti1 ) / mean( Ti1 );
   DD_Gamma1.Percept1 = skewness( Ti1 ) / ( std( Ti1 ) / mean( Ti1 ) );
else
   DD_T1.Percept1 = 0;
   DD_CV.Percept1 = 1;     % Assumes exponential distribution
   DD_Gamma1.Percept1 = 2; % Assumes exponential distribution
end

Ti2 = Pdomi(Pdomi(:,4) == -1 ,3);             % List of dominance periods for 2
if ~isempty(Ti2)
   DD_T1.Percept2 = mean( Ti2 );
   DD_CV.Percept2 = std( Ti2 ) / mean( Ti2 );
   DD_Gamma1.Percept2 = skewness( Ti2 ) / ( std( Ti2 ) / mean( Ti2 ) );
else
   DD_T1.Percept2 = 0;
   DD_CV.Percept2 = 1;     % Assumes exponential distribution
   DD_Gamma1.Percept2 = 2; % Assumes exponential distribution
end

return