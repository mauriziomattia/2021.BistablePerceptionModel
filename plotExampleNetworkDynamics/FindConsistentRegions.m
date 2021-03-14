function Pconsi = FindConsistentRegions(TS, Scrit)
%
% find candidate transitions from difference between decision levels (DD)
%

T = TS.t;
E1all = TS.X1dt;
E2all = TS.Y1dt;
D1all = TS.X2dt;   % row vector!
D2all = TS.Y2dt;


DDsign   = sign(D1all-D2all);   % sign of difference between D1 and D2 [+1,-1]
DDsize   = abs(D1all-D2all) > Scrit;   % size above criterion or not?  [0, 1]
DDsize   = DDsign .* DDsize;    % size above and sign                [+1, 0, -1]
% +1 = D1 larger
%  0 = comparable
% -1 = D2 larger

idxT = find( diff(DDsize) ~= 0 );  % list of indeces at which
% successive size/sign values are DIFFERENT
% idxT(1) indicates first difference ...
% DDsize( idxT(1) ) ~= DDsize( idxT(1)+1 )

%% assemble list of consistent regions:
%  Pconsi(n, 4) = [ start, end, duration, size of difference, Threshold Val, Z value at threshold, IDX_START, IDX_END ]

if ~isempty(idxT)
   Pconsi(1,:) = [0, T(idxT(1)), T(idxT(1)), DDsize(idxT(1)), E1all(idxT(1))-E2all(idxT(1)) , (E1all(idxT(1))+E2all(idxT(1)) )/2 , 1 , idxT(1)];  % up to and including idxT(1)
   % all DDsize are the same
   nP = 1;
   
   for i=2:length(idxT)  % go from second difference idxT(2) to last one idxT(end)
      
      if idxT(i)-idxT(i-1) == 1   % no consistent region
         continue;
      else                        % another consistent region, at least 2 longidxT(end)+1 nall
         
         nP = nP + 1;
         Pconsi(nP,:) = [T(idxT(i-1)), T(idxT(i)), T(idxT(i))-T(idxT(i-1)), DDsize(idxT(i)), E1all(idxT(i))-E2all(idxT(i)),  (E1all(idxT(i))+E2all(idxT(i)) )/2 , idxT(i-1)+1, idxT(i)];
      end
   end
   
   if numel(T) - idxT(end) ~= 1        % last difference well before last data point
      nP = nP + 1;
      Pconsi(nP,:) = [T(idxT(end)), T(end), T(end)-T(idxT(end)), DDsize(end), E1all(idxT(end))-E2all(idxT(end)),  (E1all(idxT(end))+E2all(idxT(end)) )/2 , idxT(end)+1, numel(T)];
   end
   
   
   % for i=1:nP                      % check consistency of DDsize
   %     kk = unique( [ DDsize( Pconsi(i,1):Pconsi(i,2) ), Pconsi(i,4) ] );
   %     if length(kk) ~= 1
   %       %  [num2str( i ) 'dominance period not unique!!!']
   %     end
   % end
   
else
   nP = 0;
   Pconsi = [];
end % ~isempty(idxT)
