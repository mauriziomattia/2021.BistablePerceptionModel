function TimeSeries = rebuildTimeSerie( Life, Dt , History )
%
% Rebuilds the time serie with an absolute time for convenience
%

X1 = 1; %Population indexing
X2 = 2;
Y1 = 3;
Y2 = 4;

ndxX1 = find(History.P == X1);
ndxY1 = find(History.P == Y1);
DX1 = diff(History.X(ndxX1));
DY1 = diff(History.X(ndxY1));
tDX1 = History.t(ndxX1(2:end));
tDY1 = History.t(ndxY1(2:end));
ndxX2 = find(History.P == X2);
ndxY2 = find(History.P == Y2);
DX2 = diff(History.X(ndxX2));
DY2 = diff(History.X(ndxY2));
tDX2 = History.t(ndxX2(2:end));
tDY2 = History.t(ndxY2(2:end));
t = 0:Dt:Life;

X1dt = zeros(size(t));
ndx = floor(tDX1/Dt)+1;
j = 0;
for k = ndx
   j = j + 1;
   X1dt(k) = X1dt(k) + DX1(j);
end
X1dt(1) = X1dt(1) + History.X(X1);
X1dt = cumsum(X1dt);

Y1dt = zeros(size(t));
ndx = floor(tDY1/Dt)+1;
j = 0;
for k = ndx
   j = j + 1;
   Y1dt(k) = Y1dt(k) + DY1(j);
end
Y1dt(1) = Y1dt(1) + History.X(Y1);
Y1dt = cumsum(Y1dt);


X2dt = zeros(size(t));
ndx = floor(tDX2/Dt)+1;
j = 0;
for k = ndx
   j = j + 1;
   X2dt(k) = X2dt(k) + DX2(j);
end
X2dt(1) = X2dt(1) + History.X(X2);
X2dt = cumsum(X2dt);

Y2dt = zeros(size(t));
ndx = floor(tDY2/Dt)+1;
j = 0;
for k = ndx
   j = j + 1;
   Y2dt(k) = Y2dt(k) + DY2(j);
end
Y2dt(1) = Y2dt(1) + History.X(Y2);
Y2dt = cumsum(Y2dt);


TimeSeries.t = t;
TimeSeries.X1dt = X1dt;
TimeSeries.X2dt = X2dt;
TimeSeries.Y1dt = Y1dt;
TimeSeries.Y2dt = Y2dt;



end

