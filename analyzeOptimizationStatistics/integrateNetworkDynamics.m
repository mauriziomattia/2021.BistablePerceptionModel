function [Net, History] = integrateNetworkDynamics(Net, Life)
%
%  [Net, History] = integrateNetworkDynamics(Net, Life)
%
%  Net
%  History.t
%  History.P
%  History.X
%


%% Constants and useful variables...
%
% BUFFER_DIM = round(2^16/Net.P);
BUFFER_DIM = round(2^12/Net.P);


%% Define useful empty arrays...
%
t = zeros(1,BUFFER_DIM);
P = zeros(1,BUFFER_DIM);
X = zeros(1,BUFFER_DIM);


%% Extraction of random transitions at unitary rates...
%
Durations = -log(1.-rand(Net.P,BUFFER_DIM));
Durations = Durations ./ repmat(Net.SPParam.N,1,BUFFER_DIM);
UpOrDown = rand(Net.P,BUFFER_DIM);
kDur = ones(Net.P,1);
% kUD = ones(Net.P,1);
   
   
%% Initialization of state variables...
%
Dur = zeros(Net.P,1);
for ne = 1:Net.P;
   t(ne) = Net.t;
   P(ne) = ne;
   X(ne) = Net.n(ne);
   Dur(ne) = Durations(ne,kDur(ne));
end
ne = Net.P;

SF = Net.n ./ Net.SPParam.N;


%% Numerical integration of the network dynamics...
%
while t(ne) < Life
   ne = ne + 1;
   if ne > numel(t) % Is additional memory needed to store the results?
      t = [t zeros(1,BUFFER_DIM)];
      P = [P zeros(1,BUFFER_DIM)];
      X = [X zeros(1,BUFFER_DIM)];
   end
   
   Field = Net.CParam.W*SF + Net.SPParam.Theta;
   TPup = Net.TPup(Field,Net.SPParam.Tau0);
   TPdown = Net.TPdown(Field,Net.SPParam.Tau0);
   TPRatio = TPup./TPdown;
%    NextT = Dur./(1-(1 - 1 ./ TPRatio).*SF)./TPup + Net.LastTrans;
   Den = (TPup-(TPup - TPdown).*SF);
   NextT = Inf(size(Dur));
   ndx = find(Den>0);
   NextT(ndx) = Dur(ndx)./Den(ndx) + Net.LastTrans(ndx);

   [t(ne),np] = min(NextT);
   if t(ne) <= Life
      P(ne) = np;
      if UpOrDown(np,kDur(np)) < SF(np)/((1-SF(np))*TPRatio(np)+SF(np))
         Net.n(np) = Net.n(np) - 1;
      else
         Net.n(np) = Net.n(np) + 1;
      end
      SF(np) = Net.n(np) ./ Net.SPParam.N(np);
      X(ne) = Net.n(np);
      Net.t = t(ne);
      Net.LastTrans(np) = t(ne);
      
      kDur(np) = kDur(np) + 1;
      
      if kDur(np) > BUFFER_DIM    % Are additional random numbers needed?
         Durations(np,:) = -log(1.-rand(1,BUFFER_DIM)) / Net.SPParam.N(np);
         UpOrDown(np,:) = rand(1,BUFFER_DIM);
         kDur(np) = 1;
      end
      Dur(np) = Durations(np,kDur(np));
   end
end
if t(ne) > Life
   ne = ne - 1;
end


%% Set the output of the function...
%
History.t = [t(1:ne) repmat(t(ne),1,Net.P)];
History.P = [P(1:ne) 1:Net.P];
History.X = [X(1:ne)./ Net.SPParam.N(P(1:ne))' (Net.n./Net.SPParam.N)'];
