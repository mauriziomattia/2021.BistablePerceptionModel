function out = PhiExp(mu, sigma2, tau, h, theta, tarp)
%
%  out = PhiExp(mu, sigma2, tau, h, theta, tarp)
%
%  Current-to-rate gain function of the LIF neuron under diffusion
%  approximation (Johannesma, 1968).
%
%   Copyright 2013-2019 Maurizio Mattia
%   Version: 1.1 - Jul. 19, 2019
%

%  mfpt = tau*sqrt(pi)*quad('exp(x.^2).*(1 + erf(x))', (h-mu)/sigma, (theta-mu)/sigma);

   sigma = sqrt(sigma2);
   out = zeros(size(mu));
   for n = 1:numel(mu)
      a = (h(n) - mu(n) * tau(n)) / (sigma(n) * sqrt(tau(n)));
      b = (theta(n) - mu(n) * tau(n)) / (sigma(n) * sqrt(tau(n)));
      if b > 5
         out(n) = 0;
      else
         if b < 0
%             fprintf('[PhiExp] DD regime\n')
            DDRegInt = @(w,a,b)exp(-w.^2).*(exp(2*b*w) - exp(2*a*w))./w;
            out(n) = 1 / (tau(n) * integral(@(x)DDRegInt(x,a,b), 0, Inf) + tarp(n));
         else
%             fprintf('[PhiExp] ND regime\n')
            out(n) = 1 / (tau(n) * sqrt(pi) * integral(@auxPhiExp, a, b) + tarp(n));
         end
      end
   end
end

% This function is what in LIFLibrary.nb is called Psi[w]. Actually, the
% condition for w < -3 is not exactly the same, but it should be
% equivalent.
function out = auxPhiExp(w)
   out = zeros(size(w));
   for i = 1:length(w)
      if w(i)>=-3
         out(i) = exp(w(i)^2)*(1 + erf(w(i)));
      else
         str = ['exp(-x.^2).*x./sqrt((' num2str(w(i)) ').^2 + x.^2)'];
         out(i) = 2/sqrt(pi)*quad(str, 0, 10);
      end
   end
end