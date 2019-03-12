% Maximum likelihood estimation of Weibull parameters for histogram of data.
% To evaluate shape parameter the following equation is solved using
% Newton-Raphson method:
% 1 + (1/shape)*log(shape) + (1/shape)*diGamma(1/shape) - 
% - (1/shape)* sum( |x_i^shape/sum(x_i^shape)|*
%   log(|n*x_i^shape/sum(x_i^shape)|) ) = 0,
% where n is size of input data.
% Scale parameter is equal:
% (sum(x_i^shape)/n)^(1/shape).
function [scale shape] = weibullMleHist(ax, h)

eps = 0.01; % precision of Newton-Raphson method's solution
shape = 0.1; %initial value of gamma parameter
h = h./sum(h);

shape_next = shape - weibullNewtonHist(shape, ax, h);
nIteration = 0; % number of iterations
while abs(shape_next - shape) > eps 
    % Terminating conditions when not converging
    if (isnan(shape_next) || isinf(shape_next) || shape_next > 20 || nIteration > 30)
        break;
    end
    if (shape_next <= 0)
        shape_next = 0.000001;
        break;
    end
    
    shape = shape_next;
    shape_next = shape - weibullNewtonHist(shape, ax, h);    
    
    nIteration = nIteration + 1;
end
shape = shape_next;
scale = (sum((ax.^shape).*h))^(1/shape);

function out = weibullNewtonHist(g, x, h)
% g: current value of shape parameter
% x: mean points on histograms axis
% h: histogram values
format long;

x_g = x.^g;
sum_x_g = sum(x_g.*h);
x_i = (x_g)./sum_x_g;
ln_x_i = log(x_i);

lambda = x_g.*(log(x).*sum_x_g - sum(h.*x_g.*log(x)))./(sum_x_g^2);

f = 1 + sum(ln_x_i.*h) - sum(x_i.*ln_x_i.*h);
f_prime = sum(lambda.*h.*(sum_x_g./x_g - ln_x_i - 1));

out = f./f_prime;
