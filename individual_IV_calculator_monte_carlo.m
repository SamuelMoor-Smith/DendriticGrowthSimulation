function [I] = individual_IV_calculator_monte_carlo(x, y, L, V, lambda, Rt,steps)

newL = 200e-6;

xi = x';
yi = y';
  
t2 = steps;
p = 50;

sz = size(x);
% t1 = sz(1);
% n = sz(2);

threshold_dist = 50e-6;

min_distance_to_end = min(newL - x, [], 'all');
if min_distance_to_end > threshold_dist
    % ps(i) = p;
    I = 0;
    return;
end
    
% Start with electrons at random particles in device

curr_indexes = randperm(100, p);
currx = xi(curr_indexes);
curry = yi(curr_indexes);
p = length(currx);
    
for j = 1:t2
    
    dx = repmat(xi,p,1) - currx';
    dy = repmat(yi,p,1) - curry';
    d = sqrt(dx.^2 + dy.^2);
    
    positive = dx >= 0;
    
    voltage = V * dx / newL;
    % numerator = 1 + 2000*max(voltage, 0);
    sigmoid = 1./(1+exp(-voltage));
    % sigmoid = V * ones(size(dx));    % Default to V
    % sigmoid(dx == 0) = 1;            % Set to 1 where dx == 0
    resist = Rt*exp(d/lambda);
    
    dist_weights = positive .* sigmoid ./ resist;
    dist_weights(isnan(dist_weights)) = 0;
    
    voltage_end = V.*(newL-currx')./newL;
    sigmoid_end = 1./(1+exp(-voltage_end));
    % sigmoid_end = V;
    resist_end = Rt*exp((newL-currx')/lambda);
    
    to_end = (sigmoid_end./resist_end) .* ones(p,1);
    to_end(isnan(to_end)) = 1;
    
    to_begin = zeros(p,1);

    tot_weights = [dist_weights to_begin to_end];
    
    [curr_indexes,currx,curry]=individual_get_next_index(tot_weights,xi,yi,currx,curry,curr_indexes);
    
end

I = sum(isnan(currx));