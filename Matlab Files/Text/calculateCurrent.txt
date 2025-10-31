function [I] = calculateCurrent(x, y, L, V, lambda, Rt, steps, num_e)

n = length(x); % number of particles
threshold_dist = 50;

% If the minimum distance to the end is greater than the threshold distance,
% then set current to 0 and return to avoid unnecessary calculations
min_distance_to_end = min(L - x, [], 'all');
if min_distance_to_end > threshold_dist
    I = 0;
    return;
end
    
% Start with electrons at random particles in device
% curr_indexes = randperm(100, p);

% start with num_e electrons at left electrode
LEFT_ELECTRODE = n + 1;
curr_indexes = LEFT_ELECTRODE * ones(num_e, 1);; % left electrode is index n+1

% Get the current positions of the electrons
currx = zeros(num_e, 1);; % start at x=0';
curry = nan(num_e, 1);; % start at y = nan' (y doesn't matter at electrode);

% in case num_e changes - not happening anymore
% num_e = length(currx);
    
for step = 1:steps

    dx = repmat(x',num_e,1) - currx;
    dy = repmat(y',num_e,1) - curry;

    % if dy is nan (at electrode), set to 0 since y distance doesn't matter
    dy(isnan(dy)) = 0;
    d = sqrt(dx.^2 + dy.^2);
    
    positive = dx >= 0;
    
    % voltage = V * dx / L;
    % numerator = 1 + 2000*max(voltage, 0);
    % sigmoid = 1./(1+exp(-voltage));
    % sigmoid = V * ones(size(dx));    % Default to V
    % sigmoid(dx == 0) = 1;            % Set to 1 where dx == 0
    resist = Rt*exp(d/lambda);
    
    dist_weights = positive .* 1 ./ resist;
    % dist_weights(isnan(dist_weights)) = 0;
    
    % voltage_end = V.*(L-currx')./L;
    % sigmoid_end = 1./(1+exp(-voltage_end));
    sigmoid_end = ones(num_e,1);    % Default to V
    % sigmoid_end = V;
    resist_end = Rt*exp((L-currx)/lambda);
    
    to_end = (sigmoid_end./resist_end) .* ones(num_e,1);
    to_end(curr_indexes==n+2) = 1000; % if already at end, always stay there
    
    to_begin = zeros(num_e,1);

    tot_weights = [dist_weights to_begin to_end];
    
    [curr_indexes,currx,curry]=getNextIndex(tot_weights,x',y',currx,curry,curr_indexes, L);
    
end

I = sum(currx==L);