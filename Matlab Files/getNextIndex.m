function [next_indices,nextx,nexty] = getNextIndex(tot_weights,x,y,currx,curry,curr_indexes, L)

num_e = size(tot_weights,1); % total number of electrons
n = size(tot_weights,2)-2; % number of possible particle targets (minus 2 electrodes)
LEFT_ELECTRODE = n+1;
RIGHT_ELECTRODE = n+2;

% disp(num_e + " electrons");
% disp(n + " possible targets");

% normalized hopping probs per electron with number of particles + 2 (to beginning and end) options
W_normalized = tot_weights./repmat(sum(tot_weights,2),1,n+2);
W_cdf = cumsum(W_normalized,2); % cumulative distribution function for transitions
rands = rand(num_e,1); % random numbers for each electron

C = rands<W_cdf;
[~,next_indices] = max(C,[],2); % get the index of the first column where rands < cdf, i.e. the next hopping spot

% initialize nextx and nexty
nextx = zeros(num_e,1);
nexty = zeros(num_e,1);

% set beginning positions to be at (0,0)
% the distance later will not factor in y distance so y=nan is fine
nextx(next_indices==LEFT_ELECTRODE) = 0;
nexty(next_indices==LEFT_ELECTRODE) = nan;

% set end positions to be at (L,0)
% once at end, electron is removed from simulation so y=nan is fine
nextx(next_indices==RIGHT_ELECTRODE) = L;
nexty(next_indices==RIGHT_ELECTRODE) = nan;

% set nextx and nexty to be the positions of the particles they hop to
% only do this for electrons that are hopping to a particle (not beginning or end)
not_begin_or_end = (next_indices~=n+1) & (next_indices~=n+2);
nextx(not_begin_or_end) = x(next_indices(not_begin_or_end));
nexty(not_begin_or_end) = y(next_indices(not_begin_or_end));