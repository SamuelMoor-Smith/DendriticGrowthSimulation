% function [curr_indexes,currx2,curry2] = getNextIndex(i,tot_weights,x,y,currx,curry,curr_indexes)
% 
% p = length(currx);
% n = size(tot_weights,2)-2;
% 
% W_normalized = tot_weights./repmat(sum(tot_weights,2),1,size(tot_weights,2));
% W_cdf = cumsum(W_normalized,2);
% rands = rand(size(tot_weights,1),1);
% C = rands<W_cdf;
% [~,curr_indexes] = max(C,[],2);
% 
% currx2 = nan(1,p);
% curry2 = nan(1,p);
% 
% currx2(curr_indexes==n+1) = 0;
% curry2(curr_indexes==n+1) = 0;
% 
% currx2(curr_indexes==n+3) = currx(curr_indexes==n+3);
% curry2(curr_indexes==n+3) = curry(curr_indexes==n+3);
% 
% curr_indexes(curr_indexes==n+1) = nan;
% curr_indexes(curr_indexes==n+2) = nan;
% curr_indexes(curr_indexes==n+3) = nan;
% 
% % indexes = [1:n,nan];
% % 
% % curr_indexes = nan(1,p);
% % 
% % for k = 1:p
% %     
% %     curr_indexes(k) = randsample(indexes,1,true,tot_weights(k,:));
% % 
% % end
% 
% % set currx and curry to all-nans of the correct size
% % update currx where curr_indexes is not nan
% % fill these with curr_indexes of x (excluding places where curr_indexes is nan)
% currx2(~isnan(curr_indexes)) = x(i,curr_indexes(~isnan(curr_indexes)));
% curry2(~isnan(curr_indexes)) = y(i,curr_indexes(~isnan(curr_indexes)));

function [curr_indexes,currx,curry] = getNextIndex(i,tot_weights,x,y,currx,curry,curr_indexes)

p = length(currx);
n = size(tot_weights,2)-2;

W_normalized = tot_weights./repmat(sum(tot_weights,2),1,size(tot_weights,2));
W_cdf = cumsum(W_normalized,2);
rands = rand(size(tot_weights,1),1);
C = rands<W_cdf;
[~,curr_indexes] = max(C,[],2);

currx = nan(1,p);
curry = nan(1,p);

currx(curr_indexes==n+1) = 0;
curry(curr_indexes==n+1) = 0;

curr_indexes(curr_indexes==n+1) = nan;
curr_indexes(curr_indexes==n+2) = nan;

% indexes = [1:n,nan];
% 
% curr_indexes = nan(1,p);
% 
% for k = 1:p
%     
%     curr_indexes(k) = randsample(indexes,1,true,tot_weights(k,:));
% 
% end

% set currx and curry to all-nans of the correct size
% update currx where curr_indexes is not nan
% fill these with curr_indexes of x (excluding places where curr_indexes is nan)
currx(~isnan(curr_indexes)) = x(i,curr_indexes(~isnan(curr_indexes)));
curry(~isnan(curr_indexes)) = y(i,curr_indexes(~isnan(curr_indexes)));