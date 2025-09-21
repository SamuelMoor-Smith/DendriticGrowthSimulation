function [I, pointsx, pointsy] = calculateCurrent(x, y, params)

L = params.L;
voltage = params.V;
lambda = params.lambda;
Rt = params.Rt;

p = 1000;
t2 = 1000;

sz = size(x);
t1 = sz(1);
n = sz(2);
newL = 200e-6;

pointsx = zeros(t1,n);
pointsy = zeros(t1,n);
I = zeros(t1,1);

for i = 1:t1
    
    count2 = 1;
    
    % Start with electrons at random particles
    
%     curr_indexes = randsample(1:n,p,true);
%     currx = x(i,curr_indexes);
%     curry = y(i,curr_indexes);

    currx = zeros(1,p);
    curry = zeros(1,p);
    curr_indexes = zeros(1,p);
    
    xi = x(i,:);
    yi = y(i,:);
    
    for j = 1:t2
        
        % a1 = ~ismember(currx(:),pointsx(i,:));
        % a2 = ~ismember(curry(:),pointsy(i,:));
        
        % a3 = find(a1 == 1 & a2 == 1);
        
        % if ~isempty(a3)
        %     count = length(a3);
        %     pointsx(i,count2:count2+count-1) = currx(a3);
        %     pointsy(i,count2:count2+count-1) = curry(a3);
        %     count2 = count2 + count;
        % end
            
        % Get distances to particles in +x direction
        
        dx = repmat(xi,p,1) - currx';
        dy = repmat(yi,p,1) - curry';
        
        dy(dy==repmat(yi,p,1)) = 0;
        
        d = sqrt(dx.^2 + dy.^2);
        
        positive = dx >= 0;
        
        dist_weights = positive .* 1./(Rt*exp(d/lambda));
        dist_weights(isnan(dist_weights)) = 0;
        
        to_end = (1./(Rt*exp((newL-currx')/lambda))) .* ones(p,1);
        to_end(isnan(to_end)) = 1;
        
        stay_at_beginning = currx' == 0;
        
        tot_weights = [dist_weights stay_at_beginning to_end];
        
        [curr_indexes,currx,curry]=getNextIndex(i,tot_weights,x,y,currx,curry,curr_indexes);
        
    end
    I(i) = sum(isnan(currx));
end

%     sz = size(x);
%     t1 = sz(1);
% 
%     n = params.n;
%     newL = 200e-6;
% 
%     pointsx = zeros(t1,n);
%     pointsy = zeros(t1,n);
%     I = zeros(t1,1);
% 
%     for i = 1:t1
% 
%         count2 = 1;
% 
%         xi = x(i,:);
%         yi = y(i,:);
% 
%         % Get initial electron positions
%         good = (xi > 0) .* (xi < L);
%         curr_indexes = find(good == 1);
%         w = ones(1,length(curr_indexes));
% 
%         if length(curr_indexes) > 1
%             curr_indexes = randsample(curr_indexes, num_e, true, w);
%         end
%         if length(curr_indexes) > num_e
%             curr_indexes = curr_indexes(1:num_e);
%         end
% 
%         p = length(curr_indexes);
%         if p == 0
%             break
%         end
% 
%         % Start with electrons at all particles in device
%         directions = [1,-1];
%         for dirNum=1:2
%             % Get direction value
%             direction = directions(dirNum);
% 
%             % Get initial indices
%             curr_indexes_local = curr_indexes;
%             currx_local = xi(curr_indexes_local);
%             curry_local = yi(curr_indexes_local);
% 
%             for j = 1:steps
% 
%                 % Get distances
%                 dx = repmat(xi,p,1) - currx_local';
%                 dy = repmat(yi,p,1) - curry_local';
%                 d = sqrt(dx.^2 + dy.^2);
%                 dx = dx*direction;
%                 positive = dx >= 0;
% 
%                 % Get resistances
%                 voltage = V*dx/L;
%                 sigmoid = 1./(1+exp(-voltage));
%                 resist = Rt*exp((d)/lambda);
% 
%                 dist_weights = positive .* sigmoid./resist;
%                 dist_weights(isnan(dist_weights)) = 0;
% 
%                 d = direction*(L/2 - currx_local') + L/2;
% 
%                 % Get resistance to end
%                 voltage_end = V.*(d)./L;
%                 sigmoid_end = 1./(1+exp(-voltage_end));
%                 resist_end = Rt*exp((d)/lambda);
% 
%                 to_end = (sigmoid_end./resist_end) .* ones(p,1);
%                 to_end(isnan(to_end)) = 1;
%                 to_begin = zeros(p,1);
% 
%                 tot_weights = [dist_weights to_begin to_end];
% 
%                 [curr_indexes_local,currx_local,curry_local] = getNextIndex(tot_weights, x, y, currx_local, curry_local, curr_indexes_local);
%             end
%         end
% 
%         % Add plots
%         I(i) = sum(isnan(currx_for) .* isnan(currx_back));
%         % disp(I(i))
%         % hold onx
%         % plot(i,I(i),'b.');
%         % drawnow
%     end
% end