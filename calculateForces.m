function dxdt = calculateForces(t, states, params)

    warning('off','all')

    global tindex V

    global rand_dirs_global_x rand_dirs_global_y
    if isempty(rand_dirs_global_x)
        rand_dirs_global_x = 2 * randi([0 1], params.n, 1) - 1;
        rand_dirs_global_y = 2 * randi([0 1], params.n, 1) - 1;
    end

    % ------------------------
    % Unpack some params
    % ------------------------
    n = params.n;
    eta = params.eta; % Viscosity

    % ------------------------
    % Collect Initial States
    % ------------------------
    x_p=states(1:n);
    x_v=states(n+1:2*n);
    y_p=states(2*n+1:3*n);
    y_v=states(3*n+1:4*n);
    T=states(4*n+1);


    % ------------------------
    % Calculate drag based on position
    % ------------------------
    % eta = calculate_eta(x_p, y_p, params.eta, params.xp, params.yp);
    
    if tindex <= length(params.tspan) && t > params.tspan(tindex)
        rand_dirs_global_x = 2 * randi([0 1], n, 1) - 1;
        rand_dirs_global_y = 2 * randi([0 1], n, 1) - 1;
        % % if mod(tspan_index,2) == 0
        % I_last = individual_IV_calculator_monte_carlo(x_p, y_p, 200e-6, V, 4.5e-6, 1, 120);
        % % end
        tindex = tindex + 1;
        if tindex > 1000 & V > 0
            V = 0;
            disp('Voltage set to zero after 250 iterations');
        end
        disp(num2str(tindex))
        % disp(distances(params.xp, x_p(1), params.yp, y_p(1)))
        % disp(['Computed I at t = ', num2str(tspan_index), ' | I = ', num2str(I_last)]);
        toc;
        tic;

        % % Save computed current and time
        % I_saved = [I_saved; I_last];
        % t_saved = [t_saved; t];
        
        % compliance = 1;
        % total_e = 50;
        % res = total_e*V/I_last;
        % % V = (compliance * res);
        % if I_last > compliance * total_e
        %     V = 0.8*(V);
        % end
        % disp(['Updated Voltage: ', num2str(V)]);
    end

    % ------------------------
    % Forces
    % ------------------------
    Fa_x = applied_force(n, x_p, params.alpha, V, params.Lx);
    [Fd_x, Fd_y] = drag_force(n, x_v, y_v, params.eta, params.Cd);
    FI_x = interfacial_force(n, x_p, y_p, params.wI, params.RI, params.Lx);
    % FI_x = zeros(n, 1);
    % FI_y = zeros(n, 1); % Initialize interfacial forces
    [Fp_x, Fp_y] = pinning_force(n, x_p, y_p, params.w_pin, params.R_pin, params.Lx);
    % Fp_x = zeros(n, 1); % Initialize pinning forces
    % Fp_y = zeros(n, 1); % Initialize pinning forces
    [Ft_x, Ft_y] = temperature_fluctuations(n, eta, params.T_coeff, T, rand_dirs_global_x, rand_dirs_global_y);
    % Ft_x = zeros(n, 1); % Initialize temperature fluctuation forces
    % Ft_y = zeros(n, 1); % Initialize temperature fluctuation forces

    % If the particle has reached the end, set the velocity to zero
    fin_array = finishing_array(x_p, params.Lx, params.fin);
    forces_x = Fa_x + Fd_x + FI_x + Fp_x + Ft_x; % Total force in x direction
    forces_y = 0    + Fd_y + 0    + Fp_y + Ft_y; % Total force in y direction

    %solve force balance equation in x direction
    dxdt(1:n) = states(n+1:2*n) .* fin_array;
    dxdt(n+1:2*n) = forces_x./eta .* fin_array;

    %solve force balance equation in y direction
    dxdt(2*n+1:3*n) = states(3*n+1:4*n) .* fin_array;
    dxdt(3*n+1:4*n) = forces_y./eta .* fin_array;
    
    %temperature evolution
    dxdt((4*n)+1) = (params.CT * params.Q) - params.k * (T - params.To);

    dxdt = dxdt(:);
end

% % % ------------------------
% % % Calculate eta based on position
% % % ------------------------
% function eta = calculate_eta(x_p, y_p, scale, xp, yp)
%     % Calculate the distance from each particle to the pinning sites
%     if nargin < 4 || isempty(xp) || isempty(yp)
%         xp = zeros(size(x_p)); % Default to zero if no pinning sites provided
%         yp = zeros(size(y_p)); % Default to zero if no pinning sites provided
%     end
    
%     d = sqrt((x_p - xp).^2 + (y_p - yp).^2); % Calculate distances
%     eta = scale ./ (1 + d.^2); % Calculate eta based on distance
% end
    


% % ------------------------
% % Distance Calculation
% % ------------------------

function [d,dx,dy] = distances(x_1, x_2, y_1, y_2)

    assert(all(isfinite(x_1)) && all(isfinite(x_2)), 'x input contains non-finite values');
    assert(all(isfinite(y_1)) && all(isfinite(y_2)), 'y input contains non-finite values');

    dx = x_1 - x_2;  % Calculate the difference in x-coordinates
    dy = y_1 - y_2;  % Calculate the difference in y-coordinates
    d = sqrt(dx.^2 + dy.^2);  % Calculate the Euclidean distance
end

% % ------------------------
% % Applied Force (From Electric Field)
% % ------------------------

function Fa_x = applied_force(n, x_p, alpha, V, Lx)

    % Initialize force to zero
    Fa_x = zeros(n, 1);

    % Logical index of particles within [0, Lx]
    inside = (x_p < Lx);

    Fa_x(inside) = alpha(inside) .* V / ((1/2)*Lx);
    assert(all(Fa_x == 0) | V > 0, 'Applied force should be zero if voltage is zero');
end

% % ------------------------
% % Drag Force
% % ------------------------

function [Fd_x, Fd_y] = drag_force(n, x_v, y_v, eta, Cd)
    % Fd_x=zeros(n,1);
    % Fd_y=zeros(n,1);
    Fd_x = -eta .* Cd .* x_v;
    Fd_y = -eta .* Cd .* y_v;
    % for i=1:n
    %     Fd_x(i,1)=-eta(i)*Cd*x_v(i)*abs(x_v(i));
    %     Fd_y(i,1)=-eta(i)*Cd*y_v(i)*abs(y_v(i));
    % end
end

% % ------------------------
% % Interfacial Force
% % ------------------------

function [FI_x] = interfacial_force(n, x_p, y_p, wI, RI, Lx)
    % Initialize force to zero
    FI_x = zeros(n, 1);

    % Logical index of particles within [0, Lx]
    inside = (x_p > 0) & (x_p < Lx);

    % Compute force only for those particles
    FI_x(inside) = -(2 * wI / RI^2) * ( ...
        x_p(inside) .* exp(-((x_p(inside)).^2) / RI^2) + ...
        (x_p(inside) - Lx) .* exp(-((x_p(inside) - Lx).^2) / RI^2) ...
    ); 
end

% % ------------------------
% % Pinning Force
% % ------------------------

function [Fp_x, Fp_y] = pinning_force(n, x_p, y_p, w_pin, R_pin, Lx)

    % Initialize force to zero
    Fp_x = zeros(n, 1);
    Fp_y = zeros(n, 1);

    % Logical index of particles within [0, Lx]
    inside = (x_p > 0) & (x_p < Lx);

    Fp_x(inside) = -(w_pin * pi / R_pin) * ( cos((2 * pi * (x_p(inside) + y_p(inside)) ) / R_pin) + cos((2 * pi * (x_p(inside) - y_p(inside)) ) / R_pin) );
    Fp_y(inside) = -(w_pin * pi / R_pin) * ( cos((2 * pi * (x_p(inside) + y_p(inside)) ) / R_pin) - cos((2 * pi * (x_p(inside) - y_p(inside)) ) / R_pin) );
end

% % ------------------------
% % Temperature Fluctuations
% % ------------------------

function [Ft_x, Ft_y] = temperature_fluctuations(n, eta, T_coeff, T, rand_dirs_global_x, rand_dirs_global_y)
    noise_scale = sqrt(eta * T_coeff * T);
    Ft_x = rand_dirs_global_x .* noise_scale;
    Ft_y = rand_dirs_global_y .* noise_scale;
end

% % ------------------------
% % Residual Force
% % ------------------------
% TODO: Go over this

% % Stress-strain relationship to indicate the growth of the filament

% % Stuff here... Need to figure out what states is doing. Does it encode the
% % position? Or perhaps I can try to squeeze the strain growth into the
% % dx/dt array which would encode the rate at which the x-position
% % changes... maybe?

% % Find the maximum and minimum x and y positions
% xmin_p = min([x_p; x_v]);
% xmax_p = max([x_p; x_v]);
% ymin_p = min([y_p; y_v]);
% ymax_p = max([y_p; y_v]);

% % Find the midpoint positions
% x_mid = (xmin_p + xmax_p)/2;
% y_mid = (ymin_p + ymax_p)/2;

% delta_x_pos = xmax_p - xmin_p;      % Compute the x width
% delta_y_pos = ymax_p - ymin_p;      % Compute the y width
% Fe = (kC*q^2)/(2*r_Ag)^2;             % Electrostatic force between two silver nanoparticles

% % Add a default strain amount determined by the formation pressure? Start
% % with 100 Pa, then 10 Pa and then 1000 Pa
% sigma_xx = (Fapp+Fe)/(2*delta_x_pos^2) + 1000;    % Stress experienced in x direction would be due to field
% sigma_yy = Fe/(2*delta_y_pos^2) + 1000;      % Stress experienced in y direction would be from nanoparticle interactions

% eps_xx = sigma_xx/E - nu*sigma_yy/E;    % Compute the strain of filament in x
% eps_yy = sigma_yy/E - nu*sigma_xx/E;    % Compute the strain of filament in y

% x_inc = eps_xx*delta_x_pos;     % Change in the x positions
% y_inc = eps_yy*delta_y_pos;     % Change in the y positions

% % Update x positions
% for count=1:n
%     states(count) = states(count) + x_inc*(states(count)-x_mid)^3;
% end

% % Update y positions
% for count=2*n+1:3*n
%     states(count) = states(count) + y_inc*(states(count)-y_mid)^3;
% end

% % ------------------------
% % Get Finishing Array
% % ------------------------
function fin_array = finishing_array(x_p, L, fin)
    % Check if the particles have reached the end
    if fin == 1
        fin_array = x_p < 2*L - 5e-6;
    else 
        fin_array = 1;
    end
end

 % % disp(t)
% if tspan_index <= length(tspan) && t > tspan(tspan_index)
%     tic;
%     % % if mod(tspan_index,2) == 0
%     % I_last = individual_IV_calculator_monte_carlo(x_p, y_p, 200e-6, V, 4.5e-6, 1, 120);
%     % % end
%     tspan_index = tspan_index + 1;
%     disp(num2str(tspan_index))
%     % disp(['Computed I at t = ', num2str(tspan_index), ' | I = ', num2str(I_last)]);
%     toc;

%     % % Save computed current and time
%     % I_saved = [I_saved; I_last];
%     % t_saved = [t_saved; t];
    
%     % compliance = 1;
%     % total_e = 50;
%     % res = total_e*V/I_last;
%     % % V = (compliance * res);
%     % if I_last > compliance * total_e
%     %     V = 0.8*(V);
%     % end
%     % disp(['Updated Voltage: ', num2str(V)]);
% end