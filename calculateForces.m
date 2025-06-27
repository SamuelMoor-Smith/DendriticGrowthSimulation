function dxdt = calculateForces(t, states, params)

    warning('off','all')

    global I_saved t_saved
    if isempty(I_saved)
        I_saved = []; % Initialize if empty
        t_saved = [];
    end

    % Define persistent voltage variable
    persistent V I_last tspan_index
    if isempty(V)
        V = 30; % Initial Voltage
        disp(V)
        tspan_index = 1;
        I_last = 0;
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

    % ------------------------
    % Forces
    % ------------------------
    Fa_x = applied_force(params.alpha, params.V, params.L);
    [Fd_x, Fd_y] = drag_force(n, x_v, y_v, params.eta, params.Cd);
    [FI_x, FI_y] = interfacial_force(n, x_p, y_p, params.wI, params.RI);
    [Fp_x, Fp_y] = pinning_force(n, x_p, y_p, params.wp, params.xp, params.yp, params.Rp);
    % [Ft_x, Ft_y] = temperature_fluctuations(n, eta, params.k, T);
    Ft_x = zeros(n, 1); % Initialize temperature fluctuation forces
    Ft_y = zeros(n, 1); % Initialize temperature fluctuation forces

    % If the particle has reached the end, set the velocity to zero
    fin_array = finishing_array(x_p, params.L, params.fin);
    forces_x = Fa_x + Fd_x + FI_x + Fp_x + Ft_x; % Total force in x direction
    forces_y = Fd_y + FI_y + Fp_y + Ft_y; % Total force in y direction

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
    dx = x_1 - x_2;  % Calculate the difference in x-coordinates
    dy = y_1 - y_2;  % Calculate the difference in y-coordinates
    d = sqrt(dx.^2 + dy.^2);  % Calculate the Euclidean distance
end

% % ------------------------
% % Applied Force (From Electric Field)
% % ------------------------

function Fa_x = applied_force(alpha, V, L)
    Fa_x = alpha .* V / L;
end

% % ------------------------
% % Drag Force
% % ------------------------

function [Fd_x, Fd_y] = drag_force(n, x_v, y_v, eta, Cd)
    Fd_x=zeros(n,1);
    Fd_y=zeros(n,1);
    for i=1:n
        Fd_x(i,1)=-eta(i)*Cd*x_v(i)*abs(x_v(i));
        Fd_y(i,1)=-eta(i)*Cd*y_v(i)*abs(y_v(i));
    end
end

% % ------------------------
% % Interfacial Force
% % ------------------------

function [FI_x, FI_y] = interfacial_force(n, x_p, y_p, wI, RI)
    FI_x=zeros(n,1);
    FI_y=zeros(n,1);
    for i=1:n
        [d,dx,dy] = distances(x_p, x_p(i), y_p, y_p(i));

        % Note, I played around with this function quite a bit
        % F=(2*wI).*(100e-9./d.^3).*((100e-9./d)-1);
        F=-(2*wI).*exp(-(d.^2)/RI^2);%exp(-(d-100e-9).^2/RI^2));
        
        % Sum all the forces
        FI_x = FI_x + (F.*dx);
        FI_y = FI_y + (F.*dy);
    end
end

% % ------------------------
% % Pinning Force
% % ------------------------

function [Fp_x, Fp_y] = pinning_force(n, x_p, y_p, wp, xp, yp, Rp)
    Fp_x=zeros(n,1);
    Fp_y=zeros(n,1);
    for i=1:length(wp)
        [d,dx,dy] = distances(x_p, xp(i), y_p, yp(i));
        
        % Note, I played around with this function quite a bit
        % F=2*10*wp(i).*((d.^2/Rp(i)^2).^9).*exp(-((d.^2/Rp(i)^2).^10));
        F=-(2*wp(i)).*(-d.^2/Rp^2);

        % Sum all the forces
        Fp_x = Fp_x + (F.*dx);
        Fp_y = Fp_y + (F.*dy);
    end
end

% % ------------------------
% % Temperature Fluctuations
% % ------------------------

function [Ft_x, Ft_y] = temperature_fluctuations(n, eta, k, T)
    Ft_x = (rand(n, 1)-0.5) .* sqrt(eta*k*T);
    Ft_y = (rand(n, 1)-0.5) .* sqrt(eta*k*T);
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
%     % if mod(tspan_index,2) == 0
%     I_last = individual_IV_calculator_monte_carlo(x_p, y_p, 200e-6, V, 4.5e-6, 1, 120);
%     % end
%     tspan_index = tspan_index + 1;
%     disp(num2str(tspan_index))
%     disp(['Computed I at t = ', num2str(tspan_index), ' | I = ', num2str(I_last)]);
%     toc;

%     % Save computed current and time
%     I_saved = [I_saved; I_last];
%     t_saved = [t_saved; t];
    
%     compliance = 1;
%     total_e = 50;
%     res = total_e*V/I_last;
%     % V = (compliance * res);
%     if I_last > compliance * total_e
%         V = 0.8*(V);
%     end
%     disp(['Updated Voltage: ', num2str(V)]);
% end