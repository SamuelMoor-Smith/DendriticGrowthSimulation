function dxdt = calculateForces(t, states, params)

    warning('off','all')

    global tindex V

    global I_saved t_saved I_last
    if isempty(I_saved)
        I_saved = []; % Initialize if empty
        t_saved = [];
        I_last = 0;
    end

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
    % Calculate current at intervals
    % ------------------------
    if tindex <= length(params.tspan) && t > params.tspan(tindex)
        rand_dirs_global_x = 2 * randi([0 1], n, 1) - 1;
        rand_dirs_global_y = 2 * randi([0 1], n, 1) - 1;
        if mod(tindex,2) == 0
            I_last = calculateCurrent(x_p, y_p, params.Lx, V, params.lambda, params.Rt, params.steps, params.num_e);
        end
        tindex = tindex + 1;
        if tindex > 1000 & V > 0
            V = 0;
            disp('Voltage set to zero after 1000 iterations');
        end
        disp(['Time index: ', num2str(tindex), ' / ', num2str(length(params.tspan))]);
        disp(['Computed I at t = ', num2str(tindex), ' | I = ', num2str(I_last)]);
        toc;
        tic;

        % Save computed current and time
        I_saved = [I_saved; I_last];
        t_saved = [t_saved; t];
    end

    % ------------------------
    % Forces
    % ------------------------
    % eta = calculate_eta(x_p, y_p, params.eta, params.xp, params.yp);
    Fa_x = applied_force(n, x_p, params.alpha, V, params.Lx);
    [Fd_x, Fd_y] = drag_force(n, x_v, y_v, params.eta, params.Cd);
    FI_x = interfacial_force(n, x_p, y_p, params.wI, params.RI, params.Lx);
    [Fp_x, Fp_y] = pinning_force(n, x_p, y_p, params.w_pin, params.x_pin, params.y_pin, params.R_pin, params.Lx);
    [Ft_x, Ft_y] = temperature_fluctuations(n, eta, params.T_coeff, T, rand_dirs_global_x, rand_dirs_global_y);

    % In progress 
    [Fr_x, Fr_y] = residual_force(n, x_p, y_p ,params.Lx, params.Ly);
    % [Fc_x, Fc_y] = contact_force_intrabins(n, x_p, y_p, params.Lx, params.Ly);
    % Fr_x = zeros(n, 1);
    % Fr_y = zeros(n, 1);
    Fc_x = zeros(n, 1); % Initialize contact forces
    Fc_y = zeros(n, 1); % Initialize contact forcesåå
    % In progress

    % If the particle has reached the end, set the velocity to zero
    fin_array = finishing_array(x_p, params.Lx, params.fin);

    forces_x = Fa_x + Fd_x + FI_x + Fp_x + Ft_x + Fr_x + Fc_x; % Total force in x direction
    forces_y = 0    + Fd_y + 0    + Fp_y + Ft_y + Fr_y + Fc_y; % Total force in y direction

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

function [Fp_x, Fp_y] = pinning_force(n, x_p, y_p, w_pin, x_pin, y_pin, R_pin, Lx)
    assert(isvector(x_p) && isvector(y_p), 'x_p and y_p must be vectors');
    assert(length(x_p) == length(y_p), 'x_p and y_p must be the same length');
    assert(length(x_p) == n, 'x_p must have length n');

    Fp_x=zeros(n,1);
    Fp_y=zeros(n,1);
    for i=1:length(w_pin)

        assert(isscalar(x_pin(i)) && isscalar(y_pin(i)), 'xp and yp must be scalars for each pinning site');
        [d,dx,dy] = distances(x_p, x_pin(i), y_p, y_pin(i));

        F = (2 * w_pin(i) / (R_pin(i)^2)) .* exp(-(d.^2) / (R_pin(i)^2));

        % Sum all the forces
        Fp_x = Fp_x + (F.*dx);
        Fp_y = Fp_y + (F.*dy);
    end
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
% % TODO: Residual Force
% % Fast local density gradient approximation and then force based on that
% % ------------------------

function [Fr_x, Fr_y] = residual_force(n, x_p, y_p, Lx, Ly)
    % Residual (crowding) force based on smoothed local density gradients.
    %
    % Concept: each particle feels a repulsion away from high-density regions,
    %          similar to crowd diffusion or self-exclusion.
    %
    % Parameters:
    cell_size = 5.0;    % grid resolution (µm or whatever your scale is)
    R_resid   = 1.5*cell_size; % Gaussian width
    w_resid   = 100.0;      % strength (like your pinning w)

    % ------------------------
    % Grid setup
    % ------------------------
    Nx = max(2, ceil(Lx / cell_size));
    Ny = max(2, ceil(Ly / cell_size));
    x_edges = linspace(0, Lx, Nx+1);
    y_edges = linspace(-Ly/2, Ly/2, Ny+1);
    x_centers = 0.5 * (x_edges(1:end-1) + x_edges(2:end));
    y_centers = 0.5 * (y_edges(1:end-1) + y_edges(2:end));

    % ------------------------
    % Count particles per cell (density field)
    % ------------------------
    [count, ~, ~, ~, ~] = histcounts2(y_p, x_p, y_edges, x_edges);

    % --- Accumulate forces from all occupied cells (pinning-like loop) ---
    Fr_x = zeros(n,1);
    Fr_y = zeros(n,1);
    R2   = R_resid.^2;

    for iy = 1:Ny
        for ix = 1:Nx
            m = count(iy, ix);             % occupancy of cell (iy,ix)
            if m <= 0, continue; end       % skip empty cells

            xc = x_centers(ix);
            yc = y_centers(iy);

            dx = x_p - xc;                 % vector to all particles
            dy = y_p - yc;
            d2 = dx.^2 + dy.^2;

            Fmag = (2 * w_resid * m / R2) .* exp(-d2 ./ R2);  % scalar per particle
            Fr_x = Fr_x + Fmag .* dx;
            Fr_y = Fr_y + Fmag .* dy;
        end
    end

    % column vectors
    Fr_x = Fr_x(:);
    Fr_y = Fr_y(:);
end

% % ------------------------
% % TODO: Contact Force (Fast, same-bin only)
% % ------------------------

function [Fc_x, Fc_y] = contact_force_intrabins(n, x_p, y_p, Lx, Ly)
    % Fast contact repulsion using an exponential cutoff, same-bin pairs only.
    % Particles repel with F = Kc * exp(-d / lambda) * r_hat for d <= r_cut.
    %
    % Tuning (units ~ your domain):
    cell_size = 3.0;   % bin size; choose >= r_cut for speed/coverage
    lambda    = 1.0;   % decay length (how fast force drops with distance)
    r_cut     = 3.0;   % hard cutoff distance (often ~3*lambda)
    Kc        = 2000;  % strength scale
    eps0      = 1e-12; % tiny to avoid 0/0

    % ---- grid / binning ----
    Nx = max(1, ceil(Lx / cell_size));
    Ny = max(1, ceil(Ly / cell_size));

    y0 = min(y_p);
    iy = min(max(floor((y_p - y0) ./ cell_size) + 1, 1), Ny);
    ix = min(max(floor(x_p          ./ cell_size) + 1, 1), Nx);

    bins = accumarray([iy, ix], (1:n).', [Ny, Nx], @(v){v}, {[]});
    Fc_x = zeros(n,1); Fc_y = zeros(n,1);

    r2_cut = r_cut^2;

    % ---- same-bin pairs only (upper triangle) ----
    for by = 1:Ny
        for bx = 1:Nx
            A = bins{by,bx};
            m = numel(A);
            if m < 2, continue; end

            Ax = x_p(A); Ay = y_p(A);
            dX = Ax - Ax.'; dY = Ay - Ay.';
            D2 = dX.^2 + dY.^2;

            mask = triu(true(m),1) & (D2 <= r2_cut);
            if ~any(mask,'all'), continue; end

            [I,J] = find(mask);
            dx = Ax(I) - Ax(J); dy = Ay(I) - Ay(J);
            d  = sqrt(dx.^2 + dy.^2) + eps0;

            % soft exponential repulsion
            mag = Kc * lambda^2 ./ d.^2;   % scalar per interacting pair
            fx  = mag .* (dx ./ d);
            fy  = mag .* (dy ./ d);

            % equal & opposite
            Fc_x(A(I)) = Fc_x(A(I)) + fx;  Fc_y(A(I)) = Fc_y(A(I)) + fy;
            Fc_x(A(J)) = Fc_x(A(J)) - fx;  Fc_y(A(J)) = Fc_y(A(J)) - fy;
        end
    end

    % column vectors
    Fc_x = Fc_x(:); Fc_y = Fc_y(:);
end

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