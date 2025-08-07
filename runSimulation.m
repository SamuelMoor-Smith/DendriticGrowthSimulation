function [o1,o2] = runSimulation(params)
    
    global I_saved t_saved %#ok<GVMISå>
    I_saved = [];
    t_saved = [];

    rng(1)

    % ------------------------
    % Unpack parameters
    % ------------------------
    n = params.n;
    To = params.To;
    tspan = params.tspan;

    Lx = params.Lx; % Size of material
    Ly = params.Ly; % Size of material
    electrode_width = params.electrode_width;
    electrode_height = params.electrode_height;

    % ------------------------
    % Randomize initial particle positions
    % ------------------------
    xo = 10.*rand(n,1); 
    yo = electrode_height.*rand(n,1) - 0.5*electrode_height;

    % ------------------------
    % Zero all initial particle velocities
    % ------------------------
    vox = zeros(n,1);
    voy = zeros(n,1);

    % ------------------------
    % Initial config vector
    % ------------------------
    initial = [xo.',vox.',yo.',voy.',To];

    % ------------------------
    % Run ode simulation of the differentiable force function
    % ------------------------
    tic;
        disp("Starting the solver");
        tic;
        global tindex V
        V = params.V; % Set the initial voltage
        tindex = 1; % Initialize the time index for the simulation
        [t,states]=ode45(@(t,states) calculateForces(t, states, params), tspan, initial);
    toc;

    % TODO: Go over this
    % ------------------------
    
    % A = [t,states];
    % t = A(:,1);               % First column is time
    % states = A(:,2:end);      % Remaining columns are the states

    X_pos = states(:,1:n);      % First n columns are x-positions
    Y_pos = states(:,2*n+1:3*n);  % Third group of n columns are y-positions
    
    %produce 2D plot of particle positions, save as gif
    fig = figure(4);
    h2 = tiledlayout(2,2); 
    ax1 = nexttile([1 2]);

    % Set axis limits and labels once
    xlim(ax1, [-electrode_width - 5, Lx+electrode_width + 5]);
    ylim(ax1, [-electrode_height/2 - 5, electrode_height/2 + 5]);
    xlabel(ax1, 'x-position');
    ylabel(ax1, 'y-position');

    hold(ax1, 'on')  % Keep everything on the same axes

    % Draw patch background once
    patch(ax1, ...
        [-electrode_width, 0, 0, -electrode_width], ...
        [-electrode_height/2, -electrode_height/2, electrode_height/2, electrode_height/2], ...
        [1.0, 0.8431, 0.0], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.3);  % Optional transparency
    
    patch(ax1, ...
        [Lx+electrode_width, Lx, Lx, Lx+electrode_width], ...
        [-electrode_height/2, -electrode_height/2, electrode_height/2, electrode_height/2], ...
        [1.0, 0.8431, 0.0], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.3);  % Optional transparency

    % Plot static pinning sites
    % plot(ax1, params.x_pin, params.y_pin, 'rx', 'MarkerSize', 6, 'LineWidth', 1.5);
    w_pin = params.w_pin;
    R_pin = params.R_pin;

    x = linspace(0, Lx, 100);
    y = linspace(-Ly/2, Ly/2, 100);
    [X, Y] = meshgrid(x, y);
    U = (w_pin / 2) * (sin((2 * pi * (X + Y)) / R_pin) + cos((2 * pi * (X - Y)) / R_pin));
    
    % Transparent surface at z = 0
    surf(ax1, X, Y, zeros(size(U)), U, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);  % Apply transparency

    colormap(ax1, 'pink')


    % Initialize handle for particle plot
    hParticles = plot(ax1, X_pos(1,1:n), Y_pos(1,1:n), 'b.');

    for idx=1:length(t)

        % Update particle positions
        set(hParticles, 'XData', X_pos(idx,1:n), 'YData', Y_pos(idx,1:n));

        % Update legend or annotation
        legend(ax1, num2str(t(idx)));

        drawnow

        % Save frame for GIF
        frame = getframe(fig);
        im{idx} = frame2im(frame);

    end

    filename = sprintf('dendrite_growthNEW%d.gif', i);
    for idx = 1:length(t)
        [A2,map] = rgb2ind(im{idx},256);
        if idx == 1
            imwrite(A2,map,filename,'gif','LoopCount',Inf,'DelayTime',0.1);
        else
            imwrite(A2,map,filename,'gif','WriteMode','append','DelayTime',0.1);
        end
    end

    end

