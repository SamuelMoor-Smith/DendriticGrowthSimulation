function [o1,o2] = runSimulation(params)
    
    global I_saved t_saved %#ok<GVMISå>
    I_saved = [];
    t_saved = [];

    rng(0)

    % ------------------------
    % Unpack parameters
    % ------------------------
    n = params.n;
    To = params.To;
    tspan = params.tspan;
    xoWidth = params.xoWidth;
    yoHeight = params.yoHeight;

    % ------------------------
    % Randomize initial particle positions
    % ------------------------
    xo = -xoWidth.*rand(n,1); 
    yo = linspace(-yoHeight/2,yoHeight/2,n).';

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
        [t,states]=ode45(@(t,states) calculateForces(t, states, params), tspan, initial);
    toc;

    % TODO: Go over this
    % ------------------------
    
    % A = [t,states];
    % t = A(:,1);               % First column is time
    % states = A(:,2:end);      % Remaining columns are the states

    X_pos = states(:,1:n);      % First n columns are x-positions
    Y_pos = states(:,2*n:3*n);  % Third group of n columns are y-positions
    
    %produce 2D plot of particle positions, save as gif
    fig = figure(4);
    h2 = tiledlayout(2,2); 
    ax1 = nexttile([2 1]);
    ax2 = nexttile;
    ax3 = nexttile;

    hold(ax1,'off')
    hold(ax2,'on')
    hold(ax3,'on')

    for idx=1:length(t)
        plot(ax1, X_pos(idx,1:n), Y_pos(idx,1:n),'b.')
        plot(ax2, X_pos(idx,1:n), Y_pos(idx,1:n),'b.')

        % xlim(ax1,[-100e-6,100e-6])
        xlim(ax1, [0, 100e-6])
        ylim(ax1,[-5e-6,5e-6])

        xlabel(ax1,'x-position')
        ylabel(ax1,'y-position')

        xlim(ax2,[0,100e-6])
        ylim(ax2,[-5e-6,5e-6])

        xlabel(ax2,'x-position')
        ylabel(ax2,'y-position')

        % xlim(ax3,[0,2.5e-4])
        % ylim(ax3,[0,max(I)])

        xlabel(ax3,'Time')
        ylabel(ax3,'Current')

        %pause
        legend(num2str(t(idx)))
        drawnow
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

