clear

% SET FILE OUTPUT LOCATION
output = 'D:\User\Location'; % Video Save Location

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
dt = 0.1;       % time step
terminaltime = 5000; % simulation time in "ms"
time_steps = 0:dt:terminaltime;

% Cell Parameters (Original FIB)
cellType1 = struct('a', 0.02,    'b', 2.9,   'c', 1.0, 'epsilon', .0040, 'sigma', 1.0, 'g_horz', .027, 'g_vert', .015); % Healthy
cellType2 = struct('a', 0.01,    'b', 5.3,   'c', 1.0, 'epsilon', .0040, 'sigma', 0.6, 'g_horz', .027, 'g_vert', .015); % Slow Recovery
cellType3 = struct('a', 0.02,    'b', 2.9,   'c', 1.0, 'epsilon', .0065, 'sigma', 1.6, 'g_horz', .027, 'g_vert', .015); % Fast Recovery

% FIX SCAR CHANGE G GRID MAYBE
cellType4 = struct('a', .02,    'b', 2.9,   'c', 1.0, 'epsilon', .0040, 'sigma', 0.6, 'g_horz', .027, 'g_vert', .015); % Non excitable (Scar)

% Block 2 Parameters
cellType5 = struct('a', 0.01,    'b', 5.3,   'c', 1.0, 'epsilon', .0040, 'sigma', 0.6, 'g_horz', .027, 'g_vert', .015); % Slow Recovery
cellType6 = struct('a', 0.02,    'b', 2.9,   'c', 1.0, 'epsilon', .006, 'sigma', 1.6, 'g_horz', .027, 'g_vert', .015); % Fast Recovery

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid size
num_cells_x = 41; 
num_cells_y = 41;

Block_Adjust = 0; % starting at 7 tall
Vert_Adjust = 0; % starting centered at (21,21)
Horz_Adjust = 0; % starting centered at (21,21)

% Grid cell type layout
cellTypes = ones(num_cells_y,num_cells_x);

cellTypes(18 - Block_Adjust + Vert_Adjust : 24 + Block_Adjust + Vert_Adjust, ...
    21 + Horz_Adjust : 22 + Horz_Adjust) = 2; % Bistable Long QT

cellTypes(18 - Block_Adjust + Vert_Adjust : 24 + Block_Adjust + Vert_Adjust, ...
    19 + Horz_Adjust : 20 + Horz_Adjust) = 3; % Fast Epsilon

% Line of Scar Tissue
cellTypes(17,1:24) = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters mapped to grid
a_grid       = zeros(num_cells_y,num_cells_x);
b_grid       = zeros(num_cells_y,num_cells_x);
c_grid       = zeros(num_cells_y,num_cells_x);
epsilon_grid = zeros(num_cells_y,num_cells_x);
sigma_grid   = zeros(num_cells_y,num_cells_x);
g_horz_grid  = zeros(num_cells_y,num_cells_x);
g_vert_grid  = zeros(num_cells_y,num_cells_x);

allTypes = {cellType1, cellType2, cellType3, cellType4, cellType5, cellType6};
for k = 1:6
    mask = (cellTypes == k);
    a_grid(mask)       = allTypes{k}.a;
    b_grid(mask)       = allTypes{k}.b;
    c_grid(mask)       = allTypes{k}.c;
    epsilon_grid(mask) = allTypes{k}.epsilon;
    sigma_grid(mask)   = allTypes{k}.sigma;
    g_horz_grid(mask)  = allTypes{k}.g_horz;
    g_vert_grid(mask)  = allTypes{k}.g_vert;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
v = zeros(num_cells_y, num_cells_x);
w = zeros(num_cells_y, num_cells_x);

v_history = zeros(num_cells_y, num_cells_x, length(time_steps));
w_history = zeros(num_cells_y, num_cells_x, length(time_steps));

% External stimulus
I_ext = @(t) 4 .* (t < 50);

% Locate any scar cells
scarmask = cellTypes == 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vectorized derivative function
function [dv, dw] = derivatives(v, w, ...
                                g_horz_grid, g_vert_grid, a_grid, b_grid, c_grid, ...
                                epsilon_grid, sigma_grid, I_ext_func, t, scarmask)

    % Left and Right Boundaries Insualted
    left   = [v(:,1), v(:,1:end-1)];
    right  = [v(:,2:end), v(:,end)];

    % Top and Bottom Insulated
    up     = [v(1,:); v(1:end-1,:)];
    down   = [v(2:end,:); v(end,:)];

    v_horz_neighbors = left + right; 
    v_vert_neighbors = up + down;

    % External stimulus
    stim = zeros(size(v));
    stim(1:16,1) = I_ext_func(t);
    stim(18:end,1) = I_ext_func(t);

    % Cubic reaction term
    cubic = v .* (v - a_grid) .* (c_grid - v);
    
    % Set cubic part to 0 for scar cells
    cubic(scarmask) = 0;

    dv = sigma_grid .* (cubic - w ...
        + g_horz_grid .* (v_horz_neighbors - 2*v) ...
        + g_vert_grid .* (v_vert_neighbors - 2*v)) ...
        + stim;

    dw = epsilon_grid .* (v - b_grid .* w);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main RK4 loop
for t_idx = 1:length(time_steps)
    t = time_steps(t_idx);

    v_prev = v;
    w_prev = w;

    [k1_v, k1_w] = derivatives(v_prev, w_prev, ...
                               g_horz_grid, g_vert_grid, a_grid, b_grid, c_grid, ...
                               epsilon_grid, sigma_grid, I_ext, t, scarmask);

    [k2_v, k2_w] = derivatives(v_prev + 0.5*dt*k1_v, w_prev + 0.5*dt*k1_w, ...
                               g_horz_grid, g_vert_grid, a_grid, b_grid, c_grid, ...
                               epsilon_grid, sigma_grid, I_ext, t + 0.5*dt, scarmask);

    [k3_v, k3_w] = derivatives(v_prev + 0.5*dt*k2_v, w_prev + 0.5*dt*k2_w, ...
                               g_horz_grid, g_vert_grid, a_grid, b_grid, c_grid, ...
                               epsilon_grid, sigma_grid, I_ext, t + 0.5*dt, scarmask);

    [k4_v, k4_w] = derivatives(v_prev + dt*k3_v, w_prev + dt*k3_w, ...
                               g_horz_grid, g_vert_grid, a_grid, b_grid, c_grid, ...
                               epsilon_grid, sigma_grid, I_ext, t + dt, scarmask);

    v = v_prev + (dt/6) * (k1_v + 2*k2_v + 2*k3_v + k4_v);
    w = w_prev + (dt/6) * (k1_w + 2*k2_w + 2*k3_w + k4_w);

    v_history(:,:,t_idx) = v;
    w_history(:,:,t_idx) = w;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Heatmap Video
    % Create Video Output
    videoFile = fullfile(output, 'Unnamed Cardiac Grid.mp4');

    % Create a VideoWriter object
    video = VideoWriter(videoFile, 'MPEG-4');
    video.FrameRate = 15; % adjust playback speed as desired
    open(video);

    % Create a hidden figure for off-screen rendering
    fig = figure('Visible', 'off');
    ax = axes(fig);

    % Loop through time points and generate frames
    for t_idx = 1:50:length(time_steps)
        
        % Plot the heatmap in the invisible figure
        imagesc(ax, v_history(:,:,t_idx));
        title(ax, sprintf('Membrane Potentials at t = %.1f ms', time_steps(t_idx)));
        clim(ax, [-0.2, 1.2]); % consistent color scale
        axis equal tight
        ax = gca;   % Get current axes
        ax.XColor = 'none';  % Hide x-axis line
        ax.YColor = 'none';  % Hide y-axis line
        ax.XTick = [];       % Remove x-axis ticks
        ax.YTick = [];       % Remove y-axis ticks

        hold on
        % mask annotations
            recovery_mask = (w_history(:,:,t_idx) > 0.0225) & (v_history(:,:,t_idx) < 0);
            [y, x] = find(recovery_mask);
            plot(x, y, 'x', 'Color', 'k', 'LineWidth', 1.2, 'MarkerSize', 5);
            [y, x] = find(scarmask);
            plot(x, y, 's', 'Color', 'k', 'LineWidth', 1.2, 'MarkerSize', 12);

            % Annotations (Normal FIB)
            rectangle('Position', [18.5, 17.5, 4, 7], ...
            'EdgeColor', 'k', 'LineWidth', 2);

            rectangle('Position', [18.5, 17.5, 4, 7], ...
            'EdgeColor', 'w', 'LineWidth', 2, 'LineStyle',':');

            line([20.5,20.5],[17.5,24.5],'Color','k','LineWidth',2)

            line([20.5,20.5],[17.5,24.5],'Color','w','LineWidth',2,'LineStyle',':')
        hold off
        % Force MATLAB to render internally
        drawnow;

        % Capture the current figure frame (off-screen)
        frame = getframe(fig);

        % Write the frame to the video
        writeVideo(video, frame);
    end

    % Close video file and figure
    close(video);
    close(fig);

    fprintf('Video successfully created: %s\n', videoFile);