clear

% Choose which figures to generate (1/0)

GenerateHeatMap = 1;
    RunAnimation = true;      % Toggle live animation visualization
    ExportScreenshots = false; % Toggle still frame export (no animation)
Generate3dMap = 0;
GenerateHeatMapVideo = 0;
GenerateVoltage = 0;
GenerateVDiff = 0;
GenerateVWT = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
dt = 0.1;       % time step
terminaltime = 4500; % simulation time in ms
time_steps = 0:dt:terminaltime;

gv = 0.015;
gh = 0.027;

% Cell Parameters (Original FIB)
cellType1 = struct('a', 0.02,    'b', 2.9,   'c', 1.0, 'epsilon', .0040, 'sigma', 1.0, 'g_horz', gh, 'g_vert', gv); % Healthy
cellType2 = struct('a', 0.01,    'b', 5.55,   'c', 1.0, 'epsilon', .0040, 'sigma', 0.6, 'g_horz', gh, 'g_vert', gv); % Slow Recovery
cellType3 = struct('a', 0.02,    'b', 2.9,   'c', 1.0, 'epsilon', .005, 'sigma', 1.6, 'g_horz', gh, 'g_vert', gv); % Fast Recovery

% FIX SCAR CHANGE G GRID MAYBE
cellType4 = struct('a', .02,    'b', 2.9,   'c', 1.0, 'epsilon', .0040, 'sigma', 1.0, 'g_horz', .027, 'g_vert', .015); % Non excitable (Scar)

% Block 2 Parameters
cellType5 = struct('a', 0.01,    'b', 5.3,   'c', 1.0, 'epsilon', .0040, 'sigma', 0.6, 'g_horz', .027, 'g_vert', .015); % Slow Recovery
cellType6 = struct('a', 0.02,    'b', 2.9,   'c', 1.0, 'epsilon', .006, 'sigma', 1.6, 'g_horz', .027, 'g_vert', .015); % Fast Recovery

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid size
num_cells_x = 61; 
num_cells_y = 41;

Block_Adjust = 0; % starting at 7 tall
Vert_Adjust = 0; % starting centered at (21,21)
Horz_Adjust = 30; % starting centered at (21,21)

% Grid cell type layout
cellTypes = ones(num_cells_y,num_cells_x);

cellTypes(18 - Block_Adjust + Vert_Adjust : 24 + Block_Adjust + Vert_Adjust, ...
    21 + Horz_Adjust : 22 + Horz_Adjust) = 2; % Bistable Long QT

cellTypes(18 - Block_Adjust + Vert_Adjust : 24 + Block_Adjust + Vert_Adjust, ...
    19 + Horz_Adjust : 20 + Horz_Adjust) = 3; % Fast Epsilon

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
I_ext = @(t) 4 .* (t < 50);% + 4 .* (t>180)*(t<181);

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

    % % Top and bottom edges periodic
    % up     = [v(end,:); v(1:end-1,:)];
    % down   = [v(2:end,:); v(1,:)];

    v_horz_neighbors = left + right; 
    v_vert_neighbors = up + down;

    % External stimulus
    stim = zeros(size(v));
    stim(:,1) = I_ext_func(t);

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

% Heatmap Visualization
if GenerateHeatMap == 1
    % --- User options ---
    desiredTimeSpacing  = 50; % in milliseconds (or same unit as time_steps)
    outputFolder = 'C:\Matlab\Cardiac\Cardiac Grid Outputs\PhantomFIBs';

    if ExportScreenshots && ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

    % Choose mode
    if RunAnimation
        figure;
    elseif ExportScreenshots
        % Run headless mode (no visible figure window)
        fig = figure('Visible', 'off');
    end

    % Use same loop for both modes
    for t_idx = 1:50:length(time_steps)
        % Plot the voltage heatmap
        if RunAnimation || ExportScreenshots
            clf;
            imagesc(v_history(:,:,t_idx));
            %colorbar;
            %title(sprintf('Membrane Potentials at t = %.1f ms', time_steps(t_idx)));
            %xlabel('Cell X (columns)');
            %ylabel('Cell Y (rows)');
            clim([-0.4, 1.2]);
            axis equal tight
            ax = gca;   % Get current axes
            ax.XColor = 'none';  % Hide x-axis line
            ax.YColor = 'none';  % Hide y-axis line
            ax.XTick = [];       % Remove x-axis ticks
            ax.YTick = [];       % Remove y-axis ticks
            hold on;

            % % Annotations (Normal FIB)
            rectangle('Position', [48.5, 17.5, 4, 7], ...
            'EdgeColor', 'k', 'LineWidth', 2);

            rectangle('Position', [48.5, 17.5, 4, 7], ...
            'EdgeColor', 'w', 'LineWidth', 2, 'LineStyle',':');

            line([50.5,50.5],[17.5,24.5],'Color','k','LineWidth',2)

            line([50.5,50.5],[17.5,24.5],'Color','w','LineWidth',2,'LineStyle',':')

            % Logical mask of cells where w > 0.05 & v < 0.4
            recovery_mask = (w_history(:,:,t_idx) > 0.0225) & (v_history(:,:,t_idx) < 0);
            [y, x] = find(recovery_mask);
            plot(x, y, 'x', 'Color', 'k', 'LineWidth', 1.2, 'MarkerSize', 5);
            [y, x] = find(scarmask);
            plot(x, y, 's', 'Color', 'k', 'LineWidth', 1.2, 'MarkerSize', 12);
            hold off;
        end

        % --- Export Screenshots without animation overhead ---
        if ExportScreenshots && mod(round(time_steps(t_idx)), desiredTimeSpacing) == 0
            filename = fullfile(outputFolder, sprintf('frame_%07.1fms.png', time_steps(t_idx)));
            exportgraphics(gcf, filename);
        end

        % --- Optional pause for real-time playback ---
        if RunAnimation
            pause(0.07);
        end
    end

    % Close the hidden figure if only exporting
    if ExportScreenshots && ~RunAnimation
        close(fig);
    end
end

% 3D Heatmap V-W Animation
if Generate3dMap == 1
    figure;
    for t_idx = 1:50:length(time_steps) % step every 50 frames for speed
        mesh(v_history(:,:,t_idx),w_history(:,:,t_idx),FaceColor='texturemap',FaceAlpha=.5)
        view(-37.5,60);
        c = colorbar;
        c.Label.String = 'Recovery w';
        title(sprintf('Membrane Potentials at t = %.1f ms', time_steps(t_idx)));
        xlabel('Cell X (columns)');
        ylabel('Cell Y (rows)');
        zlabel('Voltage v');
        zlim([-.4,1.2]);
        clim([-0.05, .16]);
        pause(.01)
    end
end

% Heatmap Video
if GenerateHeatMapVideo == 1
    % File path for output video
    videoFile = fullfile(pwd, 'Unnamed Cardiac Grid.mp4'); % saves in current working directory

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
        %colorbar(ax);
        title(ax, sprintf('Membrane Potentials at t = %.1f ms', time_steps(t_idx)));
        %xlabel(ax, 'Cell X (columns)');
        %ylabel(ax, 'Cell Y (rows)');
        clim(ax, [-0.2, 1.2]); % consistent color scale
        axis equal tight
        ax = gca;   % Get current axes
        ax.XColor = 'none';  % Hide x-axis line
        ax.YColor = 'none';  % Hide y-axis line
        ax.XTick = [];       % Remove x-axis ticks
        ax.YTick = [];       % Remove y-axis ticks

        hold on
        % Logical mask of cells where w > 0.05 & v < 0.4
        mask = (w_history(:,:,t_idx) > 0.0225) & (v_history(:,:,t_idx) < 0);
        [y, x] = find(mask);
        plot(x, y, 'x', 'Color', 'k', 'LineWidth', 1., 'MarkerSize', 6);
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Bistable at (21,21) and Fast Epsilon at (21,20)
if GenerateVoltage == 1
figure;
hold on
plot(time_steps, squeeze(v_history(21, 21, :)), 'b', 'LineWidth', 1.5)
plot(time_steps, squeeze(v_history(21, 20, :)), 'r', 'LineWidth', .75)
legend('Long QT', 'Fast Epsilon')
xlabel('Time (ms)')
ylabel('Membrane Potential (v)')
title('Membrane Potential at Cell (21,21)') % Add parameter details to title
hold off
end

if GenerateVWT == 1
figure
plot3(time_steps, squeeze(v_history(21, 21, :)), squeeze(w_history(21, 21, :)), 'b', 'LineWidth', 1.5);
xlabel('t - Time');
ylabel('v - Voltage')
zlabel('w - Recovery')
title('v,w,t plot')
end

if GenerateVDiff == 1
% Plot the voltage differences for 2 cells
figure;
hold on;
        plot(time_steps, squeeze(v_history(21, 21, :) - v_history(21, 20, :)), 'b', 'LineWidth', 1.5);
xlabel('Time (ms)')
ylabel('Membrane Potential (v)')
title('Voltage difference between adjacent Bistable and Fast Epsilon')
hold off;

figure;
hold on;
        plot(time_steps, squeeze(v_history(21, 22, :) - v_history(21, 19, :)), 'b', 'LineWidth', 1.5);
xlabel('Time (ms)')
ylabel('Membrane Potential (v)')
title('Voltage difference between Bistable and Fast Epsilon with a gap')
hold off; 

figure;
hold on;
        plot(time_steps, squeeze(v_history(17, 17, :) - v_history(17, 16, :)), 'b', 'LineWidth', 1.5);
xlabel('Time (ms)')
ylabel('Membrane Potential (v)')
title('Voltage difference in the path of the spiral wave')
hold off;

figure;
hold on;
        plot(time_steps, squeeze(v_history(21, 39, :) - v_history(21, 38, :)), 'b', 'LineWidth', 1.5);
xlabel('Time (ms)')
ylabel('Membrane Potential (v)')
title('Voltage difference between 2 healthy far to the right')
hold off;
end

