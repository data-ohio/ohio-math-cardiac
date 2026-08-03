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
terminaltime = 5000; % simulation time in ms
time_steps = 0:dt:terminaltime;

% Cell Parameters (Original FIB)
cellType1 = struct('a', 0.02,    'b', 2.9,   'c', 1.0, 'epsilon', .0040, 'sigma', 1.0, 'g_horz', .027, 'g_vert', .015); % Healthy
cellType2 = struct('a', 0.01,    'b', 4.95,   'c', 1.0, 'epsilon', .0040, 'sigma', 0.6, 'g_horz', .027, 'g_vert', .015); % Slow Recovery
cellType3 = struct('a', 0.02,    'b', 2.9,   'c', 1.0, 'epsilon', .0065, 'sigma', 1.6, 'g_horz', .027, 'g_vert', .015); % Fast Recovery

% FIX SCAR CHANGE G GRID MAYBE
cellType4 = struct('a', .02,    'b', 2.9,   'c', 1.0, 'epsilon', .0040, 'sigma', 0.6, 'g_horz', .027, 'g_vert', .015); % Non excitable (Scar)

% Block 2 Parameters
cellType5 = struct('a', 0.01,    'b', 5.3,   'c', 1.0, 'epsilon', .0040, 'sigma', 0.6, 'g_horz', .027, 'g_vert', .015); % Slow Recovery
cellType6 = struct('a', 0.02,    'b', 2.9,   'c', 1.0, 'epsilon', .006, 'sigma', 1.6, 'g_horz', .027, 'g_vert', .015); % Fast Recovery

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid size
num_cells_x = 51; 
num_cells_y = 51;

Block_Adjust = 0; % starting at 7 tall
Vert_Adjust = 0; % starting centered at (21,21)
Horz_Adjust = 0; % starting centered at (21,21)

% Grid cell type layout
cellTypes = ones(num_cells_y,num_cells_x);

% cellTypes(18 - Block_Adjust + Vert_Adjust : 24 + Block_Adjust + Vert_Adjust, ...
%     21 + Horz_Adjust : 22 + Horz_Adjust) = 2; % Bistable Long QT
% 
% cellTypes(18 - Block_Adjust + Vert_Adjust : 24 + Block_Adjust + Vert_Adjust, ...
%     19 + Horz_Adjust : 20 + Horz_Adjust) = 3; % Fast Epsilon

cellTypes(13:20,10:35) = 2;
cellTypes(21:28,10:35) = 3;

% % Set up 2 interracting FIBS
% % Block 1:
% cellTypes(21:27,21:22) = 3;
% cellTypes(21:27,23:24) = 2;
% % Block 2:
% cellTypes(51:57,58:59) = 6;
% cellTypes(51:57,60:61) = 5;

% % horizontal strip
% strip_depth = 3;
% half = floor(num_cells_x/2);
% cellTypes(half-strip_depth:half,:) = 2;
% cellTypes(half+1:half+strip_depth,:) = 3;

% % vert strip
% strip_depth = 3;
% half = floor(num_cells_y/2);
% cellTypes(:,half-strip_depth+1:half) = 3;
% cellTypes(:,half+1:half+strip_depth) = 2;

% % new idea
% cellTypes(1:41,half-strip_depth+1:half) = 3;
% cellTypes(1:41,half+1:half+strip_depth) = 2;

% % Wedged FIB in scar cage
% cellTypes(15:16,1:7) = 3;
% cellTypes(17:18,1:7) = 2;

% % scar block
% cellTypes(25:28,15:22) = 4;
% cellTypes(14:17,18:25) = 4;

% % Blob arrangement:
% cellTypes(18:24,23:26) = 2; % Long
% cellTypes(20:22,12:22) = 3; % Fast

% % slight asymmetry:
% cellTypes(24,19:20) = 2;

% % slight asymmetry:
% cellTypes(24,21:22) = 1;

%%%%%%%%%%%%%%%
% % New scar
% cellTypes(17,1:24) = 4;

% Cage of non-excitable for Spiral

% cellTypes(16,8:24) = 4;
% cellTypes(17:24,24) = 4;
% cellTypes(26,17:23) = 4;

% % L shape grid
% cellTypes(14:21,21:22) = 2;
% cellTypes(20:21,23:28) = 2;
% cellTypes(14:23,19:20) = 3;
% cellTypes(22:23,21:28) = 3;

%%%%%%%%%%%%%%%

% cellTypes(21:24, 19:20) = 2; % Bistable Long QT
% cellTypes(17:20, 21:22) = 2; % Bistable Long QT
% 
% cellTypes(17:20, 19:20) = 3; % Fast Epsilon
% cellTypes(21:24, 21:22) = 3; % Fast Epsilon
% 
% cellTypes(17:22, 21) = 2;
% cellTypes(17, 21:26) = 2;
% 
% cellTypes(18:22, 22) = 3;
% cellTypes(18, 22:26) = 3;

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

% % Specific Initial Conditions
% 
% % Quadrant 1 (top left)
% v(1:21,1:21) = 0;
% w(1:21,1:21) = 0;
% 
% % Quadrant 2 (top right)
% v(1:21,22:41) = 1;
% w(1:21,22:41) = 0;
% 
% % Quadrant 3 (bottom left)
% v(22:41,1:21) = 0;
% w(22:41,1:21) = 0;
% 
% % Quadrant 4 (bottom right)
% v(22:41,22:41) = 0;
% w(22:41,22:41) = .25;

% % NEW Specific Initial Conditions
% 
% % Quadrant 1 (top left)
% v(1:21,22:41) = -.211;
% w(1:21,22:41) = .072;
% 
% % Quadrant 2 (top right)
% v(1:15,1:21) = 1;
% w(1:21,1:21) = .0833;
% 
% v(16:21,1:2) = 1;
% w(16:21,1:2) = .0833;
% 
% v(17:21,1:19) = 1;
% w(17:21,1:19) = .0833;
% 
% % Quadrant 3 (bottom left)
% v(22:41,1:21) = 0.;
% w(22:41,1:21) = 0.;
% 
% % Quadrant 4 (bottom right)
% v(22:41,22:41) = -.107;
% w(22:41,22:41) = .014;

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
    desiredTimeSpacing  = 5; % in milliseconds (or same unit as time_steps)
    outputFolder = 'C:\MyProjects\Heatmap Screenshots';

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
            % rectangle('Position', [18.5, 17.5, 4, 7], ...
            % 'EdgeColor', 'k', 'LineWidth', 2);
            % 
            % rectangle('Position', [18.5, 17.5, 4, 7], ...
            % 'EdgeColor', 'w', 'LineWidth', 2, 'LineStyle',':');
            % 
            % line([20.5,20.5],[17.5,24.5],'Color','k','LineWidth',2)
            % 
            % line([20.5,20.5],[17.5,24.5],'Color','w','LineWidth',2,'LineStyle',':')

            % % Annotiations (Scar Wedge FIB)
            % rectangle('Position', [.5, 14.5, 7, 4], ...
            % 'EdgeColor', 'k', 'LineWidth', 2);
            % 
            % rectangle('Position', [.5, 14.5, 7, 4], ...
            % 'EdgeColor', 'w', 'LineWidth', 2, 'LineStyle',':');
            % 
            % line([.5,7.5],[16.5,16.5],'Color','k','LineWidth',2)
            % 
            % line([.5,7.5],[16.5,16.5],'Color','w','LineWidth',2,'LineStyle',':')

            % % Annotations (Vertical Long FIB)
            % line([37.5,37.5],[.5,81.5],'Color','k','LineWidth',2)
            % line([37.5,37.5],[.5,81.5],'Color','w','LineWidth',2,'LineStyle',':')
            % 
            % line([40.5,40.5],[.5,81.5],'Color','k','LineWidth',2)
            % line([40.5,40.5],[.5,81.5],'Color','w','LineWidth',2,'LineStyle',':')
            % 
            % line([43.5,43.5],[.5,81.5],'Color','k','LineWidth',2)
            % line([43.5,43.5],[.5,81.5],'Color','w','LineWidth',2,'LineStyle',':')

            % % Annotations (2 FIBS interracting)
            % rectangle('Position', [20.5, 20.5, 4, 7], ...
            % 'EdgeColor', 'k', 'LineWidth', 2);
            % 
            % rectangle('Position', [20.5, 20.5, 4, 7], ...
            % 'EdgeColor', 'w', 'LineWidth', 2, 'LineStyle',':');
            % 
            % line([22.5,22.5],[20.5,27.5],'Color','k','LineWidth',2)
            % 
            % line([22.5,22.5],[20.5,27.5],'Color','w','LineWidth',2,'LineStyle',':')
            % %
            % rectangle('Position', [57.5, 50.5, 4, 7], ...
            % 'EdgeColor', 'k', 'LineWidth', 2);
            % 
            % rectangle('Position', [57.5, 50.5, 4, 7], ...
            % 'EdgeColor', 'w', 'LineWidth', 2, 'LineStyle',':');
            % 
            % line([59.5,59.5],[50.5,57.5],'Color','k','LineWidth',2)
            % 
            % line([59.5,59.5],[50.5,57.5],'Color','w','LineWidth',2,'LineStyle',':')

            % % Annotations (L-FIB)
            % line([18.5,18.5],[13.5,23.5],'Color','k','LineWidth',2)
            % line([18.5,18.5],[13.5,23.5],'Color','w','LineWidth',2,'LineStyle',':')
            % 
            % line([20.5,20.5],[13.5,21.5],'Color','k','LineWidth',2)
            % line([20.5,20.5],[13.5,21.5],'Color','w','LineWidth',2,'LineStyle',':')
            % 
            % line([22.5,22.5],[13.5,19.5],'Color','k','LineWidth',2)
            % line([22.5,22.5],[13.5,19.5],'Color','w','LineWidth',2,'LineStyle',':')
            % 
            % line([28.5,28.5],[19.5,23.5],'Color','k','LineWidth',2)
            % line([28.5,28.5],[19.5,23.5],'Color','w','LineWidth',2,'LineStyle',':')
            % 
            % line([18.5,28.5],[23.5,23.5],'Color','k','LineWidth',2)
            % line([18.5,28.5],[23.5,23.5],'Color','w','LineWidth',2,'LineStyle',':')
            % 
            % line([20.5,28.5],[21.5,21.5],'Color','k','LineWidth',2)
            % line([20.5,28.5],[21.5,21.5],'Color','w','LineWidth',2,'LineStyle',':')
            % 
            % line([22.5,28.5],[19.5,19.5],'Color','k','LineWidth',2)
            % line([22.5,28.5],[19.5,19.5],'Color','w','LineWidth',2,'LineStyle',':')
            % 
            % line([18.5,22.5],[13.5,13.5],'Color','k','LineWidth',2)
            % line([18.5,22.5],[13.5,13.5],'Color','w','LineWidth',2,'LineStyle',':')
    
            % Logical mask of cells where w > 0.05 & v < 0.4
            recovery_mask = (w_history(:,:,t_idx) > 0.0225) & (v_history(:,:,t_idx) < 0);
            [y, x] = find(recovery_mask);
            plot(x, y, 'x', 'Color', 'k', 'LineWidth', 1.2, 'MarkerSize', 5);
            [y, x] = find(scarmask);
            plot(x, y, 's', 'Color', 'k', 'LineWidth', 1.2, 'MarkerSize', 12);

            fast_mask = cellTypes == 2;
            [y, x] = find(fast_mask);
            plot(x, y, 's', 'Color', 'w', 'LineWidth', 1.2, 'MarkerSize', 12);
            slow_mask = cellTypes == 3;
            [y, x] = find(slow_mask);
            plot(x, y, 's', 'Color', 'r', 'LineWidth', 1.2, 'MarkerSize', 12);
            hold off;
        end

        % --- Export Screenshots without animation overhead ---
        if ExportScreenshots && mod(round(time_steps(t_idx)), desiredTimeSpacing) == 0
            filename = fullfile(outputFolder, sprintf('frame_%07.1fms.png', time_steps(t_idx)));
            exportgraphics(gcf, filename);
        end

        % --- Optional pause for real-time playback ---
        if RunAnimation
            pause(0.1);
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
    output = 'D:\gwalt\Desktop'; % Video Save Location
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
        % mask anotations
            recovery_mask = (w_history(:,:,t_idx) > 0.0225) & (v_history(:,:,t_idx) < 0);
            [y, x] = find(recovery_mask);
            plot(x, y, 'x', 'Color', 'k', 'LineWidth', 1.2, 'MarkerSize', 5);
            [y, x] = find(scarmask);
            plot(x, y, 's', 'Color', 'k', 'LineWidth', 1.2, 'MarkerSize', 12);

        % % Annotations (Normal FIB)
            % rectangle('Position', [18.5, 17.5, 4, 7], ...
            % 'EdgeColor', 'k', 'LineWidth', 2);
            % 
            % rectangle('Position', [18.5, 17.5, 4, 7], ...
            % 'EdgeColor', 'w', 'LineWidth', 2, 'LineStyle',':');
            % 
            % line([20.5,20.5],[17.5,24.5],'Color','k','LineWidth',2)
            % 
            % line([20.5,20.5],[17.5,24.5],'Color','w','LineWidth',2,'LineStyle',':')

        % % Annotations (Vertical Long FIB)
        %     line([37.5,37.5],[.5,81.5],'Color','k','LineWidth',2)
        %     line([37.5,37.5],[.5,81.5],'Color','w','LineWidth',2,'LineStyle',':')
        % 
        %     line([40.5,40.5],[.5,81.5],'Color','k','LineWidth',2)
        %     line([40.5,40.5],[.5,81.5],'Color','w','LineWidth',2,'LineStyle',':')
        % 
        %     line([43.5,43.5],[.5,81.5],'Color','k','LineWidth',2)
        %     line([43.5,43.5],[.5,81.5],'Color','w','LineWidth',2,'LineStyle',':')

        % % Annotations (2 FIBS interracting)
        %     rectangle('Position', [20.5, 20.5, 4, 7], ...
        %     'EdgeColor', 'k', 'LineWidth', 2);
        % 
        %     rectangle('Position', [20.5, 20.5, 4, 7], ...
        %     'EdgeColor', 'w', 'LineWidth', 2, 'LineStyle',':');
        % 
        %     line([22.5,22.5],[20.5,27.5],'Color','k','LineWidth',2)
        % 
        %     line([22.5,22.5],[20.5,27.5],'Color','w','LineWidth',2,'LineStyle',':')
        %     %
        %     rectangle('Position', [57.5, 50.5, 4, 7], ...
        %     'EdgeColor', 'k', 'LineWidth', 2);
        % 
        %     rectangle('Position', [57.5, 50.5, 4, 7], ...
        %     'EdgeColor', 'w', 'LineWidth', 2, 'LineStyle',':');
        % 
        %     line([59.5,59.5],[50.5,57.5],'Color','k','LineWidth',2)
        % 
        %     line([59.5,59.5],[50.5,57.5],'Color','w','LineWidth',2,'LineStyle',':')

        % Annotations (L-FIB)
            line([18.5,18.5],[13.5,23.5],'Color','k','LineWidth',2)
            line([18.5,18.5],[13.5,23.5],'Color','w','LineWidth',2,'LineStyle',':')

            line([20.5,20.5],[13.5,21.5],'Color','k','LineWidth',2)
            line([20.5,20.5],[13.5,21.5],'Color','w','LineWidth',2,'LineStyle',':')

            line([22.5,22.5],[13.5,19.5],'Color','k','LineWidth',2)
            line([22.5,22.5],[13.5,19.5],'Color','w','LineWidth',2,'LineStyle',':')

            line([28.5,28.5],[19.5,23.5],'Color','k','LineWidth',2)
            line([28.5,28.5],[19.5,23.5],'Color','w','LineWidth',2,'LineStyle',':')

            line([18.5,28.5],[23.5,23.5],'Color','k','LineWidth',2)
            line([18.5,28.5],[23.5,23.5],'Color','w','LineWidth',2,'LineStyle',':')

            line([20.5,28.5],[21.5,21.5],'Color','k','LineWidth',2)
            line([20.5,28.5],[21.5,21.5],'Color','w','LineWidth',2,'LineStyle',':')

            line([22.5,28.5],[19.5,19.5],'Color','k','LineWidth',2)
            line([22.5,28.5],[19.5,19.5],'Color','w','LineWidth',2,'LineStyle',':')

            line([18.5,22.5],[13.5,13.5],'Color','k','LineWidth',2)
            line([18.5,22.5],[13.5,13.5],'Color','w','LineWidth',2,'LineStyle',':')
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

% figure;
% hold on
% plot(time_steps, squeeze(v_history(21, 21, :)), 'b', 'LineWidth', 1.5)
% plot(time_steps, squeeze(w_history(21, 21, :)), 'r', 'LineWidth', .75)
% legend('v', 'w')
% xlabel('Time (ms)')
% ylabel('Membrane Potential (v)')
% title('Membrane Potential at Cell (21,21)') % Add parameter details to title
% hold off
% 
% figure;
% hold on
% plot(time_steps, squeeze(v_history(15, 8, :)), 'b', 'LineWidth', 1.5)
% plot(time_steps, squeeze(w_history(15, 8, :)), 'r', 'LineWidth', .75)
% legend('v', 'w')
% xlabel('Time (ms)')
% ylabel('Membrane Potential (v)')
% title('Membrane Potential at Cell (15,8)') % Add parameter details to title
% hold off

% figure
% plot3(eavg_history, squeeze(v_history(21, 21, :)), squeeze(w_history(21,21,:)), 'b', 'LineWidth', 1.5);
% xlabel('e avg - Environment');
% ylabel('v - Voltage')
% zlabel('w - Recovery')
% title('v,w,e plot')