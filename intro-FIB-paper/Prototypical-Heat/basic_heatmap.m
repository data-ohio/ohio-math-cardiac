clear


% 1. Set up output directories BEFORE simulation
parent_outdir = fullfile('/users'); % <--- Change
timestamp = string(datetime('now','Format','MM-dd-yyyy_HHmmss'));
runfolder   = fullfile(parent_outdir, ['Run_' char(timestamp)]);
tracefolder = fullfile(runfolder, 'firing_traces');


% Make sure directories exist
if ~exist(runfolder,'dir'); mkdir(runfolder); end
if ~exist(tracefolder,'dir'); mkdir(tracefolder); end


% 2. Open logfile in output runfolder before sweep
logfile = fopen(fullfile(runfolder, 'simulation_results.txt'), 'w');
if logfile == -1
    error('Could not open log file for writing!');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters and arrays


g_horz = 0.027;
g_vert = 0.015;
dt = 0.05;
terminaltime = 2500;
time_steps = 0:dt:terminaltime;


dim = 500;
de  = 0.0065/dim;
evec = 0.003:de:.0095;
db  = 1.65/dim;
bvec = 4.3:db:5.92;


% External stimulus parameters
stim_mask = zeros(41,41);   % same size as grid
stim_mask(:,1) = 1;         % stimulate only column 1
stim_amplitude = 4;         % amplitude
stim_end_time = 0.2;        % duration (ms)


% Summary matrices
C              = zeros(length(evec), length(bvec));
avg_long_full  = zeros(length(evec), length(bvec));
avg_long_final = zeros(length(evec), length(bvec));
int_dv2_full   = zeros(length(evec), length(bvec));
int_dv2_final  = zeros(length(evec), length(bvec));
beat_count     = zeros(length(evec), length(bvec));


% Start timing
overall_tic = tic;


% Set up parallel pool
if isempty(gcp('nocreate'))
    parpool;
end


fprintf('Starting parallel computation...\n');


% Use temporary struct to store results inside parfor
tmp(length(evec)*length(bvec)) = struct();


parfor index = 1:(length(evec) * length(bvec))
    [n, m] = ind2sub([length(evec), length(bvec)], index);


    % --- Define cell types
    cellType1 = struct('a',0.02,'b',2.9,'c',1.0,'epsilon',.004,'sigma',1.0);
    cellType2 = struct('a',0.01,'b',bvec(m),'c',1.0,'epsilon',.004,'sigma',0.6);
    cellType3 = struct('a',0.02,'b',2.9,'c',1.0,'epsilon',evec(n),'sigma',1.6);
    cellType4 = struct('a', 0.00,    'b', 2.9,   'c', 1.0, 'epsilon', .0040, 'sigma', 1.0); % Non excitable (Scar)
    cellType5 = struct('a',0.02,'b',2.9,'c',1.0,'epsilon',.01, 'sigma',1.6);


    num_cells_x = 41;
    num_cells_y = 41;
    cellTypes   = ones(num_cells_y,num_cells_x);
    
    cellTypes(18:24, 21:22) = 2; % Bistable Long QT
    cellTypes(18:24, 19:20) = 3; % Fast Epsilon
    % cellTypes(17, 17:24) = 4; % top scar
    % cellTypes(25, 17:24) = 4; % bottom scar


    % --- Initial conditions
    v = -0.1 * (1 + 0.01 * rand(num_cells_y, num_cells_x));
    w = zeros(num_cells_y, num_cells_x);
    v_history = zeros(num_cells_y,num_cells_x,length(time_steps));


    % --- Precompute parameter maps
    a_map= zeros(num_cells_y,num_cells_x);
    b_map= zeros(num_cells_y,num_cells_x);
    c_map= zeros(num_cells_y,num_cells_x);
    epsilon_map= zeros(num_cells_y,num_cells_x);
    sigma_map= zeros(num_cells_y,num_cells_x);
    for i = 1:num_cells_y
        for j = 1:num_cells_x
            switch cellTypes(i,j)
                case 1, ct=cellType1;
                case 2, ct=cellType2;
                case 3, ct=cellType3;
                case 4, ct=cellType4;
                otherwise, ct=cellType5;
            end
            a_map(i,j)=ct.a;
            b_map(i,j)=ct.b;
            c_map(i,j)=ct.c;
            epsilon_map(i,j)=ct.epsilon;
            sigma_map(i,j)=ct.sigma;
        end
    end


    % --- Time Integration (RK4)
    for t_idx = 1:length(time_steps)
        t = time_steps(t_idx);


        v_prev = v; w_prev = w;


        [k1v,k1w] = compute_derivatives(v_prev, w_prev, t, ...
                                        a_map,b_map,c_map,epsilon_map,sigma_map, g_horz,g_vert, stim_mask, stim_amplitude, stim_end_time);
        [k2v,k2w] = compute_derivatives(v_prev + 0.5*dt*k1v, ...
                                        w_prev + 0.5*dt*k1w, t+0.5*dt, ...
                                        a_map,b_map,c_map,epsilon_map,sigma_map, g_horz,g_vert, stim_mask, stim_amplitude, stim_end_time);
        [k3v,k3w] = compute_derivatives(v_prev + 0.5*dt*k2v, ...
                                        w_prev + 0.5*dt*k2w, t+0.5*dt, ...
                                        a_map,b_map,c_map,epsilon_map,sigma_map, g_horz,g_vert, stim_mask, stim_amplitude, stim_end_time);
        [k4v,k4w] = compute_derivatives(v_prev + dt*k3v, ...
                                        w_prev + dt*k3w, t+dt, ...
                                        a_map,b_map,c_map,epsilon_map,sigma_map, g_horz,g_vert, stim_mask, stim_amplitude, stim_end_time);


        v = v_prev + (dt/6)*(k1v+2*k2v+2*k3v+k4v);
        w = w_prev + (dt/6)*(k1w+2*k2w+2*k3w+k4w);


        v_history(:,:,t_idx) = v;
    end


    % --- Analysis for cell (21,21) 
    fin_time_start = (terminaltime - 300) / dt;
    fin_v          = squeeze(v_history(21,21,fin_time_start+1:end));
    v_longqt_full  = squeeze(v_history(21,21,:));


    avg_fin_v = mean(fin_v);
    min_fin_v = min(fin_v);
    max_fin_v = max(fin_v);
    dv_dt     = diff(fin_v)/dt;
    int_squared_deriv = sum(dv_dt.^2)*dt;


    avg_v_longqt_full = mean(v_longqt_full);
    dv_dt_full = diff(v_longqt_full)/dt;
    int_squared_deriv_full = sum(dv_dt_full.^2)*dt;


    % Beat counting
    beat_count_number = count_beats(v_longqt_full, dt);


    % Classification

    if max(fin_v) <= 0.4
        end_behavior = 'Fully Recovered';
    elseif min(fin_v) >= 0.4 
        end_behavior = 'Stayed Depolarized';
    elseif count_beats(fin_v, dt) >= 1
        end_behavior = 'Still Firing';
    else
        end_behavior = 'Unknown';
    end


    % --- Store results in temporary struct
    tmp(index).avg_long_full  = avg_v_longqt_full;
    tmp(index).avg_long_final = avg_fin_v;
    tmp(index).int_dv2_full   = int_squared_deriv_full;
    tmp(index).int_dv2_final  = int_squared_deriv;
    tmp(index).beat_count     = beat_count_number;
    tmp(index).end_behavior   = end_behavior;


    if strcmp(end_behavior,'Still Firing')
        tmp(index).v_longqt_full = v_longqt_full;
        tmp(index).time_steps    = time_steps;
    else
        tmp(index).v_longqt_full = [];
        tmp(index).time_steps    = [];
    end
end


fprintf('Completed parallel computation.\n');


% --- Consolidate results into matrices
for index = 1:(length(evec)*length(bvec))
    [n,m] = ind2sub([length(evec), length(bvec)], index);


    avg_long_full(n,m)  = tmp(index).avg_long_full;
    avg_long_final(n,m) = tmp(index).avg_long_final;
    int_dv2_full(n,m)   = tmp(index).int_dv2_full;
    int_dv2_final(n,m)  = tmp(index).int_dv2_final;
    beat_count(n,m)     = tmp(index).beat_count;


    switch tmp(index).end_behavior
        case 'Stayed Depolarized'
            C(n,m) = 1;
        case 'Still Firing'
            C(n,m) = 2;
        case 'Fully Recovered'
            C(n,m) = 3;
        case 'Unknown'
            C(n,m) = 4;
    end

end


% --- Save summary
save(fullfile(runfolder,'summary_matrices.mat'), ...
    'avg_long_full','avg_long_final','int_dv2_full','int_dv2_final','beat_count','C','bvec','evec');


% Write labeled CSVs
write_matrix_with_labels(fullfile(runfolder,'avg_long_full.csv'), avg_long_full, bvec, evec, 'bvec','evec');
write_matrix_with_labels(fullfile(runfolder,'avg_long_final.csv'),avg_long_final,bvec,evec,'bvec','evec');
write_matrix_with_labels(fullfile(runfolder,'int_dv2_full.csv'),  int_dv2_full, bvec, evec,'bvec','evec');
write_matrix_with_labels(fullfile(runfolder,'int_dv2_final.csv'), int_dv2_final,bvec,evec,'bvec','evec');
write_matrix_with_labels(fullfile(runfolder,'end_behavior_C.csv'),C,bvec,evec,'bvec','evec');
write_matrix_with_labels(fullfile(runfolder,'beat_count.csv'),   beat_count,bvec,evec,'bvec','evec');


% --- Beat count heatmap ---
figure
imagesc(bvec, evec, beat_count);
set(gca,'YDir','normal');           % make row 1 = bottom (epsilon increasing upward)
colorbar;
xlabel('B values');
ylabel('Epsilon values');
title('Long QT Cell Beat Number');


% Save images
saveas(gcf, fullfile(runfolder, 'beat_count_heatmap.png'));
saveas(gcf, fullfile(runfolder, 'beat_count_heatmap.fig'));
print(gcf, fullfile(runfolder,'beat_count_heatmap_highres'), '-dpng', '-r300');


% --- End behavior heatmap ---
colors = [0,0,1; 1,0,0; 0,1,0]; % Blue = depolarized, Red = firing, Green = recovered


figure
imagesc(bvec, evec, C);
set(gca,'YDir','normal');  % correct orientation
clim([1 3]);
colormap(colors);
colorbar;
xlabel('B values');
ylabel('Epsilon values');
title('End Behavior of System');


% Add legend
hold on
L = line(ones(3),ones(3),'LineWidth',2);
set(L,{'color'},mat2cell(colors,ones(1,3),3));
legend('Stays Depolarized','Still Firing','Recovered');
hold off
% Enable data cursor with custom tips
dcm = datacursormode(gcf);
dcm.UpdateFcn = @(obj,event) heatmapDataTip(event, bvec, evec, C);


% Save images
saveas(gcf, fullfile(runfolder, 'end_behavior_heatmap.png'));
saveas(gcf, fullfile(runfolder, 'end_behavior_heatmap.fig'));
print(gcf, fullfile(runfolder,'end_behavior_heatmap_highres'), '-dpng', '-r300');

% Completion
fprintf(logfile, 'Simulation completed.\n');

% End timing
elapsed_time = toc(overall_tic);
fprintf('Total runtime: %.2f seconds (%.2f minutes)\n', elapsed_time, elapsed_time/60);
fprintf(logfile, 'Total runtime: %.2f seconds (%.2f minutes)\n', elapsed_time, elapsed_time/60);

fclose(logfile);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper Functions


function write_matrix_with_labels(outname,M,bvec,evec,bname,ename)
    isfile = ischar(outname) || isstring(outname);
    if isfile
        fid = fopen(outname,'w'); closeAtEnd = true;
    else
        fid = outname; closeAtEnd = false;
    end
    fprintf(fid,'%s/%s,',ename,bname);
    fprintf(fid,'%.3f,',bvec); fprintf(fid,'\n');
    for n = 1:length(evec)
        fprintf(fid,'%.3f,',evec(n));
        fprintf(fid,'%.7g,',M(n,:)); fprintf(fid,'\n');
    end
    if closeAtEnd, fclose(fid); end
end


function [dv_dt,dw_dt] = compute_derivatives(v,w,t, ...
    a_map,b_map,c_map,epsilon_map,sigma_map,g_horz,g_vert, ...
    stim_mask, stim_amplitude, stim_end_time)


    % Neighbors (insulated BC by edge copying)
    left  = [v(:,1), v(:,1:end-1)];
    right = [v(:,2:end), v(:,end)];
    up    = [v(1,:); v(1:end-1,:)];
    down  = [v(2:end,:); v(end,:)];


    v_horz_neighbors = left+right;
    v_vert_neighbors = up+down;


    dv_dt = sigma_map .* (v .* (v - a_map).*(c_map - v) ...
        - w ...
        + g_horz .* (v_horz_neighbors - 2.*v) ...
        + g_vert .* (v_vert_neighbors - 2.*v));


    % vectorized stimulus application
    if t < stim_end_time
        dv_dt = dv_dt + stim_amplitude * stim_mask;
    end


    dw_dt = epsilon_map .* (v - b_map .* w);
end


function num_beats = count_beats(trace, dt)
    % Count full action potentials: up-cross >0.4 then down-cross <=0
    num_beats=0; state=0;
    refractory = round(50/dt); 
    last_beat_idx=-Inf;
    for k=2:length(trace)
        if state==0
            if trace(k-1)<0.4 && trace(k)>=0.4
                state=1;
            end
        elseif state==1
            if trace(k-1)>0 && trace(k)<=0
                if (k-last_beat_idx) > refractory
                    num_beats=num_beats+1;
                    last_beat_idx=k;
                end
                state=0;
            end
        end
    end
end