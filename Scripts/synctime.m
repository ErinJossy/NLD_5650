% MATLAB code to calculate and plot Synchronization Time vs. Coupling d
% Specifically for Zone IV (approx. d=0.21 to d=0.79 for r=4).

clear;
close all;
clc;

fprintf('--- Script Start: Synchronization Time vs. Coupling in Zone IV ---\n');

% --- Parameters ---
r = 4.0;            % Growth parameter (fixed, r=4)
d_min_zone4 = 0.21; % Approximate START boundary of Zone IV
d_max_zone4 = 0.79; % Approximate END boundary of Zone IV
num_d = 150;        % Number of d points to test within Zone IV
d_range_zone4 = linspace(d_min_zone4, d_max_zone4, num_d); % d values within Zone IV

ic = [0.1; 0.7];    % Standard off-diagonal initial condition
sync_tol = 1e-6;    % Tolerance: abs(x-y) < sync_tol means synchronized
max_iter_sync = 30000; % Maximum number of iterations to wait

% --- Map Function ---
f_exp = @(x, r_val) x .* exp(r_val * (1 - x));

% Coupled map iteration function
function xy_next = coupled_map_exp(xy, d, r_val, f_handle)
    xn = xy(1);
    yn = xy(2);
    if isnan(xn) || isnan(yn) || isinf(xn) || isinf(yn)
         error('NaN or Inf encountered during calculation.');
    end
    fxn = f_handle(xn, r_val);
    fyn = f_handle(yn, r_val);
    x_next = d*fxn + (1-d)*fyn;
    y_next = (1-d)*fxn + d*fyn;
    xy_next = [x_next; y_next];
end

% --- Store Results ---
sync_time = NaN(1, num_d); % Initialize with NaN

% --- Calculation Loop ---
fprintf('Calculating Synchronization Time for r=%.1f across Zone IV (d=%.2f to %.2f)...\n', ...
        r, d_min_zone4, d_max_zone4);
tic; % Start timer

for id = 1:num_d
    d = d_range_zone4(id); % Use d values only from Zone IV range
    xy = ic; % Reset initial condition for each d
    n_sync = NaN;
    sim_ok = true;

    try
        for n = 1:max_iter_sync
            xy = coupled_map_exp(xy, d, r, f_exp); % Iterate

            % Check for synchronization
            if abs(xy(1) - xy(2)) < sync_tol
                n_sync = n; % Record iteration count
                break;      % Exit inner loop for this d
            end
        end % End iteration loop (n)

    catch ME
        sim_ok = false;
        % fprintf('Warning: Simulation failed for d=%.4f: %s\n', d, ME.message);
        % n_sync remains NaN
    end

    sync_time(id) = n_sync; % Store result

    % Progress indicator
    if mod(id, max(1,floor(num_d/10))) == 0
        fprintf('  d = %.3f completed (%d/%d). Sync Time: %d\n', d, id, num_d, n_sync);
    end

end % End d loop

calculation_time = toc;
fprintf('Calculation finished in %.2f seconds.\n', calculation_time);

% --- Plotting ---
fprintf('Plotting results for Zone IV...\n');
figure('Position', [100, 100, 800, 500]);
set(gcf, 'Name', 'Synchronization Time in Zone IV');

plot(d_range_zone4, sync_time, 'b.-', 'MarkerSize', 10); % Blue line with dots

grid on;
xlabel('Coupling Parameter d');
ylabel('Iterations to Synchronize (N_{sync})');
title(sprintf('Synchronization Time vs. Coupling Strength within Zone IV (r=%.1f)', r));

% --- Set plot limits EXPLICITLY to Zone IV boundaries ---
xlim([d_min_zone4, d_max_zone4]);
ylim(bottom=0); % Ensure y-axis starts at 0

% Add text for points that didn't sync (should ideally not happen in Zone IV)
nan_indices = find(isnan(sync_time));
if ~isempty(nan_indices)
    hold on;
    plot(d_range_zone4(nan_indices), zeros(size(nan_indices)), 'rx', 'MarkerSize', 8, 'LineWidth', 1.5);
    legend({'Sync Time', sprintf('No Sync within %d iter', max_iter_sync)}, 'Location', 'best');
    hold off;
    fprintf('Warning: %d points within Zone IV did not synchronize within %d iterations.\n', length(nan_indices), max_iter_sync);
else
    legend({'Sync Time'}, 'Location', 'best');
end

fprintf('--- Script End ---\n');