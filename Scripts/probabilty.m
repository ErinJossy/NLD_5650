% MATLAB code to recreate Figure 5 from
% Udwadia & Raju, Physica D 111 (1998) 16-26
% Plots probability of periodic/synchronous behavior vs. r
% by scanning coupling parameter d for each r.

clear;
close all;
clc;

fprintf('--- Script Start: Recreating Figure 5 ---\n');

% --- Parameters ---
% Simulation Quality (Increase for final plot, decrease for testing)
num_r = 100;        % Number of r points (Paper likely used ~200-300)
num_d_scan = 500;   % Number of d values scanned (Paper used 5000)
num_ic = 3;         % Number of initial conditions (Paper used 3)

% Simulation Details
r_min = 0.07;       % Start r value (approximate from Fig 5)
r_max = 4.9;        % End r value
r_values = linspace(r_min, r_max, num_r);
d_scan_values = linspace(0.001, 0.999, num_d_scan); % Avoid exact 0, 1

N_transient = 1500; % Iterations to discard
N_analyze = 1000;   % Iterations to analyze for sync/periodicity
P_check = 100;      % Window size (last P points) for checking behavior

% Detection Tolerances
sync_tol = 1e-5;    % Tolerance for average |x-y| < tol for sync
period_tol = 1e-6;  % Tolerance for std(x) and std(y) < tol for periodic

% Initial Conditions (Off-diagonal)
ic_list = {[0.1; 0.7], [0.2; 0.9], [0.8; 0.3]}; % List of ICs

% --- Map Functions ---
f_exp = @(x, r_val) x .* exp(r_val .* (1 - x));

% Coupled map iteration function
function xy_next = coupled_map_exp(xy, d, r_val, f_handle)
    xn = xy(1);
    yn = xy(2);
    fxn = f_handle(xn, r_val);
    fyn = f_handle(yn, r_val);
    x_next = d*fxn + (1-d)*fyn;
    y_next = (1-d)*fxn + d*fyn;
    xy_next = [x_next; y_next];
end

% --- Initialization ---
prob_periodic = zeros(1, num_r);
prob_sync = zeros(1, num_r);
prob_total = zeros(1, num_r);

% --- Main Calculation Loop ---
fprintf('Starting calculation (%d r values, %d d scans, %d ICs)...\n', num_r, num_d_scan, num_ic);
total_sims_per_r = num_d_scan * num_ic;
tic;

% Use PARFOR if Parallel Computing Toolbox is available
% parfor i_r = 1:num_r % Uncomment for parallel execution
for i_r = 1:num_r % Comment out if using parfor
    r = r_values(i_r);
    count_periodic_current_r = 0;
    count_sync_current_r = 0;
    count_either_current_r = 0;
    sims_ok_current_r = 0; % Count successful simulations

    for k_ic = 1:num_ic
        xy0 = ic_list{k_ic};

        for j_d = 1:num_d_scan
            d = d_scan_values(j_d);
            xy = xy0; % Reset IC

            % Store last P points for analysis
            history_x = zeros(1, P_check);
            history_y = zeros(1, P_check);
            history_diff = zeros(1, P_check);
            transient_ok = true;
            analysis_ok = true;

            try
                % Transient
                for n = 1:N_transient
                    xy = coupled_map_exp(xy, d, r, f_exp);
                    if any(isnan(xy)) || any(isinf(xy)) || max(abs(xy))>1e7; error('Diverged'); end
                end

                % Analysis phase - store history
                idx_hist = 1; % Circular buffer index (simple approach)
                for n = 1:N_analyze
                    xy = coupled_map_exp(xy, d, r, f_exp);
                    if any(isnan(xy)) || any(isinf(xy)) || max(abs(xy))>1e7; error('Diverged'); end
                    % Store only the last P points
                    history_x(idx_hist) = xy(1);
                    history_y(idx_hist) = xy(2);
                    history_diff(idx_hist) = abs(xy(1) - xy(2));
                    idx_hist = mod(idx_hist, P_check) + 1; % Move index, wrap around
                end

            catch ME
                transient_ok = false; % Mark simulation as failed
                analysis_ok = false;
                % fprintf('Warning: Sim failed r=%.3f, d=%.3f, ic=%d: %s\n', r, d, k_ic, ME.message);
            end

            % --- Classify Behavior (only if simulation finished) ---
            if transient_ok && analysis_ok
                sims_ok_current_r = sims_ok_current_r + 1; % Count successful sim
                is_sync = false;
                is_periodic = false;

                % Check for Synchronization
                avg_diff = mean(history_diff);
                if avg_diff < sync_tol
                    is_sync = true;
                    count_sync_current_r = count_sync_current_r + 1;
                    count_either_current_r = count_either_current_r + 1;
                end

                % Check for Periodicity (only if NOT synchronized)
                if ~is_sync
                    std_x = std(history_x);
                    std_y = std(history_y);
                    % Check if BOTH x and y settled (fixed point or periodic)
                    if std_x < period_tol && std_y < period_tol
                        is_periodic = true;
                        count_periodic_current_r = count_periodic_current_r + 1;
                        count_either_current_r = count_either_current_r + 1;
                    end
                end
            end % end if sim ok
        end % end d loop
    end % end ic loop

    % Calculate probabilities for this r
    % Denominator could be total attempts or only successful sims.
    % Using total attempts (total_sims_per_r) matches the idea of
    % picking a random d/ic and seeing the outcome, including divergence.
    if total_sims_per_r > 0
        prob_sync(i_r) = count_sync_current_r / total_sims_per_r;
        prob_periodic(i_r) = count_periodic_current_r / total_sims_per_r;
        prob_total(i_r) = count_either_current_r / total_sims_per_r; % Sync OR Periodic
    end

    % Progress update
    if mod(i_r, max(1, floor(num_r/10))) == 0 || i_r == num_r
        fprintf('  r = %.3f (%d/%d) completed. Sync: %.3f, Per: %.3f, Tot: %.3f (Elapsed: %.1f s)\n', ...
                r, i_r, num_r, prob_sync(i_r), prob_periodic(i_r), prob_total(i_r), toc);
    end

end % end r loop

calculation_time = toc;
fprintf('Calculation finished in %.2f seconds.\n', calculation_time);

% --- Plotting ---
fprintf('Plotting results...\n');
figure('Position', [100, 100, 800, 550]);
hold on;

% Plot data - use styles similar to paper
plot(r_values, prob_total, 'k-', 'LineWidth', 1.5', 'DisplayName', 'Total (Sync or Periodic)');
plot(r_values, prob_sync, 'k-.', 'LineWidth', 1.2, 'DisplayName', 'Synchronous orbit'); % Dash-dot
plot(r_values, prob_periodic, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Periodic orbit'); % Dashed

hold off;
grid on;
xlabel('Growth Parameter r');
ylabel('Probability');
title({'Probability of Periodic/Synchronous Behavior vs. r', ...
       sprintf('(Based on %d d-values and %d ICs per r)', num_d_scan, num_ic)});
legend('show', 'Location', 'best');
ylim([0, 1.05]); % Probability range
xlim([r_min, r_max]);

fprintf('--- Script End ---\n');