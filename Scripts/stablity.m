% MATLAB code to recreate Figure 3 from
% Udwadia & Raju, Physica D 111 (1998) 16-26
% Shows theoretical and numerical range of d for stable synchronization vs. r.

clear;
close all;
clc;

fprintf('--- Script Start: Recreating Figure 3 ---\n');

% --- Parameters ---
r_min = 2.5;        % Start r value (approximate from Fig 3)
r_max = 5.0;        % End r value
num_r = 250;        % Number of r points
r_values = linspace(r_min, r_max, num_r);

% For LE calculation (Panel a)
x0_single = 0.5;        % Initial condition for single map
N_trans_single = 2000;  % Transient iterations for single map
N_iter_single = 5000;   % Iterations for LE averaging

% For synchronization check (Panel b)
xy0_coupled = [0.1; 0.7]; % Initial condition for coupled map (off-diagonal)
N_trans_coupled = 2000;   % Transient iterations for coupled map
N_check_coupled = 500;    % Iterations to check sync at the end
sync_tol = 1e-5;          % Tolerance for checking |x - y| < tol
num_d_scan = 401;         % Resolution for scanning d values (0 to 1)
d_scan_values = linspace(0, 1, num_d_scan);

% --- Helper Functions ---

% Function to calculate LE for the single exponential map
function lambda = calculate_LE_single_map(r, x0, N_trans, N_iter)
    f = @(x, r_val) x .* exp(r_val .* (1 - x));
    f_prime = @(x, r_val) (1 - r_val .* x) .* exp(r_val .* (1 - x));
    x = x0;
    lambda = NaN; % Default to NaN
    try
        % Transient
        for n = 1:N_trans
            x = f(x, r);
            if isnan(x) || isinf(x) || abs(x)>1e6; error('Diverged'); end
        end
        % LE calculation
        le_sum = 0;
        valid_iter = 0;
        for n = 1:N_iter
            deriv = f_prime(x, r);
             if abs(deriv) < eps % Handle superstable points
                 le_term = -Inf;
            elseif isnan(deriv) || isinf(deriv)
                 error('Invalid derivative'); % Derivative exploded
            else
                 le_term = log(abs(deriv));
            end
            if ~isinf(le_term) % Only add finite terms to sum
                le_sum = le_sum + le_term;
                valid_iter = valid_iter + 1;
            end

            x = f(x, r); % Iterate map
            if isnan(x) || isinf(x) || abs(x)>1e6; error('Diverged'); end
        end
        if valid_iter > 0
            lambda = le_sum / valid_iter;
        elseif le_sum == -Inf % Handle case where it hit superstable immediately
            lambda = -Inf;
        end
    catch ME
      % fprintf('LE Calc Warning for r=%.4f: %s\n', r, ME.message);
      lambda = NaN; % Keep as NaN if error
    end
end

% Coupled map iteration function (redefined for clarity, same as before)
function xy_next = coupled_map(xy, d, r, f_handle)
    xn = xy(1); yn = xy(2);
    fxn = f_handle(xn, r); fyn = f_handle(yn, r);
    x_next = d*fxn + (1-d)*fyn;
    y_next = (1-d)*fxn + d*fyn;
    xy_next = [x_next; y_next];
end

% --- Calculation Part 1: Theoretical Boundaries (Panel a) ---
fprintf('Calculating theoretical boundaries (Panel a)...\n');
tic;
lambda_values = zeros(1, num_r);
d_lower_theory = NaN(1, num_r);
d_upper_theory = NaN(1, num_r);

for i = 1:num_r
    r = r_values(i);
    lambda = calculate_LE_single_map(r, x0_single, N_trans_single, N_iter_single);
    lambda_values(i) = lambda; % Store LE for potential debugging

    % Calculate boundaries only if LE is positive (chaotic uncoupled map)
    % and finite (not NaN or Inf)
    if ~isnan(lambda) && ~isinf(lambda) && lambda > 0
        exp_neg_lambda = exp(-lambda);
        d_lower_theory(i) = (1 - exp_neg_lambda) / 2;
        d_upper_theory(i) = (1 + exp_neg_lambda) / 2;
    end
     if mod(i, 25) == 0; fprintf('  a: r = %.2f done\n', r); end
end
toc;
fprintf('Theoretical calculation complete.\n');

% --- Calculation Part 2: Numerical Simulation Boundaries (Panel b) ---
fprintf('Calculating numerical boundaries (Panel b)...\n');
tic;
d_lower_sim = NaN(1, num_r);
d_upper_sim = NaN(1, num_r);
f_map = @(x, r_val) x .* exp(r_val .* (1 - x)); % Define map handle once

for i = 1:num_r
    r = r_values(i);
    is_sync = false(1, num_d_scan); % Track synchronization for each d

    for j = 1:num_d_scan
        d = d_scan_values(j);
        xy = xy0_coupled; % Reset IC
        sim_ok = true;
        try
            % Transient
            for n = 1:N_trans_coupled
                xy = coupled_map(xy, d, r, f_map);
                if any(isnan(xy)) || any(isinf(xy)) || max(abs(xy))>1e7; error('Diverged'); end
            end
            % Check sync at the end
            diff_sum = 0;
            count = 0;
            for n = 1:N_check_coupled
                 xy = coupled_map(xy, d, r, f_map);
                 if any(isnan(xy)) || any(isinf(xy)) || max(abs(xy))>1e7; error('Diverged'); end
                 % Check sync only in the last few iterations
                 if n > N_check_coupled - 50
                    diff_sum = diff_sum + abs(xy(1) - xy(2));
                    count = count + 1;
                 end
            end
            avg_diff = diff_sum / count;
            if avg_diff < sync_tol
                is_sync(j) = true;
            end
        catch ME
            sim_ok = false; % Mark simulation as failed for this (r,d)
            % fprintf('Sim Warning r=%.3f, d=%.3f: %s\n', r, d, ME.message);
        end
    end % End d loop

    % Find first and last d where synchronization occurred
    sync_indices = find(is_sync);
    if ~isempty(sync_indices)
        d_lower_sim(i) = d_scan_values(sync_indices(1));
        d_upper_sim(i) = d_scan_values(sync_indices(end));
    end

    if mod(i, 25) == 0; fprintf('  b: r = %.2f done\n', r); end
end % End r loop
toc;
fprintf('Numerical simulation complete.\n');


% --- Plotting ---
fprintf('Plotting results...\n');
figure('Position', [100, 100, 700, 800]); % Make figure taller

% Panel (a) - Theoretical
ax_a = subplot(2, 1, 1);
plot(ax_a, r_values, d_upper_theory, 'k-', 'LineWidth', 1);
hold(ax_a, 'on');
plot(ax_a, r_values, d_lower_theory, 'k-', 'LineWidth', 1);
hold(ax_a, 'off');
grid(ax_a, 'on');
xlabel(ax_a, 'r');
ylabel(ax_a, 'd');
title(ax_a, 'Theoretical Synchronization Boundaries');
ylim(ax_a, [0, 1]);
xlim(ax_a, [r_min, r_max]);
text(ax_a, 0.05, 0.9, '(a)', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold');


% Panel (b) - Numerical
ax_b = subplot(2, 1, 2);
plot(ax_b, r_values, d_upper_sim, 'k-', 'LineWidth', 1);
hold(ax_b, 'on');
plot(ax_b, r_values, d_lower_sim, 'k-', 'LineWidth', 1);
hold(ax_b, 'off');
grid(ax_b, 'on');
xlabel(ax_b, 'r');
ylabel(ax_b, 'd');
title(ax_b, 'Numerical Simulation Synchronization Boundaries');
ylim(ax_b, [0, 1]);
xlim(ax_b, [r_min, r_max]);
text(ax_b, 0.05, 0.9, '(b)', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold');

sgtitle('Figure 3 Recreation: Synchronization Range vs. r', 'FontWeight', 'bold');
fprintf('--- Script End ---\n');