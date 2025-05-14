% MATLAB code to calculate and plot a 2D map of the
% Largest Lyapunov Exponent (LLE, lambda_1) vs. (r, d)
% for the coupled exponential map system.
% Corresponds to "Idea 1" based on the Udwadia & Raju paper context.

clear;
close all;
clc;

fprintf('--- Script Start: 2D LLE Map Calculation ---\n');

% --- Parameters ---
% Scan Resolution (Increase for higher quality, decrease for speed)
num_r = 100;  % Number of points along r-axis
num_d = 120;  % Number of points along d-axis

% Parameter Ranges
r_min = 2.5; r_max = 5.0;
d_min = 0;   d_max = 1;

r_range = linspace(r_min, r_max, num_r);
d_range = linspace(d_min, d_max, num_d);

% Simulation Details for LE Calculation
N_trans = 1500;     % Transient iterations before LE calc
N_iter_LE = 3000;   % Iterations for LE averaging
ic = [0.1; 0.7];    % Standard off-diagonal initial condition

% --- Map Functions ---
f_exp = @(x, r_val) x .* exp(r_val * (1 - x));
df_exp = @(x, r_val) (1 - r_val * x) .* exp(r_val * (1 - x)); % Derivative f'(x)

% --- Store Results ---
LE1_map = NaN(num_r, num_d); % Initialize map with NaN

% --- Main Calculation Loop ---
% Consider using parfor here if Parallel Computing Toolbox is available
% parfor ir = 1:num_r % Uncomment for parallel execution
fprintf('Calculating 2D LLE Map (%d x %d grid)...\n', num_r, num_d); tic;
for ir = 1:num_r     % Comment out if using parfor
    r = r_range(ir);
    % Create temporary row for parallel execution if used
    % temp_row = NaN(1, num_d); % Uncomment for parfor

    for id = 1:num_d
        d = d_range(id);

        % --- Inner LE calculation (QR method) ---
        xy = ic;      % Reset IC
        Q = eye(2);   % Reset basis
        le_sum = zeros(1,2);
        sim_ok = true;

        try
            % Transient iterations
            for n=1:N_trans
                xn=xy(1); yn=xy(2);
                fxn = f_exp(xn, r); fyn = f_exp(yn, r);
                xy = [d*fxn + (1-d)*fyn; (1-d)*fxn + d*fyn];
                if any(isnan(xy))||any(isinf(xy))||max(abs(xy))>1e7; error('Diverged'); end
            end

            % LE calculation iterations
            for n=1:N_iter_LE
                xn=xy(1); yn=xy(2);

                % Jacobian
                dfxn = df_exp(xn, r); dfyn = df_exp(yn, r);
                if isnan(dfxn) || isinf(dfxn) || isnan(dfyn) || isinf(dfyn); error('Derivative unstable'); end
                J = [d*dfxn, (1-d)*dfyn; (1-d)*dfxn, d*dfyn];

                % QR step
                Z = J * Q;
                [Q_new, R] = qr(Z);
                Q = Q_new;

                % Accumulate log stretching factors
                diagR = diag(R);
                if any(abs(diagR)<eps)
                     le_term = [-Inf, -Inf];
                else
                     le_term = log(abs(diagR'));
                     if any(isnan(le_term)) || any(isinf(le_term)); error('Invalid log term'); end
                end
                le_sum = le_sum + le_term(isfinite(le_term)); % Sum only finite parts

                % Iterate map
                fxn = f_exp(xn, r); fyn = f_exp(yn, r);
                xy = [d*fxn + (1-d)*fyn; (1-d)*fxn + d*fyn];
                if any(isnan(xy))||any(isinf(xy))||max(abs(xy))>1e7; error('Diverged'); end
            end % end LE iteration loop

            % Store largest LE
            LEs = sort(le_sum / N_iter_LE, 'descend');
            current_LLE = LEs(1);
            % Check for very large negative numbers, often indicates instability before NaN
            if current_LLE < -20 % Threshold for likely instability/error
                 current_LLE = NaN;
            end
             % Assign to map (or temp row for parfor)
            LE1_map(ir, id) = current_LLE; % Direct assignment for regular for loop
            % temp_row(id) = current_LLE; % Uncomment for parfor

        catch ME
            % LE remains NaN if any error occurred
            % fprintf('Warning: Sim failed r=%.3f, d=%.3f: %s\n', r, d, ME.message);
        end
        % --- End Inner LE ---
    end % end d loop

    % Assign temp row to map if using parfor
    % LE1_map(ir, :) = temp_row; % Uncomment for parfor

    % Progress Indicator
     if mod(ir, max(1, floor(num_r/10))) == 0
        fprintf('  r row %d/%d completed (Elapsed: %.1f s)\n', ir, num_r, toc);
     end

end % end r loop

calculation_time = toc;
fprintf('Calculation finished in %.2f seconds.\n', calculation_time);

% --- Plotting ---
fprintf('Plotting results...\n');
figure('Position', [100, 100, 800, 600]);
set(gcf, 'Name', '2D LLE Map (r, d)');

imagesc(d_range, r_range, LE1_map);
axis xy; % Put d on x-axis, r on y-axis, origin at bottom-left
colorbar; % Show color scale for LLE values
xlabel('Coupling Parameter d');
ylabel('Growth Parameter r');
title('Largest Lyapunov Exponent (\lambda_1) vs. (r, d)');

% Improve colormap for visualization
% Use a diverging map centered around zero
max_abs_val = max(abs(LE1_map(:)), [], 'omitnan'); % Find max absolute value for scaling
caxis([-max_abs_val/2, max_abs_val]); % Center color axis roughly
cmap = bluewhitered(256); % Or use other diverging maps like coolwarm, bwr
colormap(cmap);

% Optional: Add contour line for LLE = 0 (boundary between chaos and stability)
hold on;
contour(d_range, r_range, LE1_map, [0 0], 'k-', 'LineWidth', 1.2, 'DisplayName', '\lambda_1 = 0');
hold off;
% legend('show'); % Show legend for contour if desired

fprintf('--- Script End ---\n');


% Helper function for bluewhitered colormap (if not built-in)
function cmap = bluewhitered(n)
    if nargin < 1, n = size(get(gcf,'colormap'),1); end
    % Define key colors: Blue -> White -> Red
    colors = [0 0 1; 1 1 1; 1 0 0];
    % Interpolate between key colors
    xi = linspace(0, 1, size(colors, 1));
    xq = linspace(0, 1, n);
    cmap = interp1(xi, colors, xq, 'linear');
end