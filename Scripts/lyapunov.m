% MATLAB code to recreate Figure 6 from
% Udwadia & Raju, Physica D 111 (1998) 16-26
% Calculates and plots Lyapunov exponents vs. d for coupled exponential map.

clear;
close all;
clc;

fprintf('--- Script Start: Recreating Figure 6 ---\n');

% --- Parameters ---
r = 4.0;            % Growth parameter (specified in Fig 6 caption)
num_d = 401;        % Number of d points (adjust resolution as needed)
d_values = linspace(0, 1, num_d); % Range of coupling parameter d

N_transient = 2000; % Iterations to discard (let trajectory settle)
N_iterate = 10000;  % Iterations for LE calculation (longer is better)

% Initial conditions
ic_a = [0.2; 0.8];  % Off-diagonal for Fig 6(a) (e.g., paper's x0 != y0)
ic_b = [0.5; 0.5];  % On-diagonal for Fig 6(b) (x0 = y0)

% --- Map function and its derivative ---
f_exp = @(x, r_val) x .* exp(r_val * (1 - x));
df_exp = @(x, r_val) (1 - r_val * x) .* exp(r_val * (1 - x)); % Derivative f'(x)

% --- Store results ---
LEs_a = NaN(num_d, 2); % For off-diagonal ICs (initialize with NaN)
LEs_b = NaN(num_d, 2); % For on-diagonal ICs

% --- Main Calculation Loop ---
fprintf('Calculating Lyapunov Exponents (r=%.1f)...\n', r);
tic; % Start timer

for i_d = 1:num_d
    d = d_values(i_d);

    % --- Case (a): Off-diagonal IC ---
    xy = ic_a;
    Q = eye(2); % Initial orthonormal basis for tangent space
    le_sum_a = zeros(1, 2);
    sim_ok_a = true;

    try
        % Transient iterations
        for n = 1:N_transient
            xn = xy(1); yn = xy(2);
            fxn = f_exp(xn, r); fyn = f_exp(yn, r);
            xy = [d*fxn + (1-d)*fyn; (1-d)*fxn + d*fyn];
            if any(isnan(xy)) || any(isinf(xy)) || max(abs(xy))>1e7; error('Diverged'); end
        end

        % Iterations for LE calculation
        for n = 1:N_iterate
            xn = xy(1); yn = xy(2);

            % Calculate Jacobian J = [J11 J12; J21 J22]
            dfxn = df_exp(xn, r); dfyn = df_exp(yn, r);
            if isnan(dfxn) || isinf(dfxn) || isnan(dfyn) || isinf(dfyn); error('Derivative unstable'); end
            J11 = d * dfxn;     J12 = (1 - d) * dfyn;
            J21 = (1 - d) * dfxn; J22 = d * dfyn;
            J = [J11, J12; J21, J22];

            % Evolve tangent vectors and QR decomposition
            Z = J * Q;
            [Q_new, R] = qr(Z); % Q_new is the new basis
            Q = Q_new;         % Update basis for next step

            % Check R diagonal elements before taking log
            diagR = diag(R);
            if any(abs(diagR) < eps) % Check for values too close to zero
                 le_term = [-Inf, -Inf]; % Assign -Inf if singular
            else
                 le_term = log(abs(diagR')); % Sum log of stretching factors
                 if any(isnan(le_term)) || any(isinf(le_term)) % Check result of log
                     error('Invalid log term');
                 end
            end
            % Only add finite parts to the sum (handles -Inf case)
            le_sum_a = le_sum_a + le_term(isfinite(le_term));


            % Iterate the map
            fxn = f_exp(xn, r); fyn = f_exp(yn, r);
            xy = [d*fxn + (1-d)*fyn; (1-d)*fxn + d*fyn];
            if any(isnan(xy)) || any(isinf(xy)) || max(abs(xy))>1e7; error('Diverged'); end
        end
        % Calculate final average LEs
        LEs_a(i_d, :) = sort(le_sum_a / N_iterate, 'descend');

    catch ME
        % fprintf('Warning: Problem during calculation (a) for d=%.4f: %s\n', d, ME.message);
        sim_ok_a = false; % LEs for this d remain NaN
    end

    % --- Case (b): On-diagonal IC ---
    xy = ic_b;
    Q = eye(2); % Reset basis
    le_sum_b = zeros(1, 2);
    sim_ok_b = true;

    try
        % Transient iterations (trajectory stays on diagonal)
        for n = 1:N_transient
            xn = xy(1); % xn=yn
            xy = [f_exp(xn, r); f_exp(xn, r)]; % Simplified iteration
            if any(isnan(xy)) || any(isinf(xy)) || max(abs(xy))>1e7; error('Diverged'); end
        end

        % Iterations for LE calculation
        for n = 1:N_iterate
            xn = xy(1); % xn=yn

            % Calculate Jacobian J (using xn for both derivatives)
            dfxn = df_exp(xn, r);
            if isnan(dfxn) || isinf(dfxn); error('Derivative unstable'); end
            J11 = d * dfxn;     J12 = (1 - d) * dfxn;
            J21 = (1 - d) * dfxn; J22 = d * dfxn;
            J = [J11, J12; J21, J22];

            % Evolve tangent vectors and QR decomposition
            Z = J * Q;
            [Q_new, R] = qr(Z);
            Q = Q_new;

            % Check R diagonal elements
            diagR = diag(R);
            if any(abs(diagR) < eps)
                 le_term = [-Inf, -Inf];
            else
                 le_term = log(abs(diagR'));
                 if any(isnan(le_term)) || any(isinf(le_term))
                     error('Invalid log term');
                 end
            end
            le_sum_b = le_sum_b + le_term(isfinite(le_term));

            % Iterate the map (staying on diagonal)
            xy = [f_exp(xn, r); f_exp(xn, r)];
            if any(isnan(xy)) || any(isinf(xy)) || max(abs(xy))>1e7; error('Diverged'); end
        end
        LEs_b(i_d, :) = sort(le_sum_b / N_iterate, 'descend');

    catch ME
        % fprintf('Warning: Problem during calculation (b) for d=%.4f: %s\n', d, ME.message);
        sim_ok_b = false; % LEs remain NaN
    end

    % Progress indicator
    if mod(i_d, max(1,floor(num_d/10))) == 0
        fprintf('  d = %.3f completed (%d/%d)\n', d, i_d, num_d);
    end

end % End d loop

calculation_time = toc;
fprintf('Calculation finished in %.2f seconds.\n', calculation_time);

% --- Plotting ---
fprintf('Plotting results...\n');
figure('Position', [100, 100, 650, 800]); % Make figure taller
set(gcf, 'Name', 'Figure 6 Recreation');

% --- Subplot (a) - Off-diagonal IC ---
ax_a = subplot(2, 1, 1);
plot(ax_a, d_values, LEs_a(:, 1), 'k-', 'LineWidth', 1); % LE1
hold(ax_a, 'on');
plot(ax_a, d_values, LEs_a(:, 2), 'k-', 'LineWidth', 1); % LE2
plot(ax_a, d_values, zeros(size(d_values)), 'k:'); % Zero line for reference
hold(ax_a, 'off');
grid(ax_a, 'on');
ylim(ax_a, [-5, 1.1]); % Match paper's y-limits, slight padding
xlim(ax_a, [0, 1]);
ylabel(ax_a, '\lambda');
title(ax_a, sprintf('r = %.1f', r));
text(ax_a, 0.05, 0.1, ' ', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom');
set(ax_a, 'XTickLabel', []); % Remove x-axis labels for top plot

% --- Subplot (b) - On-diagonal IC ---
ax_b = subplot(2, 1, 2);
plot(ax_b, d_values, LEs_b(:, 1), 'k-', 'LineWidth', 1); % LE1
hold(ax_b, 'on');
plot(ax_b, d_values, LEs_b(:, 2), 'k-', 'LineWidth', 1); % LE2
plot(ax_b, d_values, zeros(size(d_values)), 'k:'); % Zero line
hold(ax_b, 'off');
grid(ax_b, 'on');
ylim(ax_b, [-5, 1.1]); % Match paper's y-limits, slight padding
xlim(ax_b, [0, 1]);
xlabel(ax_b, 'Coupling Parameter d');
ylabel(ax_b, '\lambda');
text(ax_b, 0.05, 0.1, ' ', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom');

% Optional overall title
% sgtitle('Figure 6 Recreation: Lyapunov Exponents vs. d (Coupled Exponential Map)');

fprintf('--- Script End ---\n');