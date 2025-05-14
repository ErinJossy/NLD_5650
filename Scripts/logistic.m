% MATLAB code to recreate Figure 4 from
% Udwadia & Raju, Physica D 111 (1998) 16-26
% Shows bifurcation diagrams (x, y, x-y) vs. d for the coupled logistic map.

clear;
close all;
clc;

fprintf('--- Script Start: Recreating Figure 4 ---\n');

% --- Parameters ---
a = 3.7;            % Growth rate for logistic map (specified in paper)
x0 = 0.1;           % Initial condition x0
y0 = 0.75;          % Initial condition y0
xy0 = [x0; y0];     % Combined initial condition vector

d_min = 0;          % Start d value
d_max = 1;          % End d value
num_d = 500;        % Number of d points to simulate
d_values = linspace(d_min, d_max, num_d); % Range of coupling parameter d

N_transient = 2000; % Iterations to discard for transient
N_plot_bif  = 300;  % Iterations to plot for bifurcation diagram per d

% --- Logistic Map and Coupled Map Functions ---
f_logistic = @(z, a_val) a_val .* z .* (1 - z);

% Coupled map iteration function (using logistic map)
function xy_next = coupled_map_logistic(xy, d, a_val, f_handle)
    xn = xy(1);
    yn = xy(2);
    fxn = f_handle(xn, a_val); % Use logistic map handle
    fyn = f_handle(yn, a_val); % Use logistic map handle
    x_next = d*fxn + (1-d)*fyn;
    y_next = (1-d)*fxn + d*fyn;
    xy_next = [x_next; y_next];
end

% --- Setup Figure ---
figure('Position', [100, 100, 700, 850]); % Taller figure for 3 plots
set(gcf, 'Name', 'Figure 4 Recreation');

% Get axes handles - allows plotting into specific subplots from the loop
ax_a = subplot(3, 1, 1); % Top plot for x vs d
ax_b = subplot(3, 1, 2); % Middle plot for y vs d
ax_c = subplot(3, 1, 3); % Bottom plot for x-y vs d

hold(ax_a, 'on'); grid(ax_a, 'on');
hold(ax_b, 'on'); grid(ax_b, 'on');
hold(ax_c, 'on'); grid(ax_c, 'on');

% --- Main Calculation and Plotting Loop ---
fprintf('Calculating Bifurcation Diagrams...\n');
tic;

for i_d = 1:num_d
    d = d_values(i_d);
    xy = xy0; % Reset IC for each d

    % Transient calculation
    transient_ok = true;
    try
        for n = 1:N_transient
            xy = coupled_map_logistic(xy, d, a, f_logistic);
            % Check for divergence (logistic map usually stays bounded for a=3.7)
            if any(isnan(xy)) || any(isinf(xy)) || any(xy < -0.1) || any(xy > 1.1)
                error('Diverged during transient');
            end
        end
    catch ME
        % fprintf('Warning: Transient failed for d=%.4f: %s\n', d, ME.message);
        transient_ok = false;
    end

    % Plotting iterations (only if transient was okay)
    if transient_ok
        % Pre-allocate temporary storage (optional, can plot directly)
        x_points = NaN(1, N_plot_bif);
        y_points = NaN(1, N_plot_bif);
        diff_points = NaN(1, N_plot_bif);
        plot_ok = true;

        for n = 1:N_plot_bif
            try
                xy = coupled_map_logistic(xy, d, a, f_logistic);
                if any(isnan(xy)) || any(isinf(xy)) || any(xy < -0.1) || any(xy > 1.1)
                     error('Diverged during plotting');
                end
                 % Store points *after* iteration
                x_points(n) = xy(1);
                y_points(n) = xy(2);
                diff_points(n) = xy(1) - xy(2);
            catch ME
                % fprintf('Warning: Plotting failed for d=%.4f: %s\n', d, ME.message);
                plot_ok = false;
                break; % Stop plotting for this d
            end
        end

        % Plot the collected points if plotting went okay
        if plot_ok
            plot(ax_a, d * ones(1, N_plot_bif), x_points, 'k.', 'MarkerSize', 2);
            plot(ax_b, d * ones(1, N_plot_bif), y_points, 'k.', 'MarkerSize', 2);
            plot(ax_c, d * ones(1, N_plot_bif), diff_points, 'k.', 'MarkerSize', 2);
        end
    end % end if transient_ok

    % Progress Indicator
    if mod(i_d, 50) == 0
        fprintf('  d = %.3f completed (%d/%d)\n', d, i_d, num_d);
    end

end % end main d loop

toc;
fprintf('Calculation complete. Finalizing plots...\n');

% --- Finalize Plots ---
hold(ax_a, 'off');
hold(ax_b, 'off');
hold(ax_c, 'off');

% Panel (a) settings
ylabel(ax_a, 'x');
title(ax_a, sprintf('Coupled Logistic Map (r = %.1f)', a));
ylim(ax_a, [0 1]); % Logistic map range
xlim(ax_a, [d_min d_max]);
set(ax_a, 'XTickLabel', []); % Remove x-axis numbers for top plot

% Panel (b) settings
ylabel(ax_b, 'y');
ylim(ax_b, [0 1]);
xlim(ax_b, [d_min d_max]);
set(ax_b, 'XTickLabel', []); % Remove x-axis numbers for middle plot

% Panel (c) settings
ylabel(ax_c, 'x - y');
ylim_c = [-1 1]; % Expected range for difference
ylim(ax_c, ylim_c);
xlim(ax_c, [d_min d_max]);
xlabel(ax_c, 'Coupling Parameter d');

% Add Zone Indicators to Plot (c) - Using same boundaries as Fig 2 for consistency
zone_boundaries = [0, 0.03, 0.13, 0.21, 0.79, 0.87, 0.97, 1.0];
zone_labels = {'1', '2', '3', '4', '5', '6', '7'};
zone_centers = diff(zone_boundaries)/2 + zone_boundaries(1:end-1);
label_y_pos = ylim_c(2) * 0.80; % Adjust vertical position for this plot's scale
label_y_pos_main = ylim_c(2) * 0.90;
hold(ax_c, 'on');
for k = 1:length(zone_labels)
    text(ax_c, zone_centers(k), label_y_pos, zone_labels{k}, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
         'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'none');
end
text(ax_c, mean(zone_centers), label_y_pos_main, 'ZONES', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
     'FontSize', 11, 'FontWeight', 'bold');

line_color = [0 0.5 0.8]; % A shade of blue
line_style = '-';
line_width = 1.0;
internal_boundaries = zone_boundaries(2:end-1); % Boundaries excluding 0 and 1
for d_bnd = internal_boundaries
    plot(ax_c, [d_bnd d_bnd], ylim_c, ... % Draw line from y_min to y_max at d_bnd
        'Color', line_color, ...
        'LineStyle', line_style, ...
        'LineWidth', line_width, ...
        'HandleVisibility', 'off'); % Don't show in legend if one is added later
end

for k = 1:length(zone_labels)
    text(ax_c, zone_centers(k), label_y_pos, zone_labels{k}, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
         'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'none');
end
text(ax_c, mean(zone_centers), label_y_pos_main, 'ZONES', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
     'FontSize', 11, 'FontWeight', 'bold');
hold(ax_c, 'off');

% Optional super title
% sgtitle(gcf, 'Figure 4 Recreation: Coupled Logistic Map Dynamics');

fprintf('--- Script End ---\n');