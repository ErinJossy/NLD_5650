clear;
close all;
clc;

fprintf('--- Script Start ---\n');

% --- Parameters ---
r = 4.0;            % Growth parameter (fixed as per Fig 2 caption/text)
N_transient = 2000; % Iterations to discard
N_plot_phase= 5000; % Iterations to plot for phase plots (a, b)
N_plot_hopf = 10000;% Iterations to plot for Hopf loops (d)
N_plot_bif  = 200;  % Iterations to plot for bifurcation diagram (c) per d

% --- Map Functions ---
f = @(x, r) x .* exp(r * (1 - x));

% Coupled map iteration function (define once)
function xy_next = coupled_map(xy, d, r, f_handle)
    xn = xy(1);
    yn = xy(2);
    fxn = f_handle(xn, r);
    fyn = f_handle(yn, r);
    x_next = d*fxn + (1-d)*fyn;
    y_next = (1-d)*fxn + d*fyn;
    xy_next = [x_next; y_next];
end

% --- Initial Condition ---
xy0 = [0.1; 0.7]; % Generic, off-diagonal

% --- Generate Figure 1: Plots (a), (b), (c) ---
fprintf('\n--- Generating Figure 1: Plots (a, b, c) ---\n');
figure(1); % Explicitly create/select Figure 1
set(gcf, 'Position', [50, 50, 800, 750], 'Name', 'Figure 2(a,b,c)'); % Set size and name

% --- Plot (a): d = 0.16 ---
ax1 = subplot(2, 2, 1); % Top-left
d_a = 0.16;
xy = xy0;
fprintf('Calculating plot (a) for d = %.2f...\n', d_a);
% Transient and Plotting (Combined logic for brevity)
xy_transient = xy;
transient_ok = true;
try
    for n = 1:N_transient
        xy_transient = coupled_map(xy_transient, d_a, r, f);
        if any(isnan(xy_transient)) || any(isinf(xy_transient)); error('Diverged'); end
    end
catch
    transient_ok = false; warning('Transient failed for plot (a)');
end
x_a = NaN(1, N_plot_phase); y_a = NaN(1, N_plot_phase);
if transient_ok
    xy = xy_transient; % Start plotting from end of transient
    for n = 1:N_plot_phase
        try
            xy = coupled_map(xy, d_a, r, f);
            if any(isnan(xy)) || any(isinf(xy)); error('Diverged'); end
            x_a(n) = xy(1); y_a(n) = xy(2);
        catch
            warning('Plotting failed for plot (a)'); break;
        end
    end
end
plot(ax1, x_a, y_a, 'b.', 'MarkerSize', 1);
hold(ax1, 'on');
plot(ax1, [0 5], [0 5], 'k--', 'LineWidth', 0.5);
hold(ax1, 'off');
axis(ax1, [0 5 0 5]); axis(ax1, 'square');
xlabel(ax1, 'x'); ylabel(ax1, 'y'); title(ax1, sprintf('d = %.2f', d_a));
text(ax1, 0.05, 0.9, '(a)', 'Units', 'normalized', 'FontSize', 11, 'FontWeight', 'bold');
set(ax1, 'XTick', 0:1:5, 'YTick', 0:1:5);

% --- Plot (b): d = 0.21 ---
ax2 = subplot(2, 2, 2); % Top-right
d_b = 0.21;
xy = xy0; % Reset IC
fprintf('Calculating plot (b) for d = %.2f...\n', d_b);
% Transient and Plotting
xy_transient = xy; transient_ok = true;
try
    for n = 1:N_transient
        xy_transient = coupled_map(xy_transient, d_b, r, f);
        if any(isnan(xy_transient)) || any(isinf(xy_transient)); error('Diverged'); end
    end
catch
    transient_ok = false; warning('Transient failed for plot (b)');
end
x_b = NaN(1, N_plot_phase); y_b = NaN(1, N_plot_phase);
if transient_ok
    xy = xy_transient;
    for n = 1:N_plot_phase
        try
            xy = coupled_map(xy, d_b, r, f);
            if any(isnan(xy)) || any(isinf(xy)); error('Diverged'); end
            x_b(n) = xy(1); y_b(n) = xy(2);
        catch
            warning('Plotting failed for plot (b)'); break;
        end
    end
end
plot(ax2, x_b, y_b, 'b.', 'MarkerSize', 1);
hold(ax2, 'on');
plot(ax2, [0 5], [0 5], 'k--', 'LineWidth', 0.5);
hold(ax2, 'off');
axis(ax2, [0 5 0 5]); axis(ax2, 'square');
xlabel(ax2, 'x'); ylabel(ax2, 'y'); title(ax2, sprintf('d = %.2f', d_b));
text(ax2, 0.05, 0.9, '(b)', 'Units', 'normalized', 'FontSize', 11, 'FontWeight', 'bold');
set(ax2, 'XTick', 0:1:5, 'YTick', 0:1:5);


% --- Plot (c): x-y vs d with Zones ---
ax3 = subplot(2, 1, 2); % Spans bottom row
d_values_c = linspace(0.001, 0.999, 500); % Use more points, avoid exact 0/1
num_d = length(d_values_c);
fprintf('Calculating Bifurcation Diagram (c)...\n');
tic;
hold(ax3, 'on'); % Hold on for plotting points for each d

for i_d = 1:num_d
    d = d_values_c(i_d);
    xy = xy0; % Reset IC for each d
    % Transient
    transient_ok = true;
    try
        for n = 1:N_transient
            xy = coupled_map(xy, d, r, f);
            if any(isnan(xy)) || any(isinf(xy)) || max(abs(xy))>1e6; error('Diverged'); end
        end
    catch
        transient_ok = false;
    end
    % Store points for plotting difference only if transient was okay
    if transient_ok
        xy_diff = NaN(1, N_plot_bif); % Pre-allocate with NaN
        plot_ok = true;
        for n = 1:N_plot_bif
            try
                xy = coupled_map(xy, d, r, f);
                if any(isnan(xy)) || any(isinf(xy)) || max(abs(xy))>1e6; error('Diverged'); end
                xy_diff(n) = xy(1) - xy(2);
            catch
                plot_ok = false; break; % Stop plotting for this d
            end
        end
        if plot_ok
             plot(ax3, d * ones(1, N_plot_bif), xy_diff, 'b.', 'MarkerSize', 2);
        end
    end
    if mod(i_d, 50) == 0; fprintf('  c: d = %.3f done\n', d); end
end
hold(ax3, 'off');
toc;
xlabel(ax3, 'Coupling Parameter d'); ylabel(ax3, 'x - y');
ylim_c = [-5 5]; ylim(ax3, ylim_c); xlim(ax3, [0 1]); grid(ax3, 'on');
title(ax3, sprintf('Synchronicity Diagram (r = %.1f)', r));
text(ax3, 0.01, 0.95, 'c)', 'Units', 'normalized', 'FontSize', 11, 'FontWeight', 'bold');

% Add Zone Indicators to Plot (c)
zone_boundaries = [0, 0.03, 0.13, 0.21, 0.79, 0.87, 0.97, 1.0];
zone_labels = {'1', '2', '3', '4', '5', '6', '7'};
zone_centers = diff(zone_boundaries)/2 + zone_boundaries(1:end-1);
label_y_pos = ylim_c(2) * 0.85; label_y_pos_main = ylim_c(2) * 0.95;
hold(ax3, 'on');
for k = 1:length(zone_labels)
    text(ax3, zone_centers(k), label_y_pos, zone_labels{k}, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
         'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'none');
end
text(ax3, mean(zone_centers), label_y_pos_main, 'ZONES', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
     'FontSize', 11, 'FontWeight', 'bold');

line_color = [0 0.5 0.8]; % A shade of blue
line_style = '-';
line_width = 1.0;
internal_boundaries = zone_boundaries(2:end-1); % Boundaries excluding 0 and 1
for d_bnd = internal_boundaries
    plot(ax3, [d_bnd d_bnd], ylim_c, ... % Draw line from y_min to y_max at d_bnd
        'Color', line_color, ...
        'LineStyle', line_style, ...
        'LineWidth', line_width, ...
        'HandleVisibility', 'off'); % Don't show in legend if one is added later
end

hold(ax3, 'off');
fprintf('--- Figure 1 Complete ---\n');


% Parameters
r = 4; % Growth parameter for exponential map
d1 = 0.016; % Coupling parameter for left plot (zone 2)
d2 = 0.984; % Coupling parameter for right plot (zone 6)
n = 10000; % Number of iterations
transient = 5000; % Transient iterations to discard

% Define the coupled exponential map
f = @(x) x.*exp(r*(1-x));

% Initialize
x1 = zeros(1, n); y1 = zeros(1, n);
x2 = zeros(1, n); y2 = zeros(1, n);

% Initial conditions (not on diagonal)
x1(1) = 0.1; y1(1) = 0.75;
x2(1) = 0.1; y2(1) = 0.75;

% Iterate the maps
for i = 1:n-1
    % Left plot (d = 0.16)
    fx = f(x1(i));
    fy = f(y1(i));
    x1(i+1) = d1*fx + (1-d1)*fy;
    y1(i+1) = (1-d1)*fx + d1*fy;
    
    % Right plot (d = 0.84)
    fx = f(x2(i));
    fy = f(y2(i));
    x2(i+1) = d2*fx + (1-d2)*fy;
    y2(i+1) = (1-d2)*fx + d2*fy;
end

% Discard transient
x1 = x1(transient:end);
y1 = y1(transient:end);
x2 = x2(transient:end);
y2 = y2(transient:end);

% Create figure
figure('Position', [100, 100, 1200, 500]);

% Left subplot (d = 0.16)
subplot(1,2,1);
plot(x1, y1, '.', 'MarkerSize', 1, 'Color', [0 0.5 0]);
hold on;
plot([0 4], [0 4], 'r--'); % Diagonal line
xlabel('x_n');
ylabel('y_n');
title(sprintf('d = %.2f (Zone 2 to Zone 1)', d1));
axis([0 4 0 4]);
grid on;

% Right subplot (d = 0.84)
subplot(1,2,2);
% Plot alternating points to show the two loops
plot(x2(1:2:end), y2(1:2:end), '.', 'MarkerSize', 1, 'Color', [0 0 0.7]);
hold on;
plot(x2(2:2:end), y2(2:2:end), '.', 'MarkerSize', 1, 'Color', [0.7 0 0]);
plot([0 4], [0 4], 'r--'); % Diagonal line
xlabel('x_n');
ylabel('y_n');
title(sprintf('d = %.2f (Zone 6 to Zone 7)', d2));
axis([0 4 0 4]);
grid on;
sgtitle('d): Hopf Bifurcations in Coupled Exponential Map (r=4)');