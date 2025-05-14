% Parameters
r_min = 0;       % Minimum r value
r_max = 5;       % Maximum r value
r_step = 0.01;   % Step size for r
n_transient = 500; % Number of transient iterations
n_keep = 100;    % Number of iterations to keep after transients
x0 = 0.5;        % Initial condition

% Initialize
r_values = r_min:r_step:r_max;
x_values = zeros(length(r_values), n_keep);

% Main computation loop
for i = 1:length(r_values)
    r = r_values(i);
    x = x0;
    
    % Transient iterations
    for j = 1:n_transient
        x = x * exp(r*(1 - x));
    end
    
    % Store kept iterations
    for j = 1:n_keep
        x = x * exp(r*(1 - x));
        x_values(i, j) = x;
    end
end

% Plotting
figure;
hold on;
for i = 1:length(r_values)
    plot(r_values(i)*ones(1, n_keep), x_values(i, :), 'b.', 'MarkerSize', 1);
end
hold off;

% Formatting
xlabel('Parameter r');
ylabel('Population x_n');
title('Bifurcation Diagram for Ricker Map: x_{n+1} = x_n e^{r(1-x_n)}');
xlim([r_min r_max]);
grid on;

%%
% Parameters
r = 4;               % Chaotic regime for Ricker map
n_iter = 1000;       % Total iterations
n_transient = 500;   % Transient iterations to discard
d_values = 0:0.005:1; % Coupling parameter range (fine resolution for smooth plot)
x0 = 0.2;            % Initial condition for x
y0 = 0.8;            % Initial condition for y (asymmetric to avoid diagonal)

% Initialize arrays to store results
x_final = zeros(length(d_values), 100); % Store last 100 x-values for each d
y_final = zeros(length(d_values), 100); % Store last 100 y-values for each d

% Define Ricker map
f = @(x) x .* exp(r * (1 - x));

% Main loop over coupling strengths d
for i = 1:length(d_values)
    d = d_values(i);
    x = x0;
    y = y0;
    
    % Transient iterations (discarded)
    for k = 1:n_transient
        x_new = d * f(x) + (1 - d) * f(y);
        y_new = (1 - d) * f(x) + d * f(y);
        x = x_new;
        y = y_new;
    end
    
    % Store post-transient iterations
    for k = 1:100
        x_new = d * f(x) + (1 - d) * f(y);
        y_new = (1 - d) * f(x) + d * f(y);
        x = x_new;
        y = y_new;
        x_final(i, k) = x;
        y_final(i, k) = y;
    end
end

% Plot bifurcation diagram for x_n vs. d
figure;
subplot(2, 1, 1);
hold on;
for i = 1:length(d_values)
    plot(d_values(i) * ones(1, 100), x_final(i, :), 'b.', 'MarkerSize', 1);
end
xlabel('Coupling strength (d)');
ylabel('x_n');
title('Bifurcation Diagram for Coupled Ricker Maps (r=4)');
xlim([0 1]);
ylim([0 5]);
grid on;

% Plot bifurcation diagram for y_n vs. d
subplot(2, 1, 2);
hold on;
for i = 1:length(d_values)
    plot(d_values(i) * ones(1, 100), y_final(i, :), 'b.', 'MarkerSize', 1);
end
xlabel('Coupling strength (d)');
ylabel('y_n');
xlim([0 1]);
ylim([0 5]);
grid on;

