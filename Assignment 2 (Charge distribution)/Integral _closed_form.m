clear;  % Clear workspace variables
clc;    % Clear command window

% Define the lower and upper integration limits
a = -50:10:50;  % Lower limits
b = -30:10:70;  % Upper limits

% Define the x values for which the integral is computed
x_values = -20:5:30;  % x values from -20 to 30 in steps of 5

% Define the integrand as a function of x and y
integrand = @(x, y) log(abs(x - y));  % Example: log(|x - y|)

% Preallocate matrices for numerical and closed-form results
num_limits = length(a);  % Number of integration limit pairs
num_x = length(x_values);  % Number of x values
result = zeros(num_x, num_limits);  % Numerical integration results
result_closed_form = zeros(num_x, num_limits);  % Closed-form results

%% Loop over x values
for x_idx = 1:num_x
    x = x_values(x_idx);  % Current x value
    fun = @(y) integrand(x, y);  % Define the function to integrate

    % Loop over integration limits
    for i = 1:num_limits
        % Compute the numerical integral
        result(x_idx, i) = integral(fun, a(i), b(i));

        % Compute the closed-form solution (if available)
        result_closed_form(x_idx, i) = (b(i) - x) * log(abs(x - b(i)) + 1e-10) ...
                                       + a(i) - b(i) ...
                                       + (x - a(i)) * log(abs(x - a(i)) + 1e-10);
    end
end

% Check for NaN or Inf values in the results
if any(isnan(result(:))) || any(isinf(result(:)))
    warning('Numerical integration produced NaN or Inf values.');
end
if any(isnan(result_closed_form(:))) || any(isinf(result_closed_form(:)))
    warning('Closed-form calculation produced NaN or Inf values.');
end

%% Plotting the results
figure;
hold on;
for row = 1:num_x
    % Plot numerical results
    plot(result(row, :), 'b-o', 'DisplayName', ['Numerical (x = ' num2str(x_values(row)) ')']);
    % Plot closed-form results
    plot(result_closed_form(row, :), 'r--x', 'DisplayName', ['Closed-Form (x = ' num2str(x_values(row)) ')']);
end
hold off;

% Create custom x-axis labels for the (a(i), b(i)) pairs
x_labels = cell(1, num_limits);
for i = 1:num_limits
    x_labels{i} = sprintf('(%d, %d)', a(i), b(i));
end

% Set x-axis labels and adjust plot appearance
set(gca, 'XTick', 1:num_limits, 'XTickLabel', x_labels);
xlabel('Integration Limits (a(i), b(i))');
ylabel('Integral Value');
title('Numerical vs Closed-Form Integration Results');
legend show;
grid on;
xtickangle(45);  % Rotate x-axis labels for better readability

% Adjust y-axis limits to ensure all data is visible
ylim([min([result(:); result_closed_form(:)]), max([result(:); result_closed_form(:)])]);