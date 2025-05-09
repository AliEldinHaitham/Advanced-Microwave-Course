clc;
clear all;
close all;
%% Main Script: Generate Points and Call Charge Distribution Function

% Step 1: Input Validation for Number of Points
while true
    numPoints = input('Enter the number of points (must be an integer greater than 1): ');
    if isnumeric(numPoints) && numPoints > 1 && mod(numPoints, 1) == 0
        break; % Valid input, exit the loop
    else
        disp('Invalid input. Please enter an integer greater than 1.');
    end
end

% Step 2: Initialize and Input Coordinates
points = zeros(numPoints, 2); % Initialize matrix to store points
for i = 1:numPoints
    while true
        fprintf('Enter coordinates for point %d (x, y): ', i);
        coords = input('', 's'); % Input as a string
        coords = str2num(coords); % Convert string to numeric array
        if numel(coords) == 2
            points(i, :) = coords; % Store valid coordinates
            break; % Exit the loop
        else
            disp('Invalid input. Please enter two numbers separated by a space (e.g., "1 2").');
        end
    end
end

% Step 3: Check if the Curve is Closed
if isequal(points(1, :), points(end, :))
    disp('Closed curve detected (first and last points are the same).');
else
    disp('Open curve detected (first and last points are not the same).');
end

% Step 4: Ask User for Connection Type
while true
    connectionType = input('Do you want to connect the points with a straight line or a curved line? \n (Enter "straight" or "curved"): ', 's');
    if strcmpi(connectionType, 'straight') || strcmpi(connectionType, 'curved')
        break; % Valid input, exit the loop
    else
        disp('Invalid input. Please enter "straight" or "curved".');
    end
end

% Step 5: Predefined Set of Generated Points
generatedPointOptions = [11, 21, 51, 81, 121]; % Set of values for generated points
generatedPoints = cell(length(generatedPointOptions), 1); % Cell array to store generated points

% Step 6: Generate Points for All Values in the Set
for idx = 1:length(generatedPointOptions)
    numGeneratedPoints = generatedPointOptions(idx); % Current number of points to generate

    % Initialize array to store generated points
    tempPoints = [];

    if strcmpi(connectionType, 'straight')
        % Total number of segments
        numSegments = numPoints - 1 + isequal(points(1, :), points(end, :)); % Closed curve adds a segment
        pointsPerSegment = floor(numGeneratedPoints / numSegments); % Evenly divide points per segment
        remainingPoints = numGeneratedPoints - pointsPerSegment * numSegments; % Extra points to distribute

        for i = 1:numPoints-1
            % Calculate the number of points for this segment
            nPoints = pointsPerSegment + (remainingPoints > 0); % Add extra point if needed
            if remainingPoints > 0, remainingPoints = remainingPoints - 1; end

            % Linearly interpolate between consecutive points
            t = linspace(0, 1, nPoints + 1); % `nPoints + 1` ensures no repeated points
            segmentPoints = (1 - t') .* points(i, :) + t' .* points(i+1, :); % Interpolated points on segment

            % Append all points except the last point (to avoid duplication at connections)
            tempPoints = [tempPoints; segmentPoints(1:end-1, :)];
        end

        % If closed curve, connect the last point back to the first
        if isequal(points(1, :), points(end, :))
            nPoints = pointsPerSegment + (remainingPoints > 0); % Handle extra points
            if remainingPoints > 0, remainingPoints = remainingPoints - 1; end
            t = linspace(0, 1, nPoints + 1);
            lastSegment = (1 - t') .* points(end, :) + t' .* points(1, :);
            tempPoints = [tempPoints; lastSegment];
        else
            % Add the final point explicitly for open curves
            tempPoints = [tempPoints; points(end, :)];
        end

    elseif strcmpi(connectionType, 'curved')
        % Use spline interpolation to generate points on a smooth curve
        t = 1:numPoints; % Parameter for interpolation
        tFine = linspace(1, numPoints, numGeneratedPoints + 1); % Exact number of points

        % Interpolate x and y coordinates separately
        xSmooth = spline(t, points(:, 1), tFine);
        ySmooth = spline(t, points(:, 2), tFine);
        tempPoints = [xSmooth', ySmooth'];
    end


    % Store the generated points in the cell array
    generatedPoints{idx} = tempPoints;
end

figure(1); % Create a new figure for each set of generated points
    hold on;
    scatter(points(:, 1), points(:, 2), 'filled'); % Plot the input points
    plot(tempPoints(:, 1), tempPoints(:, 2), '-','linewidth',2, 'DisplayName', 'Generated Curve/Line');
    title('Generated Curve/Line');
    hold off;
% Step 7: Call the Charge Distribution Function
calculateChargeDistribution(generatedPoints);


%% Function: Calculate Charge Distribution and Generate Plots
function calculateChargeDistribution(generatedPoints)
    % Constants
    start_time = cputime; % Get the initial CPU time
    V0 = 1; % Applied potential (Volts)
    eps0 = 1; % Permittivity of free space (normalized to 1)
    N_vals = [11, 21, 51, 81, 121]; % Increasing number of unknowns for convergence study
    markers = {'-s', '-o', '-d', '-*', '-.', '-x'};

    % Initialize total charge array
    total_charge = zeros(length(generatedPoints), 1);

    % Loop through each set of generated points
    for idx = 1:length(generatedPoints)
        points = generatedPoints{idx}; % Current set of points
        x_n = points(:, 1); % x-values of the generated points
        y_n = points(:, 2); % y-values of the generated points

        % Calculate midpoints of the segments
        x_mid = (x_n(1:end-1) + x_n(2:end)) / 2; % Midpoints of x-values
        y_mid = (y_n(1:end-1) + y_n(2:end)) / 2; % Midpoints of y-values

        % Fill the MoM matrix using the modified integral
        S = zeros(length(x_mid), length(x_mid)); % Initialize matrix
        for row = 1:length(x_mid)
            x_obs = x_mid(row); % Observation point (x)
            y_obs = y_mid(row); % Observation point (y)
            for col = 1:length(x_mid)
                % Integration limits for the col-th segment
                a = x_n(col); % Start of the segment (x)
                b = x_n(col+1); % End of the segment (x)
                y_a = y_n(col); % Start of the segment (y)
                y_b = y_n(col+1); % End of the segment (y)

                if a == b
                    % Vertical segment: integrate over y_prime
                    integrand = @(y_prime) log(sqrt((x_obs - a+1e-10)^2 + (y_obs - y_prime+1e-10).^2));
                    S(row, col) = -(1 / (2 * pi * eps0)) * integral(integrand, y_a, y_b, 'AbsTol', 1e-10);
                else
                    % Non-vertical segment: integrate over x_prime
                    slope = (y_b - y_a) / (b - a); % Slope of the segment
                    c = y_a - slope * a; % Intercept of the segment
                    integrand = @(x_prime) log(sqrt((x_obs - x_prime+1e-10).^2 + (y_obs+1e-10 - (slope * x_prime + c)).^2)) * sqrt(1 + slope^2);
                    S(row, col) = -(1 / (2 * pi * eps0)) * integral(integrand, a, b, 'AbsTol', 1e-10);
                end
            end
        end

        % Solve for charge densities (q is N x 1)
        q = S \ (V0 * ones(length(x_mid), 1)); % Ensure V0 is a column vector

        % Compute total charge by summing q * segment lengths
        segment_lengths = sqrt(diff(x_n).^2 + diff(y_n).^2); % Lengths of all segments
        total_charge(idx) = sum(q .* segment_lengths); % Sum of charge densities multiplied by segment lengths

        % Plot charge distribution
        axis_plot = linspace(0, 1, N_vals(idx));
        figure(2);
        hold on;
        plot(axis_plot, q, markers{idx}, 'DisplayName', ['N = ' num2str(length(x_mid))]);
    end

    % Finalize charge distribution plot
    figure(2);
    xlabel('Position on Strip (m)');
    ylabel('Charge Density (C/m)');
    title('Charge Distribution on PEC Strip');
    legend show;
    grid on;

    % Compute change in total charge
    deltaQ = abs(diff(total_charge));

    % Plot Convergence Study of Total Charge
    figure(3);
    plot(N_vals(2:end), deltaQ, '-o');
    xlim([0, 130]);
    xlabel('Number of Unknowns (N)');
    ylabel('Change in Total Charge (C)');
    title('Convergence of Total Charge');
    grid on;

    % Normalize total charge: Q * L / eps0
    L = sqrt((x_n(end) - x_n(1))^2 + (y_n(end) - y_n(1))^2); % Length of the strip
    normalized_charge = total_charge * L / eps0;

    % Compute 1/N
    inverse_N = 1 ./ N_vals;

    % Plot QL/eps0 vs 1/N
    figure(4);
    plot(inverse_N, normalized_charge, '-o', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('1/N');
    ylabel('QL / \epsilon_0');
    title('Normalized Convergence Plot');
    grid on;

    % Approximate Asymptotic Total Charge
    asymp_charge = total_charge(end);
    disp(['Asymptotic Total Charge: ', num2str(asymp_charge), ' C']);

    % Compute the potential V_strip(x) along the strip
    X = linspace(min(x_n), max(x_n), 5 * length(x_mid)); % Fine grid for potential evaluation
    Y = linspace(min(y_n), max(y_n), 5 * length(y_mid)); % Fine grid for potential evaluation
    V_strip = zeros(size(X)); % Initialize potential array
    for i = 1:length(X)
        for j = 1:length(q)
            % Integration limits for the j-th segment
            a = x_n(j);
            b = x_n(j+1);
            y_a = y_n(j);
            y_b = y_n(j+1);

            if a == b
                % Vertical segment: integrate over y_prime
                integrand = @(y_prime) log(sqrt((Y(i) - y_prime).^2));
                V_strip(i) = V_strip(i) - (q(j) / (2 * pi * eps0)) * integral(integrand, y_a, y_b, 'AbsTol', 1e-10);
            else
                % Non-vertical segment: integrate over x_prime
                slope = (y_b - y_a) / (b - a); % Slope of the segment
                c = y_a - slope * a; % Intercept of the segment
                integrand = @(x_prime) log(sqrt((X(i) - x_prime).^2 + (Y(i) - (slope * x_prime + c)).^2)) * sqrt(1 + slope^2);
                V_strip(i) = V_strip(i) - (q(j) / (2 * pi * eps0)) * integral(integrand, a, b, 'AbsTol', 1e-10);
            end
        end
    end

    % Compute the error in potential
    error_potential = V_strip - V0; % Error = Computed Potential - Applied Potential

    % Plot the error in potential along the strip
    figure(5);
    axis_plot = linspace(0, 1, length(error_potential));
    plot(axis_plot, error_potential, 'LineWidth', 1);
    xlabel('Position on Strip (m)');
    ylabel('Error in Potential (V)');
    title('Error in Potential Distribution Along the Strip');
    grid on;

    % Plot the potential along the strip
    figure(6);
    axis_plot = linspace(0, 1, length(V_strip));
    plot(axis_plot, V_strip, 'r-');
    xlabel('x (along the strip)', 'Interpreter', 'latex');
    ylabel('V(x)', 'Interpreter', 'latex');
    title('Electric Potential Along the Strip', 'Interpreter', 'latex');
    grid on;

    % Compute and plot electric potential in the surrounding region
    [X_grid, Y_grid] = meshgrid(linspace(min(x_n) - 1, max(x_n) + 1, 50), linspace(min(y_n) - 1, max(y_n) + 1, 50));
    V_grid = zeros(size(X_grid)); % Initialize potential array

    % Compute the potential at each grid point
    for n = 1:length(x_mid)
        a = x_n(n); % Start of the segment
        b = x_n(n+1); % End of the segment
        y_a = y_n(n); % Start of the segment (y)
        y_b = y_n(n+1); % End of the segment (y)

        for i = 1:numel(X_grid)
            if a == b
                % Vertical segment: integrate over y_prime
                integrand = @(y_prime) log(sqrt((X_grid(i) - a)^2 + (Y_grid(i) - y_prime).^2));
                V_grid(i) = V_grid(i) - (q(n) / (2 * pi * eps0)) * integral(integrand, y_a, y_b, 'AbsTol', 1e-10);
            else
                % Non-vertical segment: integrate over x_prime
                slope = (y_b - y_a) / (b - a); % Slope of the segment
                c = y_a - slope * a; % Intercept of the segment
                integrand = @(x_prime) log(sqrt((X_grid(i) - x_prime).^2 + (Y_grid(i) - (slope * x_prime + c)).^2)) * sqrt(1 + slope^2);
                V_grid(i) = V_grid(i) - (q(n) / (2 * pi * eps0)) * integral(integrand, a, b, 'AbsTol', 1e-10);
            end
        end
    end

    % Plot the potential distribution
    figure(7);
    contourf(X_grid, Y_grid, V_grid, 15); % Contour plot of the potential
    colorbar;
    hold on;
    plot(x_n, y_n, 'k-', 'LineWidth', 2); % Plot the strip
    contour(X_grid, Y_grid, V_grid, [0, 0], 'LineColor', 'k', 'LineWidth', 2); % Zero potential contour
    hold off;
    xlabel('x (m)');
    ylabel('y (m)');
    title('Potential Distribution around PEC Strip');

    % Plot the potential distribution in 3D
    figure(8);
    surf(X_grid, Y_grid, V_grid); % 3D surface plot of the potential
    shading interp;
    colormap(jet);
    colorbar;
    hold on;
    % Define the zero potential plane
    x_plane = [min(X_grid(:)), max(X_grid(:)), max(X_grid(:)), min(X_grid(:))]; % x-coordinates of the plane
    y_plane = [min(Y_grid(:)), min(Y_grid(:)), max(Y_grid(:)), max(Y_grid(:))]; % y-coordinates of the plane
    z_plane = [0, 0, 0, 0]; % z-coordinates (V = 0)
    patch(x_plane, y_plane, z_plane, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Red, semi-transparent plane
    contour3(X_grid, Y_grid, V_grid, [0, 0], 'LineColor', 'k', 'LineWidth', 2); % Black, thick line
    hold off;
    xlabel('x (m)');
    ylabel('y (m)');
    zlabel('Potential (V)');
    title('3D Potential Distribution around PEC Strip');

    % Compute and plot electric field in the surrounding region
    [Ex, Ey] = gradient(-V_grid);
    figure(9);
    quiver(X_grid, Y_grid, Ex, Ey, 2, 'LineWidth', 1);
    xlabel('x (m)');
    ylabel('y (m)');
    title('Electric Field Distribution around PEC Strip');
    grid on;
    axis equal;

    elapsed_time = cputime - start_time; % Calculate the elapsed CPU time
    fprintf('Elapsed CPU time: %.4f seconds\n', elapsed_time);
end