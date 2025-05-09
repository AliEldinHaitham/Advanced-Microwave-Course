clc;
clear all;
close all;

%% Main Script: Generate Points and Call Current Distribution Function

% Constants
mu0 = pi*4e-7;          % Permeability of free space [H/m]
eps0 = 8.8541878128e-12; % Permittivity of free space [F/m]
c = 1/sqrt(mu0*eps0);    % Speed of light [m/s]
eta0 = sqrt(mu0/eps0);   % Intrinsic impedance [Ohm]

% Simulation Parameters
freq = 300e6;           % Frequency [Hz]
lambda = c/freq;         % Wavelength [m]
k0 = 2*pi/lambda;        % Wavenumber [rad/m]
omega = 2*pi*freq;       % Angular frequency [rad/s]
E0 = 1;                  % Incident field amplitude [V/m]
phi_inc = 0;            % Incident angle (180 deg, +x direction)

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

% Step 7: Call the Current Distribution Function
calculateCurrentDistribution(generatedPoints, k0, omega, mu0, E0, phi_inc, eta0,lambda);

%% Function: Calculate Current Distribution and Generate Plots
function calculateCurrentDistribution(generatedPoints, k0, omega, mu0, E0, phi_inc, eta0,lambda)
    start_time = cputime; % Get the initial CPU time
    N_vals = [11, 21, 51, 81, 121]; % Increasing number of unknowns for convergence study
    markers = {'-s', '-o', '-d', '-*', '-.', '-x'};

    % Loop through each set of generated points
    for idx = 1:length(generatedPoints)
        points = generatedPoints{idx}; % Current set of points
        x_n = points(:, 1); % x-values of the generated points
        y_n = points(:, 2); % y-values of the generated points

        % Calculate midpoints of the segments
        x_mid = (x_n(1:end-1) + x_n(2:end)) / 2; % Midpoints of x-values
        y_mid = (y_n(1:end-1) + y_n(2:end)) / 2; % Midpoints of y-values

        % Calculate segment lengths
        segment_lengths = sqrt(diff(x_n).^2 + diff(y_n).^2);

        % Compute incident field at midpoints (TM polarization - Ez component)
        E_inc = zeros(length(x_mid), 1);
        for i = 1:length(x_mid)
            E_inc(i) = E0 * exp(1j*k0*(cos(phi_inc)*x_mid(i) + sin(phi_inc)*y_mid(i)));
        end

        % Fill the MoM matrix using Hankel function
        Z = zeros(length(x_mid), length(x_mid)); % Initialize impedance matrix
        for row = 1:length(x_mid)
            x_obs = x_mid(row); % Observation point (x)
            y_obs = y_mid(row); % Observation point (y)
            for col = 1:length(x_mid)
                % Integration limits for the col-th segment
                a = x_n(col); % Start of the segment (x)
                b = x_n(col+1); % End of the segment (x)
                y_a = y_n(col); % Start of the segment (y)
                y_b = y_n(col+1); % End of the segment (y)

                % Distance threshold for singularity extraction
                R = sqrt((x_obs - (a+b)/2)^2 + (y_obs - (y_a+y_b)/2)^2);
                
                if R < 3*segment_lengths(col) % Near or self term
                    if row == col % Self term
                        % Small argument approximation for Hankel function
                        Z(row,col) = (omega*mu0/4)*segment_lengths(col)*...
                            (1 - 1j*(2/pi)*(log(k0*segment_lengths(col)/4) + 0.577215665 - 1));
                    else % Near term
                        if a == b
                            % Vertical segment: integrate over y_prime
                            integrand = @(y_prime) besselh(0,2,k0*sqrt((x_obs-a)^2 + (y_obs-y_prime).^2));
                            Z(row,col) = (omega*mu0/4)*integral(integrand, y_a, y_b, 'AbsTol', 1e-10);
                        else
                            % Non-vertical segment: integrate over x_prime
                            slope = (y_b - y_a) / (b - a); % Slope of the segment
                            c = y_a - slope * a; % Intercept of the segment
                            integrand = @(x_prime) besselh(0,2,k0*sqrt((x_obs-x_prime).^2 + (y_obs-(slope*x_prime+c)).^2));
                            Z(row,col) = (omega*mu0/4)*integral(integrand, a, b, 'AbsTol', 1e-10)*sqrt(1 + slope^2);
                        end
                    end
                else % Far term
                    if a == b
                        % Vertical segment: integrate over y_prime
                        integrand = @(y_prime) besselh(0,2,k0*sqrt((x_obs-a)^2 + (y_obs-y_prime).^2));
                        Z(row,col) = (omega*mu0/4)*integral(integrand, y_a, y_b, 'AbsTol', 1e-10);
                    else
                        % Non-vertical segment: integrate over x_prime
                        slope = (y_b - y_a) / (b - a); % Slope of the segment
                        c = y_a - slope * a; % Intercept of the segment
                        integrand = @(x_prime) besselh(0,2,k0*sqrt((x_obs-x_prime).^2 + (y_obs-(slope*x_prime+c)).^2));
                        Z(row,col) = (omega*mu0/4)*integral(integrand, a, b, 'AbsTol', 1e-10)*sqrt(1 + slope^2);
                    end
                end
            end
        end

        % Solve for current densities (J is N x 1)
        reg_param = 1e-12 * max(abs(diag(Z)));
        Z_reg = Z + reg_param*eye(size(Z));
        J = Z_reg \ (-E_inc); % Solve for current density
        
        % Normalize current by incident magnetic field (H_inc = E0/eta0)
        J_normalized = J * eta0;

        % Create normalized position vector (0 to 1)
        normalized_position = linspace(0, 1, length(x_mid))';

        % Plot current density magnitude and phase
        figure(1 + idx);
        subplot(2,1,1);
        plot(normalized_position, abs(J_normalized), 'b-o', 'LineWidth', 1.5);
        xlabel('Normalized Position', 'FontSize', 12);
        ylabel('|J_z / H_{inc}|', 'FontSize', 12);
        title(['Normalized Surface Current Density Magnitude (N = ' num2str(length(x_mid)) ')'], 'FontSize', 14);
        grid on;

        subplot(2,1,2);
        plot(normalized_position, unwrap(angle(J_normalized))*180/pi, 'r-o', 'LineWidth', 1.5);
        xlabel('Normalized Position', 'FontSize', 12);
        ylabel('Phase [deg]', 'FontSize', 12);
        title(['Surface Current Density Phase (N = ' num2str(length(x_mid)) ')'], 'FontSize', 14);
        grid on;
        hold off;
    end
        %% Compute and Plot Total Electric Field in the Surrounding Region
        % Create a grid around the scatterer
        x_min = min(x_n) - 2*lambda;
        x_max = max(x_n) + 2*lambda;
        y_min = min(y_n) - 2*lambda;
        y_max = max(y_n) + 2*lambda;
        
        [X_grid, Y_grid] = meshgrid(linspace(x_min, x_max, 50), linspace(y_min, y_max, 50));
        
        % Initialize total field arrays
        Ez_total = zeros(size(X_grid));
        Hx_total = zeros(size(X_grid));
        Hy_total = zeros(size(X_grid));
        
        % Compute incident field on grid
        Ez_inc = -E0 * exp(1j*k0*(cos(phi_inc)*X_grid + sin(phi_inc)*Y_grid));
        
        % Compute scattered field from current distribution
        for n = 1:length(x_mid)
            % Segment coordinates
            x1 = x_n(n);
            x2 = x_n(n+1);
            y1 = y_n(n);
            y2 = y_n(n+1);
            
            % Current on this segment
            J_seg = J(n);
            
            if x1 == x2 % Vertical segment
                % Integrate over y_prime
                for i = 1:numel(X_grid)
                    integrand = @(y_prime) besselh(0,2,k0*sqrt((X_grid(i)-x1).^2 + (Y_grid(i)-y_prime).^2));
                    Ez_total(i) = Ez_total(i) - (omega*mu0/4)*J_seg*integral(integrand, y1, y2, 'AbsTol', 1e-6);
                end
            else % Non-vertical segment
                % Slope and intercept
                m = (y2 - y1)/(x2 - x1);
                c = y1 - m*x1;
                
                % Integrate over x_prime
                for i = 1:numel(X_grid)
                    integrand = @(x_prime) besselh(0,2,k0*sqrt((X_grid(i)-x_prime).^2 + (Y_grid(i)-(m*x_prime+c)).^2)) * sqrt(1 + m^2);
                    Ez_total(i) = Ez_total(i) - (omega*mu0/4)*J_seg*integral(integrand, x1, x2, 'AbsTol', 1e-6);
                end
            end
        end
        
        % Total field is sum of incident and scattered fields
        Ez_total = Ez_total + Ez_inc;
        
        % Compute magnetic field components from Ez
        [Ey, Ex] = gradient(-Ez_total);
        Hx_total = (1j/(omega*mu0)) * Ey;
        Hy_total = (-1j/(omega*mu0)) * Ex;
        
        % Plot total electric field magnitude
        figure(2 + length(N_vals));
        pcolor(X_grid, Y_grid, abs(Ez_total));
        shading interp;
        colorbar;
        hold on;
        plot(x_n, y_n, 'k-', 'LineWidth', 2);
        title(['Total Electric Field |Ez| (N = ' num2str(length(x_mid)) ')']);
        xlabel('x [m]');
        ylabel('y [m]');
        axis equal;
        hold off;
       
        % Plot magnetic field vectors
        figure(3 + length(N_vals));
        quiver(X_grid, Y_grid, real(Hx_total), real(Hy_total), 2);
        hold on;
        plot(x_n, y_n, 'k-', 'LineWidth', 2);
        title('Magnetic Field H (Real Part)');
        xlabel('x [m]');
        ylabel('y [m]');
        axis equal;
        hold off;
        % Plot Poynting vectors (time-average power flow)
        figure(4 + length(N_vals));
        Px = 0.5*real(Ez_total .* conj(Hy_total));
        Py = -0.5*real(Ez_total .* conj(Hx_total));
        quiver(X_grid, Y_grid, Px, Py, 2);
        hold on;
        plot(x_n, y_n, 'k-', 'LineWidth', 2);
        title('Time-Average Poynting Vector');
        xlabel('x [m]');
        ylabel('y [m]');
        axis equal;
        hold off;

    elapsed_time = cputime - start_time; % Calculate the elapsed CPU time
    fprintf('Elapsed CPU time: %.4f seconds\n', elapsed_time);
end