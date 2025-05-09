clc;        % Clear the command window
clear;      % Clear the workspace
close all;  % Close all open figures

%% Step 1: Input Validation for Number of Points
% Ask the user for the number of points and validate the input
while true
    numPoints = input('Enter the number of points (must be an integer greater than 1): ');
    if isnumeric(numPoints) && numPoints > 1 && mod(numPoints, 1) == 0
        break; % Valid input, exit the loop
    else
        disp('Invalid input. Please enter an integer greater than 1.');
    end
end

%% Step 2: Initialize and Input Coordinates
% Initialize an empty matrix to store the points
points = zeros(numPoints, 2);

% Ask the user to input the coordinates of each point
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

%% Step 3: Check if the Curve is Closed
% Check if the first and last points are the same (closed curve)
if isequal(points(1, :), points(end, :))
    disp('Closed curve detected (first and last points are the same).');
else
    disp('Open curve detected (first and last points are not the same).');
end

%% Step 4: Ask User for Connection Type
% Ask the user if they want to connect the points with a straight line or a curved line
while true
    connectionType = input('Do you want to connect the points with a straight line or a curved line? \n (Enter "straight" or "curved"): ', 's');
    if strcmpi(connectionType, 'straight') || strcmpi(connectionType, 'curved')
        break; % Valid input, exit the loop
    else
        disp('Invalid input. Please enter "straight" or "curved".');
    end
end

%% Step 5: Plot the Points
figure(1);
hold on;
scatter(points(:, 1), points(:, 2), 'filled'); % Plot the points as markers
title('Connecting Points Based on User Input');

%% Step 6: Connect Points Based on User's Choice
if strcmpi(connectionType, 'straight')
    % Connect with straight lines
    plot(points(:, 1), points(:, 2), '-', 'DisplayName', 'Straight Line');
    disp('Points connected with straight lines.');

    % Check if the points form a straight line
    if ~isStraightLine(points)
        disp('Sorry, the points do not form a straight line. Charge distribution calculation aborted.');
    else
        hold off;
        % Calculate and plot the charge distribution
        calculateChargeDistribution(points);
    end
elseif strcmpi(connectionType, 'curved')
    % Use spline interpolation to create a smooth curve
    n = size(points, 1); % Number of points
    t = 1:n; % Parameter for interpolation
    tFine = linspace(1, n, 1000); % Fine grid for smooth curve

    % Interpolate x and y coordinates separately
    xSmooth = spline(t, points(:, 1), tFine);
    ySmooth = spline(t, points(:, 2), tFine);

    % Plot the smooth curve
    plot(xSmooth, ySmooth, '-', 'DisplayName', 'Smooth Curve');
    disp('Points connected with a smooth curve.');
    hold off;
end


%%  Check if Points Form a Straight Line
function isStraight = isStraightLine(points)
    % Check if all y-coordinates are the same (for horizontal line)
    if all(points(:, 2) == points(1, 2))
        isStraight = true;
    % Check if all x-coordinates are the same (for vertical line)
    elseif all(points(:, 1) == points(1, 1))
        isStraight = true;
    % Check if the slope between all consecutive points is the same
    else
        slopes = diff(points(:, 2)) ./ diff(points(:, 1));
        isStraight = all(abs(slopes - slopes(1)) < 1e-6); % Tolerance for floating-point errors
    end
end

%% Calculate and Plot Charge Distribution
function calculateChargeDistribution(points)
    % Constants
    % Extract the width w from the points
    w_x = abs(points(end, 1) - points(1, 1)) ; % Half the distance between endpoints
    w_y = abs(points(end, 2) - points(1, 2)) ; % Half the distance between endpoints along y-axis
    w = sqrt(w_x^2 + w_y^2)/2;
    start_time = cputime; % Get the initial CPU time
    L = 2*w;                % Length of the strip (meters)
    V0 = 1;                 % Applied potential (Volts)
    eps0 = 1;               % Permittivity of free space (normalized to 1)
    N_vals = [11,20,50,80,120]; % Increasing number of unknowns for convergence study
    markers = {'-s', '-o', '-d', '-*','-.','-x'};

    
    % Initialize total charge array
    total_charge = zeros(size(N_vals));

    % Loop over different values of N
    for idx = 1:length(N_vals)
        N = N_vals(idx); % Number of unknowns
        Delta = L/N;     % Segment length
        x_n = linspace(-w, w, N+1); % Segment edges
        x_mid = (x_n(1:end-1) + x_n(2:end))/2; % Midpoints

        % Fill the MoM matrix using the analytical integral
        S = zeros(N, N);
        for m = 1:N
            x_m = x_mid(m); % Observation point
            for n = 1:N
                % Integration limits for the n-th segment
                a = x_n(n);
                b = x_n(n+1);
                % Analytical integral for the Green's function
                S(m, n) = -(1/(2*pi*eps0)) * ((b - x_m) * log(abs(x_m - b) + 1e-10) + a - b ...
                                           + (x_m - a) * log(abs(x_m - a) + 1e-10));
            end
        end

        % Solve for charge densities (q is N x 1)
        q =S\ (V0*ones(N,1)) ; % Ensure V0 is a column vector of size N x 1

        % Compute total charge using the basis functions
        total_charge(idx) = sum(q) * Delta; % Sum of charge densities multiplied by segment length

        % Plot charge distribution
        figure(2); hold on;
        plot(x_mid, q,markers{idx}, 'DisplayName', ['N = ' num2str(N)]);
    end

    % Finalize charge distribution plot
    figure(2);
    xlabel('Position on Strip (m)'); ylabel('Charge Density (C/m)');
    xlim([-0.55,0.55]);
    title('Charge Distribution on PEC Strip');
    legend show;
    grid on;

    % Compute change in total charge
    deltaQ = abs(diff(total_charge));

    % Plot Convergence Study of Total Charge
    figure(3);
    plot(N_vals(2:end), deltaQ, '-o');
    xlim([0,130]);
    xlabel('Number of Unknowns (N)');
    ylabel('Change in Total Charge (C)');
    title('Convergence of Total Charge');
    grid on;

    % Normalize total charge: Q * L / eps0
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
    X = linspace(-w, w, 20*N); % Fine grid for potential evaluation
    V_strip = zeros(size(X)); % Initialize potential array
    for i = 1:length(X)
        for j = 1:length(q)
            % Integration limits for the j-th segment
            a = x_n(j);
            b = x_n(j+1);
            % Analytical formula for the integral
            V_strip(i) = V_strip(i) + (-q(j) / (2 * pi * eps0)) * ...
                ((b - X(i)) * log(abs(X(i) - b) + 1e-10) + a - b ...
                + (X(i) - a) * log(abs(X(i) - a) + 1e-10));
        end
    end

    % Compute the error in potential
    error_potential = V_strip - V0; % Error = Computed Potential - Applied Potential

    % Plot the error in potential along the strip
    figure(5);
    plot(X, error_potential,'LineWidth',1);
    xlabel('Position on Strip (m)');
    ylabel('Error in Potential (V)');
    title('Error in Potential Distribution Along the Strip');
    grid on;

    % Plot the potential along the strip
    figure(6);
    plot(X, V_strip, 'r-');
    xlabel('x (along the strip)', 'Interpreter', 'latex');
    ylabel('V(x)', 'Interpreter', 'latex');
    title('Electric Potential Along the Strip', 'Interpreter', 'latex');
    grid on;

% ============================================================
% Step 7: Compute and plot electric potential in the surrounding region
% ============================================================
[X, Y] = meshgrid(linspace(-2*L, 2*L, 50), linspace(-2*L, 2*L, 50)); % Define a 2D grid
V = zeros(size(X)); % Initialize potential array

% Compute the potential at each grid point
for n = 1:N
    a = x_n(n); % Start of the segment
    b = x_n(n+1); % End of the segment
    for i = 1:numel(X)
        % Compute the contribution of the n-th segment to the potential at (X(i), Y(i))
        V(i) = V(i) - (1/(2*pi*eps0)) * integral(@(x) log(sqrt((X(i)-x).^2 + Y(i).^2)), a, b, 'AbsTol', 1e-10) * q(n);
    end
end

% Plot the potential distribution
figure(7);
contourf(X, Y, V, 15); % Contour plot of the potential
colorbar;
hold on; % Hold the plot to add the zero contour line

% Overlay the zero potential contour line
contour(X, Y, V, [0, 0], 'LineColor', 'k', 'LineWidth', 2); % Black, thick line
hold off;
xlabel('x (m)');
ylabel('y (m)');
title('Potential Distribution around PEC Strip');

% Plot the potential distribution in 3D
figure(8);
surf(X, Y, V); % 3D surface plot of the potential
shading interp; % Smooth shading
colormap(jet); % Use a color map 
colorbar; % Add a color bar
hold on; % Hold the plot to add the zero contour line
% Define the zero potential plane
x_plane = [min(X(:)), max(X(:)), max(X(:)), min(X(:))]; % x-coordinates of the plane
y_plane = [min(Y(:)), min(Y(:)), max(Y(:)), max(Y(:))]; % y-coordinates of the plane
z_plane = [0, 0, 0, 0]; % z-coordinates (V = 0)

% Add the zero potential plane using patch
patch(x_plane, y_plane, z_plane, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Red, semi-transparent planehold off;
contour3(X, Y, V, [0, 0], 'LineColor', 'k', 'LineWidth', 2); % Black, thick line
xlabel('x (m)');
ylabel('y (m)');
zlabel('Potential (V)');
title('3D Potential Distribution around PEC Strip');

% ============================================================
% Step 8: Compute and plot electric field in the surrounding region
% ============================================================
% Compute the electric field as the negative gradient of the potential
[Ex, Ey] = gradient(-V);
% Plot the electric field distribution
figure(10);
quiver(X, Y, Ex, Ey, 2, 'LineWidth', 1); % Scale arrows by 1.5xlabel('x (m)');
ylabel('y (m)');
title('Electric Field Distribution around PEC Strip');
grid on;
axis equal;
elapsed_time = cputime - start_time; % Calculate the elapsed CPU time
fprintf('Elapsed CPU time: %.4f seconds\n', elapsed_time);
end