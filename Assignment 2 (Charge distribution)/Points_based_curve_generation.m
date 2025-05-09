clc; clear; close all;
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

%% Step 5: Predefined Set of Generated Points
% Predefined set of values for the number of generated points
generatedPointOptions = [11, 21, 51, 81, 121];

% Initialize a cell array to store the generated points for each option
generatedPoints = cell(length(generatedPointOptions), 1);

%% Step 6: Generate Points for All Values in the Set
for idx = 1:length(generatedPointOptions)
    numGeneratedPoints = generatedPointOptions(idx); % Current number of points to generate

    % Generate points on the curve or line
    if strcmpi(connectionType, 'straight')
        % Generate points along straight lines connecting the input points
        tempPoints = zeros(numGeneratedPoints, 2);
        for i = 1:numPoints-1
            % Linearly interpolate between consecutive points
            t = linspace(0, 1, ceil(numGeneratedPoints / (numPoints-1)));
            tempPoints((i-1)*length(t)+1:i*length(t), :) = ...
                (1 - t') .* points(i, :) + t' .* points(i+1, :);
        end
        % Trim excess points if necessary
        tempPoints = tempPoints(1:numGeneratedPoints, :);
    elseif strcmpi(connectionType, 'curved')
        % Use spline interpolation to generate points on a smooth curve
        t = 1:numPoints; % Parameter for interpolation
        tFine = linspace(1, numPoints, numGeneratedPoints); % Fine grid for smooth curve

        % Interpolate x and y coordinates separately
        xSmooth = spline(t, points(:, 1), tFine);
        ySmooth = spline(t, points(:, 2), tFine);
        tempPoints = [xSmooth', ySmooth'];
    end

    % Store the generated points in the cell array
    generatedPoints{idx} = tempPoints;

    %% Step 7: Plot the Points and Generated Curve/Line
    figure(idx); % Create a new figure for each set of generated points
    hold on;
    scatter(points(:, 1), points(:, 2), 'filled'); % Plot the input points
    plot(tempPoints(:, 1), tempPoints(:, 2), '-', 'DisplayName', 'Generated Curve/Line');
    title(sprintf('Generated Points on Curve/Line (N = %d)', numGeneratedPoints));
    hold off;

    % %% Step 8: Call the Charge Distribution Function
    % calculateChargeDistribution(tempPoints);
end

