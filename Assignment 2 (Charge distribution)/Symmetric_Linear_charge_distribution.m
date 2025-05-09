clc; close all; clear;
start_time = cputime; % Get the initial CPU time
    w=0.5;
    L = 2*w;                % Length of the strip (meters)
    V0 = 1;                 % Applied potential (Volts)
    eps0 = 1;               % Permittivity of free space (normalized to 1)
    N_vals = [10,20,50,80,120]; % Increasing number of unknowns for convergence study
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
        S_initial=zeros(N,N/2);
        for m = 1:N
            x_m = x_mid(m); % Observation point
            for n = 1:N/2
                % Integration limits for the n-th segment
                a = x_n(n);
                b = x_n(n+1);
                % Analytical integral for the Green's function
                S_initial(m, n) = -(1/(2*pi*eps0)) * ((b - x_m) * log(abs(x_m - b) + 1e-10) + a - b ...
                                           + (x_m - a) * log(abs(x_m - a) + 1e-10));
            end
        end

        S=horzcat(S_initial,flipud(fliplr(S_initial)));

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
    V_strip_initial = zeros(1,length(X)/2); % Initialize potential array
    for i = 1:length(X)/2
        for j = 1:length(q)
            % Integration limits for the j-th segment
            a = x_n(j);
            b = x_n(j+1);
            % Analytical formula for the integral
            V_strip_initial(i) = V_strip_initial(i) + (-q(j) / (2 * pi * eps0)) * ...
                ((b - X(i)) * log(abs(X(i) - b) + 1e-10) + a - b ...
                + (X(i) - a) * log(abs(X(i) - a) + 1e-10));
        end
    end
  
V_strip = horzcat(V_strip_initial,fliplr(V_strip_initial));

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
V_initial = zeros(length(X),length(X)/2); % Initialize potential array
% Compute the potential at each grid point
for n = 1:N
    a = x_n(n); % Start of the segment
    b = x_n(n+1); % End of the segment
    for i = 1:numel(X)/2
        % Compute the contribution of the n-th segment to the potential at (X(i), Y(i))
        V_initial(i) = V_initial(i) - (1/(2*pi*eps0)) * integral(@(x) log(sqrt((X(i)-x).^2 + Y(i).^2)), a, b, 'AbsTol', 1e-10) * q(n);
    end
end

V = horzcat(V_initial,flipud(fliplr(V_initial)));
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
