%% Complete Quadrature Coupler Analysis with 3D Surface Plot
clc; clear; close all;

%% Parameters (91 points ensures exact 90° phase)
deltaA_dB = linspace(0, 8, 181);       % Amplitude imbalance (0 to 5 dB)
deltaPhi_deg = linspace(0, 180, 181);    % Phase imbalance (0° to 180°)
[DeltaA, DeltaPhi] = meshgrid(deltaA_dB, deltaPhi_deg);

%% Calculations
Ex = 1;                                 % Reference amplitude
Ey = 10.^(-DeltaA/20);                  % Imbalanced amplitude
DeltaPhi_rad = deg2rad(DeltaPhi);       % Convert to radians

% Axial Ratio Calculation
AR_num = 1 + sqrt(1 - 4*(Ex.*Ey.*sin(DeltaPhi_rad)).^2 ./ (Ex.^2 + Ey.^2).^2);
AR_den = 1 - sqrt(1 - 4*(Ex.*Ey.*sin(DeltaPhi_rad)).^2 ./ (Ex.^2 + Ey.^2).^2);
AR_dB = 10*log10(AR_num ./ AR_den);

%% Plotting
% 1. Magnitude Imbalance vs AR at exactly 90° phase
figure;
plot(deltaA_dB, AR_dB(90,:), 'b', 'LineWidth', 2); % 90° at index 46
xlabel('Magnitude Imbalance (dB)');
ylabel('Axial Ratio (dB)');
title('AR vs Magnitude at 90° Phase');
grid on;
set(gca, 'FontSize', 12);

% 2. Phase Imbalance vs AR at exactly 0 dB magnitude
figure;
plot(deltaPhi_deg, AR_dB(:,1), 'r', 'LineWidth', 2); % 0 dB at index 1
xlabel('Phase Imbalance (degrees)');
ylabel('Axial Ratio (dB)');
xlim([0,180]);
title('AR vs Phase at 0 dB Magnitude difference');
grid on;
set(gca, 'FontSize', 12);

% 3. 3D Surface Plot
figure;
surf(DeltaA, DeltaPhi, AR_dB, 'EdgeColor', 'none');
shading interp;
colormap(jet);
colorbar;
xlabel('Magnitude Imbalance (dB)');
ylim([0,180]);
ylabel('Phase Imbalance (degrees)');
zlabel('Axial Ratio (dB)');
title('3D Surface: AR vs Magnitude and Phase Imbalance');
view(45, 30); % Adjust viewing angle
set(gca, 'FontSize', 12);
