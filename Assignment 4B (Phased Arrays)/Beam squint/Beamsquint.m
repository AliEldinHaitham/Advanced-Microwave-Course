% ======================================================================
% Comprehensive Beam Squint Analysis in Phased Arrays
% This script analyzes beam squint effects vs array size, frequency,
% bandwidth, and scan angle, comparing phase shifters and true-time delay
% ======================================================================

clc; clear; close all;

%% ====================== SYSTEM PARAMETERS ======================
fc = 50e9;               % Center frequency (50 GHz)
c = 3e8;                 % Speed of light (m/s)
lambda = c/fc;           % Wavelength
d = lambda/2;            % Half-wavelength spacing
theta0 = 30;             % Scan angle (degrees)
N = 64;                  % Number of array elements

%% ====================== ARRAY SIZE ANALYSIS ======================
arraySizes = [15, 30, 50, 80, 120]; % Different array sizes
BW_percent = 20;         % Fixed bandwidth percentage

% Calculate upper band edge frequency
fu = fc * (1 + BW_percent/200);

% Calculate beam squint angle (independent of array size)
delta_theta = abs(real(asind((fc/fu)*...
                        sind(theta0))) - theta0);

%% Plot Beam Squint vs Array Size
figure('Position', [100 100 900 500]);
plot(arraySizes, delta_theta*ones(size(arraySizes)), '-o', ...
    'Color', [0 0.4470 0.7410], 'LineWidth', 2.5, 'MarkerSize', 8, ...
    'MarkerFaceColor', [0 0.4470 0.7410]);

xlabel('Number of Array Elements (N)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Beam Squint Angle (degrees)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Beam Squint vs Array Size\n(f_c = %.1f GHz, BW = %d%%, \\theta_0 = %d°)', ...
    fc/1e9, BW_percent, theta0), ...
    'FontSize', 14, 'FontWeight', 'bold');
legend('Beam Squint Angle', 'Location', 'northeast');
grid on;
set(gca, 'FontSize', 11, 'LineWidth', 1.5, 'XTick', arraySizes);
xlim([min(arraySizes)-5 max(arraySizes)+5]);

%% Calculate beamwidths for reference
beamwidths_deg = zeros(size(arraySizes));
for idx = 1:length(arraySizes)
    beamwidths_deg(idx) = 102/(arraySizes(idx));
end

% Create second y-axis for beamwidth
figure('Position', [100 100 900 500]);
plot(arraySizes, beamwidths_deg,'-o', ...
    'Color', [0 0.4470 0.7410], 'LineWidth', 2.5, 'MarkerSize', 8, ...
    'MarkerFaceColor', [0 0.4470 0.7410]);
ylabel('Beamwidth (degrees)','FontSize', 12, 'FontWeight', 'bold');
xlabel('Number of Array Elements (N)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Bandwidth vs Array Size\n(f_c = %.1f GHz, BW = %d%%, \\theta_0 = %d°)', ...
    fc/1e9, BW_percent, theta0), ...
    'FontSize', 14, 'FontWeight', 'bold');
legend( 'Beamwidth', 'Location', 'northeast');
grid on;
set(gca, 'FontSize', 11, 'LineWidth', 1.5, 'XTick', arraySizes);
%% ====================== FREQUENCY AND BANDWIDTH ANALYSIS ======================
% Define frequencies for analysis (relative to center frequency)
freqs = [0.5*fc, fc, 1.5*fc];

% Define scan angles to analyze
scanAngles = [-45, -30, -15, 15, 30, 45]; % Different scan angles (degrees)

% Define full frequency range for sweep
freqRange = linspace(1e9, 100e9, 500); % 1-100 GHz range

% Define bandwidth percentages for analysis
BW_percentages = linspace(1, 40, 50); % 1% to 40% bandwidth

%% Calculate beam squint effects
% 1. Frequency sweep calculation
delta_theta_freq = real(asind((fc./freqRange)*sind(theta0))) - theta0;

% 2. Beam squint vs bandwidth for all analysis frequencies
delta_theta_BW = zeros(length(freqs), length(BW_percentages));
for fIdx = 1:length(freqs)
    f = freqs(fIdx);
    for bwIdx = 1:length(BW_percentages)
        % Calculate band edges
        fl = f * (1 - BW_percentages(bwIdx)/200);
        fu = f * (1 + BW_percentages(bwIdx)/200);
        % Calculate squint at band edges (keep sign)
        delta_lower = real(asind((f/fl)*sind(theta0))) - theta0;
        delta_upper = real(asind((f/fu)*sind(theta0))) - theta0;
        % Store maximum squint (keep sign)
        tempArray = [delta_lower, delta_upper];
        [~, maxIdx] = max(abs(tempArray));
        delta_theta_BW(fIdx, bwIdx) = tempArray(maxIdx);
    end
end

% 3. Beam squint vs bandwidth at fixed frequency for different scan angles
delta_theta_BW_angles = zeros(length(scanAngles), length(BW_percentages));
for aIdx = 1:length(scanAngles)
    theta = scanAngles(aIdx);
    for bwIdx = 1:length(BW_percentages)
        % Calculate band edges
        fl = fc * (1 - BW_percentages(bwIdx)/200);
        fu = fc * (1 + BW_percentages(bwIdx)/200);
        % Calculate squint at band edges (keep sign)
        delta_lower = real(asind((fc/fl)*sind(theta))) - theta;
        delta_upper = real(asind((fc/fu)*sind(theta))) - theta;
        % Store maximum squint (keep sign)
        tempArray = [delta_lower, delta_upper];
        [~, maxIdx] = max(abs(tempArray));
        delta_theta_BW_angles(aIdx, bwIdx) = tempArray(maxIdx);
    end
end

%% ====================== PLOT RESULTS ======================
% 1. Plot beam squint vs frequency for each analysis frequency
for fIdx = 1:length(freqs)
    figure('Position', [100 100 900 500]);
    hold on;
    colors = lines(length(scanAngles));
    
    for aIdx = 1:length(scanAngles)
        % Calculate beam squint for current angle (keep sign)
        delta = real(asind((freqs(fIdx)./freqRange).*sind(scanAngles(aIdx)))) - scanAngles(aIdx);
        
        plot(freqRange/1e9, delta, ...
            'Color', colors(aIdx,:), 'LineWidth', 2.5, ...
            'DisplayName', sprintf('\\theta_0 = %d°', scanAngles(aIdx)));
    end
    
    xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Beam Squint Angle (degrees)', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('Beam Squint vs Frequency (F = %.1f GHz)', ...
        freqs(fIdx)/1e9), 'FontSize', 14, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 11);
    grid on;
    xlim([1 100]);
    set(gca, 'FontSize', 11, 'LineWidth', 1.5);
end
% 2. Plot beam squint vs bandwidth for all analysis frequencies
figure('Position', [100 100 900 500]);
hold on;
colors = lines(length(freqs));

for fIdx = 1:length(freqs)
    plot(BW_percentages, delta_theta_BW(fIdx,:), ...
        'Color', colors(fIdx,:), 'LineWidth', 2.5, ...
        'DisplayName', sprintf('f = %.1f GHz', freqs(fIdx)/1e9));
end

xlabel('Bandwidth (% of frequency)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Beam Squint Angle (degrees)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Beam Squint vs Bandwidth for Different Frequencies\n(N = %d, \\theta_0 = %d°)', ...
    N, theta0), 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 11);
grid on;
set(gca, 'FontSize', 11, 'LineWidth', 1.5);


% 3. Plot beam squint vs bandwidth at fixed frequency for different scan angles
figure('Position', [100 100 900 500]);
hold on;
colors = lines(length(scanAngles));
lineStyles = {'-', '--', ':','-', '--', ':'};

for aIdx = 1:length(scanAngles)
    plot(BW_percentages, delta_theta_BW_angles(aIdx,:), ...
        'Color', colors(aIdx,:), 'LineStyle', lineStyles{aIdx}, ...
        'LineWidth', 2.5, 'DisplayName', sprintf('\\theta_0 = %d°', scanAngles(aIdx)));
end

xlabel('Bandwidth (% of frequency)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Beam Squint Angle (degrees)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Beam Squint vs Bandwidth at Different Scan Angles\n(Frequency = %.1f GHz, N = %d)', ...
    fc/1e9, N), 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 11);
grid on;
set(gca, 'FontSize', 11, 'LineWidth', 1.5);

%% ====================== RADIATION PATTERN ANALYSIS ======================
% Parameters for pattern analysis
centerFrequency_Hz = fc;
frequencyRange_Hz = linspace(0.9*fc, 1.1*fc, 100); % ±10% bandwidth
scanAngle_deg = theta0;
arraySizeElements = N;
elementSpacing_m = d;
wavelength_m = lambda;
speedOfLight_mps = c;

% Select frequencies near center frequency for pattern visualization
patternFrequencies_GHz = [centerFrequency_Hz/1e9-8, centerFrequency_Hz/1e9, centerFrequency_Hz/1e9+8];
[~, freqIndices] = min(abs(frequencyRange_Hz/1e9 - patternFrequencies_GHz'), [], 2);

% Define angle range centered on scan angle
angleRange_deg = linspace(scanAngle_deg-30, scanAngle_deg+30, 1000);

% Initialize gain matrices
phaseShifterGain = zeros(length(frequencyRange_Hz), length(angleRange_deg));
trueTimeDelayGain = zeros(length(frequencyRange_Hz), length(angleRange_deg));

% Calculate patterns for all frequencies
for freqIdx = 1:length(frequencyRange_Hz)
    currentWavelength_m = speedOfLight_mps/frequencyRange_Hz(freqIdx);
    
    % Phase shifter pattern (exhibits squint)
    spatialFrequency = sind(angleRange_deg)/currentWavelength_m - sind(scanAngle_deg)/wavelength_m;
    phaseShifterGain(freqIdx,:) = abs(sinc(arraySizeElements*elementSpacing_m*spatialFrequency)).^2;
    
    % True-time delay pattern (no squint)
    spatialFrequencyTTD = (sind(angleRange_deg) - sind(scanAngle_deg))/currentWavelength_m;
    trueTimeDelayGain(freqIdx,:) = abs(sinc(arraySizeElements*elementSpacing_m*spatialFrequencyTTD)).^2;
end

%% Plot radiation patterns
figure('Position', [100 100 900 700]);

% Phase shifters pattern
subplot(2,1,1); hold on;
phaseShifterLegendEntries = cell(1,length(patternFrequencies_GHz));

for freqIdx = 1:length(patternFrequencies_GHz)
    % Convert to dB scale
    gain_dB = 10*log10(phaseShifterGain(freqIndices(freqIdx),:));
    plot(angleRange_deg, gain_dB, 'LineWidth', 2.5);
    
    % Find and display peak angle
    [~,peakIdx] = max(gain_dB);
    peakAngle_deg = angleRange_deg(peakIdx);
    phaseShifterLegendEntries{freqIdx} = sprintf('%d GHz (%.1f°)', ...
        patternFrequencies_GHz(freqIdx), peakAngle_deg);
end

xlabel('Angle (degrees)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Gain (dB)', 'FontSize', 12, 'FontWeight', 'bold');
title('Phase Shifters: Beam Squint Effect', 'FontSize', 14, 'FontWeight', 'bold');
legend(phaseShifterLegendEntries, 'Location', 'best', 'FontSize', 11);
grid on; ylim([-40 0]);
set(gca, 'FontSize', 11, 'LineWidth', 1.5);

% True-time delay pattern
subplot(2,1,2); hold on;
trueTimeDelayLegendEntries = cell(1,length(patternFrequencies_GHz));

for freqIdx = 1:length(patternFrequencies_GHz)
    % Convert to dB scale
    gain_dB = 10*log10(trueTimeDelayGain(freqIndices(freqIdx),:));
    plot(angleRange_deg, gain_dB, 'LineWidth', 2.5);
    
    % Find and display peak angle
    [~,peakIdx] = max(gain_dB);
    peakAngle_deg = angleRange_deg(peakIdx);
    trueTimeDelayLegendEntries{freqIdx} = sprintf('%d GHz (%.1f°)', ...
        patternFrequencies_GHz(freqIdx), peakAngle_deg);
end

xlabel('Angle (degrees)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Gain (dB)', 'FontSize', 12, 'FontWeight', 'bold');
title('True-Time Delay Units: No Beam Squint', 'FontSize', 14, 'FontWeight', 'bold');
legend(trueTimeDelayLegendEntries, 'Location', 'best', 'FontSize', 11);
grid on; ylim([-40 0]);
set(gca, 'FontSize', 11, 'LineWidth', 1.5);