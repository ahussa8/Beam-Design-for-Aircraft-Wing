clear;clc;close all;
% Constants
g = 9.81; % acceleration due to gravity (m/s^2)
L = 4; % length of each wing in meters
max_takeoff_weight = 1113; % maximum takeoff weight in kg
max_load_factor = 4.5; % maximum load factor
wing_load = 24560; % maximum load on one wing N
N = 10000;

%Part 1
% (a)
% Calculate w0 for distributed load function
w0 = 4 * wing_load / (pi * L); % calculate w0 in Newtons per meter
fprintf("w0(N/m): %f\n", 4 * wing_load / (pi * 4));

% Distributed load function as a handle
w_x = @(x) w0 * sqrt(1 - (x/L).^2);

% (b) Maximum internal shear force and bending moment using numerical integration
PloadsArray = zeros(1, N);
xAxis = linspace(0, L, N + 1);
deltaX = L/N;
for i = 1 : N
   PloadsArray(i) = (w_x(xAxis(i)) + w_x(xAxis(i+1))) / 2 * deltaX;
end
V_x = zeros(1, N + 1);
M_x = zeros(1, N + 1);
for i = 1:N
    for j = i:N
        V_x(i) = V_x(i) + PloadsArray(j);
        M_x(i) = M_x(i) + PloadsArray(j) * (xAxis(j) + deltaX / 2 - xAxis(i));
    end
end

quiver(xAxis(1:end-1), 0*zeros(1,N), 0*PloadsArray, PloadsArray, 0);
hold on
plot( [0, L], [0, 0], 'k' )
xlabel('Load Position (m)')
ylabel('Load Magnitude (N)')
title('Point Loads')

% Plot the internal shear force diagram
figure;
area(xAxis, V_x, 'LineWidth', 2);
title('Internal Shear Force in the Beam');
xlabel('Position (m)');
ylabel('Shear Force (kN)');

% Plot the internal bending moment diagram
figure;
area(xAxis, M_x, 'LineWidth', 2);
title('Internal Bending Moment in the Beam');
xlabel('Position (m)');
ylabel('Bending Moment (kNm)')

% Calculate maximum shear force and bending moment at the root of the wing (x=L)
max_shear_force = V_x(1);
max_bending_moment = M_x(1);

fprintf('Maximum shear force (N): %f\n', round(max_shear_force));
fprintf('Maximum bending moment (Nm): %f\n', round(max_bending_moment));

% (c) Spar Design for strength
material_yield_strength = 414e6; % 2014-T6 Aluminum yield strength in Pascals
safety_factor = 1.5;

% Constants
bh = 0.2;% Total height of the beam in m
ww = 0.002; % Web thickness in m
fw = 0.1; % Flange width in m
%graph stress vs Falnge thickness
tfVals = linspace(0, 10e-3, N);
stress = linspace(0, 10e-3, N);

for i = 1:length(tfVals)
    fh = tfVals(i);
    wh = bh - (2 * fh);
    A_flange = fw * fh;
    A_web = ww * wh;
    y_NA = ((fh + wh + 0.5 * fh) * A_flange + (fh + 0.5 * wh) * A_web + 0.5 * fh * A_flange) / (2 * A_flange + A_web);
    d_flange = abs(y_NA - (fh + wh + 0.5 * fh));
    d_web = abs(y_NA - (fh + 0.5 * wh));
    I = 2 * ((1/12) * fw * fh^3 + A_flange * d_flange^2) + (1/12) * ww * wh^3 + A_web * d_web^2;
    stress(i) = ((max_bending_moment*(bh/2)) / (I));
end


allowable_stress = material_yield_strength / safety_factor;
figure;
plot(tfVals, stress, 'LineWidth', 2);
title('Stress vs Flange Thickness');
xlabel('Flange Thickness (m)');
ylabel('Stress (MPa)');

hold on;
plot(tfVals, tfVals*0 + allowable_stress)

fh = 0;

%find where thickness and stress intersect on the two lines
idx = find(stress <= allowable_stress, 1, 'first');
if ~isempty(idx)
    fh = tfVals(idx);
    fprintf('Minimum safe flange thickness (m): %f\n', fh);
else
    fprintf('No safe flange thickness found within the range.\n');
end

%Von Mises stress
%For Web
wh = bh - (2 * fh);
A_flange = fw * fh;
A_web = ww * wh;
y_NA = ((fh + wh + 0.5 * fh) * A_flange + (fh + 0.5 * wh) * A_web + 0.5 * fh * A_flange) / (2 * A_flange + A_web);
d_flange = abs(y_NA - (fh + wh + 0.5 * fh));
d_web = abs(y_NA - (fh + 0.5 * wh));
I = 2 * ((1/12) * fw * fh^3 + A_flange * d_flange^2) + (1/12) * ww * wh^3 + A_web * d_web^2;
Q = ((fh/2)+(wh/2)) * fw * fh;
tau = (max_shear_force * Q) / (I * ww); 
stress = ((max_bending_moment*(bh/2)) / (I));
vonMisesStressFlange = (stress)*1e-6; % Von Mises Stress
vonMisesStressWeb = sqrt(3*(tau)^2)*1e-6;
fprintf("Maximum von-Mises stress for Flange(MPa) %f\n" , vonMisesStressFlange);
fprintf("Maximum von-Mises stress for Web(MPa): %f\n" , vonMisesStressWeb);

if((vonMisesStressFlange < allowable_stress) && (vonMisesStressWeb < allowable_stress))
    fprintf('Selected flange thickness is sufficient to prevent yielding of the spar with the required safety factor because %f < %f and %f < %f.\n', vonMisesStressFlange, allowable_stress*1e-6, vonMisesStressWeb, allowable_stress*1e-6);
end

% (e) Weight calculation of the spar
density_aluminum = 2789; % density of 2014-T6 Aluminum in kg/m^3
volume_flange = fw * fh * L * 2; % Volume of both flanges
volume_web = ww * (bh - (2 * fh)) * L; % Volume of the web
total_volume = volume_flange + volume_web;
total_weight = total_volume * density_aluminum;

fprintf('Total weight of the spar (kg): %f\n', total_weight);

%Part 2
%(b)
% Constants for stress calculation
I = 2 * ((1/12) * fw * fh^3 + A_flange * d_flange^2) + (1/12) * ww * wh^3 + A_web * d_web^2; % Moment of inertia for a rectangular section minus web
c = bh / 2; % Distance from the neutral axis to the top or bottom fiber of the beam

% Calculate bending stress and safety factor at each point
bending_stress = (abs(M_x) * c) / I;
safety_factor = material_yield_strength./bending_stress;

% Plotting the safety factor along the spar
figure;
plot(xAxis(1:end-1), safety_factor(1:end-1), 'b', 'LineWidth', 2); % Adjust the x-axis to match the size of safety_factor array
title('Safety Factor Along the Wing Spar');
xlabel('Position along the spar (m)');
ylabel('Safety Factor');
ylim([0, 12]); % Set y-limits to visually enhance the plot
xlim([0, 4]); % Set y-limits to visually enhance the plot

%Part 3
% Constants
bh = 0.2; % Total height of the beam in meters
ww = 0.002; % Web thickness in meters 
wh = bh - (2 * fh);
allowable_stress = material_yield_strength / 1.5;
total_volume_flange = 0;
total_volume_web = 0;

% Initialize variables for minimum flange width calculation
min_flange_width = zeros(size(xAxis));

% Calculate minimum flange width required to maintain the desired safety factor
for i = 1:length(xAxis)
    fw = 0.002; % Start with the minimum possible flange width
    while fw <= 1 % Maximum flange width
        % Calculate the height of the web
        
        % Calculate areas of flange and web
        A_flange = fw * fh;
        A_web = ww * wh;
        
        % Calculate the neutral axis (y_NA)
        y_NA = ((fh + wh + 0.5 * fh) * A_flange + (fh + 0.5 * wh) * A_web + 0.5 * fh * A_flange) / (2 * A_flange + A_web);
        
        % Calculate distances from the neutral axis
        d_flange = abs(y_NA - (fh + wh + 0.5 * fh));
        d_web = abs(y_NA - (fh + 0.5 * wh));
        
        % Calculate the moment of inertia (I)
        I = 2 * ((1/12) * fw * fh^3 + A_flange * d_flange^2) + (1/12) * ww * wh^3 + A_web * d_web^2;
        
        % Calculate bending stress
        bending_stress = (abs(M_x(i)) * (bh / 2)) / I;
        
        % Check if bending stress is within allowable limits
        if bending_stress <= (allowable_stress)
            min_flange_width(i) = fw;
            total_volume_flange = total_volume_flange + A_flange * deltaX;
            total_volume_web = total_volume_web + A_web * deltaX;
            break;
        end
        
        % Increment flange width
        fw = fw + 0.0001;
    end
end

% Plot minimum flange width required along the spar
figure;
plot(xAxis, min_flange_width, 'LineWidth', 2);
title('Minimum Flange Width Along the Wing Spar');
xlabel('Position along the spar (m)');
ylabel('Flange Width (m)');
xlim([0, 4]); 

%Weight
total_volume = (2*total_volume_flange) + total_volume_web;
total_weight = total_volume * density_aluminum;
fprintf('Total weight of the spar with changing flange width(kg): %f\n', total_weight);

equidistant_points = linspace(0, L, 11);
equidistant_points = equidistant_points(2:end); % Ignore the first point to get exactly 10 points
sampled_widths = interp1(xAxis, min_flange_width, equidistant_points);

% Display the flange widths at the 10 points
disp('Flange widths at 10 equidistant points along the beam:');
for i = 1:length(sampled_widths)
    fprintf('At x = %.2f m, Flange width = %.4f m\n', equidistant_points(i), sampled_widths(i));
end

PloadsArray1 = 0;
PloadsArray2 = 0;
PloadsArray3 = 0;
PloadsArray4 = 0;
PloadsArray5 = 0;
PloadsArray6 = 0;
PloadsArray7 = 0;
PloadsArray8 = 0;
PloadsArray9 = 0;
PloadsArray10 = 0;


for i = 1:length(PloadsArray)
    if(xAxis(i) >= 0 && xAxis(i) <= 0.4)
        PloadsArray1 = PloadsArray1 + PloadsArray(i);
    elseif(xAxis(i) > 0.4 && xAxis(i) <= 0.8)
        PloadsArray2 = PloadsArray2 + PloadsArray(i);
    elseif(xAxis(i) > 0.8 && xAxis(i) <= 1.2)
        PloadsArray3 = PloadsArray3 + PloadsArray(i);
    elseif(xAxis(i) > 1.2 && xAxis(i) <= 1.6)
        PloadsArray4 = PloadsArray4 + PloadsArray(i);
    elseif(xAxis(i) > 1.6 && xAxis(i) <= 2)
        PloadsArray5 = PloadsArray5 + PloadsArray(i);
    elseif(xAxis(i) > 2 && xAxis(i) <= 2.4)
        PloadsArray6 = PloadsArray6 + PloadsArray(i);
    elseif(xAxis(i) > 2.4 && xAxis(i) <= 2.8)
        PloadsArray7 = PloadsArray7 + PloadsArray(i);
    elseif(xAxis(i) > 2.8 && xAxis(i) <= 3.2)
        PloadsArray8 = PloadsArray8 + PloadsArray(i);
    elseif(xAxis(i) > 3.2 && xAxis(i) <= 3.6)
        PloadsArray9 = PloadsArray9 + PloadsArray(i);
    elseif(xAxis(i) > 3.6 && xAxis(i) <= 4)
        PloadsArray10 = PloadsArray10 + PloadsArray(i);
    end
end

fprintf("load 0-0.4 (N):%f\n", PloadsArray1); 
fprintf("load 0.4-0.8 (N):%f\n", PloadsArray2); 
fprintf("load 0.8-1.2 (N):%f\n", PloadsArray3); 
fprintf("load 1.2-1.6 (N):%f\n", PloadsArray4); 
fprintf("load 1.6-2 (N):%f\n", PloadsArray5); 
fprintf("load 2-2.4 (N):%f\n", PloadsArray6); 
fprintf("load 2.4-2.8 (N):%f\n", PloadsArray7); 
fprintf("load 2.8-3.2 (N):%f\n", PloadsArray8); 
fprintf("load 3.2-3.6 (N):%f\n", PloadsArray9); 
fprintf("load 3.6-4 (N):%f\n", PloadsArray10); 


