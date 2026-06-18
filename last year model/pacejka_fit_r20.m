clear all
close all
%% extract and inspect data
% data access: load the Hoosier 43075 16x7.5-10 R20 (8 inch rim) cornering run
% directly instead of using the interactive file picker.
% resolve the data path relative to THIS script's location so it works
% no matter what MATLAB's current folder is (script is in Tire data/temp).
thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir); thisDir = pwd; end   % fallback when running a selection
datafile = fullfile(thisDir, '..', 'Data', 'RunData_Cornering_Matlab_SI_Round9', 'B2356run7.mat');  % R20, 8 inch rim
t = importdata(datafile);  % import data
s = importdata(datafile);  % import data
% extract fields
ET = extractfield(t, 'ET');  % elapsed time
FZ = extractfield(t, 'FZ');  % normal load
IA = extractfield(t, 'IA');  % camber angle
SA = extractfield(t, 'SA');  % slip angle
SR = extractfield(t, 'SR');  % slip ratio
FY = extractfield(t, 'FY');  % lateral force
FX = extractfield(t, 'FX');  % longitudinal force
MZ = extractfield(t, 'MZ');  % aligning moment
ETL = extractfield(s, 'ET');  % elapsed time
FZL = extractfield(s, 'FZ');  % normal load
IAL = extractfield(s, 'IA');  % camber angle
SAL = extractfield(s, 'SA');  % slip angle
SRL = extractfield(s, 'SR');  % slip ratio
FYL = extractfield(s, 'FY');  % lateral force
FXL = extractfield(s, 'FX');  % longitudinal force
MZL = extractfield(s, 'MZ');  % aligning moment
% create function to plot inclination angle, slip angle, and vertical load vs. time
tiledlayout(3,1);
nexttile
plot(ET, FZ);
title('Normal Load (N) vs. Elapsed Time (s)')
nexttile
plot(ET, IA);
title('Inclination Angle (deg) vs. Elapsed Time (s)')
nexttile
plot(ET, SA);
title('Slip Angle (deg) vs. Elapsed Time (s)')
%nexttile
%plot(SA, MZ);
%title('Slip Angle (deg) vs. Aligning Moment (Nm)')
% R20s
ET_range = ET > 320 & ET < 445  % IA = 2 deg
% LCOs
ET_rangeL = ETL > 490 & ETL < 620  %  IA = 2 deg
% original
% ET_range = ET > 290 & ET < 450;  % get range of ET values that sweeps SA and FZ values w/ IA =
% get relevant SA, FZ, and FY values
SA_range = SA(ET_range); FZ_range = FZ(ET_range); FY_range = FY(ET_range); IA_range = IA(ET_range); FX_range = FX(ET_range); SR_range = SR(ET_range); MZ_range = MZ(ET_range);
SA_rangeL = SA(ET_rangeL); FZ_rangeL = FZ(ET_rangeL); FY_rangeL = FY(ET_rangeL); IA_rangeL = IA(ET_rangeL); FX_rangeL = FX(ET_rangeL); SR_rangeL = SR(ET_rangeL); MZ_rangeL = MZ(ET_rangeL);
%% fitting pacejka function to the data
x_data = [SA_range; FZ_range];  % [2 x 5578] matrix
y_data = FY_range;  % [1 x 5578] array
x_dataL = [SA_rangeL; FZ_rangeL];
y_dataL = FY_rangeL;
x0 = [1 1 1 1];  % starting point for coefficient guesses
% define pacejka function
% P1 = B
% P2 = C
% P3 = D1
% P4 = D2
% create inputs P (pacejka coefficients) and x_data (slip angle and normal load data) to define function
% input function into lsqcurvefit()
% find P
Pacejka4_Model = @(P, x_data) ((P(1) + P(2)/1000*x_data(2,:)).*x_data(2,:)).*sin(P(4)*atan(P(3).*x_data(1,:)));
[P, f_resnorm, f_residual] = lsqcurvefit(Pacejka4_Model, x0, x_data, y_data, [], []);
Pacejka4_ModelL = @(PL, x_data) ((PL(1) + PL(2)/1000*x_data(2,:)).*x_data(2,:)).*sin(PL(4)*atan(PL(3).*x_data(1,:)));
[PL, f_resnorm, f_residual] = lsqcurvefit(Pacejka4_ModelL, x0, x_dataL, y_dataL, [], []);
%% creating a surface based on our pacejka coefficients/model
% define a linear space for SA
SA_space = linspace(-12, 12, 100);
% define a linear space for FZ
FZ_space = linspace(-1200, -200, 100);
% determine values of FY w/ swept SA and FZ inputs
FY_mf = zeros(length(SA_space), length(FZ_space));
for i=1:length(SA_space)  % sweep SA
    for j=1:length(FZ_space)  % sweep FZ
        FY_mf(j,i) = Pacejka4_Model(P, [SA_space(i); FZ_space(j)]);
    end
end
CF_mf = FY_mf.*cosd(SA_space);  % pacejka
FY_mfL = zeros(length(SA_space), length(FZ_space));
for i=1:length(SA_space)  % sweep SA
    for j=1:length(FZ_space)  % sweep FZ
        FY_mfL(j,i) = Pacejka4_ModelL(PL, [SA_space(i); FZ_space(j)]);
    end
end
% space = pacejka
% data = range
figure
% r20 = surf(SA_space, FZ_space, FY_mf);
% set(r20, 'FaceColor', [1,0,0])
% hold on
surf(SA_space, FZ_space, FY_mfL);
% %% plot raw data agains pacejka surface
% % raw data plot
% figure
% title('Pacejka Tyre Model for TTC Round 8 Run 15 LCO Tyres (IA = 2 deg)')
% xlabel('Slip Angle (deg)')
% ylabel('Normal Load (N)')
% zlabel('Lateral Force (N)')
% %plot3(SA_range, FZ_range, FY_range.*cosd(SA_range), '-')  % data
% plot3(SA_range, FZ_range, FY_range, '-')
% hold on
% %surf(SA_space, FZ_space, CF_mf)  % pacejka
% surf(SA_space, FZ_space, FY_mf)
% hold on
%
% % plot LCOs
% plot3(SA_rangeL, FZ_rangeL, FY_rangeL, '-')
% hold on
% surf(SA_spaceL, FZ_spaceL, FY_mfL)
% figure
% plot(FZ_space, FY_mf,'-');
% title('R20 Normal Load vs. Lateral Force')
%view(180,0);
%colormap jet
%figure
% FY_diff is a matrix of [100 x 100], rows the difference with changing normal force
% columns the difference with changing slip angle
% FY_diff = abs(FY_mfL - FY_mf);
% surf(SA_space, FZ_space, FY_diff);
% colorbar
% plot FY_mfL
% get rid of color
% new color according to the values stored in the matrix FY_diff
xlabel('Slip Angle (deg)')
ylabel('Normal Load (N)')
zlabel('Lateral Force (N)')
%coeff_L = FY_mfL./FZ_rangeL;
%coeff = FY_mf./FZ_range;
%figure
%plot(FZ_range, coeff);
%plot(FZ_rangeL, coeffL);
%figure
%plot3(SA_range, FZ_range, FX_range - (FY_range.*sind(SA_range)));
%plot3(SA_range, FZ_range, FX_range - abs(FY_range.*sind(SA_range)));
%title('Slip Angle and Vertical Load vs. Longitudinal Force')
%xlabel('Slip Angle (deg)')
%ylabel('Normal Load (N)')
%zlabel('Longitudinal Force (N)')
%figure
%plot3(SA_range, FZ_range, MZ_range);
%title('Slip Angle (deg) vs. Aligning Moment (Nm)');
%xlabel('Slip Angle (deg)');
%ylabel('Normal Load (N)');
%zlabel('Aligning Moment (Nm)');
