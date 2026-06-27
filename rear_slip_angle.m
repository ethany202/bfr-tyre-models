% rear_slip_angle.m
% =========================================================================
% Rear-axle slip angle for the Hoosier 43075 16x7.5-10 R20, computed by
% inverting the fitted Pacejka lateral model against the lateral force the
% rear axle must generate in steady-state cornering.
%
% WHY THIS METHOD: the MoTeC IMU lateral-accel / yaw-rate channels on
% Driving Day 06 logged all zeros (sensor not plugged in), so slip angle
% cannot be read directly. Instead we use the operating point (radius +
% speed -> lateral acceleration), split the lateral force onto the rear
% tyres with the team's Milliken load-transfer numbers, and solve the tyre
% model for the slip angle that produces that force.
%
%   alpha_rear  such that  Fy_outer(alpha) + Fy_inner(alpha) = m_rear * Ay
%
% Inputs come from camber_gain_optimization.m (vehicle) and
% hoosier_r20_tire_params.mat (tyre).
% =========================================================================
clear; clc;

%% ---- tyre model ----
S = load('hoosier_r20_tire_params.mat');   % gives p_fit / tireParams.coeffs
p = S.p_fit(:)';

%% ---- vehicle (from "2026 SUS Calculations.xlsx" Milliken sheet) ----
statF = 674.9;   statR = 688.7;            % static load per tyre [N]
trF   = 462.9/1.2; trR = 496.1/1.2;        % lateral load transfer per tyre per g [N/g]
camber_beneficial = 1.5;                   % |rear static camber| used as beneficial IA [deg]
Wf = 2*statF;  Wr = 2*statR;               % axle static loads [N]
L  = 1.535;                                % wheelbase [m]

%% ---- operating point ----
% Requested condition: R = 5.3-5.5 m at ~28 km/h.
R  = 5.4;                                   % corner radius [m]
V  = 28/3.6;                                % speed [m/s]
Ay = V^2 / (R*9.81);                        % lateral acceleration [g]
fprintf('Operating point: R = %.2f m, V = %.1f km/h  ->  Ay = %.2f g\n\n', R, V*3.6, Ay);

aR = axleSlip(Ay, statR, trR, Wr, camber_beneficial, p);
aF = axleSlip(Ay, statF, trF, Wf, camber_beneficial, p);
fprintf('  alpha_REAR  = %.2f deg\n', aR);
fprintf('  alpha_FRONT = %.2f deg   (balance: %+.2f deg, %s)\n\n', ...
        aF, aF-aR, ternary(aF>aR,'understeer','oversteer'));

%% ---- slip angle vs lateral acceleration (reusable table) ----
fprintf('%6s %9s %9s %9s\n','Ay[g]','V@5.4m','alpha_R','alpha_F');
for Ay = [0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.12]
    Vk = sqrt(Ay*9.81*5.4)*3.6;
    aR = axleSlip(Ay,statR,trR,Wr,camber_beneficial,p);
    aF = axleSlip(Ay,statF,trF,Wf,camber_beneficial,p);
    fprintf('%6.2f %8.1f  %8.2f %9.2f\n',Ay,Vk,aR,aF);
end

% ===================== local functions ===================================
function a = axleSlip(Ay, stat, tr, W, gamma, p)
    Fzo = stat + tr*Ay;                 % outer tyre load [N]
    Fzi = max(stat - tr*Ay, 1);         % inner tyre load [N]
    Freq = W * Ay;                      % required axle lateral force [N]
    f = @(al) Fcorner(al,Fzo,gamma,p) + Fcorner(al,Fzi,gamma,p) - Freq;
    a = fzero(f, 2.0);
end

function F = Fcorner(alpha, Fz, gamma, p)
    % cornering-force magnitude at slip-angle magnitude alpha (+gamma = beneficial)
    F = -mf_lateral(p, alpha, Fz, gamma);
end

function FY = mf_lateral(p, alpha, Fz, gamma)
    a0=p(1);a1=p(2);a2=p(3);a3=p(4);a4=p(5);a5=p(6);a6=p(7);a7=p(8);
    a8=p(9);a9=p(10);a10=p(11);a11=p(12);a12=p(13);a13=p(14);
    C=a0; D=a1.*Fz.^2+a2.*Fz;
    BCD=a3.*sin(2.*atan(Fz./a4)).*(1-a5.*abs(gamma));
    B=BCD./(C.*D); E=a6.*Fz+a7;
    SH=a8.*gamma+a9.*Fz+a10; SV=a11.*Fz.*gamma+a12.*Fz+a13;
    phi=alpha+SH;
    FY=D.*sin(C.*atan(B.*phi-E.*(B.*phi-atan(B.*phi))))+SV;
end

function out = ternary(c,a,b); if c, out=a; else, out=b; end; end
