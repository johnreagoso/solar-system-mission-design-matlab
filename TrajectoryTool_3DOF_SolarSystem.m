%% 3-DOF Trajectory Tool
clear all ;
fclose all;
close all;
clc;

addpath('C:\Users\jreagoso\Desktop\DeepSpace_AstroModeling');

%% Initial Conditions(ICs): SUN MJ200 (Earth-Centered-Earth-Fixed)
global Julian_base_epoch;
global AU;

AU = 1.495978707e8;

YearDep = 2030;     MonthDep = 8;      DayDep = 9;
YearArr1 = 2031;    MonthArr1 = 8;    DayArr1 = 8;
YearArr2 = 2032;    MonthArr2 = 1;     DayArr2 = 1;

[DepartMJD,  Julian_cent_epoch_dep_initial] = Julian_datecalc(YearDep, MonthDep, DayDep);
[Arrive1MJD1, Julian_cent_epoch_arr1_initial] = Julian_datecalc(YearArr1, MonthArr1, DayArr1);
[Arrive1MJD2, Julian_cent_epoch_arr2_initial] = Julian_datecalc(YearArr2, MonthArr2, DayArr2);

%% Execute
[DEPART_POSVECT, DEPART_VELVECT, POSTFLYBY_POSVECT, POSTFLYBY_VELVECT] = ...
    SolarSystemMission_tool_ExecSingle('Mars', 'Earth', 'Mars', DepartMJD, Arrive1MJD1, Arrive1MJD2);

%% Execute Runge-Kutta for trajectory propagation for 1st leg:
Julian_base_epoch = DepartMJD;

tol     = 1.00e-7;       
TOF1    = Arrive1MJD1-DepartMJD;
maxTOF1 = TOF1*24*60*60;   % days -> sec
to1     = 0.00:60*60*24:maxTOF1; % sec

x01 = [DEPART_POSVECT(1)+7000; DEPART_POSVECT(2)+7000; DEPART_POSVECT(3);...
       DEPART_VELVECT(1); DEPART_VELVECT(2); DEPART_VELVECT(3)];
   
options = odeset('RelTol',tol);    
[time1, state_vector1] = ode23(@Function_3DOF_SolarSystem, to1, x01, options);

X1 = state_vector1(:,1);     Y1 = state_vector1(:,2);

%% Execute Runge-Kutta for trajectory propagation for 2nd leg:
Julian_base_epoch = Arrive1MJD1;

tol     = 1.00e-7;         
TOF2    = Arrive1MJD2-Arrive1MJD1;
maxTOF2 = TOF2*24*60*60;   % days -> sec
to2     = 0.00:60*60*24:maxTOF2; % sec

x02  = [POSTFLYBY_POSVECT(1)+3500; POSTFLYBY_POSVECT(2)+3500; POSTFLYBY_POSVECT(3);...
       POSTFLYBY_VELVECT(1); POSTFLYBY_VELVECT(2); POSTFLYBY_VELVECT(3)];
   
options = odeset('RelTol',tol);    
[time2, state_vector2] = ode23(@Function_3DOF_SolarSystem, to2, x02, options);

X2 = state_vector2(:,1);     Y2 = state_vector2(:,2);

%% Planet Orbit Plotting:
earth_day_incr = 0.0:0.250:365;
for kk = 1:numel(earth_day_incr)

    JulianCentEpochArrival = ((Julian_base_epoch + earth_day_incr(kk) - 2451545.0)/ 36525);
    [earthPosVector, ~, ~, ~] = planet_orbit_parameters('Earth',...
        JulianCentEpochArrival);

    earth_X(kk) = earthPosVector(1)/AU; 
    earth_Y(kk) = earthPosVector(2)/AU;
end

mars_day_incr = 0.0:0.250:782;
for kk = 1:numel(mars_day_incr)

    JulianCentEpochArrival = ((Julian_base_epoch + mars_day_incr(kk) - 2451545.0)/ 36525);
    [marsPosVector, ~, ~, ~] = planet_orbit_parameters('Mars', ...
        JulianCentEpochArrival);

    mars_X(kk) = marsPosVector(1)/AU; 
    mars_Y(kk) = marsPosVector(2)/AU;
end

venus_day_incr = 0.0:0.250:270;
for kk = 1:numel(venus_day_incr)

    JulianCentEpochArrival = ((Julian_base_epoch +venus_day_incr(kk) - 2451545.0)/ 36525);
    [venusPosVector, ~, ~, ~] = planet_orbit_parameters('Venus', ...
        JulianCentEpochArrival);

    venus_X(kk) = venusPosVector(1)/AU; 
    venus_Y(kk) = venusPosVector(2)/AU;
end

figure;
a1 = plot(earth_X, earth_Y,'--b', ...
     mars_X, mars_Y,'--r', venus_X, venus_Y,'--m'   ,'LineWidth', 0.01); 

axis equal;
%% Planet Position Plotting
   
[earthPosVectorAtDep, ~, ~, ~] = planet_orbit_parameters('Earth',  Julian_cent_epoch_dep_initial);
[MarsPosVectorAtDep, ~, ~, ~] = planet_orbit_parameters('Mars',   Julian_cent_epoch_dep_initial);        
[VenusPosVectorAtDep, ~, ~, ~] = planet_orbit_parameters('Venus',   Julian_cent_epoch_dep_initial);     

[MarsPosVectorAtArr1, ~, ~, ~] = planet_orbit_parameters('Mars',   Julian_cent_epoch_arr1_initial);
[EarthPosVectorAtArr1, ~, ~, ~] = planet_orbit_parameters('Earth',   Julian_cent_epoch_arr1_initial);
[VenusPosVectorAtArr1, ~, ~, ~] = planet_orbit_parameters('Venus',   Julian_cent_epoch_arr1_initial);

[earthPosVectorAtArr2, ~, ~, ~] = planet_orbit_parameters('Earth', Julian_cent_epoch_arr2_initial);

hold on;
% Depart:
p1 = plot(earthPosVectorAtDep(1)/AU,  earthPosVectorAtDep(2)/AU, 'ob', 'MarkerFaceColor', 'blue', 'MarkerSize', 4.0);
text(earthPosVectorAtDep(1)/AU,  earthPosVectorAtDep(2)/AU, 'Earth Depart:  23 Feb 2023', 'FontSize', 9, 'VerticalAlignment','bottom','HorizontalAlignment','left');

%p1a = plot(VenusPosVectorAtDep(1)/AU,  VenusPosVectorAtDep(2)/AU, 'or', 'MarkerFaceColor', 'red', 'MarkerSize', 4.0);
%text(MarsPosVectorAtDep(1)/AU,  MarsPosVectorAtDep(2)/AU, 'Mars at Earth Depart', 'FontSize', 9, 'VerticalAlignment','bottom','HorizontalAlignment','left');
%text(VenusPosVectorAtDep(1)/AU,  VenusPosVectorAtDep(2)/AU, 'Venus at Earth Depart', 'FontSize', 9, 'VerticalAlignment','bottom','HorizontalAlignment','left');

%Flyby:
p2 = plot(VenusPosVectorAtArr1(1)/AU,  VenusPosVectorAtArr1(2)/AU, 'or', 'MarkerFaceColor', 'red', 'MarkerSize',  4.0);
%text(MarsPosVectorAtArr1(1)/AU,  MarsPosVectorAtArr1(2)/AU, 'Mars Flyby', 'FontSize', 9, 'VerticalAlignment','top','HorizontalAlignment','left');
text(VenusPosVectorAtArr1(1)/AU,  VenusPosVectorAtArr1(2)/AU, 'Venus Flyby:  21 Sep 2023 ', 'FontSize', 9, 'VerticalAlignment','top','HorizontalAlignment','left');

p2a = plot(EarthPosVectorAtArr1(1)/AU,  EarthPosVectorAtArr1(2)/AU, 'ob', 'MarkerFaceColor', 'blue', 'MarkerSize',  4.0);
text(EarthPosVectorAtArr1(1)/AU,  EarthPosVectorAtArr1(2)/AU, 'Earth at Flyby', 'FontSize', 9, 'VerticalAlignment','top','HorizontalAlignment','left');

% Arrive2:
p3 = plot(earthPosVectorAtArr2(1)/AU, earthPosVectorAtArr2(2)/AU, 'ob', 'MarkerFaceColor', 'blue', 'MarkerSize',4.0);
text(earthPosVectorAtArr2(1)/AU, earthPosVectorAtArr2(2)/AU, 'Earth Return:  12 Jun 2024', 'FontSize', 9, 'VerticalAlignment','top','HorizontalAlignment','left');

%% Trajectory Plotting:
a1 = plot(X1/AU, Y1/AU,'-g', X2/AU, Y2/AU,'-g');

xlim([-2.5, 2.5]); 
ylim([-2.50, 2.50]);
hold on;
plot(0,0,'+k');
xlabel('AU'); ylabel('AU');
title('Earth-Venus-Earth Free Return Trajectory (Launch: 23 Feb 2023)', 'FontSize', 16);
set(0,'defaultAxesFontSize',14);
ax = gca;
ax.FontSize = 12;
ax.FontWeight = 'bold';


