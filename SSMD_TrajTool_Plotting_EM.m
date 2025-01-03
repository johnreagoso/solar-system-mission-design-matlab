%% Solar System Trajectory Plotting Tool
% This plotting utility provides the ability to show an interplanetary
% trajectory through the solar system and orbit profile of planets. 

%init: JReagoso Aug 2023

clear all ;
fclose all;
close all;
clc;

%% Initial Conditions(ICs): SUN MJ200 (Earth-Centered-Earth-Fixed)
global Julian_base_epoch;
global AU;

AU = 1.495978707e8;

% YearDep = 2030;     MonthDep = 8;      DayDep = 9;
% YearArr1 = 2031;    MonthArr1 = 8;    DayArr1 = 8;
% YearArr2 = 2032;    MonthArr2 = 1;     DayArr2 = 1;

YearDep = 2030;     MonthDep = 12;      DayDep = 11;
YearArr1 = 2031;    MonthArr1 = 7;    DayArr1 = 5;
YearArr2 = 2032;    MonthArr2 = 1;     DayArr2 = 1;

[DepartMJD,  Julian_cent_epoch_dep_initial] = Julian_datecalc(YearDep, MonthDep, DayDep);
[Arrive1MJD1, Julian_cent_epoch_arr1_initial] = Julian_datecalc(YearArr1, MonthArr1, DayArr1);

%% Execute
[DEPART_POSVECT, DEPART_VELVECT] = SolarSystemMission_tool_ExecSingleNoFlyby('Mars', 'Earth', DepartMJD, Arrive1MJD1);

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

%% Planet Orbit Plotting:
earth_day_incr = 0.0:0.250:365;
for kk = 1:numel(earth_day_incr)

    earthJulianCentEpochArrival = ((Julian_base_epoch + earth_day_incr(kk) - 2451545.0)/ 36525);
    [earthPosVector, ~, ~, ~] = planet_orbit_parameters('Earth',...
        earthJulianCentEpochArrival);

    earth_X(kk) = earthPosVector(1)/AU; 
    earth_Y(kk) = earthPosVector(2)/AU;
end

mars_day_incr = 0.0:0.250:782;
for kk = 1:numel(mars_day_incr)

    marsJulianCentEpochArrival = ((Julian_base_epoch + mars_day_incr(kk) - 2451545.0)/ 36525);
    [marsPosVector, ~, ~, ~] = planet_orbit_parameters('Mars', ...
        marsJulianCentEpochArrival);

    mars_X(kk) = marsPosVector(1)/AU; 
    mars_Y(kk) = marsPosVector(2)/AU;
end

venus_day_incr = 0.0:0.250:270;
for kk = 1:numel(venus_day_incr)

    venusJulianCentEpochArrival = ((Julian_base_epoch +venus_day_incr(kk) - 2451545.0)/ 36525);
    [venusPosVector, ~, ~, ~] = planet_orbit_parameters('Venus', ...
        venusJulianCentEpochArrival);

    venus_X(kk) = venusPosVector(1)/AU; 
    venus_Y(kk) = venusPosVector(2)/AU;
end

mercury_day_incr = 0.0:0.250:88;
for kk = 1:numel(mercury_day_incr)

    mercuryJulianCentEpochArrival = ((Julian_base_epoch +mercury_day_incr(kk) - 2451545.0)/ 36525);
    [mercuryPosVector, ~, ~, ~] = planet_orbit_parameters('Mercury', ...
        mercuryJulianCentEpochArrival);

    mercury_X(kk) = mercuryPosVector(1)/AU; 
    mercury_Y(kk) = mercuryPosVector(2)/AU;
end

hold on;
a1 = plot(earth_X, earth_Y,'--b', mars_X, mars_Y,'--r', venus_X, venus_Y,'--m','LineWidth', 0.01); 

hold on;
a1a = plot(mercury_X, mercury_Y,'--', 'color',[0.8500 0.3250 0.0980],'LineWidth', 0.01); 
xlim([-1.75,1.75]); ylim([-1.75,1.75]);
 
axis equal;
%% Planet Position Plotting
[earthPosVectorAtDep, ~, ~, ~]    = planet_orbit_parameters('Earth',    Julian_cent_epoch_dep_initial);
[MarsPosVectorAtDep, ~, ~, ~]     = planet_orbit_parameters('Mars',     Julian_cent_epoch_dep_initial);        
[VenusPosVectorAtDep, ~, ~, ~]    = planet_orbit_parameters('Venus',    Julian_cent_epoch_dep_initial);   
[MercuryPosVectorAtDep, ~, ~, ~]  = planet_orbit_parameters('Mercury',  Julian_cent_epoch_dep_initial);   

[EarthPosVectorAtArr1, ~, ~, ~]   = planet_orbit_parameters('Earth',    Julian_cent_epoch_arr1_initial);
[MarsPosVectorAtArr1, ~, ~, ~]    = planet_orbit_parameters('Mars',     Julian_cent_epoch_arr1_initial);
[VenusPosVectorAtArr1, ~, ~, ~]   = planet_orbit_parameters('Venus',    Julian_cent_epoch_arr1_initial);
[MercuryPosVectorAtArr1, ~, ~, ~] = planet_orbit_parameters('Mercury',  Julian_cent_epoch_arr1_initial);

%% Depart:
% Earth:
p1dep = plot(earthPosVectorAtDep(1)/AU,  earthPosVectorAtDep(2)/AU, 'ob', 'MarkerFaceColor', 'blue', 'MarkerSize', 5.0, 'MarkerEdgeColor','k');
text(earthPosVectorAtDep(1)/AU,  earthPosVectorAtDep(2)/AU, 'Earth at ERO Depart', 'FontSize', 9, 'VerticalAlignment','top','HorizontalAlignment','left');

% Mars:
p2dep = plot(MarsPosVectorAtDep(1)/AU,  MarsPosVectorAtDep(2)/AU, 'ob', 'MarkerFaceColor', 'red', 'MarkerSize', 5.0);
text(MarsPosVectorAtDep(1)/AU,  MarsPosVectorAtDep(2)/AU, 'Mars ERO Depart: 11-Dec-30', 'FontSize', 9, 'VerticalAlignment','bottom','HorizontalAlignment','left');

%% Arrival:
% Earth:
p2arr = plot(EarthPosVectorAtArr1(1)/AU,  EarthPosVectorAtArr1(2)/AU, 'ob', 'MarkerFaceColor', 'blue', 'MarkerSize',  5.0);
text(EarthPosVectorAtArr1(1)/AU,  EarthPosVectorAtArr1(2)/AU, 'ERO Earth Arrival: 5-Jul-31', 'FontSize', 9, 'VerticalAlignment','top','HorizontalAlignment','left');

% Mars:
p3arr = plot(MarsPosVectorAtArr1(1)/AU, MarsPosVectorAtArr1(2)/AU, 'ob', 'MarkerFaceColor', 'red', 'MarkerSize',5.0);
text(MarsPosVectorAtArr1(1)/AU, MarsPosVectorAtArr1(2)/AU, 'Mars at ERO Arrival', 'FontSize', 9, 'VerticalAlignment','top','HorizontalAlignment','left');

% Venus:
p2arr = plot(VenusPosVectorAtArr1(1)/AU,  VenusPosVectorAtArr1(2)/AU, 'or', 'MarkerFaceColor', 'magenta', 'MarkerSize', 5.0, 'MarkerEdgeColor','k');
text(VenusPosVectorAtArr1(1)/AU,  VenusPosVectorAtArr1(2)/AU, 'Venus at Earth Arrival', 'FontSize', 9, 'VerticalAlignment','top','HorizontalAlignment','right');

% Mercury:
p2arr = plot(MercuryPosVectorAtArr1(1)/AU,  MercuryPosVectorAtArr1(2)/AU, 'or', 'MarkerFaceColor', [0.8500 0.3250 0.0980], 'MarkerSize',  5.0, 'MarkerEdgeColor','k');
text(MercuryPosVectorAtArr1(1)/AU,  MercuryPosVectorAtArr1(2)/AU, 'Mercury at Earth Arrival', 'FontSize', 9, 'VerticalAlignment','top','HorizontalAlignment','right');

%% Trajectory Plotting:
traj1 = plot(X1/AU, Y1/AU,'-g', 'LineWidth', 1.0);
hold off;

xlim([-2.0, 2.0]); 
ylim([-2.0, 2.0]);
hold on;
plot(0,0,'+k');
xlabel('AU'); ylabel('AU');
%title('Mars-Earth ERO Ballistic Transfer (TEI: 09 Aug 2030)', 'FontSize', 14);
title('Mars-Earth ERO Ballistic Transfer (TEI: 11 Dec 2030)', 'FontSize', 14);

set(0,'defaultAxesFontSize',14);
ax = gca;
ax.FontSize = 12;
ax.FontWeight = 'bold';
disp('plotting complete');

