%% 3-DOF Trajectory Tool
clear all ;
fclose all;
close all;
clc;

global Julian_base_epoch;
global AU;

AU = 1.495978707e8;

YearDep = 2030;     MonthDep = 1;      DayDep = 1;
YearArr1 = 2031;    MonthArr1 = 4;     DayArr1 = 6;
YearArr2 = 2039;    MonthArr2 = 10;     DayArr2 = 16;

[DepartMJD,  Julian_cent_epoch_dep_initial] = Julian_datecalc(YearDep, MonthDep, DayDep);
[Arrive1MJD1, Julian_cent_epoch_arr1_initial] = Julian_datecalc(YearArr1, MonthArr1, DayArr1);
[Arrive1MJD2, Julian_cent_epoch_arr2_initial] = Julian_datecalc(YearArr2, MonthArr2, DayArr2);

%% Execute
[DEPART_POSVECT, DEPART_VELVECT, POSTFLYBY_POSVECT, POSTFLYBY_VELVECT] = ...
    SolarSystemMission_tool_ExecSingle('Earth', 'Jupiter', 'Pluto', DepartMJD, Arrive1MJD1, Arrive1MJD2);

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
earth_day_incr = 0.0:0.25:365;
for kk = 1:numel(earth_day_incr)

    JulianCentEpochArrival = ((Julian_base_epoch + earth_day_incr(kk) - 2451545.0)/ 36525);
    [earthPosVector, ~, ~, ~] = planet_orbit_parameters('Earth',...
        JulianCentEpochArrival);

    earth_X(kk) = earthPosVector(1)/AU; 
    earth_Y(kk) = earthPosVector(2)/AU;
end
clearvars JulianCentEpochArrival;

jupiter_day_incr = 0.0:0.25:4333;
for kk = 1:numel(jupiter_day_incr)

    JulianCentEpochArrival = ((Julian_base_epoch + jupiter_day_incr(kk) - 2451545.0)/ 36525);
    [JupiterPosVector, ~, ~, ~] = planet_orbit_parameters('Jupiter', ...
        JulianCentEpochArrival);

    Jupiter_X(kk) = JupiterPosVector(1)/AU; 
    Jupiter_Y(kk) = JupiterPosVector(2)/AU;
end
clearvars JulianCentEpochArrival;

saturn_day_incr = 0.0:0.25:10756;
for kk = 1:numel(saturn_day_incr)

    JulianCentEpochArrival = ((Julian_base_epoch + saturn_day_incr(kk) - 2451545.0)/ 36525);
    [SaturnPosVector, ~, ~, ~] = planet_orbit_parameters('Saturn', ...
        JulianCentEpochArrival);

    Saturn_X(kk) = SaturnPosVector(1)/AU; 
    Saturn_Y(kk) = SaturnPosVector(2)/AU;
end

clearvars JulianCentEpochArrival;

uranus_day_incr = 0.0:0.25:30687;
for kk = 1:numel(uranus_day_incr)

    JulianCentEpochArrival = ((Julian_base_epoch + uranus_day_incr(kk) - 2451545.0)/ 36525);
    [UranusPosVector, ~, ~, ~] = planet_orbit_parameters('Uranus', ...
        JulianCentEpochArrival);

    Uranus_X(kk) = UranusPosVector(1)/AU; 
    Uranus_Y(kk) = UranusPosVector(2)/AU;
end
clearvars JulianCentEpochArrival;

neptune_day_incr = 0.0:0.25:60182;
for kk = 1:numel(neptune_day_incr)

    JulianCentEpochArrival = ((Julian_base_epoch + neptune_day_incr(kk) - 2451545.0)/ 36525);
    [NeptunePosVector, ~, ~, ~] = planet_orbit_parameters('Neptune', ...
        JulianCentEpochArrival);

    Neptune_X(kk) = NeptunePosVector(1)/AU; 
    Neptune_Y(kk) = NeptunePosVector(2)/AU;
end
clearvars JulianCentEpochArrival;

pluto_day_incr = 0.0:0.25:90560;
for kk = 1:numel(pluto_day_incr)

    JulianCentEpochArrival = ((Julian_base_epoch +pluto_day_incr(kk) - 2451545.0)/ 36525);
    [PlutoPosVector, ~, ~, ~] = planet_orbit_parameters('Pluto', ...
        JulianCentEpochArrival);

    Pluto_X(kk) = PlutoPosVector(1)/AU; 
    Pluto_Y(kk) = PlutoPosVector(2)/AU;
end

figure;
a1 = plot(earth_X, earth_Y,'--b');
hold on;
a2 = plot(Jupiter_X, Jupiter_Y, '--', 'Color', [0.8500 0.3250 0.0980]);
hold on;
a3 = plot(Saturn_X, Saturn_Y, '--', 'Color', [0.6 0.4 0.2]);
hold on;
a4 = plot(Uranus_X, Uranus_Y, '--', 'Color', [0.6350 0.0780 0.1840]);
hold on;
a5 = plot(Neptune_X, Neptune_Y,'--', 'Color', [0 0.4470 0.7410]); %, Pluto_X, Pluto_Y,'--c','LineWidth', 0.01); 
hold on;
a6 = plot(Pluto_X, Pluto_Y,'--', 'Color', [0 0 0.6]); %, Pluto_X, Pluto_Y,'--c','LineWidth', 0.01); 


axis equal;

%% Planet Position Plotting   
[earthPosVectorAtDep, ~, ~, ~] = planet_orbit_parameters('Earth',  Julian_cent_epoch_dep_initial);
%[JupiterPosVectorAtDep, ~, ~, ~] = planet_orbit_parameters('Jupiter',   Julian_cent_epoch_dep_initial);        
%[PlutoPosVectorAtDep, ~, ~, ~] = planet_orbit_parameters('Pluto',   Julian_cent_epoch_dep_initial);     

%% Flyby Locations:
[JupiterPosVectorAtArr1, ~, ~, ~] = planet_orbit_parameters('Jupiter',   Julian_cent_epoch_arr1_initial);
[SaturnPosVectorAtArr1, ~, ~, ~] = planet_orbit_parameters('Saturn',   Julian_cent_epoch_arr1_initial);
%[EarthPosVectorAtArr1, ~, ~, ~] = planet_orbit_parameters('Earth',   Julian_cent_epoch_arr1_initial);
[PlutoPosVectorAtArr1, ~, ~, ~] = planet_orbit_parameters('Pluto',   Julian_cent_epoch_arr1_initial);

%% Arrival locations:
%[EarthPosVectorAtArr1, ~, ~, ~] = planet_orbit_parameters('Earth',   Julian_cent_epoch_arr1_initial);
[JupiterPosVectorAtArr2, ~, ~, ~] = planet_orbit_parameters('Jupiter', Julian_cent_epoch_arr2_initial);
[SaturnPosVectorAtArr2, ~, ~, ~]  = planet_orbit_parameters('Saturn',  Julian_cent_epoch_arr2_initial);
[UranusPosVectorAtArr2, ~, ~, ~]  = planet_orbit_parameters('Uranus',  Julian_cent_epoch_arr2_initial);
[NeptunePosVectorAtArr2, ~, ~, ~] = planet_orbit_parameters('Neptune', Julian_cent_epoch_arr2_initial);
[PlutoPosVectorAtArr2, ~, ~, ~]   = planet_orbit_parameters('Pluto',   Julian_cent_epoch_arr2_initial);

hold on;
% Depart:
p1 = plot(earthPosVectorAtDep(1)/AU,  earthPosVectorAtDep(2)/AU, 'ob', 'MarkerFaceColor', 'blue', 'MarkerEdgeColor','k','MarkerSize', 5.0);
text(earthPosVectorAtDep(1)/AU,  earthPosVectorAtDep(2)/AU, 'Earth Depart: 1-Jan-2030', 'FontSize', 9, 'VerticalAlignment','bottom','HorizontalAlignment','left');

% Jupiter Flyby:0
p2 = plot(JupiterPosVectorAtArr1(1)/AU,  JupiterPosVectorAtArr1(2)/AU, 'or', 'MarkerFaceColor', [0.8500 0.3250 0.0980],'MarkerEdgeColor', 'k','MarkerSize',  5.0);
text(JupiterPosVectorAtArr1(1)/AU,  JupiterPosVectorAtArr1(2)/AU, 'Jupiter Gravity Assist Flyby:  6-Apr-2031 ', 'FontSize', 9, 'VerticalAlignment','top','HorizontalAlignment','left');

p2b = plot(SaturnPosVectorAtArr1(1)/AU,  SaturnPosVectorAtArr1(2)/AU, 'ob', 'MarkerFaceColor', [0.6 0.4 0.2], 'MarkerSize',  5.0);
text(SaturnPosVectorAtArr1(1)/AU,  SaturnPosVectorAtArr1(2)/AU, 'Saturn at Jupiter Flyby', 'FontSize', 9, 'VerticalAlignment','top','HorizontalAlignment','right');

p2 = plot(PlutoPosVectorAtArr1(1)/AU,  PlutoPosVectorAtArr1(2)/AU, 'or', 'MarkerFaceColor', [0 0 0.6], 'MarkerEdgeColor', 'k','MarkerSize',  5.0);
text(PlutoPosVectorAtArr1(1)/AU,  PlutoPosVectorAtArr1(2)/AU, 'Pluto at Jupiter Flyby', 'FontSize', 9, 'VerticalAlignment','top','HorizontalAlignment','left');

%% Pluto Arrival:

%p2a = plot(JupiterPosVectorAtArr2(1)/AU,  JupiterPosVectorAtArr2(2)/AU, 'ob', 'MarkerFaceColor', 'blue', 'MarkerSize',  4.0);
%text(JupiterPosVectorAtArr2(1)/AU,  JupiterPosVectorAtArr1(2)/AU, 'Jupiter at Flyby', 'FontSize', 9, 'VerticalAlignment','top','HorizontalAlignment','left');

p2c = plot(UranusPosVectorAtArr2(1)/AU,  UranusPosVectorAtArr2(2)/AU, 'ob', 'MarkerFaceColor', [0.6350 0.0780 0.1840], 'MarkerSize',  5.0);
text(UranusPosVectorAtArr2(1)/AU,  UranusPosVectorAtArr2(2)/AU, 'Uranus at Pluto Arrival', 'FontSize', 9, 'VerticalAlignment','top','HorizontalAlignment','left');

p2d = plot(NeptunePosVectorAtArr2(1)/AU,  NeptunePosVectorAtArr2(2)/AU, 'ob', 'MarkerEdgeColor','k','MarkerFaceColor', [0 0.4470 0.7410], 'MarkerSize',  5.0);
text(NeptunePosVectorAtArr2(1)/AU,  NeptunePosVectorAtArr2(2)/AU, 'Neptune at Pluto Arrival', 'FontSize', 9, 'VerticalAlignment','top','HorizontalAlignment','left');


p2e = plot(PlutoPosVectorAtArr2(1)/AU,  PlutoPosVectorAtArr2(2)/AU, 'ob', 'MarkerFaceColor', [0 0 0.6], 'MarkerEdgeColor', 'k','MarkerSize',  5.0);
text(PlutoPosVectorAtArr2(1)/AU,  PlutoPosVectorAtArr2(2)/AU, 'Pluto Arrival 16-Oct-2039', 'FontSize', 9, 'VerticalAlignment','top','HorizontalAlignment','left');

% Arrive2:
% p3 = plot(earthPosVectorAtArr2(1)/AU, earthPosVectorAtArr2(2)/AU, 'ob', 'MarkerFaceColor', 'blue', 'MarkerSize',4.0);
% text(earthPosVectorAtArr2(1)/AU, earthPosVectorAtArr2(2)/AU, 'Earth Return:  12 Jun 2024', 'FontSize', 9, 'VerticalAlignment','top','HorizontalAlignment','left');

%% Trajectory Plotting:
hold on;
a1 = plot(X1/AU, Y1/AU,'-g', X2/AU, Y2/AU,'Color', [0.2 0.4 0.4]);

xlim([-40, 52]); 
ylim([-35, 47]);
hold on;
plot(0,0,'+k', 'MarkerSize',  1.0);
xlabel('AU'); ylabel('AU');
title('Newer Horizons  Earth-Jupiter-Pluto Trajectory   (Launch: 1-Jan-2030)', 'FontSize', 16);
%% 
set(0,'defaultAxesFontSize',14);
ax = gca;
ax.FontSize = 12;
ax.FontWeight = 'bold';
disp('complete');


