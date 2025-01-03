%{
Function Name:  SolarSystemMission_tool_TourPlanner.m
Purpose:        This functionThis program computes the basic interplanetary mission 
                trajectory components- positions of the planets at departure and 
                arrival, and flyby, all velocity vectors (heliocentric/planetary(inf)) 
                of the spacecraft at departure and arrival, and the elements of the 
                 transfer orbits (a, e etc.)

Function Inputs:
Parameter                                   Description
-----------------           --------------------------------------------
N/A - All variable defined inside of this script

Function Outputs:
Parameter                                  Description
-----------------           ---------------------------------------------


External References:
Function Name                                Purpose
-------------               ----------------------------------------------
N/A

Development History:
    Name               Date                  Description of Change
------------        ----------          -----------------------------
John Reagoso        12/20/2010          Created.
                     
                    09/09/2012          Completed adding final config.m
                                        elements/ Updated hdf5 output
                                        parameters.

                    09/19/2012          Added FreeFlyer executable
                                        location string.
  
                    01/09/2013          Integrated with TourPlanner
                                        Configuration m-file.
 
                    02/15/2013          Separated Flyby_tool functionality
                                        as a separate m-file (for future
                                        increased modularity supporting
                                        MGAs.
 
                    04/18/2013          Updated .mat file data x-fer output
 
                    05/28/2015          Modified Planet-B flyby
                                        check-criteria
                                        
 %--------------------------------------------------------------------------
 %}
clc
fclose all
clear all

%% Start the clock.. 
vClockTimeStart = clock;

%% Planetary parameters required
global mu_sun;
global mu_earth;
global mu_venus;
global mu_mars;
global mu_mercury;
global mu_jupiter;
global mu_saturn;
global mu_neptune;
global mu_uranus;
global mu_pluto;
% global mu_asteroid_select;

mu_sun     = 1.327124e11;
mu_earth   = 3.9865e5;

mu_mercury = 0.0553* mu_earth; 
mu_venus   = 0.8149* mu_earth;

mu_mars    = 0.1074* mu_earth; 
mu_jupiter = 317.938* mu_earth;
mu_saturn  = 95.181* mu_earth;
mu_uranus  = 14.531* mu_earth;
mu_neptune = 17.135* mu_earth;
mu_pluto   = 0.0022* mu_earth;

%% Mathematical Parameters
global pi;
pi = 3.14159265;

%%  Configuration Inputs

[ConfigParameters] = SolarSystemMissionTourPlanner_ConfigFile();

%% Departure window(s):user-defined mission elements provided by configuration file:
%initial:
YearDepInitial   = ConfigParameters{5,1};
MonthDepInitial  = ConfigParameters{6,1};
DayDepInitial    = ConfigParameters{7,1};

%final:
YearDepFinal  = ConfigParameters{8,1};
MonthDepFinal = ConfigParameters{9,1};
DayDepFinal   = ConfigParameters{10,1};

Code_obj    = 2 ;  % 1 for planet-1 to planet-2, with no fly-by to a planet-3.
                    % 2 for planet-1 to planet-2, with fly-by to a planet-3.

sCode_obj= 'No_Flyby';
%sCode_obj = '1_Flyby';
%sCode_obj = '2_Flyby';
%sCode_obj = '3_Flyby';
%sCode_obj = '4_Flyby';

FlyBy_counter = 0;
      
%% Sensitivity Settings:
%   Here we set the basic mission parameter(s)- min/max flight time outbound, min/max flight 
%   time inbound, max delta-V, max v_in_vector delta:

    MinTOF         = ConfigParameters{11,1}; 
    MaxTOF         = ConfigParameters{12,1};  

    process_count   = 0     ;
    DeltaV_limit   = ConfigParameters{4,1};
        
    MinFlyTOF      = ConfigParameters{13,1};     % min return segment TOF post-flyby to planet "3".
    MaxFlybyDays   = ConfigParameters{14,1};     % max return segment TOF post-flyby to planet "3".
    
    DeltaV_infTol  = ConfigParameters{17,1};     % might need to re-evaluate this one..
    DecadeWindow   = 1     ; 

    vDeltaV_limit = ConfigParameters{18,1};
    
%% time travel not possible with this code:
if (YearDepFinal < YearDepInitial)
    
    disp('error- final departure input listed pre initial departure date. Time travel not possible yet');
    
end
%% Julian day conversion
[JulianDayDepInitial,  Julian_cent_epoch_dep_initial]    = Julian_datecalc(YearDepInitial,  MonthDepInitial,  DayDepInitial);
[JulianDayDepartFinal, Julian_cent_epoch_dep_final]      = Julian_datecalc(YearDepFinal, MonthDepFinal, DayDepFinal);

%The following function utilizes the Julian Date computation m-file to determine the Julian date
%of the desired departure timeframe.
PlanetDepart      = ConfigParameters{1,1};
PlanetArrival1    = ConfigParameters{2,1};
PlanetArrival2    = ConfigParameters{3,1};

%This function below utlizes the JPL ephemeris date for the requested departure planet and computes
%the (e) and the eccentric anomaly for the arrival planet
%% Departure planet and arrival planet heliocentric position/velocity elements plus delta-angle/ TOF

%% Recommend utilizing a double FOR - loop process for both the departure and arrival times as the
% software "travels" through every potential mission possibility based on the mission template/
%time_span = abs(Julian_day_dep_date - Julian_day_arr_date) ;

JulianDayDepartListVector = JulianDayDepInitial:1.0:JulianDayDepartFinal;
JJ_end = length(JulianDayDepartListVector) ;

aJulianDayArrivalListVector = JulianDayDepartListVector(1):1.0:(JulianDayDepartListVector(JJ_end) + MaxTOF) ;

vPlot_count = 1 ;
element_entries = JJ_end* (MaxTOF - MinTOF + 1) ;

M_count1 = 0;

%porkchop_plot_matrix = zeros(181,1);
% mission_specs = zeros(30000, 8) ;  % believe that 65,536 rows are the max for an Excel spreadsheet.

vFlyByCounter = 0;
vPkNum = 0;
MissionSpecs = {};
MissionSpecsExSum = {};


%% for-loops that cycle through all combinations of departure and arrival date(s):
for jj = 1:1.0:JJ_end   %how can datewe cycle through al  l these computations using a vectorization scheme, instead a for-loop?

    JulianDayDep = JulianDayDepartListVector(jj) ;
    [GregDepartDay, GregDepartMonth, GregDepartYear] = JulianDay_to_Greg_cal(JulianDayDep) ;
    disp('======================================================================================');
    disp('departure date is:') ;    disp(GregDepartDay) ;
    disp('======================================================================================');   
    
    vDate_num                 = double(GregDepartDay) ;
    vJulianCentEpochDepart   = ((JulianDayDep - 2451545.0)/ 36525);
    
    % Date arrival loop, max_TOF iteration for every departure date starting from the earliest allowed arrival date (min TOF) to the
    % latest possible arrival date (max_TOF):
    
        kk_init = jj + MinTOF ;
        kk_end  = jj + MaxTOF ;
        
        traj_counter = 0 ;
        
    for kk = kk_init:1.0:kk_end    %for kk = min_TOF: 1.0: max_TOF
        
        % Arrival planet heliocentric orbtial parameters/ position/velocity determination:
        JulianDayArrival = aJulianDayArrivalListVector(kk) ;
        [GregArrivalDate, ~, GregArrivalDay, GregArrivalYear] = JulianDay_to_Greg_cal(JulianDayArrival) ; %#ok<*NASGU,*ASGLU>

        fprintf('Depart date:   %s  || Arrive-1 date:   %s\n', GregDepartDay,  GregArrivalDate)
        
        TOF = JulianDayArrival - JulianDayDep;
        
        JulianCentEpochArrival = ((JulianDayArrival - 2451545.0)/ 36525);
      
        [InfVelVector1, InfVelVector2, PosVectorDepart, VelVectorDepart, PosVectorAtArrival, VelVectorAtArrival, TrueLongAtDepart, TrueLongAtArrival, ThetaSweptRad, Eccen, Incl, SMA, Energy1, TrueAnomDep, TrueAnomArr]...
        = patched_sc_vector_elements_compute(PlanetDepart, vJulianCentEpochDepart, PlanetArrival1, JulianCentEpochArrival, JulianDayDep, JulianDayArrival) ;

        %%  Mission spacecraft energy/delta-v/ C3 energy/ delta-v requirements:
        C3Mgntd           = (norm(InfVelVector1))^2 ;
        vV_sc_p            = sqrt(C3Mgntd + 2*mu_earth/6578);
        scVelEarthOrbit   = sqrt(mu_earth/6578);

        vV_sc_p            = sqrt(C3Mgntd + 2*42828.37/3776);
        scVelEarthOrbit   = sqrt(42828.37/3776);



        delta_v           = vV_sc_p - scVelEarthOrbit;
        
        Vel_inf_Norm = norm(InfVelVector1);
        
        vPkNum = vPkNum + 1;
        PlanetSeqStickData{vPkNum,1} = {JulianDayDep, JulianDayArrival, SMA, Eccen, Incl, InfVelVector1(1),InfVelVector1(2),InfVelVector1(3), InfVelVector2(1),InfVelVector2(2),InfVelVector2(3),delta_v};
        
        %PorkChopPlot{vPkNum,1} = {GregDepartDay, GregArrivalDate, JulianDayDep, JulianDayArrival, C3Mgntd};
        %PorkChopPlot(vPkNum,:)  = [JulianDayDep, JulianDayArrival, C3Mgntd];
        
 %%     Flyby(s) Mission Architecture computed here...      
 vFlagCheck = 0;
 
 if (delta_v <= vDeltaV_limit) && strcmp(sCode_obj, 'No_Flyby')== 0      
% if (Vel_inf_Norm <= 8.00 && strcmp(sCode_obj, 'No_Flyby')== 0)       
     
   %  disp(sGregArrivalDate);
            
            [vFlybyComplete, caFlyBySpecsAdd, vFlagCheck] = Flyby_tool(TOF, PosVectorAtArrival, VelVectorAtArrival, InfVelVector2, ...
                MinFlyTOF, JulianDayArrival,MaxFlybyDays, JulianCentEpochArrival, PlanetArrival1, PlanetArrival2);

            vFlyByCounter = vFlyByCounter + 1;
            
    if (vFlagCheck == 1)
        
        PreCheckSize = size(MissionSpecs,1);        
        vSpecRowCount = size(caFlyBySpecsAdd,1);
    
            for vRows2Add = 1:1:vSpecRowCount 
                
                    MissionSpecs{PreCheckSize +vRows2Add, 1}  = JulianDayDep;
                    MissionSpecs{PreCheckSize +vRows2Add, 2}  = JulianDayArrival ;
                    MissionSpecs{PreCheckSize +vRows2Add, 3}  = GregDepartDay; 
                    MissionSpecs{PreCheckSize +vRows2Add, 4}  = GregArrivalDate;
                    MissionSpecs{PreCheckSize +vRows2Add, 5}  = TOF;     
                    MissionSpecs{PreCheckSize +vRows2Add, 6}  = SMA;    
                    MissionSpecs{PreCheckSize +vRows2Add, 7}  = Eccen;    
                    MissionSpecs{PreCheckSize +vRows2Add, 8}  = Incl *180/pi;
                    MissionSpecs{PreCheckSize +vRows2Add, 9}  = Energy1;        
                    MissionSpecs{PreCheckSize +vRows2Add, 10} = TrueAnomDep *180/pi;
                    MissionSpecs{PreCheckSize +vRows2Add, 11} = TrueAnomArr *180/pi;
                    MissionSpecs{PreCheckSize +vRows2Add, 12} = ThetaSweptRad *180/pi;
                    MissionSpecs{PreCheckSize +vRows2Add, 13} = TrueLongAtDepart *180/pi;
                    MissionSpecs{PreCheckSize +vRows2Add, 14} = TrueLongAtArrival *180/pi;
                    MissionSpecs{PreCheckSize +vRows2Add, 15} = norm(PosVectorDepart);
                    MissionSpecs{PreCheckSize +vRows2Add, 16} = norm(PosVectorAtArrival);
                    MissionSpecs{PreCheckSize +vRows2Add, 17} = norm(VelVectorDepart); 
                    MissionSpecs{PreCheckSize +vRows2Add, 18} = norm(VelVectorAtArrival); 
                    MissionSpecs{PreCheckSize +vRows2Add, 19} = norm(InfVelVector1); 
                    MissionSpecs{PreCheckSize +vRows2Add, 20} = norm(InfVelVector2);
                    MissionSpecs{PreCheckSize +vRows2Add, 21} = C3Mgntd;
                    MissionSpecs{PreCheckSize +vRows2Add, 22} = delta_v;
                    
                    MissionSpecs{PreCheckSize +vRows2Add, 23} = caFlyBySpecsAdd{vRows2Add, 1};    % PreFlybyVector... must match with Row-20..
                    MissionSpecs{PreCheckSize +vRows2Add, 24} = caFlyBySpecsAdd{vRows2Add, 21};   % JulianDayArriveAtB
                    MissionSpecs{PreCheckSize +vRows2Add, 25} = caFlyBySpecsAdd{vRows2Add, 3};    % sGregDateArrival
                    MissionSpecs{PreCheckSize +vRows2Add, 26} = caFlyBySpecsAdd{vRows2Add, 4};    % TOFpastFlyby
                    MissionSpecs{PreCheckSize +vRows2Add, 27} = caFlyBySpecsAdd{vRows2Add, 22};   % vSMAarc2
                    MissionSpecs{PreCheckSize +vRows2Add, 28} = caFlyBySpecsAdd{vRows2Add, 5};    % Eccen Arc2
                    MissionSpecs{PreCheckSize +vRows2Add, 29} = caFlyBySpecsAdd{vRows2Add, 6};    % InclArc 2
                    
                    MissionSpecs{PreCheckSize +vRows2Add, 30} = caFlyBySpecsAdd{vRows2Add, 23};   % vEnergyArc
                    MissionSpecs{PreCheckSize +vRows2Add, 31} = caFlyBySpecsAdd{vRows2Add, 19};   % TrueAnomAtC
                    MissionSpecs{PreCheckSize +vRows2Add, 32} = caFlyBySpecsAdd{vRows2Add, 7};    %  %% LOOK HERE !!!
                    MissionSpecs{PreCheckSize +vRows2Add, 33} = caFlyBySpecsAdd{vRows2Add, 15};   % TruLongAtC
                    MissionSpecs{PreCheckSize +vRows2Add, 34} = caFlyBySpecsAdd{vRows2Add, 8};    % norm of PostFlybyInfVec
                    
                    MissionSpecs{PreCheckSize +vRows2Add, 35} = caFlyBySpecsAdd{vRows2Add, 9};    % vDeltaVInf
                    MissionSpecs{PreCheckSize +vRows2Add, 36} = caFlyBySpecsAdd{vRows2Add, 10};   % TurnAngleRad
                    MissionSpecs{PreCheckSize +vRows2Add, 37} = caFlyBySpecsAdd{vRows2Add, 11};   % EccenTrunRad
                    MissionSpecs{PreCheckSize +vRows2Add, 38} = caFlyBySpecsAdd{vRows2Add, 12};   % vRpTrialRad
                    MissionSpecs{PreCheckSize +vRows2Add, 39} = caFlyBySpecsAdd{vRows2Add, 13};   % VTurnComp
                    MissionSpecs{PreCheckSize +vRows2Add, 40} = caFlyBySpecsAdd{vRows2Add, 1};    % norm of PreFlybyInfVec
                    MissionSpecs{PreCheckSize +vRows2Add, 41} = caFlyBySpecsAdd{vRows2Add, 26};   % norm of Inv Velocity Vector at Planet-C
                    MissionSpecs{PreCheckSize +vRows2Add, 42} = caFlyBySpecsAdd{vRows2Add, 28};   % VinfFly Delta Percentage 

                    MissionSpecsExSum{PreCheckSize +vRows2Add, 1}  = GregDepartDay; 
                    MissionSpecsExSum{PreCheckSize +vRows2Add, 2}  = GregArrivalDate;
                    MissionSpecsExSum{PreCheckSize +vRows2Add, 3}  = TOF;     
                    MissionSpecsExSum{PreCheckSize +vRows2Add, 4}  = Energy1;        
                    MissionSpecsExSum{PreCheckSize +vRows2Add, 5}  = ThetaSweptRad *180/pi;
                    MissionSpecsExSum{PreCheckSize +vRows2Add, 6}  = TrueLongAtDepart  *180/pi;
                    MissionSpecsExSum{PreCheckSize +vRows2Add, 7}  = TrueLongAtArrival *180/pi;
                    MissionSpecsExSum{PreCheckSize +vRows2Add, 8}  = delta_v;
                    MissionSpecsExSum{PreCheckSize +vRows2Add, 9}  = caFlyBySpecsAdd{vRows2Add,3};   % sGregDateArrival
                    MissionSpecsExSum{PreCheckSize +vRows2Add, 10} = caFlyBySpecsAdd{vRows2Add,4};   % TOFpastFlyby
                    MissionSpecsExSum{PreCheckSize +vRows2Add, 11} = caFlyBySpecsAdd{vRows2Add,2};   % Energy2 
                    MissionSpecsExSum{PreCheckSize +vRows2Add, 12} = TOF + MissionSpecsExSum{PreCheckSize +vRows2Add, 10};  %Total TOF
                    MissionSpecsExSum{PreCheckSize +vRows2Add, 13} = caFlyBySpecsAdd{vRows2Add,26};  % norm(Inv VelVec at Planet-C)
                    MissionSpecsExSum{PreCheckSize +vRows2Add, 14} = caFlyBySpecsAdd{vRows2Add, 12}; % Flyby Distance
                    
                    if MissionSpecs{PreCheckSize +vRows2Add, 38} < 8000 && MissionSpecsExSum{PreCheckSize +vRows2Add, 12} <= 400
                       disp(MissionSpecs{PreCheckSize +vRows2Add, 38});
                       disp(MissionSpecsExSum{PreCheckSize +vRows2Add, 12});
                    %  pause; 
                    end

                    %break; 
                    
            end    
    end
    
 end
%     if (vFlagCheck == 1)
%         
%         PreCheckSize = size(caMissionSpecs,1);        
%         vSpecRowCount = size(caFlyBySpecsAdd,1);
%     
%         for vRows2Add = 1:1:vSpecRowCount 
%             caMissionSpecs{PreCheckSize +vRows2Add,1} = caFlyBySpecsAdd(vRows2Add,:);
%         end      
%     end
%   axis([xmin,xmax,ymin,ymax])
    
    end
    
end

sMatFileNameSpecs = ['C:\Users\JREAGOSO\Desktop/', 'StickTime2a.mat'];
save(sMatFileNameSpecs, 'PlanetSeqStickData');

vClockTimeEnd = clock;
vTour_complete = 1;

sSpecsName = [PlanetDepart,'_',PlanetArrival1,'_',PlanetArrival2, '_', num2str(vClockTimeEnd(1,1)),'_',num2str(vClockTimeEnd(1,2)),...
    '_',num2str(vClockTimeEnd(1,3)),'_',num2str(vClockTimeEnd(1,4)), num2str(vClockTimeEnd(1,5)),'_RunMissionSpecs'];

sMatFileNameSpecs = ['C:\Users\JREAGOSO\Desktop/', sSpecsName, '.mat'];
sMatFileNameSpecsExecSum = ['C:\Users\JREAGOSO\Desktop/', sSpecsName, 'ExecSum.mat'];

save(sMatFileNameSpecs, 'MissionSpecs');
save(sMatFileNameSpecsExecSum, 'MissionSpecsExSum');

%save('C:\Users\JREAGOSO\Desktop\DSWG SVN\Deep Space\DSWG Projects\Manned Venus Flyby\Matlab Software/EVE2018v223apr.mat'          ,'caMissionSpecs') ;
%save('C:\Users\JREAGOSO\Desktop\DSWG SVN\Deep Space\DSWG Projects\Manned Venus Flyby\Matlab Software/EVE2018ExecSumv223Apr.mat'   ,'caMissionSpecsExSum') ;

xls_file_string1 = ['check_Earth_to_', PlanetArrival1, sSpecsName, '.xls'];
xls_file_string2 = ['check_Earth_to_', PlanetArrival1, sSpecsName, 'ExecSum.xls'];
    
try
    
xlswrite(xls_file_string1, MissionSpecs, 'Specs_', 'A1') ;
xlswrite(xls_file_string2, MissionSpecsExSum, 'SpecsSum_', 'A1') 
catch %#ok<CTCH>
    disp('                                                                 ');
    disp('=================================================================');
    disp('was unable to successfully transfer .mat(s) to Excel spreadsheets');
    disp('=================================================================');
    
end
    
%% Report Clock Start/End Time(s):
disp('                                                                 ');
disp('SSMD code started:');
disp(vClockTimeStart);
disp('SSMD code finished:');
disp(vClockTimeEnd);

%% Plotting Results:

% MJD Departure vs. Total TOF:
vSpecsSize = size(MissionSpecs,1);

for ii = 1:1:vSpecsSize

    mTOF(ii,1)          = MissionSpecs{ii,5} + MissionSpecs{ii,26};
    mDepartMJD(ii,1)    = MissionSpecs{ii,1};    
    mArriveMJD(ii,1)    = mDepartMJD(ii,1) + mTOF(ii,1);
    mC3(ii,1)           = MissionSpecs{ii,21};
    mDeltaV(ii,1)       = MissionSpecs{ii,22};
    mPlCInfVel(ii,1)    = MissionSpecs{ii,41};
    aa(ii,1)            = ii;
    flyby_rp(ii,1)      = MissionSpecs{ii,38};

end

% jj = 1;
% for ii = 1:1:vSpecsSize
% 
%     Vel_InfFlyby_in_norm  = MissionSpecs{ii, 40};
%     Vel_InfFlyby_out_norm = MissionSpecs{ii, 34};
% 
%     VelInf_diff = abs(Vel_InfFlyby_out_norm - Vel_InfFlyby_in_norm);
%     VelInf_halfsum_avg = abs((Vel_InfFlyby_out_norm + Vel_InfFlyby_in_norm)/2);
% 
%     vInfVel_delta_percent(ii) = VelInf_diff/VelInf_halfsum_avg;
%     
%     if MissionSpecs{ii,38} > 5000.0
%     %if vInfVel_delta_percent(ii) < 0.01
% 
%         mTOF2(jj,1)          = MissionSpecs{ii,5} + MissionSpecs{ii,26};
%         mDepartMJD2(jj,1)    = MissionSpecs{ii,1};    
%         mArriveMJD2(jj,1)    = mDepartMJD(ii,1) + mTOF(ii,1);
%         mC32(jj,1)           = MissionSpecs{ii,21};
%         jj = jj + 1;
%     end
% 
% end

try 
    figure;
    mTOFcorr = mTOF/365.0;

    pointsize = 12;
    % scatter(mDepartMJD, mArriveMJD, pointsize, mC3,'o', 'filled'); 
    scatter(mDepartMJD, mTOFcorr, pointsize, mC3,'o', 'filled'); 
    xlabel('MJD Departure'); 
    ylabel('Total TOF(days)'); 
    %axis([caPorkChopPlot{1}{3}, caPorkChopPlot{size(caPorkChopPlot,1)}, 300, 550])   %    axis square; 
    %axis square;
    grid on; 
    %ylim([1.0, 4.0])
    title 'Earth-Mars-Earth Free Return Trajectories (Jan 1995 - Dec 2004)'; 
    legend 'C3 (km^2/sec^2)'; 
    colorbar

catch
    disp('Code unable to provide plotting graphics.');
    
end

try
   pointsize = 20;
   mPorkChopDepart    = PorkChopPlot(:,1);
   mPorkkChopArrival  = PorkChopPlot(:,2);
   mPorkChopC3        = PorkChopPlot(:,3);
   
   scatter(mPorkChopDepart, mPorkkChopArrival, pointsize, mPorkChopC3, 'x');
catch
    disp('here');
end

%% Here is the porkchop plot generator. 
% vArr_date_timespan = (vMaxTOF - vMinTOF)+ jj;
% 
% vMissionPlot_length = size(aMission_plot,1);
% 
% vFirstDepDate = aMission_plot{1,1};
% vLastDepDate  = aMission_plot{vMissionPlot_length,1};
% 
% vFirstArrDate = aMission_plot{1,2};
% vLastArrDate  = aMission_plot{vMissionPlot_length,2};
% 
% vMission_dep_span = vLastDepDate - vFirstDepDate + 1;
% vMission_arr_span = vLastArrDate - vFirstArrDate + 1;
% 
% caPorkChopMJDmap = cell(vMission_dep_span, vMission_arr_span);
% vPorkChopMJDmap_size = numel(caPorkChopMJDmap);
% 
% 
% aPorkChopSpecs = {sPlanetDepart; sPlanetArrival1; vMonthDepInitial; vYearDepInitial; vMonthDepFinal; vYearDepFinal};
% 
% PorkChopPlotGenerator(aMission_plot, vArr_date_timespan, vMinTOF);

% save('C:\Users\JREAGOSO\Desktop\DSWG SVN\Deep Space\DSWG Projects\Manned Venus Flyby\Matlab Software/mission_specs_2way_2020_2030.mat'  ,'mission_specs_flyby') ;
% date_max = Julian_day_arr ;
% 
% % xls_file_string1 = ['check_Earth_to_',    planet_arrival,    '.xls' ];
% % xlswrite(xls_file_string1, mission_plot, 'Earth_Venus_specs_', 'A1') 
% 
% % xls_file_string2 = ['Earth_to_', planet_arrival, 'w_flyby_to', planet_arrival_3, '.xls' ];
% % xlswrite(xls_file_string2, mission_specs, 'Earth_Venus_specs_', 'A1') 
% 
% delta_flight = max_TOF- min_TOF ;
% disp(mission_plot) ; 
% mission_plot = mission_plot' ;
% 
% %multi_end = length(Julian_day_dep_vector)/(max_TOF-min_TOF) ;
% multi_end = 10 ;
% 
% 
%        plot_kk = plot_kk + 1;
%         vCount = vCount + 1;
%         if traj_counter == 1
%             break ; 
%         end
% 
% for ii = 1:1:length(Julian_day_dep_vector)
% 
%     multiple = ii* delta_flight ;
%     beginner = multiple - (delta_flight - 1) ;
%     
%     pork_chop_plot                  = mission_plot(beginner:multiple) ;
%     pork_chop_plot_matrix(:,ii)     = flipud(pork_chop_plot) ;
%         
%     
% end
% 
% xls_file_string1 = ['check_Earth_to_', planet_arrival, '.xls' ];
% xlswrite(xls_file_string1, mission_specs_flyby, 'Earth_Venus_specs_', 'A1') 
%  
% save('C:\Users\JREAGOSO\Desktop\DSWG SVN\Deep Space\DSWG Projects\Manned Venus Flyby\Matlab Software/mission_specs_aligned_1way_2015_2020','mission_specs_aligned') ;
% 
% xls_file_string3 = ['Earth_to_', planet_arrival, '_aligned.xls' ];
% xlswrite(xls_file_string3, mission_specs_aligned, 'Earth_Venus_specs_aligned', 'A1') ;
% contour(mission_specs)   ;
% contour(pork_chop_plot_matrix) ;
% contour(mission_specs_aligned)  ;





%%      Pork-chop plot generator:                   
%             if (strcmp(sCode_obj, 'No_Flyby')== 0 )
%         
%                 aMission_plot{vPlot_count, 1} = vJulianDayDep;        %#ok<*AGROW>
%                 aMission_plot{vPlot_count, 2} = aJulianDayArrival;      
%                 aMission_plot{vPlot_count, 3} = sGregDepartDay;
%                 aMission_plot{vPlot_count, 4} = sGregArrivalDate;
%                 %   mission_plot{plot_count,5}  = C3_mgntd ;   
%         
%                     if vC3Mgntd <= 30   
%                         aMission_plot{vPlot_count,5} = vC3Mgntd;
%                 
%                     else
%                         aMission_plot{vPlot_count,5} = 0.000;
%             
%                     end
% 
%             vPlot_count = vPlot_count + 1 ;      
% 
%             end
%         
%             
%     if (vDelta_v <= vDeltaV_limit) && strcmp(sCode_obj, 'No_Flyby')== 0     
%         
%                     vDecadeWindow = 1 ;
%                     mMission_specs{vM_count1, 1,  vDecadeWindow} = sGregDepartDay;              
%                     mMission_specs{vM_count1, 2,  vDecadeWindow} = sGregArrivalDate;            
%                     mMission_specs{vM_count1, 3,  vDecadeWindow} = vSMA;                         
%                     mMission_specs{vM_count1, 4,  vDecadeWindow} = vEnergy1;                  
%                     mMission_specs{vM_count1, 5,  vDecadeWindow} = vTOF;                       
%                     mMission_specs{vM_count1, 6,  vDecadeWindow} = vEccen;                         
%                     mMission_specs{vM_count1, 7,  vDecadeWindow} = vIncl *180/pi;              
%                     mMission_specs{vM_count1, 8,  vDecadeWindow} = vThetaSweptRad* 180/pi;          
%                     mMission_specs{vM_count1, 9,  vDecadeWindow} = vC  3Mgntd;                 
%                     mMission_specs{vM_count1, 10, vDecadeWindow} = vDelta_v;                    
%                     mMission_specs{vM_count1, 11, vDecadeWindow} = aSCInfVelVector1;            
%                     mMission_specs{vM_count1, 12, vDecadeWindow} = aSCInfVelVector2;            
%                     mMission_specs{vM_count1, 13, vDecadeWindow} = vJulianDayDep;              
%                     mMission_specs{vM_count1, 14, vDecadeWindow} = aJulianDayArrival;
%                     
%                     vM_count1 = vM_count1 + 1;
%     end
 






