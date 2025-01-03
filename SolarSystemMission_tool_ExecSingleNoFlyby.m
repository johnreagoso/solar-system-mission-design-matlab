
function [PosVectorDepart, VelVectorDepart] = SolarSystemMission_tool_ExecSingleNoFlyby(...
    PlanetDepart, PlanetArrival1, JulianDayDep, JulianDayArrival1)

%{
Function Name:  SolarSystemMission_tool_ExecSingle.m

Purpose:        This functionThis program computes the basic interplanetary mission 
                trajectory components- positions of the planets at a SINGLE departure and 
                arrival, and flyby, transfer trajectory data (pos/vel vectors (heliocentric/planetary(inf)) 
                of the spacecraft at departure and arrival(s).
                Main purpose of this script is to assist with trajectory
                plotting output as well as provide framing for heuristic
                tool future development. 

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
John Reagoso        02/20/2024          Created.
                     
 %--------------------------------------------------------------------------
 %}

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
   
    vFlyByCounter = 0;
    vPkNum = 0;
    MissionSpecs = {};
    MissionSpecsExSum = {};
    
    
    %% Main Body:
    [GregDepartDay, ~, ~] = JulianDay_to_Greg_cal(JulianDayDep) ;
    disp('======================================================================================');
    disp('departure date is:') ;    disp(GregDepartDay) ;
    disp('======================================================================================');   
    
    vJulianCentEpochDepart   = ((JulianDayDep - 2451545.0)/ 36525);
    
    % Arrival planet heliocentric orbtial parameters/ position/velocity determination:
    %[GregArrival1Date, ~, ~, ~] = JulianDay_to_Greg_cal(JulianDayArrival1);

    %[GregArrival2Date, ~, GregArrival2Day, GregArrival2Year] = JulianDay_to_Greg_cal(JulianDayArrival2);

    %fprintf('Depart date:   %s  || Arrive-1 date:   %s  || Arrive-2 date:    %s\n', GregDepartDay,  GregArrival1Date, GregArrival2Date);
    %TOF1 = JulianDayArrival1 - JulianDayDep;
    %TOF2 = JulianDayArrival2 - JulianDayArrival1;

    JulianCentEpochArrival1 = ((JulianDayArrival1 - 2451545.0)/ 36525);
    %JulianCentEpochArrival2 = ((JulianDayArrival2 - 2451545.0)/ 36525);

    [InfVelVector1, ~, PosVectorDepart, VelVectorDepart, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~]...
    = patched_sc_vector_elements_compute(PlanetDepart, vJulianCentEpochDepart, PlanetArrival1, JulianCentEpochArrival1, JulianDayDep, JulianDayArrival1) ;
    
    %%  Mission spacecraft energy/delta-v/ C3 energy/ delta-v requirements:
    C3Mgntd           = (norm(InfVelVector1))^2 ;
    vV_sc_p            = sqrt(C3Mgntd + 2* mu_earth/6578); % assumes 200km alt orbit.. 
    
    scVelEarthOrbit   = sqrt(mu_earth/6578);
    delta_v           = vV_sc_p - scVelEarthOrbit;
    
    Vel_inf_Norm = norm(InfVelVector1);
    
     %%     Flyby(s) Mission Architecture computed here...      
     
%     if (Vel_inf_Norm >= 8.00)
%         disp('please check mission plan as it requires over 64 km2/sec2 in C3')
%     end
%     
%     [~, caFlyBySpecsAdd, ~] = Flyby_tool(TOF1, PosVectorAtArrival, VelVectorAtArrival, InfVelVector2, ...
%                     TOF2, JulianDayArrival1, TOF2, JulianCentEpochArrival1, PlanetArrival1, PlanetArrival2);
%     
%     MissionSpecs{1, 1}  = JulianDayDep;
%     MissionSpecs{1, 2}  = JulianDayArrival1 ;
%     MissionSpecs{1, 3}  = GregDepartDay; 
%     MissionSpecs{1, 4}  = GregArrival1Date;
%     MissionSpecs{1, 5}  = TOF1;     
%     MissionSpecs{1, 6}  = SMA;    
%     MissionSpecs{1, 7}  = Eccen;    
%     MissionSpecs{1, 8}  = Incl *180/pi;
%     MissionSpecs{1, 9}  = Energy1;        
%     MissionSpecs{1, 10} = TrueAnomDep *180/pi;
%     MissionSpecs{1, 11} = TrueAnomArr *180/pi;
%     MissionSpecs{1, 12} = ThetaSweptRad *180/pi;
%     MissionSpecs{1, 13} = TrueLongAtDepart *180/pi;
%     MissionSpecs{1, 14} = TrueLongAtArrival *180/pi;
%     MissionSpecs{1, 15} = norm(PosVectorDepart);
%     MissionSpecs{1, 16} = norm(PosVectorAtArrival);
%     MissionSpecs{1, 17} = norm(VelVectorDepart); 
%     MissionSpecs{1, 18} = norm(VelVectorAtArrival); 
%     MissionSpecs{1, 19} = norm(InfVelVector1); 
%     MissionSpecs{1, 20} = norm(InfVelVector2);
%     MissionSpecs{1, 21} = C3Mgntd;
%     MissionSpecs{1, 22} = delta_v;
%     
%     MissionSpecs{1, 23} = caFlyBySpecsAdd{1, 1};    % PreFlybyVector... must match with Row-20..
%     MissionSpecs{1, 24} = caFlyBySpecsAdd{1, 21};   % JulianDayArriveAtB
%     MissionSpecs{1, 25} = caFlyBySpecsAdd{1, 3};    % sGregDateArrival
%     MissionSpecs{1, 26} = caFlyBySpecsAdd{1, 4};    % TOFpastFlyby
%     MissionSpecs{1, 27} = caFlyBySpecsAdd{1, 22};   % vSMAarc2
%     MissionSpecs{1, 28} = caFlyBySpecsAdd{1, 5};    % Eccen Arc2
%     MissionSpecs{1, 29} = caFlyBySpecsAdd{1, 6};    % InclArc 2
%     
%     MissionSpecs{1, 30} = caFlyBySpecsAdd{1, 23};   % vEnergyArc
%     MissionSpecs{1, 31} = caFlyBySpecsAdd{1, 19};   % TrueAnomAtC
%     MissionSpecs{1, 32} = caFlyBySpecsAdd{1, 7};    %  %% LOOK HERE !!!
%     MissionSpecs{1, 33} = caFlyBySpecsAdd{1, 15};   % TruLongAtC
%     MissionSpecs{1, 34} = caFlyBySpecsAdd{1, 8};    % norm of PostFlybyInfVec
%     
%     MissionSpecs{1, 35} = caFlyBySpecsAdd{1, 9};    % vDeltaVInf
%     MissionSpecs{1, 36} = caFlyBySpecsAdd{1, 10};   % TurnAngleRad
%     MissionSpecs{1, 37} = caFlyBySpecsAdd{1, 11};   % EccenTrunRad
%     MissionSpecs{1, 38} = caFlyBySpecsAdd{1, 12};   % vRpTrialRad
%     MissionSpecs{1, 39} = caFlyBySpecsAdd{1, 13};   % VTurnComp
%     MissionSpecs{1, 40} = caFlyBySpecsAdd{1, 1};    % norm of PreFlybyInfVec
%     MissionSpecs{1, 41} = caFlyBySpecsAdd{1, 26};   % norm of Inv Velocity Vector at Planet-C    
%      MissionSpecs{1, 42} = caFlyBySpecsAdd{1, 27};   % Post flyby Heliocentric vel vector at Planet-B
% %                  
%      sc_Heliocen_postFlyby =  MissionSpecs{1, 42};
% % 
% % 
%     MissionSpecsExSum{1, 1}  = GregDepartDay; 
%     MissionSpecsExSum{1, 2}  = GregArrival1Date;
%     MissionSpecsExSum{1, 3}  = TOF1;     
%     MissionSpecsExSum{1, 4}  = Energy1;        
%     MissionSpecsExSum{1, 5}  = ThetaSweptRad *180/pi;
%     MissionSpecsExSum{1, 6}  = TrueLongAtDepart  *180/pi;
%     MissionSpecsExSum{1, 7}  = TrueLongAtArrival *180/pi;
%     MissionSpecsExSum{1, 8}  = delta_v;
%     MissionSpecsExSum{1, 9}  = caFlyBySpecsAdd{1,3};   % sGregDateArrival
%     MissionSpecsExSum{1, 10} = caFlyBySpecsAdd{1,4};   % TOFpastFlyby
%     MissionSpecsExSum{1, 11} = caFlyBySpecsAdd{1,2};   % Energy2 
%     MissionSpecsExSum{1, 12} = TOF1 + MissionSpecsExSum{1, 10};  %Total TOF
%     MissionSpecsExSum{1, 13} = caFlyBySpecsAdd{1,26};  % norm(Inv VelVec at Planet-C)
%     MissionSpecsExSum{1, 14} = caFlyBySpecsAdd{1, 12}; % Flyby Distance


end

  
