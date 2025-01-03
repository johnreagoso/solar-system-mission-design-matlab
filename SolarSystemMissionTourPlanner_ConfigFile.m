%{
    ##########################__79_Character_Header__###########################
    PROCEDURE NAME: SolarSystemMissionTourPlanner_Config.m

    PURPOSE:    The following code below allow the operator to trigger the Solary
                System Tour Planner, and set up specific configuration settings, 
                based on the required desired analysis (number of fly-by(s), 
                departure/arrival/flyby(s) windows, number of flybys, output file 
                location etc.

	External References:
 	Name                                Purpose
	----                                -------
    n/a                                  n/a

    INPUTS:
		1.  User-selection input          -
	
	OUTPUTS:
		1. 	Cell Array (17x1)           Transfer User-Defined inputs to main 
          ('caConfigParameters')        functional script (SolarSystemMission_tool_TourPlanner.m)
 
  ##########################__79_Character_Header__###########################
   External References:
   NAME								PURPOSE
   ----                             -------  (none currently)
 
   Development History:
   NAME				DATE			DESCRIPTION OF CHANGE
   Reagoso, J.D.  	2012/1/9        Script initialized  
 ---------------------------------------------------------------------------------------------------
%}

function [caConfigParameters] = SolarSystemMissionTourPlanner_ConfigFile()
    caConfigParameters = cell(17,1);
    
%%  Departure Planet:
    sDeparture_planet = 'Mars';

%%  Arrival-1 Planet (or 1st Fly-by planet):
    sArrivalPlanet1   = 'Earth';
    
%%  Arrival-2 Planet (2nd Fly-by, or return  planet):
    sArrivalPlanet2   = 'Earth';
    
%%  Arrival-2 Planet (2nd Fly-by, or return  planet):
    sArrivalPlanet3   = 'Saturn';
    
%%  Enter the max-allowable delta-v permitted, 
%   This values represents the maximum impulsive-burn propulsion parameter that will be used to power 
%   the first ballistic trajectory-leg. The code will not analyze follow-on trajectory analysis that 
%   required delta-v(s) above this value:
    vMaxDeltaVAllowed = 6.0; % in km/sec
   
%%  Enter departure window data in following Gregorian format (yyyy, m(m))
%   Provide the departure window for the first planet, i.e the earliest possible departure date, 
%   vDepartureInitialYear(Month, day) etc., followed by the last possible departure data for analysis,
%   vDepartureFinalYear(Month, day) etc.:

%   (Earliest date):
    vDepartureInitialYear   = 2030;
    vDepartureInitialMonth  = 8;
    vDepartureInitialDay    = 9;
    
%   (Final date):
    vDepartureFinalYear     = 2031;
    vDepartureFinalMonth    = 1;
    vDepartureFinalDay      = 12;
       
%   ("Remember- 30 days has September, April, May.. and  November.")
%% Enter the first trajectory min/max flight time(s)(days);
%  recommend keeping this day no less 60 days to ensure continuous stability in the Lamberts solver 
%  embedded within this code.
    vMinimumTOF         = 364;
    vMaxiumumTOF        = 371;
    
%% Enter the second, post-flyby trajectory-leg min/max flight time(s)(days);
    vCodeObjOneFlyby    = 1;    % 1- yes, 0- no
    
    vMinimumFlybyTOF    = 60; 
    vMaxiumumFlybyTOF   = 365*1;
%%    
    vFlybyMinRp         = 0.0;  % km
    
%% Enter the third, post-2nd flyby trajectory-leg min/max flight time(s)(days);
    vCodeObjTwoFlyby    = 0;   % 1- yes, 0- no
    vMinimumFlyby2TOF   = 60; 
    vMaxiumumFlyby2TOF  = 365*1;         
    vFlyby2MinRp        = 0; % km
 
%% Enter the fourth, post-3nd flyby trajectory-leg min/max flight time(s)(days);
%     vCodeObjThreeFlyby  = 0;  % 1- yes, 0- no
%     vMinimumFlyby3TOF   = 90; 
%     vMaxiumumFlyby3TOF  = 270;    
%     vFlyby3MinRp        = 500; % km
   
    %%  Flyby Sensitivity Parameters:
%   Enter the percentage of permitted 'slop' in patching the pre-flyby V_infinity(in) magnitude, 
%   with the post-flyby V_infinity(out) magnitude. The larger the percentage, the increased chance 
%   that a particular trajectory patch will be matched for a mission profile.
%   Likewise, there is an increased risk of patching (2) trajectories that are not realistic physically.
    vFlybyPatchTolerance  = 0.025;  % Both vectors have to match within (100 percent, minus this value)
    

%   Cell Array passed to main functional script below:    
    caConfigParameters{1,1}  = sDeparture_planet;
    caConfigParameters{2,1}  = sArrivalPlanet1;
    caConfigParameters{3,1}  = sArrivalPlanet2;
    
    caConfigParameters{4,1}  = vMaxDeltaVAllowed;
    caConfigParameters{5,1}  = vDepartureInitialYear;
    caConfigParameters{6,1}  = vDepartureInitialMonth;
    caConfigParameters{7,1}  = vDepartureInitialDay;
    
    caConfigParameters{8,1}  = vDepartureFinalYear;
    caConfigParameters{9,1}  = vDepartureFinalMonth;
    caConfigParameters{10,1} = vDepartureFinalDay;     
    
    caConfigParameters{11,1} = vMinimumTOF;
    caConfigParameters{12,1} = vMaxiumumTOF;
    
    caConfigParameters{13,1} = vMinimumFlybyTOF;
    caConfigParameters{14,1} = vMaxiumumFlybyTOF;
    
    caConfigParameters{15,1} = vMinimumFlyby2TOF;
    caConfigParameters{16,1} = vMaxiumumFlyby2TOF;
    
    caConfigParameters{17,1} = vFlybyPatchTolerance;
    caConfigParameters{18,1} = vMaxDeltaVAllowed;


end




