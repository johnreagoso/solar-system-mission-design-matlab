

%%  Solar System Flyby Tool
%   This tool will combine all the original flyby code functionality of the
%   original Solar System Tour software into (1) specific m-file, that is
%   placed within an appropriate FOR look structure in order to compute
%   multiple flyby(s).
%{
    ##########################__79_Character_Header__###########################
    PROCEDURE NAME: Flyby_tool.m

    PURPOSE:    The following code below...

	External References:
 	Name                                Purpose
	----                                -------
    n/a                                  n/a

    INPUTS:
		1.  1st trajectory arc details    -
	
	OUTPUTS:
		1. 	Cell Array (17x1)           Transfer User-Defined inputs to main 
          ('caConfigParameters')        functional script (SolarSystemMission_tool_TourPlanner.m)

  ##########################__79_Character_Header__###########################
   External References:
   NAME								PURPOSE
   ----                             -------  (none currently)
 
   Development History:
   NAME				DATE			DESCRIPTION OF CHANGE
   Reagoso, J.D.  	2012/1/09        Script initialized
                    
                    2013/4/30       vDeltaVinfTol algorithm modified

                    2023/2/13       Improved flyby altitude screening logic
 ---------------------------------------------------------------------------------------------------
%}
function[vFlybyComplete, caFlyBySpecs, vFlagCheck] = Flyby_tool(vEarlierLegTOF, aSCHelioPosVecAtFlyby, aFlybyPlanetHelioVelVec, aV2_sc_inf_vector, ...
                vMin_fly_TOF, vJulianDayArrAtPlanetA, vMax_flyby_days, vJulianCentEpochPlanetA, sPlanetB, sPlanetC)
    
    caFlyBySpecs = {1,36}; 
    vFlagCheck = 0;
            
    vDeltaVinfTol = 0.025;
    % vDeltaVinfTol = 0.0750;

    % here is the spacecraft position at fly-by at planet-2 utilizes the same position vector as the plan based on the patched-conic theory.
    aSC_pos_vec_arr_IncomingFlyby   = aSCHelioPosVecAtFlyby;   

    % here we utilize the flyby planet's heliocentric velocity [km/sec]
    aFlybyPlanetVelVecArr           = aFlybyPlanetHelioVelVec ; 
    
   % aFlybyPlanetVelVecArrEQ = ECI2EQ(aFlybyPlanetVelVecArr);
              
    % here we initilize the s/c v inf vector (pre-flyby)
    aSCPreFlyby_InfVelVector       = aV2_sc_inf_vector ;   
    
   % aSCPreFlybyV_InfVelVectorEQ = ECI2EQ(aSCPreFlybyV_InfVelVector);
     
    % here we initialize and set the first 2nd leg TOF iteration to the minimum value.
    vTOFdaysPostFlyby               = vMin_fly_TOF ;         
      
    % here we initialize the s/c flyby date 
    vFlybyJulianDayArrival= vJulianDayArrAtPlanetA ;             
    
    FlybyPlanet = sPlanetB;

    [~, PlanetRadius, MinFlybyDist, MuPlanet] = planet_specs(FlybyPlanet);
    
%%  This while loop will compute all potential trajectories possible from planet-2 to planet-3, from the shortest possible TOF 
%   to the longest possible TOF for the one specific fly-by date:

    vArcDataRowCount = 0;
            
            while (vTOFdaysPostFlyby <= vMax_flyby_days) % && traj_counter ==0 )
                                
                JulianDayArrAtPlanetB        = vFlybyJulianDayArrival + vTOFdaysPostFlyby ;
                
                JulianCentEpochPlanetB       = ((JulianDayArrAtPlanetB - 2451545.0)/ 36525);
                
                [GregDateArrivalB, ~, ~, ~]  = JulianDay_to_Greg_cal(JulianDayArrAtPlanetB) ;
                 
                [aSCPostFlyby_InfVelVector, aSCPlanetCInfVelVector, aPosVectorDepart, aVelVectorDepart, aPosVectorAtArrival2, aVelVectorAtArrival2, ...
                            vTrueLongAtPlanetB, vTrueLongAtPlanetC, vThetaSweptRad, vEccenArc2, vInclArc2, vSMAarc2, vEnergyArc2, vTrueAnomAtB, vTrueAnomAtC]...
                                                            = patched_sc_vector_elements_compute(sPlanetB, vJulianCentEpochPlanetA, sPlanetC, ...
                                                                                JulianCentEpochPlanetB, vJulianDayArrAtPlanetA, JulianDayArrAtPlanetB) ; %#ok<*ASGLU>
                                                                            
                                                                           
%                     if (norm(aPosVecArrA) ~= norm(aSCHelioPosVecAtFlyby)) || (vTrue_long_arr ~= vTrueLongAtPlanetA) || (norm(aVelVecAtA) ~= norm(aFlybyPlanetHelioVelVec)) || (vTOF ~= vTOFdaysPostFlyby)
%                         disp('mis-match in patch scheme, fix code') ; 
%                         vError_match  = 1 ;
%                     end  

         %       aSCPostFlybyInfVelVectorEQ = ECI2EQ(aSCPostFlybyInfVelVector);         
                
%                 vDeltaVinf  =  1 - abs(norm(aSCPostFlyby_InfVelVector)/norm(aSCPreFlyby_InfVelVector));
%                 vDeltaVinf2 =  1 - abs(norm(aSCPreFlyby_InfVelVector)/norm(aSCPostFlyby_InfVelVector));
%                
%                 vDeltaVinfAvg = 0.5*(abs(vDeltaVinf) + abs(vDeltaVinf2));
%                 
                Vel_InfFlyby_in_norm  = norm(aSCPreFlyby_InfVelVector);
                Vel_InfFlyby_out_norm = norm(aSCPostFlyby_InfVelVector);

                VelInf_diff = abs(Vel_InfFlyby_out_norm - Vel_InfFlyby_in_norm);
                VelInf_halfsum_avg = abs((Vel_InfFlyby_out_norm + Vel_InfFlyby_in_norm)/2);

                vInfVel_delta_percent = VelInf_diff/VelInf_halfsum_avg;

                TurnAngleRad  =  acos(dot(aSCPreFlyby_InfVelVector,aSCPostFlyby_InfVelVector )/(norm(aSCPreFlyby_InfVelVector)*norm(aSCPostFlyby_InfVelVector))) ;
                
                EccenTurnRad    =  (sin(TurnAngleRad/2))^-1 ;
                
                % because we are using a simply approximation of the v_inf vector, including an arbitrary tolerance value, we use
                % the average of the pre/pos v_inf vectors: 
                vSCVInfVelVectorAvg = 0.5*norm(aSCPreFlyby_InfVelVector) + 0.5*norm(aSCPostFlyby_InfVelVector) ;
                
                maxSCTurn   = 2*asin(1/(1 + (MinFlybyDist* vSCVInfVelVectorAvg^2)/MuPlanet)) ; 
                
                vRpFlybyDistRad  = ((EccenTurnRad - 1)* MuPlanet)/vSCVInfVelVectorAvg^2 ;
                
               % vSCTurnCompute  = 2*asin(1/(1 + vRpFlybyDistRad* vSCVInfVelVectorAvg^2 /vMuPlanet)) ; 
                vSCTurnCompute = 2*asin(1/EccenTurnRad);  %check complete on 25 June 2015 JDR
                
                Total_TOF = vEarlierLegTOF + vTOFdaysPostFlyby ;
              
                check = 0;
                %% =================================================================================================
                %% Output                
                %if (abs(vDeltaVinfAvg) <= vDeltaVinfTol) && (vSCTurnAngleRad > vMaxSCTurn) 
                if (abs(vInfVel_delta_percent) <= vDeltaVinfTol) && (TurnAngleRad > maxSCTurn)
                    %     disp('case1');
                       check = check + 1;

                    if (abs(vInfVel_delta_percent)<= vDeltaVinfTol) && (MinFlybyDist > vRpFlybyDistRad)
                 %   if (abs(vDeltaVinfAvg)<= vDeltaVinfTol) && (vMinFlybyDist > vRpFlybyDistRad)
                 %      disp('case2');
                       check = check + 1;
                    end
       
                    if check == 1
                       disp('problem here'); 
                       pause
                    end
                    
                end
                             
                %% =================================================================================================
                %if (abs(vDeltaVinfAvg) <= vDeltaVinfTol) && (vSCTurnAngleRad <= vMaxSCTurn) && (vMinFlybyDist < vRpFlybyDistRad)
                
                %if (abs(vDeltaVinfAvg) <= vDeltaVinfTol) && (vMinFlybyDist < vRpFlybyDistRad)
                % ONLY USE FOR EXAMPLE:
                %if (vInfVel_delta_percent <= vDeltaVinfTol) && (0.0 < vRpFlybyDistRad)
                
                if (vInfVel_delta_percent <= vDeltaVinfTol) && (MinFlybyDist< vRpFlybyDistRad)
                %if (vInfVel_delta_percent <= vDeltaVinfTol) && (0.0 < vRpFlybyDistRad)
               
                    vArcDataRowCount = vArcDataRowCount + 1;
                    
                    caFlyBySpecs{vArcDataRowCount, 1}  = norm(aSCPreFlyby_InfVelVector);
                    caFlyBySpecs{vArcDataRowCount, 2}  = vEnergyArc2;
                    caFlyBySpecs{vArcDataRowCount, 3}  = GregDateArrivalB; 
                    caFlyBySpecs{vArcDataRowCount, 4}  = vTOFdaysPostFlyby;
                    caFlyBySpecs{vArcDataRowCount, 5}  = vEccenArc2;     
                    caFlyBySpecs{vArcDataRowCount, 6}  = vInclArc2 *180/pi;     
                    caFlyBySpecs{vArcDataRowCount, 7}  = vThetaSweptRad *180/pi;    
                    caFlyBySpecs{vArcDataRowCount, 8}  = norm(aSCPostFlyby_InfVelVector);     
                    %caFlyBySpecs{vArcDataRowCount, 9}  = vDeltaVinf;       
                    caFlyBySpecs{vArcDataRowCount, 9}  = VelInf_diff;       
                    caFlyBySpecs{vArcDataRowCount, 10} = TurnAngleRad;
                    caFlyBySpecs{vArcDataRowCount, 11} = EccenTurnRad;
                    caFlyBySpecs{vArcDataRowCount, 12} = vRpFlybyDistRad;
                    caFlyBySpecs{vArcDataRowCount, 13} = vSCTurnCompute;
                    caFlyBySpecs{vArcDataRowCount, 14} = vTrueLongAtPlanetB *180/pi;
                    caFlyBySpecs{vArcDataRowCount, 15} = vTrueLongAtPlanetC *180/pi;
                    caFlyBySpecs{vArcDataRowCount, 16} = vSMAarc2;
                    caFlyBySpecs{vArcDataRowCount, 17} = Total_TOF; 
                    caFlyBySpecs{vArcDataRowCount, 18} = vTrueAnomAtB; 
                    caFlyBySpecs{vArcDataRowCount, 19} = vTrueAnomAtC; 
                        
                    caFlyBySpecs{vArcDataRowCount, 20} = vTrueLongAtPlanetC; 
                    caFlyBySpecs{vArcDataRowCount, 21} = JulianDayArrAtPlanetB;   
                    caFlyBySpecs{vArcDataRowCount, 22} = vSMAarc2;
                    caFlyBySpecs{vArcDataRowCount, 23} = vEnergyArc2;
                    caFlyBySpecs{vArcDataRowCount, 24} = norm(aPosVectorAtArrival2);
                    caFlyBySpecs{vArcDataRowCount, 25} = norm(aVelVectorAtArrival2);
                    caFlyBySpecs{vArcDataRowCount, 26} = norm(aSCPlanetCInfVelVector);
                    caFlyBySpecs{vArcDataRowCount, 27} = aVelVectorDepart;
                    caFlyBySpecs{vArcDataRowCount, 28} = vInfVel_delta_percent;
                    
                    vFlagCheck = 1 ;
                                        
                end
                      
                vTOFdaysPostFlyby = vTOFdaysPostFlyby + 1.0;
              
            end
vFlybyComplete = 1;
         