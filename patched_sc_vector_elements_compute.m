function[aV1scInfVector, aV2scInfVector, aPosVector1,aSCvector1, aPosVector2, aSCvector2, vTrueLong1, vTrueLong2, vThetaSwept, e, vIncl, vSMA, vEnergy, vTrueAnom1, vTrueAnom2]...
    = patched_sc_vector_elements_compute(sPlanet1, vJulianCentEpoch1, sPlanet2, vJulianCentEpoch2, vJulianDay1, vJulianDay2)

vMuSun      = 1.327124e11  ;

[aPosVector1, ~, aVelVector1, ~, vTrueLong1] = planet_orbit_parameters(sPlanet1, vJulianCentEpoch1) ;

[aPosVector2, ~, aVelVector2, ~, vTrueLong2] = planet_orbit_parameters(sPlanet2, vJulianCentEpoch2) ;

vPosVectorMagntd1 = norm(aPosVector1) ; 
vPosVectorMagntd2 = norm(aPosVector2) ;

%this following portion of code assumes the s/c will depart earth in a direct trajectory:
vCross = aPosVector1(1)*aPosVector2(2) - aPosVector1(2)*aPosVector2(1) ; 

% Spacecraft angle swept (theta_swept) computation:
    if (vCross >= 0 )
    
        vThetaSwept = acos((dot(aPosVector1, aPosVector2)/(vPosVectorMagntd1*vPosVectorMagntd2))) ;
    
    else
        vThetaSwept = 2*pi - acos((dot(aPosVector1, aPosVector2)/(vPosVectorMagntd1*vPosVectorMagntd2))) ;
    
    end

% Time of Flight (TOF) computation:
vTOF = vJulianDay2 - vJulianDay1 ;

%% Heliocentric transfer trajectory computations:
% Determining x-fer trajectory orbital parameters, plus the initial and final velocity vectors of the spacecraft.
% Lambert solver utilized below (based on the Prussing/Conway trigonometric derivation and utilizing the MATLAB fzero functionality:

  vTOF = vTOF*24*60*60;

  [aSCvector1, aSCvector2] = LambertsOrbitDetermine(aPosVector1, aPosVector2, vThetaSwept, vTOF, vMuSun) ;

%   [aSCvector1, aSCvector2] = LambertsSolver_UniversalVarBisect(aPosVector1, aPosVector2, vThetaSwept, vTOF, vMuSun);
    
    vSMA = (2/norm(aPosVector1) - (norm(aSCvector1)^2)/vMuSun)^-1 ; 

    vEnergy = -vMuSun/(2*vSMA) ;

[e, vIncl, vTrueAnom1]       = state_elements(vMuSun, aPosVector1, aSCvector1) ;

[~, ~, vTrueAnom2]          = state_elements(vMuSun, aPosVector2, aSCvector2) ;

aV1scInfVector            = aSCvector1 - aVelVector1 ;
aV2scInfVector            = aSCvector2 - aVelVector2 ;  % here we compute the V_inf (pre_flyby)

end