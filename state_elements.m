%{
    ##########################__79_Character_Header__###########################
    PROCEDURE NAME:    state_elements.m

    PURPOSE:    This m-file computes the orbital elements of a spacecraft utilizing only
    its position/velocity vectors required to be synched at a common epoch. 
    (Common epochs are assumed within this script).

	External References:
 	Name                                                        Purpose
	----                                                        -------
    n/a                                                         n/a

    INPUTS:
		1.  (1) position vector          -
        2.  (1) velocity vector          -

	OUTPUTS:
		1. 	eccentricity                         (e)              -
        2.  inclination (rad)                    (vIncl)          -
        3.  rt ascention of Ascending Node(rad)  (vRAAN)          -
        4.  argument of periapse (rad)           (vArgPer)        -
        5.  true anomaly (rad)                   (vTrueAnom)      -
        6.  true longitude (rad)                 (vTrueLong)      -

  ##########################__79_Character_Header__###########################
   External References:
   NAME								PURPOSE
   ----                             -------  (none currently)
 
   Development History:
   NAME				DATE			DESCRIPTION OF CHANGE

   Reagoso, J.D.  	2010/12/15      Script initialized  

                    2013/2/22       Updated with vRAAN, vArgPer, vTrueLong
 ---------------------------------------------------------------------------------------------------
%}

%% Main Body 
function[e, vIncl, vTrueAnom] = state_elements(vMu, aPosVector, aVelVector)

kVector = [0;0;1];          % kHat unit vector

vR = norm(aPosVector);      % range Magnitude

vV = norm(aVelVector);      % velocity magnitude

aHvector    = cross(aPosVector, aVelVector) ;       % Spec angular momentum vector

vH          = norm(aHvector);                       % Spec angular momentum magnitude

vIncl       = acos(aHvector(3)/norm(aHvector)) ;    % inclination (rad)

aNodeVector = cross(kVector, aHvector);             % node vector

eVector     = (1/vMu)*(( vV^2 - vMu/vR)*aPosVector - dot(aPosVector, aVelVector)*aVelVector); % eccen vector                                                                      );

e           = norm(eVector) ;                       % eccentricity (unitless)

%% Rigth Ascension of the Ascending Node(RAAN) (rad):

if (aNodeVector(2) < 0.00)
    vRAAN = 2*pi - acos(aNodeVector(1)/norm(aNodeVector));    

else
    vRAAN = acos(aNodeVector(1)/norm(aNodeVector));

end

%% Argument of Periapse (rad):

if (eVector(3) < 0.0)
    vArgPer = 2*pi - acos(dot(aNodeVector, eVector)/(norm(aNodeVector)*norm(eVector))); 
    
else
    vArgPer = acos(dot(aNodeVector, eVector)/(norm(aNodeVector)*norm(eVector))); 
    
end

%% True Anomaly (rad)
if (dot(aPosVector, aVelVector) < 0.0)
    vTrueAnom = 2*pi - acos(dot(eVector, aPosVector)/(e*vR)) ;

else
    vTrueAnom = acos(dot(eVector, aPosVector)/(e*vR)) ;
    
end

%% True Longitude (rad)

vTrueLong = vRAAN + vArgPer + vTrueAnom;

end