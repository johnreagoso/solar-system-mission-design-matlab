function [aState] = Function_3DOF_SolarSystem(elapsed_sec, aStateInput)

global Julian_base_epoch;

%  State vector         --> [X   Y   Z   Vx  Vy  Vz]
% (d/dt) State vector   --> [Vx  Vy  Vz  AccelX  AccelY  AccelZ] 

% Planetary Gravitation Constants (km^3/sec^2):
EarthMu   = 3.98605e5;
SunMu     = 1.32712e11; 
MarsMu    = 4.28284e4;
JupiterMu = 1.26713e8;

Rvector = [aStateInput(1); aStateInput(2); aStateInput(3)];     % km
Vvecfor = [aStateInput(4); aStateInput(5); aStateInput(6)];     % km/sec

dx  = Vvecfor(1);   dy  = Vvecfor(2);   dz  = Vvecfor(3);

%% Epoch Time determine:
JulianCentEpochArrival = ((Julian_base_epoch + elapsed_sec - 2451545.0)/ 36525);

%% Planet Position Vectors
% earth:
[earthPosVector, ~, ~, ~] = planet_orbit_parameters('Earth', JulianCentEpochArrival) ;
% mars:
[marsPosVector, ~, ~, ~] = planet_orbit_parameters('Mars', JulianCentEpochArrival) ;
% jupiter:
[jupiterPosVector, ~, ~, ~] = planet_orbit_parameters('Jupiter', JulianCentEpochArrival) ;

%% Gravitational Acceleration Compute: 
% km/sec^2
ddx_sun = (-1* SunMu*Rvector(1)/(norm(Rvector))^3);       
ddy_sun = (-1* SunMu*Rvector(2)/(norm(Rvector))^3);      
ddz_sun = (-1* SunMu*Rvector(3)/(norm(Rvector))^3);

sc_earthvector = Rvector - earthPosVector; 
%disp('distance from Earth:')
%disp(norm(sc_earthvector));

ddx_earth = (-1* EarthMu*sc_earthvector(1)/(norm(sc_earthvector))^3);           
ddy_earth = (-1* EarthMu*sc_earthvector(2)/(norm(sc_earthvector))^3);     
ddz_earth = (-1* EarthMu*sc_earthvector(3)/(norm(sc_earthvector))^3);
%  
sc_marsvector = Rvector - marsPosVector;
ddx_mars = (-1* MarsMu*sc_marsvector(1)/(norm(sc_marsvector))^3);           
ddy_mars = (-1* MarsMu*sc_marsvector(2)/(norm(sc_marsvector))^3);     
ddz_mars = (-1* MarsMu*sc_marsvector(3)/(norm(sc_marsvector))^3);
% 
sc_jupitervector = Rvector - jupiterPosVector;
ddx_jupiter = (-1* JupiterMu*sc_jupitervector (1)/(norm(sc_jupitervector ))^3);           
ddy_jupiter = (-1* JupiterMu*sc_jupitervector (2)/(norm(sc_jupitervector ))^3);     
ddz_jupiter = (-1* JupiterMu*sc_jupitervector (3)/(norm(sc_jupitervector ))^3);

%% Final Acceleration Compute:
ddx = ddx_sun + ddx_earth + ddx_mars + ddx_jupiter;  
ddy = ddy_sun + ddy_earth + ddy_mars + ddy_jupiter;   
ddz = ddz_sun + ddz_earth + ddz_mars + ddz_jupiter;    

aState = [dx; dy; dz; ...
    ddx; ddy; ddz];

end

