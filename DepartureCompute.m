%% Departure Hyperbolic Trajectory Compute

% Purpose- this script provides for required parking orbit placement in
% order to conduct a planetary departure maneveuver. All the user has to
% provide is the outbound cartesian hyperbolic asymptote velocity
% coordinates (this data typically comes from a patched conic
% solution) and this code will compute the DEC, RA and C3 of the
% departure asymptote, as well as the parking orbit TA, RAAN for maneuver
% start position.

% Algorithm comes from: https://www.thebackshed.org/forum/uploads/cdeagle/2017-07-20_185149_hyper1_matlab.pdf

% vel_inf_vector = [1.59174874145; 1.99711032827; -1.23556580364];
% vel_inf_vector = [-0.0195357052812; 2.88081260692; 0.0812706402012];
vel_inf_vector = [2.41033888468; 1.37646618474; 0.591410094825];


Vinf_mag       = norm(vel_inf_vector);

alt  = 380.0;       % circular parking orbit altitude (km)
incl = 27.0000;     % parking orbit inclination (km)

%% Departure Hyperbola Parameters:

C3 = Vinf_mag^2 ;
DECL = asind(vel_inf_vector(3)/norm(vel_inf_vector));
RLA  = atand(vel_inf_vector(2)/vel_inf_vector(1));

%% Compute true anomaly and RAAN at parking orbit for departure"
denom = 1 + (Mars_GenPhysCons.RE_EQ + alt)*C3 /Mars_GenPhysCons.GM_KM;

n = asind(1/denom);

beta = 90.0-DECL;
sma = ((Mars_GenPhysCons.RE_EQ + alt) + (Mars_GenPhysCons.RE_EQ + alt))/2.0 ;

trueAnom1 = acosd(cosd(beta)/sind(incl)) - n;
trueAnom2 = -acosd(cosd(beta)/sind(incl)) - n;

raan1 = 180 + RLA + asind(cotd(beta)/tand(incl));
raan2 = 360 + RLA - asind(cotd(beta)/tand(incl));

if raan1 > 360; raan1 = 360.0 - raan1; end
if raan2 > 360; raan2 = 360.0 - raan2; end



