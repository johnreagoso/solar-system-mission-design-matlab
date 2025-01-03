%This m-file will compute the specific planetry ephemeris data for a
%specific Julian Date

%function[a,e,i,omega,w,w_bar,L,M] = ephemeris_compute(Julian_day_mumber)

function[vSMAmean, vPlanetRadius, vMinFlybyDist, vMuPlanetConst] = planet_specs(planet_arrival)

switch planet_arrival
    
    case 'Mercury'
        vSMAmean  = .38709927;
        vPlanetMassPercntg = .0553;
        vRadiusPercntg = .382;
        vMinRpDist = 0.0;
        
    case 'Venus'
        vSMAmean = .72333566; 
        vPlanetMassPercntg = .8149;
        vRadiusPercntg = .949;
        vMinRpDist = 300.00;
        
    case 'Earth'
        vSMAmean  = 1.00000261;
        vPlanetMassPercntg = 1.0;
        vRadiusPercntg = 1.0;
        vMinRpDist = 100.0;    
        
    case 'Mars'
        vSMAmean = 1.52371034;
        vPlanetMassPercntg = .1074;
        vRadiusPercntg = .532;
        vMinRpDist = 100.0;
        
    case 'Jupiter'
        vSMAmean = 5.20288700;
        vPlanetMassPercntg = 317.938;
        vRadiusPercntg = 11.209;
        vMinRpDist = 100.0;
        
    case 'Saturn'
        vSMAmean = 9.53667594;
        vPlanetMassPercntg = 95.181;
        vRadiusPercntg = 9.49;
        vMinRpDist = 100.0;        
        
    case 'Uranus'
        vSMAmean = 19.18916464;
        vPlanetMassPercntg = 14.531;
        vRadiusPercntg = 4.007;
        vMinRpDist = 100.0;
        
    case 'Neptune'
        vSMAmean = 30.06992276;
        vPlanetMassPercntg = 17.135;
        vRadiusPercntg = 3.83;
        vMinRpDist = 100.0;
        
    case 'Pluto'
        vSMAmean = 39.48211675;
        vPlanetMassPercntg = .0022;
        vRadiusPercntg = 0.18;
        vMinRpDist = 0.0;
        
    otherwise
        disp('Please enter a Solar System planet only');
        
end

AU = 1.495978e8;
vSMAmean = vSMAmean * AU;

vEarthRadius = 6378.12; %in km
vEarthMuConst = 3.9865e5; % km^3/s^2

vPlanetRadius  = vRadiusPercntg * vEarthRadius;
vMuPlanetConst = vPlanetMassPercntg * vEarthMuConst;
vMinFlybyDist  = vPlanetRadius + vMinRpDist;

end
