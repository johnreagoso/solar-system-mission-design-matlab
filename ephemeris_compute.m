%{
=====================================================================
  These data are to be used as described in the related document
  titled "Keplerian Elements for Approximate Positions of the
  Major Planets" by E.M. Standish (JPL/Caltech) available from
  the JPL Solar System Dynamics web site (http://ssd.jpl.nasa.gov/).
=====================================================================

Keplerian elements and their rates, with respect to the mean ecliptic
and equinox of J2000, valid for the time-interval 1800 AD - 2050 AD
%}


%This m-file will compute the specific planetry ephemeris data for a
%specific Julian Date

function[a_mean, planet_mu, radius_planet, a, e, i, omega, argmt_periapsis, M] = ephemeris_compute(planet, Julian_cent_epoch)
                                             
earth_radius = 6378.12; %km
earth_mu = 3.986e5;     %km^3/ sec^2

switch planet
    
    case 'Earth'
        a_mean          = 1.0;
        planet_mu       = 1* earth_mu;
        radius_planet   = 1* earth_radius;
        a0 = 1.00000261; 
        a1 = 0.00000562;
        e0 = 0.01671123;  
        e1 = -0.00004392;
        i0 = -0.00001531; 
        i1 = -0.01294668;
        L0 = 100.46457166; 
        L1 = 35999.37344981;
        w_bar0 = 102.93768193; 
        w_bar1 = 0.32327364;
        omega0 = 0.00; 
        omega1 = 0.00;
        
    case 'Mercury'
        
        a_mean          = .38709927;
        planet_mu       = .0553* earth_mu;
        radius_planet   = .382* earth_radius;
        a0 = .38709927;
        a1 = .00000037;
        e0 = .20563593;
        e1 = .00001906;
        i0 =  7.00497902;
        i1 = -.00594749;
        L0 =  252.25032350;
        L1 = 149472.67411175;
        w_bar0 = 77.45779628; 
        w_bar1 = .16047689;
        omega0 = 48.33076593; 
        omega1 = -.12534081;
        
    case 'Venus'
        
        a_mean          = .72333566; 
        planet_mu       = .8149* earth_mu;
        radius_planet   = .949* earth_radius;
        a0 = .72333566;
        a1 = .00000390;
        e0 = .00677672;  
        e1 = -.00004107;
        i0 = 3.39467605; 
        i1 = -.00078890;
        L0 = 181.97909950; 
        L1 = 58517.81538729;
        w_bar0 = 131.60246718; 
        w_bar1 = .00268329;
        omega0 = 76.67984255; 
        omega1 = -.27769418;
        
    case 'Mars'
        
        a_mean          = 1.52371034;
        planet_mu       = .1074* earth_mu;
        radius_planet   = .532* earth_radius;
        a0 = 1.52371034; 
        a1 = .00001847;
        e0 = .09339410;
        e1 = .00007882;
        i0 = 1.84969142;
        i1 = -.00813131;
        L0 = -4.553343205;
        L1 = 19140.30268499;
        w_bar0 = -23.94362959; 
        w_bar1 = .44441088;
        omega0 = 49.55953891;
        omega1 = -.29257343;
        
    case 'Jupiter'
        
        a_mean          = 5.20288700;
        planet_mu       = 317.938* earth_mu;
        radius_planet   = 11.209* earth_radius;
        a0 = 5.20288700;
        a1 = -.00011607;
        e0 =  .04838624;
        e1 = -.00013253;
        i0 = 1.30439695;
        i1 = -.00183714;
        L0 = 34.39644051;
        L1 = 3034.74612775;
        w_bar0 = 14.72847983;
        w_bar1 = .21252668;
        omega0 = 100.47390909;
        omega1 = .20469106;
    
    case 'Saturn'
        
        a_mean          = 9.53667594;
        planet_mu       = 95.181 * earth_mu;
        radius_planet   = 9.49 * earth_radius;
        a0 = 9.53667594;
        a1 = -.00125060;
        e0 = .05386179;
        e1 = -.00050991;
        i0 = 2.48599187;
        i1 = .00193609;
        L0 = 49.95424423;
        L1 = 1222.49362201;
        w_bar0 = 92.59887831;
        w_bar1 = -.41897216;
        omega0 = 113.66242448;
        omega1 = -.28867794;
    
    case 'Uranus'
        
        a_mean          = 19.18916464;
        planet_mu       = 14.531* earth_mu;
        radius_planet   = 4.007* earth_radius;
        a0 =  19.18916464;
        a1 =  -.00196176;
        e0 = .04725744;
        e1 = -.00004397;
        i0 = .77263783;
        i1 = -.00242939;
        L0 = 313.23810451;
        L1 = 428.48202785;
        w_bar0 = 170.95427630;
        w_bar1 = .40805281;
        omega0 = 74.01692503;
        omega1 = .04240589;
    
    case 'Neptune'
        
        a_mean          = 30.06992276;
        planet_mu       = 17.135* earth_mu;
        radius_planet   = 3.83* earth_radius;
        a0 = 30.06992276; 
        a1 = .00026291;
        e0 = .00859048;
        e1 = .00005105;
        i0 = 1.77004347;
        i1 = .00035372;
        L0 = -55.12002969;
        L1 = 218.45945325;
        w_bar0 = 44.96476227;
        w_bar1 = -.32241464;
        omega0 = 131.78422574;
        omega1 = -.00508664;
    
    case 'Pluto'
        a_mean          = 39.5294;
        planet_mu       = .0022* earth_mu;
        radius_planet   = 0.18* earth_radius;
        a0 = 39.48211675; 
        a1 = -.00031596;
        e0 = .24882730;
        e1 = .00005170;
        i0 = 17.14001206;
        i1 = .00004818;
        L0 = 238.92903833;
        L1 = 145.20780515;
        w_bar0 = 224.06891629;
        w_bar1 = -.04062942;
        omega0 = 110.30393684;
        omega1 = -.01183482;
        
    otherwise
        disp('Please re-enter a Solar System planet with the correct spelling and/or the first letter capitalized only');
end
%--------------------------------------------------------------------------
AU = 1.495978e8;                        

a_mean = a_mean * AU;

a = a0 + a1* (Julian_cent_epoch);       

a = a* AU;

e = e0 + e1* (Julian_cent_epoch);

%JPL planetary ephemeris data are provided in degrees/degrees per century
%convert to radians as well:

i = (i0 + i1* (Julian_cent_epoch)) *pi/180 ;

omega = (omega0 + omega1* (Julian_cent_epoch)) *pi/180 ;

            [omega] = angle_minimize(omega);

w_bar = (w_bar0 + w_bar1* (Julian_cent_epoch)) *pi/180 ;

            [w_bar] = angle_minimize(w_bar) ;

argmt_periapsis = (w_bar - omega) ;

            [argmt_periapsis] = angle_minimize(argmt_periapsis);

L = (L0 + L1* (Julian_cent_epoch))*pi/180 ;
    
            [L] = angle_minimize(L);% angle_minimize function not required here, but left anyway for clarity for 
            %subsequent code review. 

M = L - w_bar;

            [M] = angle_minimize(M); 

end
