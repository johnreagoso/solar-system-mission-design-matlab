%% planetary (arrival or departure) heliocentric orbital parameters
function[pos_vec_planet, pos_vec_magntd_planet, vel_vec_planet, vel_magntd_planet, true_long_planet] = planet_orbit_parameters(sPlanet, Julian_cent_epoch)

[a_planet_mean, u_planet, r_planet, a_planet, e_planet, i_planet, omega_planet, argmt_periapsis_planet, M_planet] = ephemeris_compute(sPlanet,...
                Julian_cent_epoch);

%[E_planet, nu_planet] = kepler_solver_vector(e_planet, M_planet); %this function returns the E, nu as requested, utilizing the e, M provided earlier

[E_planet, nu_planet] = kepler_solver(e_planet, M_planet); %this function returns the E, nu as requested, utilizing the e, M provided earlier


true_long_planet = nu_planet + omega_planet + argmt_periapsis_planet; %computes the true longitude of the departure planet at spacecraft departure. 

[pos_vec_magntd_planet, pos_vec_planet, vel_magntd_planet, vel_vec_planet, h_magntd_planet, flt_angle_planet] = state_compute(a_planet, e_planet, i_planet,...
                        argmt_periapsis_planet, omega_planet, nu_planet) ;
                    
end