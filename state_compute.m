%   1. This m-file will compute a planet's/spacecraft's position and velocity vector
%   utilizing provided orbital elements. 

%   2. ****************All values are for heliocentric orbits. If this code it to be utilized
%   for earth-centered trajectories, the mu below must me
%   modified*********************

function [pos_vector_magntd, pos_vector, vel_vector_magntd, vel_vector, h_spc_ang_magntd, flight_path_angle] = state_compute(a, e, i, w, omega, nu)

mu  = 1.327124e11 ; %obviously, for heliocentric orbits only.
h = (a* (1 - e^2)* mu)^ 0.5 ;
% first we compute the s/c's or planet's pos_vector, vel_vector in
% perifocal coordinates, then rotate to ecliptic. 

r_perifocal  = (h^2/mu)* (1/(1 + e* cos(nu)))* [cos(nu); sin(nu); 0] ;
v_perifocal  = (mu/h)* [-1* sin(nu); e + cos(nu); 0] ;

R_w     = [cos(w) sin(w) 0; -1* sin(w) cos(w) 0; 0  0  1];

R_i     = [1  0  0; 0 cos(i) sin(i); 0 -1* sin(i) cos(i)] ;

R_omega = [cos(omega) sin(omega) 0; -1* sin(omega) cos(omega) 0; 0  0  1] ;

Rot_ecliptic = transpose(R_w* R_i* R_omega) ; %rotates the perifocal pos/vel vectors into the 
%3-D perifocal frame of coordinates, then takes the transpose. 

%Rot_ecliptic = transpose(Rot_ecliptic) ; 

pos_vector = Rot_ecliptic* r_perifocal ;
pos_vector_magntd = norm(pos_vector) ; 

vel_vector = Rot_ecliptic* v_perifocal ;
vel_vector_magntd = norm(vel_vector) ;

h_spc_ang_magntd = norm(Fcross(pos_vector, vel_vector)) ;

flight_path_angle = acos(h_spc_ang_magntd/(pos_vector_magntd* vel_vector_magntd)) ;

return 

 function[aVectorCrossOut] = Fcross(aVector1, aVector2)

     aVectorCrossOut(1,1) = aVector1(2)*aVector2(3)-aVector1(3)*aVector2(2);
     
     aVectorCrossOut(2,1) = aVector1(3)*aVector2(1)-aVector1(1)*aVector2(3); 
     
     aVectorCrossOut(3,1) = aVector1(1)*aVector2(2)-aVector1(2)*aVector2(1);

return 



 
 
 
 
 
 
