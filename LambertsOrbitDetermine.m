%{
    ##########################__79_Character_Header__###########################
    PROCEDURE NAME: LambertsOrbitDetermine.m

    PURPOSE:    The following code computes via (2) potential separate Lamberts
                algorithms the transfer trajectory velocity vectors at r1 and r2. 
                
                For elliptical trajectories, the code first will attemp to utilize the
                Matlab Fzero functionality, which solves for the root (SMA) of Lambert's equations
                based on the trigonometric Prussing/Conway Lamberts algorithm.

                If this fails, the module (via a Try/Catch statement) calls 
                a Lamberts Bisection routine(which uses a Universal Variables approach.)
 
                For hyperbolic trajectories, the code utlizes the Curtis
                Universal Variables algorithm. 

NOTE: Bisection algorithm is used as a back-up 

	External References:
 	Name                                Purpose
	----                                -------
    Prussing/Conway Text                Reference
    'Orbital Mechanics' (1993)

    Curtis 'Orbital Mechanics           Reference
    for Engineering Students'(2005)

    INPUTS:
		1.  Position Vector at t1       aPosVector1
        2.  Position Vector at t2       aPosVector2
        3.  Angle swept                 vThetaSwept
        4.  Mu                          vMu_

	OUTPUTS:
		1. 	(2) Velocity Vectors            -

  ##########################__79_Character_Header__###########################
   External Code Utilized:
   NAME                                 PURPOSE
 
   f_g_series_compute_Lambert           Compute orbital elements

   DEVELOPMENTAL HISTORY:
   NAME				DATE			DESCRIPTION OF CHANGE
   Reagoso, J.D.  	2012/1/9        Script initialized  

                    2013/2/15       Curtis Universal Variables
                                    functionality added

                    2013/2/27       Try/Catch statement added (with Bisection
                                    Lamberts algorithm provided)
 ---------------------------------------------------------------------------------------------------
%}


function [aVelVector1, aVelVector2] = LambertsOrbitDetermine(aPosVector1, aPosVector2, vThetaSwept, TOF, vMu_)

%% global parameters
global r1  r2  d  s  theta  vMu   pi   A   
global parabolic_tflight  vTOF  tflight_min
global a_min a_max

vMu = vMu_*1 ;

pi  = 3.141592653589 ;
%% Basic user-defined mission geometry computations:

%Provided will be the time of flight (TOF), distances at two distinct
%times (t1, t2) and the estimated angle (theta) between the two sightings.

r1 = norm (aPosVector1) ;
r2 = norm (aPosVector2) ;

theta = vThetaSwept*1;

d = sqrt(r1^2 + r2^2 - 2*r1*r2* cos(theta)) ;
s = 0.50*(r1 + r2 + d) ;

%tflight = TOF*8.640e4 ; %(ToF x 24*60*60) % this assumes days input, then converts to seconds
vTOF = TOF;

%% minimum energy orbit parameters/tflight_min calculations
a_min = (r1 + r2 + d)/4;

alpha_min = pi ;

beta_min= 2*asin(sqrt((s - d)/s)) ;

    if (theta >= pi)   
    
        beta_min= -2*asin(sqrt((s - d)/s)) ;  

    end
%-------------------------------------------------------------------------
%% Orbit-type determination (elliptical/hyperbolic/parabolic:

parabolic_tflight = ((r1 + r2 + d)^1.5)/(6* sqrt(vMu))  -  ((r1 + r2 - d)^1.5)/(6* sqrt(vMu))  ;

a_max = a_min* 100 ;

tflight_min = ((a_min^1.5)/sqrt(vMu))* ((alpha_min - sin(alpha_min)) - (beta_min - sin(beta_min)))  ;

%% Here, the code will either return a fail-message to the main script and transition to the universal variables methodology 
%  or continue with computing the sma using the Matlab fzero functionality.  
 
if parabolic_tflight < vTOF      % this is an elliptical orbit or hyperbolic orbit
    
    try     % Try/Catch error statement here:

   %   ====================================================================================================================
   %   First attempt will be to utilizing the Matlab Fzero functionality to find the root (SMA) for this specific geometry. 
   %   This should work for all cases of ellipitcal orbits. Fzero does sometimes fail near parabolic orbits (e ~ 0.99)
   %   ====================================================================================================================
        
        vSMA = fzero(@lamberts_e, [a_min, a_max], 1.0e-9) ;
        
        vP = 4*vSMA*((s - r1)*(s - r2))/d^2  *sin(0.5*(alpha_e(vSMA) + beta_e(vSMA)))^2 ;
        
        v_e = sqrt(1 - vP/vSMA) ;
        
        [aVelVector1, aVelVector2] = f_g_series_compute_Lambert(vMu, vSMA, v_e, vThetaSwept, aPosVector1, aPosVector2) ;
        
    catch sErrorMsg
     %   disp(sErrorMsg);

%   ====================================================================================================================
%   Assuming that the Matlab Fzero root finder method did not work,the script will transition to utilizng the a universal
%   variable approach, incorporating a Newton's method to find the root(z). 
%   ====================================================================================================================
   
                        A = sin(vThetaSwept)* sqrt((r1* r2)/(1 - cos(vThetaSwept))) ;
                        vZ = 0.010;

                        vFz  = F(vZ, vTOF);

                %%  Here we calculate an initial guess for (z):
                %   (identify where z, vFz crosses the x-axis)...
                      ii = 0;
                            while vFz < 0

                                vZ = vZ + .20 ;
                                vFz = F(vZ,vTOF);

                                ii = ii + 1;

                            end

                %%  Main Body: Newton Solver to determine (z):
                    vTol = 1.00e-8 ;     vNmax = 500 ;
                    vRatio = 1;

                    jj = 0 ;

                    while (abs(vRatio) > vTol) && (jj<= vNmax)

                        jj    = jj + 1 ;
                        vRatio = F(vZ,vTOF)/dFdz(vZ) ;

                        vZ     = vZ - vRatio ;
                    end

                %% Error Warning:
                     if (jj == vNmax)
                        fprintf('\n\n **Number of iterations exceeds')
                        fprintf(' %g \n\n', vNmax)
                     end

                %% Velocity Vector Computation (f, g Series)

                    f = 1 - y(vZ)/r1 ;

                    g = A* sqrt(y(vZ)/vMu) ;

                    gdot = 1 - y(vZ)/ r2 ;

                    aVelVector1 = aPosVector2/g - f*aPosVector1/g ;
                    aVelVector2 = (gdot* aPosVector2 - aPosVector1)/g ;
                    
    end
    
else
   
    %% Hyperbolic orbit
    %  Orbit is hyperbolic (below uses the universal variables algorithm):

        A = sin(vThetaSwept)* sqrt((r1* r2)/(1 - cos(vThetaSwept))) ;
        
        vZ = -100.00;  % better pick for a hyperbolic orbit.. 
 
        vFz  = F(vZ, vTOF);

%%  Here we calculate an initial guess for (z):
%   (identify where z, vFz crosses the x-axis)...
      ii = 0;
            while vFz < 0
    
                vZ = vZ + .20 ;
                vFz = F(vZ,vTOF);

                ii = ii + 1;

            end
            
%%  Main Body: Newton Solver to determine (z):
    vTol = 1.00e-8 ;     vNmax = 500 ;
    vRatio = 1;

    jj = 0 ;

    while (abs(vRatio) > vTol) && (jj<= vNmax)
    
        jj    = jj + 1 ;
        vRatio = F(vZ,vTOF)/dFdz(vZ) ;
    
        vZ     = vZ - vRatio ;
    end

%% Error Warning:
     if (jj == vNmax)
        fprintf('\n\n **Number of iterations exceeds')
        fprintf(' %g \n\n', vNmax)
     end

%% Hyperbolic Velocity Vector Computation (f, g Series)

    f = 1 - y(vZ)/r1 ;

    g = A* sqrt(y(vZ)/vMu) ;

    gdot = 1 - y(vZ)/ r2 ;

    aVelVector1 = aPosVector2/g - f*aPosVector1/g ;
    aVelVector2 = (gdot* aPosVector2 - aPosVector1)/g ;

end
    
return

%% Subroutines used by main body for elliptical Lamberts Solver (Prussing/Conway) and Hyperbolic (Univ Variables)

function[dummy] = lamberts_e(a)
%% sub-routine that returns the elliptical lamberts equation for the Matlab fzero functionality 
global vTOF  vMu 
      
    dummy = a^1.5/sqrt(vMu)* (alpha_e(a) - sin(alpha_e(a))- (beta_e(a) - sin(beta_e(a)))) - vTOF;

return

% function[dummy] = multi_lamberts_e(a)
% %% sub-routine that returns the elliptical lamberts equation for the Matlab fzero functionality 
% global vTOF  vMu   vN      
%     dummy = a^1.5/sqrt(vMu)* (2*pi*vN +   alpha_e(a) - sin(alpha_e(a))- (beta_e(a) - sin(beta_e(a)))) - vTOF;
% 
% return


function[dummy] = lamberts_h(a)
%% sub-routine that returns the hyperbolic lamberts equation for the Matlab fzero functionality 

    global vTOF vMu 
    
    dummy = -a^1.5/sqrt(vMu)* (sinh(alpha_h(a)) - alpha_h(a) - (sinh(beta_h) - beta_h(a))) - vTOF;

return 



function dummy = alpha_e(a)
%% sub-routine that returns alpha(a) elliptical components for the Lamberts iteration
    global r1   r2   d   pi
    
    global vTOF  tflight_min
    
        if vTOF > tflight_min 

             dummy = 2*pi - 2* asin( sqrt((r1 + r2 + d)/(4* a))) ;   %here tflight < tflight_min
         
        else 
        
             dummy = 2* asin( sqrt((r1 + r2 + d)/(4* a))) ;          %here tflight > tflight_min
            
        end

return 


function dummy = beta_e(a)
%% sub-routine that returns elliptical beta(a) components for the Lamberts iteration
    global r1   r2   d  theta   pi
   
        if theta <= pi  
           
                dummy = 2* asin( sqrt((r1 + r2 - d)/(4* a))) ; %(this is for when 0 <= theta < pi)
            
        else            

                dummy = -2* asin( sqrt((r1 + r2 - d)/(4* a))) ; %(this is for when pi<= theta < 2pi)

        end
           
return 

function dummy = alpha_h(a)
%% sub-routine that returns hyperbolic alpha(a) components for the Lamberts iteration depending on the input geometry:
    global r1   r2   d  
    
    global tflight  tflight_min 
   
        if tflight > tflight_min  
           
                dummy =  2*asinh( sqrt((r1 + r2 + d)/(4* -a))) ; 
                
        else            
                dummy =  2* asinh( sqrt((r1 + r2 + d)/(4* -a))) ;

        end
           
return 


function dummy = beta_h(a)
%% sub-routine that returns hyberbolic beta(a) components for the Lamberts iteration
    global r1   r2  d  theta
   
        if theta <= pi  %(this is for when 0 <= theta < pi)
           
                dummy =  2* asinh( sqrt((r1 + r2 - d)/(4* -a))) ;
            
        else            %(this is for when pi<= theta < 2pi)

                dummy =  -2* asinh( sqrt((r1 + r2 - d)/(4* -a))) ;

        end
           
return

function dum = y(z)
       
    global r1  r2   A
    
    dum = r1  + r2 + A* (z*S(z) - 1)/ sqrt(C(z)) ;      
    
return
    
function dum = F(z,t)
            
    global vMu A
    dum = (y(z)/C(z))^1.5 *S(z) + A*sqrt(y(z)) - sqrt(vMu)*t ;

return
            
function dum = dFdz(z)
global A 
        
        if z == 0
            dum = sqrt(2)/40*y(0)^1.5 + (A/8)*(sqrt(y(0))...
                + A* sqrt(1/2/y(0))) ;
        
        else 
            dum = (y(z)/C(z))^1.5*(1/2/z*(C(z) - 3*S(z)/2/C(z))...
                + 3*S(z)^2/4/C(z))...
                + A/8*(3*S(z)/C(z)*sqrt(y(z))...
                + A*sqrt(C(z)/y(z)));
        end
      
return

function dum = S(z)

   if z > 0
        dum  = (sqrt(z) - sin(sqrt(z)))/(sqrt(z))^3;
   elseif z < 0
        dum = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z))^3;
   else    
        dum = 1/6 ;
   end
   
return

function dum = C(z)

if z > 0
    dum = (1 - cos(sqrt(z)))/z;
elseif z < 0
    dum = (cosh(sqrt(-z)) - 1)/(-z);
else 
    dum = 0.50;
    
end

return
        



