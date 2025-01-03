 function[sDate_string, Greg_month, Greg_day, Greg_year] = JulianDay_to_Greg_cal(Julian_day)

% utilizes algorithm from http://quasar.as.utexas.edu/BillInfo/JulianDatesG.html
% assumes that the Julian day hr/min-sec is 0000

%Julian_day = 2443589.9900 ;

z = Julian_day ; %+ 0.50 ;

w_top = z - 1867216.25 ;

w = fix(w_top/ 36524.25) ; %,'round') ;

x = fix(w/4) ;

a = z + 1 + w - x ;

b = a + 1524 ;

c_top = b - 122.1 ;

c = fix(c_top/ 365.25) ;

d = fix(365.25*c) ;

e_top = b - d ;

e = fix(e_top/ 30.6001) ;

f = fix(30.6001*e) ;

Greg_day = b-d-f ;

month1 = e-1 ; 
month2 = e - 13 ;

if (month1 < 1 || month1 > 12) 
        month = month2 ;
    
else
        month = month1 ;
    
end    
    

if (month == 1 || month == 2)
    Greg_year = c - 4715 ;
    
else
    Greg_year = c - 4716 ;
    
end


if (month == 1)

    Greg_month = 'Jan ' ;

elseif (month == 2)
    
    Greg_month = 'Feb ' ;
    
elseif (month == 3)
    
    Greg_month = 'Mar ' ;
    
elseif (month == 4)
    
    Greg_month = 'Apr ' ;
    
elseif (month == 5)
    
    Greg_month = 'May ' ;
    
elseif (month == 6)
    
    Greg_month = 'Jun ' ;
    
elseif (month == 7)
    
    Greg_month = 'July ' ;

elseif (month == 8)
    
    Greg_month = 'Aug ' ;
    
elseif (month == 9)
    
    Greg_month = 'Sep ' ;
    
elseif (month == 10)
    
    Greg_month = 'Oct ' ;

elseif (month == 11)
    
    Greg_month = 'Nov ' ;
    
elseif (month == 12)
    
    Greg_month = 'Dec ' ;
   
end

    sGreg_day_string = num2str(Greg_day + 0.50) ;
                   
    sGreg_year_string = num2str(Greg_year ) ;
                                 
    sDate_string = [Greg_month,  sGreg_day_string,' ',  sGreg_year_string ] ;

end


