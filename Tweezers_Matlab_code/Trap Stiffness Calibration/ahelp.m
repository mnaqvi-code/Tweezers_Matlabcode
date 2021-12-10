% ***********************Help file for program tweezercalib 2.0****************************
% 
% 1) Start Matlab and make current directory the same as the path where you saved the program folder e.g.
% type  cd('c:\.....\perfectfit') in the Matlab window.
% 
% 2) Type 'start_fit' in the Matlab window: >>start_fit
% 
% 3 An interface will appear, click the buttom 'load_time_series' 
%  
% 4) Tha data file you load should time series for x and/or y and possibly the z signal recorded from a photodiode. The
% time series should be ordered as coloumns: x is first coloumn y is second and z is third coloumn.
% 
% 5) Fill in the relevant numbers and click 'view data' 
% 
% 7) If you trap with IR light and your data are recorded by a photodiode which is not manufactured for IR-light, choose 
%     'position detector treated as a filter with two parameters'
%     
% 8) If you have lowpass filters in the detector system choose the relevant numbers in 'filter corrections'.
% 
% 6) Finally a fit can be done by clicking 'fit_powerspectrum'
% 
% PS-If the number of data points per bin is smaller than 200 the algorithm changes to a slower more precise 
% fitting algorithm
