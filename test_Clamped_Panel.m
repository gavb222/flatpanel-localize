
% By convention:
%
% - We will use a sampling frequency of 44100 Hz for the filters
%
% - We will look at the panel responses from [20 Hz 3000 Hz] with a
%       resolution of 10 Hz, unit of arrays is Hz
%
% - Magnitude responses are given in linear scale, not in dB scale, so that
%       they are compatible with the F_i(w) term in the governing equation
%
% - Driver locations (xi, yi) are given as percentages of Lx and Ly


% Example with a 4-driver panel:

% Panel information:
Lx = 0.3;
Ly = 0.5;
driver_locations = [0.25 0.25 ; 0.75 0.75; 0.25 0.75 ; 0.75 0.25];

% Example designing two filters: 
fs = 44100;
coeffs_driver1 = [1 1 1 1 0.5 0.5];
coeffs_driver2 = [1 0.7 0.7 1 0.5 0.5];
coeffs_driver3 = [1 1 0.7 1 0.5 0.5];
coeffs_driver4 = [1 0.7 1 1 0.5 0.5];



% call helper funtion get_biquad_response
[response1, frequencies] = get_biquad_response(coeffs_driver1,44100);
[response2, ~]           = get_biquad_response(coeffs_driver2,44100);
[response3, ~]           = get_biquad_response(coeffs_driver2,44100);
[response4, ~]           = get_biquad_response(coeffs_driver2,44100);

% package together:
%driver_responses = [response1 ; response2; response3 ; response4];

% ^Under normal circumstances, we would use the filter coefficients spit out
% from Gavin's model to filter the drivers on the panel. For testing
% purposes, lets just use an all-pass filter, or all ones, as shown below.
% A gain of 3 is applied. Drivers 1 and 4 are driven out of phase in order
% to drive the (2,1) mode for testing purposes 


driver_responses = ones(4,length(frequencies))*3;
driver_responses(1,:) = -3;
driver_responses(4,:) = -3;


% instantiate the panel:
panel_under_test = Clamped_Panel(driver_locations,driver_responses,frequencies,Lx,Ly);

% plot the response of a given frequency:
panel_under_test.view_total_scan(200,1);

