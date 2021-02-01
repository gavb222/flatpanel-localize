# flatpanel-localize
Clamped_Panel - main class, includes methods for viewing each drivers scan individually or summed to see the total response at each frequency

Single_Driver_Scan - helper class, saves the displacement matrix (u) of each driver's scan separately. Do not call this class on it's own, it's meant to be contained within Clamped_Panel

get_biquad_response - helper function for getting the magnitude response of a biquad filter based on filter coefficients. This function also enforces that the resposne is limited to the frequencies of interest

test_Clamped_Panel - shows an example of instantiating a Clamped_Panel. This may be the most helpful script to peruse, as when your model spits out filter coefficients, this script shows the next steps
