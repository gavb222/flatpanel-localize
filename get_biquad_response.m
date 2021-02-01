function [mag_resp,freqs] = get_biquad_response(coeffs,fs)

    % returns the magnitude repsonse at the frequencies contained in
    % 'freqs' of the biquad filter given by the input coefficients and
    % sampling frequency
    %
    % 10 Hz resolution is enforced, and frequency limits are contained
    % between [20Hz 1000Hz] (modifyable by changing min_frequency and
    % max frequency)
    min_frequency = 20;
    max_frequency = 1000;


    % make biquad filter object from coefficientw
    biquadfilter = dsp.BiquadFilter(coeffs,1);

    % get magnitude and frequency response arrays
    [resp,freqs] = freqz(biquadfilter,round(fs/10),'whole',fs);
    
    % take the absolute value of the frequency response to get the
    % magnitude resposne in linear scale
    mag_resp = (abs(resp));
    
    % shorten to appropriate length:
    [~ , minidx] = min(abs(min_frequency - freqs));
    [~ , maxidx] = min(abs(max_frequency - freqs));
    
    % format output
    freqs = freqs(minidx:maxidx)';
    mag_resp = mag_resp(minidx:maxidx)';
    
    
    
end

