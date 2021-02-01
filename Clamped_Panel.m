classdef Clamped_Panel < handle
    
    %
    
    properties

        
        % Panel physical parameters
        Lx;     % Length of Panel in X dimension (meters)
        Ly;     % Length of Panel in Y dimension (meters)
        h;      % thickness of panel (meters)
        E;      % Young's modulus of the panel (pa)
        v;      % Poissons ratio of the panel (unitless)
        rho;    % Density of the panel (kg/meters^3)
        mass;   % Mass of the panel (kg)
        D;      % Bending stiffness of plate E*h^3 / 12(1-v^2)
        
        % Driver specific paramaters:
        driver_locations;
        driver_responses;
        num_drivers;
        
        % Paramater to save the frequencies at which the driver_responses
        % are provided for
        frequencies;
        
        % Store the number of modes in each dimension:
        m_modes;
        n_modes;
        
        % Create a matrix to store the quality factor of each mode
        Qs;
        
        % Driver scans
        Driver_Scans;
        
        
        % Panel grid parameters
        Lx_samp;        % spatial sampling period in x dimension (meters)
        Ly_samp;        % spatial sampling period in y dimension (meters)
        Lx_pixels;      % number of pixels in the x dimension
        Ly_pixels;      % number of pixels in the y dimension
        grid_xs;        % locations (in the x dimension) of the pixels
        grid_ys;        % locations (in the y dimension) of the pixels
        grid_XV;        % meshgrid Xs from grid_xs and grid_ys
        grid_YV;        % meshgrid Ys from grid_xs and grid_ys
        grid_XV_flat;   % flattened grid_XV to a row vector
        grid_YV_flat;   % flattened grid_XV to a row vector
        
   
    end
    
    methods
        
% -------------------------------------------------------------------------
% CONSTRUCTOR
        function self = Clamped_Panel(driver_locations,driver_responses,frequencies,Lx,Ly,h,E,v,rho,mode_limits,varargin)
            
            
            
            % ------------------------------------------------------------
            % CHECK FOR NUMBER OF INPUTS AND ASSIGN DEFAULT VALUE
            
            % TODO: delete disp statements after testing
            
            % set defaults for all parameters, which is a 1m x 1m panel
            % made out of 1 mm thick aluminum material
            default_h = .001;       % meters 
            default_E = 68.9E9;     % pascals
            default_v = 0.334;      
            default_rho = 2700;     % kg / m^2
            
            % default number of modes in each dimension
            default_mode_limits = [15 15];
            
            if nargin < 10
                self.m_modes = default_mode_limits(1);
                self.n_modes = default_mode_limits(2);
                disp('Default mode limits used');
            else
                self.rho = rho;
            end
            
            if nargin < 9
                self.rho = default_rho;
                disp('Default rho used');
            else
                self.rho = rho;
            end
            
            if nargin < 8
                self.v = default_v;
                disp('Default v used');
            else
                self.v = v;
            end
                
            if nargin < 7
                self.E = default_E;
                disp('Default E used');
            else
                self.E = E;
            end
            
            if nargin < 6
                self.h = default_h;
                disp('Default h used');
            else
                self.h = h;
            end
            
            if nargin < 5 && nargin >= 3
                msg = 'Must assign length Lx and width Ly of the panel';
                error(msg);
            end
            
            if nargin < 3
                msg = 'Must assign all driver parameters';
                error(msg);
            end
            
            
            
            
            % UNPACK THE DRIVER INPUT PARAMETERS:
            
            % Check to see how many drivers there are
            driver_locations_size = size(driver_locations);
            self.num_drivers = driver_locations_size(1);
            
            % save the input driver parameters
            self.driver_locations = [driver_locations(:,1)*Lx driver_locations(:,2)*Ly];
            self.driver_responses = driver_responses;
            
            % the frequencies input variable has to be the same length as
            % the number of columns in the driver_responses array
            self.frequencies = frequencies;
            
            
            % UNPACK THE PHYSICAL PANEL PARAMETERS:
            self.Lx = Lx;
            self.Ly = Ly;
            self.mass = self.rho * self.h * self.Lx * self.Ly;
            self.D = (self.E * (self.h^3)) / (12 * (1 - (self.v^2)));
            
            
            % ------------------------------------------------------------
            % PARSE OPTIONAL INPUT ARGUMENTS
            
            % add a Q matrix, add max and min mode numbers, write a
            % function to read in vibrometer data and extract the Q values,
            % or see if the vibrometer can do a full Q analysis. Ask bocko
            % about this. Ability to save those Q values.
            

            % 'SpatialSampling' -> Defaults to 0.01, must be a number
            %                      smaller than min(Lx, Ly)
            %
            % 'Qs'              -> The (i,j)th element of the Qs matrix gives
            %                      the quality factor of the (i,j)th mode
            %                      on the panel. If not specified, it
            %                      defaults to 10 for each mode.
            %
            % 
            %
            %
            % Handle the optional input paramaters
            p = inputParser;                
            default_spatial_sampling = .005; % self.Lx_samp and self.Ly_samp
            default_Qs = ones(self.m_modes,self.n_modes);%*10;
            
            % Parse the input for optional parameters
            addParameter(p,'SpatialSampling',default_spatial_sampling);
            addParameter(p,'Qs',default_Qs);
            parse(p,varargin{:});
            
            % Populate optional params into class properties/control variables
            spatial_sampling = p.Results.SpatialSampling;
            self.Qs = p.Results.Qs;
            
            
            % ------------------------------------------------------------
            % INPUT VALIDATION
            
            % Check to make sure the input driver data is the right sizes:
            driver_responses_size = size(driver_responses);
            frequencies_size = size(frequencies);
            % Checks:
            % 1. frequencies must be a 1 x n matrix
            % 2. the length of frequencies must equal the number
            %    of columns in driver_responses 
            % 3. the number of rows of driver_locations and
            %    driver_responses must be equal
            % 4. driver locations must ber a 2 x n matrix
            
            if frequencies_size(1) ~= 1
                msg = 'The frequencies array must be a singular row vector';
                error(msg);
            end
            
            if frequencies_size(2) ~= driver_responses_size(2)
                msg = 'The number of frequencies in the frequencies vector must match the number of frequencies in the driver responses matrix';
                error(msg);
            end
            
            if driver_locations_size(1) ~= driver_responses_size(1)
                msg = 'The number of rows in the driver locations matrix must equal the number of rows in the driver responses matrix, giving a 1-1 pairing between each driver location and its response';
                error(msg)
            end
            
            if driver_locations_size(2) ~= 2
                msg = 'The number of columns in the driver locations matrix must be 2, such that (xi,yi) data can be stored along the columns';
                error(msg);
            end
            
            % ------------------------------------------------------------
            % MAKE THE GRID DOMAIN:
            
            % Define the grid domain
            self.Lx_samp = spatial_sampling;
            self.Ly_samp = spatial_sampling;
            self.Lx_pixels = floor(Lx / self.Lx_samp);
            self.Ly_pixels = floor(Ly / self.Ly_samp);
            % offset the pixels so that the bottom left pixel (origin pixel)
            % isn't at (0,0), but at (Lx_samp/2, Ly_samp/2) 
            self.grid_xs = ([0:self.Lx_pixels-1] * self.Lx_samp) + (self.Lx_samp/2);
            self.grid_ys = ([0:self.Ly_pixels-1] * self.Ly_samp) + (self.Ly_samp/2);
            
            % define the meshgrid from the xs and ys, and flatten
            [self.grid_XV, self.grid_YV] = meshgrid(self.grid_xs,self.grid_ys);
            self.grid_XV_flat = self.grid_XV(:)';
            self.grid_YV_flat = self.grid_YV(:)';
            

            % Make an empty array, which will contain Single_Driver_Scan
            % objects.
            self.Driver_Scans = [];
        
        
            % populate the self.Driver_Scans array by calling the
            % Single_Driver_Scan constructor function
            for driver_number = 1:self.num_drivers
                
                disp(['Processing Driver ' num2str(driver_number)]);
            
                curr_driver_location = self.driver_locations(driver_number,:);
                curr_driver_response = self.driver_responses(driver_number,:);
            
                self.Driver_Scans = [self.Driver_Scans   Single_Driver_Scan(curr_driver_location,curr_driver_response,self.frequencies,self.Lx,self.Ly,self.h,self.E,self.v,self.rho,self.m_modes,self.n_modes,self.grid_XV,self.grid_YV,self.Qs,self.D)];
            
            end
            
        end
    
            
            
            
           

   
            
       
% -------------------------------------------------------------------------
% get_mode_shape

        function [mode_shape] = get_mode_shape(self,m,n)
            % This function returns an array containing the shape of the
            % (m,n)th mode
            
            mode_shape = sin( (m * pi * self.grid_XV)/self.Lx) .* sin( (n * pi * self.grid_YV)/self.Ly);
            
        end
        
        
% -------------------------------------------------------------------------
% view_total_scan

        function [scan_shape] = view_total_scan(self,freq,show)
            % This function outputs the vibration shape shown on the Panel
            % (with physical data, these are called "vibrometer scans",
            % hence the name of the function). 
            %
            % freq is the desired frequency at which to view the scan,
            % the closest frequency in Hz that is populated in this class
            % will be shown
            %
            % the 'show' input controls if imagesc is called on the
            % scan_shape output parmater. show == 1 enables this function
            % to call imagesc natively
            
            
            % find the cloest frequency
            [~, minidx] = min(abs(freq - self.frequencies));
            
            % initialize the output to zeros, then populate by summing the
            % Single_Driver_Scans at the desired frequency:
            scan_shape = zeros(size(self.grid_XV));
            for idx = 1:self.num_drivers
                scan_shape = scan_shape + self.Driver_Scans(idx).u(:,:,minidx);
            end
            
            % show if enabled
            if show == 1
                imagesc((abs(scan_shape)));
                title(['Frequency = ' num2str(self.frequencies(minidx))])
            end
            
        end

        
% -------------------------------------------------------------------------
% view_total_scan

        function [scan_shape] = view_single_driver_scan(self,freq,driver_number,show)
            % This function outputs the vibration shape excited on the
            % panel from one of the drivers on the panel
            %
            % freq is the desired frequency at which to view the scan,
            % the closest frequency in Hz that is populated in this class
            % will be shown
            %
            % driver_number is the driver who's scan we'll show (this is
            % the Signal_Driver_Scan object in column driver_number in self.Driver_Scans)
            %
            % the 'show' input controls if imagesc is called on the
            % scan_shape output parmater. show == 1 enables this function
            % to call imagesc natively
            
            
            % find the cloest frequency
            [~, minidx] = min(abs(freq - self.frequencies));
            
            % find the closest driver number:
            driver_number = round(driver_number);
            if driver_number < 1
                driver_number = 1;
            elseif driver_number > self.num_drivers
                driver_number = self.num_drivers
            end
            
            
            scan_shape = self.Driver_Scans(driver_number).u(:,:,minidx);

            
            % show if enabled
            if show == 1
                imagesc((abs(scan_shape)));
                title(['Frequency = ' num2str(self.frequencies(minidx))])
            end
            
        end
        
        
        
        
        
        
        
    end
    
end

        
       
