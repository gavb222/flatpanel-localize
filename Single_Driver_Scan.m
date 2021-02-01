classdef Single_Driver_Scan < handle
    
    properties
        
        
        % Panel physical parameters
        Lx;     % Length of Panel in X dimension (meters)
        Ly;     % Length of Panel in Y dimension (meters)
        h;      % thickness of panel (meters)
        E;      % Young's modulus of the panel (pa)
        v;      % Poissons ratio of the panel (unitless)
        rho;    % Density of the panel (kg/meters^3)
        D;      % Bending stiffness of plate E*h^3 / 12(1-v^2)
        
        % Driver specific paramaters, there is only one driver
        driver_location;
        driver_response
        
        % Paramater to save the frequencies at which the driver_responses
        % are provided for
        frequencies;
        
        % Store the number of modes in each dimension:
        m_modes;
        n_modes;
        
        % Create a matrix to store the quality factor of each mode
        Qs;
        
        
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
        
        
        % Response matrix and parameters
        u;
        xi;
        yi;
     
          
   
    end
    
    
    
    
    
    methods
        
        
        function self = Single_Driver_Scan(driver_location,driver_response,frequencies,Lx,Ly,h,E,v,rho,m_modes,n_modes,grid_XV,grid_YV,Qs,D)
            
            % attach input arguments. Recall that this class will not be
            % instantiated by hand, but handled by the Clamped_Panel class
            self.driver_location = driver_location;
            self.driver_response = driver_response;
            self.frequencies = frequencies;
            self.Lx = Lx;
            self.Ly = Ly;
            self.h = h;
            self.E = E;
            self.rho = rho;
            self.v = v;
            self.m_modes = m_modes;
            self.n_modes = n_modes;
            self.grid_XV = grid_XV;
            self.grid_YV = grid_YV;
            self.Qs = Qs;
            self.D = D;
            
            
            
            
            % pick driver location
            self.xi = driver_location(1);
            self.yi = driver_location(2);
           
            % initialize the panel grid, and thus the response matrix (u
            % for displacement)
            grid_size = size(self.grid_XV);
            self.u = zeros(grid_size(1),grid_size(2),length(self.frequencies));
            
            % compute the governing equation for the response of the panel
            % due to a point force at (xi,yi)
            coeff1 = 4 / (self.rho * self.Lx * self.Ly * self.h);
            
            for f_idx = 1:length(self.frequencies)
                f = self.frequencies(f_idx);  % expecting radians
                w = 2 * pi * f;
                
                running_total_at_w = zeros(grid_size(1),grid_size(2));
                
                for m = 1:self.m_modes
                    for n = 1:self.n_modes
                           
                        w_mn = self.get_mode_frequency(m,n);
                        curr_Q = self.Qs(m,n);
                       
                        curr_shape = self.get_mode_shape(m,n);
                        curr_shape = curr_shape * coeff1;
                        
               
                        coeff2 = 1 / ( (w^2) - (w_mn^2) + ( (1i * w * w_mn) / (curr_Q)));
                        
                        curr_shape = curr_shape * coeff2;
                        
                        coeff3 = self.driver_response(f_idx) * sin((m*pi*self.xi)/self.Lx) * sin((n*pi*self.yi)/self.Ly);
                        
                        curr_shape = curr_shape * coeff3;
                        
                        running_total_at_w = running_total_at_w + curr_shape;
                        
                        
                    end
                    
                    
                    self.u(:,:,f_idx) = running_total_at_w;
                    
                end
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
 

% -------------------------------------------------------------------------
% get_mode_frequency
        function [w] = get_mode_frequency(self,m,n)
            % This function the RADIAL (w) resonant frequency of the (m,n)th
            % mode. This uses 'edge effect factors' proposed by mitchell
            % and hazel

            dm = 1 /((((n * self.Lx) / (m * self.Ly))^2) + 2);
            dn = 1/((((m * self.Ly) / (n * self.Lx))^2) + 2);
            
            coeff1 = (pi^4 * self.D) / (self.rho * self.h);
            km_div_pi = (m + dm) / self.Lx;
            kn_div_pi = (n + dn) / self.Ly;
            coeff2 = (km_div_pi^2) + (kn_div_pi^2);
            coeff2 = coeff2^2;
            w_squared = coeff1*coeff2;
            w = sqrt(w_squared);
            
            
        end
% -------------------------------------------------------------------------
 


        

        
        
        
        
    end
end

