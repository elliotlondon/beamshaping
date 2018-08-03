% This program simulates the intensity distribution of a Gaussian laser
% beam passing through a focusing lens using Fresnel Optics. This program
% then calculates the 1/e^2 value of the beam at various points along the
% focus, using a loop with specified parameters.

[fresn_left_2D,image_2D_xl,image_2D_yl] = fresnel(0.1,0.4885);
[fresn_centre_2D,image_2D_xc,image_2D_yc] = fresnel(0.1,0.5);
[fresn_right_2D,image_2D_xr,image_2D_yr] = fresnel(0.1,0.5115);

%% Functions %%
% function [Fresn,image_x,image_y,e2Width,z] = fresnel(start_z,last_z)
    % Initial parameters
    N = 2.^10;
    f = 0.4; 
    points = 51;
    aperture = 25.4e-3;
    wavelength = 633e-9;
    z = linspace(start_z,last_z,points);
    
    % Make a linear array num_points long from -20 to 20 times aperture width
    % and then calculate length of object plane in metres
    lens_x = linspace(-20*aperture,20*aperture,N);
    lens_y = transpose(linspace(-20*aperture,20*aperture,N));
    Lx = abs(lens_x(1)) + abs(lens_x(end));
    Ly = abs(lens_y(1)) + abs(lens_y(end));
    
    % Create a Gaussian input beam (width dependent upon measured beam size)
    beamradius = 6.4e-3/2;      % Value obtained from measured beam 1/e diameter
    g = exp(-((lens_x.^2 + lens_y.^2)/beamradius.^2));
    
    % Create an image mask
    mask = zeros(N);
    Im = zeros(N);
    I = 1:N;
    Ux = I - N/2;
    Uy = N/2 - I;
    [xmask,ymask] = meshgrid(Ux,Uy);
    aperture = (xmask.^2 + ymask.^2) <= (aperture/beamradius).^2;
    mask(aperture) = 1;
    mask = complex(mask, Im);
    
    % Multiply input beam by lens phase term and mask
    F = exp(1i*pi*(lens_x.^2 + lens_y.^2)/(wavelength*f));       % Calculate phase shift due to lens
    g = g.*F.*mask;
    
    e2Width = 1:points;
    gouy = 1:points;
    for i = 1:points
        % Fresnel transform a distance from initial plane
        LX = N * wavelength * z(i) / Lx;              % Calculate length of image plane in x
        LY = N * wavelength * z(i) / Ly;              % Calculate length of image plane in y
        image_x = linspace(-LX/2,LX/2,N);        % Make array of image plane points spanning +/-LX/2
        image_y = transpose(linspace(-LY/2,LY/2,N));        % Make array of image plane points spanning +/-LY/2
        
        % Fresnel equation
        h = exp(-1i*pi*(lens_x.^2 + lens_y.^2)/(wavelength*z(i)));    % Calculate free space propagation term
        H = exp((1j*pi*(image_x.^2 + image_y.^2))./(wavelength*z(i)));    % Calculate prefactor
        Fresn = H.*(fftshift(fft2(fftshift(h.*g))));
        
        gouy(i) = -angle(Fresn(end/2,end/2));   % Save angle(Fresn) as the Gouy phase
        Fresn = abs(Fresn).^2./max(abs(Fresn(:)).^2);     % Make Fresn absolute
        
        % Find the 1/e^2 of the beam and store this in e2Width
        for j = 1:N
            if Fresn(N/2 + 1,j) >= exp(-2)
                e2Width(i) = abs(image_x(j));
                break
            end
        end
    end
    
    %%% Calculate waist evolution and Gouy phase? %%%
    z = z - (last_z + start_z)/2;       % Scale z so that 0 is in the centre
    
    % Create a Gaussian waist curve using the min(e2Width)
    simfit = 1:points;
    for i = 1:points
        simfit(i) = min(e2Width).*sqrt(1 + ((wavelength.*z(i))./(pi.*(min(e2Width).^2))).^2);
    end
    
    % Calculate the Rayleigh range
    simrayleigh = pi*min(simfit).^2/wavelength;
    
    % Plot the simulated curve
    subplot(2,1,1)
    hold on, grid on
    plot(z,e2Width,'x',z,simfit)
    plot([-min(simfit),-min(simfit)],[0,12e-5],'k-')
    title('Simulated 1/e^2 of a Gaussian beam passing through a focusing lens (-Z_R to Z_R) [2D]')
    xlabel 'Distance (m)', ylabel 'Width (m)'
    set(gca,'XLim',[z(1) z(end)],'Ylim',[4e-5 12e-5])
    simstring = sprintf('Waist size = %.3f \\mum \nRayleigh range = %.3f mm',min(simfit)*1e5,simrayleigh*1e3);
    annotation('textbox',[.375 .5 .3 .394],'String',simstring,'FitBoxToText','on','BackgroundColor','w');
    
    % Calculate the theoretical Gouy phase shift
    realgouy = 1:points;
    for i = 1:points
        realgouy(i) = -atan(z(i)./simrayleigh);
    end
    
    % Plot the width evolution and Gouy phase 
    subplot(2,1,2)
    hold on, grid on
    plot(z,gouy,z,realgouy)
    plot([-min(simfit),-min(simfit)],[-pi/2,pi/2],'k-')
    title('Gouy phase shift of a Gaussian beam passing through a focusing lens')
    xlabel 'Width (m)', ylabel 'Angle (rads)'
    set(gca,'Xlim',[z(1) z(end)],'Ylim',[-pi/2 pi/2],'YTick',(-pi/2:pi/4:pi/2),...,
            'YTickLabels',({'-\pi/2','-\pi/4','0','\pi/4','\pi/2'}))
    
    figure, surf(image_x,image_y,Fresn,'LineStyle','None'), view(2)
    
    % Apply scaling for use in 3x3 subplots (curvefitter50.m)
    scaling_x = -200/image_x(1);
    image_x = image_x*scaling_x;
    scaling_y = -200/image_y(1);
    image_y = image_y*scaling_y;
% end
