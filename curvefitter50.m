% 1. This program extracts data from a .txt file created by Labview, then plots a curve 
% from the equation w(z) = w0*sqrt(1+((\lambda*z)/pi*w0^2)^2) to the extracted data, 
% along with simulated data and finds a curve minimum.
% 2. Plots 3D meshes of the raw data
% 3. Plots 3D meshes of w(z) = w0*sqrt(1+((\lambda*z)/pi*w0^2)^2) with w0
% derived from the raw data
% 4. Create subplots comparing the widths and lineouts at -z_R, 0 and z_R.

%% PART 1 %%
% Load the data and create initial parameters
load 2018-06-18.scan001-data.txt
fname = X2018_06_18_scan001_data;           % Store the filename in a var.
n = 209;                                    % Number of rows
p = 50;                                     % Number of data points
pixel_size = 2.4*(-99.5:99.5);              % Define pixel size

% Extract the required data and store in suitable variables
Xdata = fname(1:4:end-8,:);
Xgauss = fname(2:4:end-4,:);
Ydata = fname(3:4:end-4,:);
Ygauss = fname(4:4:end-4,:);
e2x = fname(n-4,1:p+1)*2.4;
e2y = fname(n-3,1:p+1)*2.4;
xerror = fname(n-2,1:p+1)*1e-6;
yerror = fname(n-1,1:p+1)*1e-6;

X = fname(n,1:(p+1));               % Grab the X values
X = X.*1e-3;                        % Put X in meters (starts with mm)
X = X - max(X)/2;                   % Shift X so that the central point is zero

% Raw data is FWHM, change to 1/e^2 and metres
e2x = e2x*1e-6;
e2y = e2y*1e-6;

% Calculate the function to be solved, a(1) == w_0. a(2) is an x-axis shift
funx = @(a,X) a(1).*sqrt(1 + ((633e-9.*X)./(pi.*(a(1).^2)) + a(2)).^2);
funy = @(a,X) a(1).*sqrt(1 + ((633e-9.*X)./(pi.*(a(1).^2)) + a(2)).^2);

% Generate the nonlinear fits
a0x = min(e2x);         % x initial guess
a0y = min(e2y);         % y initial guess
w0x = lsqcurvefit(funx,[a0x,5e-6],X,e2x);         % Coefficient for fitting x
w0y = lsqcurvefit(funy,[a0y,6e-6],X,e2y);         % Coefficient for fitting y

% Find minima and Rayleigh ranges, and change to appropriate units for plotting
curveMinX = min(funx(w0x,X));
curveMinY = min(funy(w0y,X));
xRayleigh = pi*curveMinX.^2/633e-9*1e3;
yRayleigh = pi*curveMinY.^2/633e-9*1e3;
curveMinX = curveMinX*1e6;
curveMinY = curveMinY*1e6;

% Calculate the M2 values for the fwhm data against the fitted curves
% Requisite: If the ratio is smaller than one, take the reciprocal
xM2 = (1:(p+1));
yM2 = (1:(p+1));
for i = 1:(p+1)
    if funx(w0x,X(i)) >= e2x(i)
        xM2(i) = funx(w0x,X(i))./e2x(i);
    else
        xM2(i) = e2x(i)./funx(w0x,X(i));
    end
    if funy(w0y,X(i)) >= e2y(i)
        yM2(i) = funx(w0y,X(i))./e2y(i);
    else
        yM2(i) = e2y(i)./funx(w0y,X(i));
    end
end
xM2 = sum(xM2)./(p+1);
yM2 = sum(yM2)./(p+1);

% Plot the data
figure
errorbar(X,e2x,xerror,'x')              % x data with errors
hold on, grid on
plot(X,funx(w0x,X))                     % x fit
errorbar(X,e2y,yerror,'x')              % y data with errors
plot(X,funy(w0y,X))                     % y fit
% plot(X,simfit)                          % Simulated data
legend('x-axis 1/e^2','x Best fit','y-axis 1/e^2','y Best Fit', ...,
       'Simulated data','Location','NorthWest')
title('Measured data for an aberrated beam, post-correction');
xlabel 'Distance (m)', ylabel 'Width (m)'
set(gca,'Xlim',[-12.5e-3 12.5e-3],'XTick',(-12.5e-3:2.5e-3:12.5e-3),'Ylim',[1e-5 15e-5])

% Add the Rayleigh ranges and w_0 values to a texbox on the graph
textboxStringX = sprintf('w_0x = %.3f \\mum \n \t z_Rx = %.3f mm \n',curveMinX,xRayleigh);
textboxStringY = sprintf('\nw_0y = %.3f \\mum \n \t z_Ry = %.3f mm \n',curveMinY,yRayleigh);
% textboxStringsim = sprintf('\nSimulated min = %.3f \\mum \n Simulated z_R = %.3f mm', ...,
%                          min(simfit)*1e6,simrayleigh*1e3);
textboxStringM2 = sprintf('\nM^2x = %.3f \t \n M^2y = %.3f',xM2,yM2);
textboxTotalString = strcat(textboxStringX,textboxStringY,textboxStringM2);
annotation('textbox',[.375 .5 .3 .394],'String',textboxTotalString,'FitBoxToText','on','BackgroundColor','w');

%% PART 2 %%
% Plot the raw x data
figure
surf(pixel_size,X,Xdata,'LineStyle','None'), view(90,90)
title('Width/Distance for raw x data, no mask/aberrations');
xlabel 'Width (\mum)', ylabel 'Distance (m)', zlabel 'Measured Intensity'
set(gca,'XLim',[-2.4*99.5 2.4*99.5],'YLim',[-12.5e-3 12.5e-3])
% Plot the raw y data
figure, surf(pixel_size,X,Ydata,'LineStyle','None'), view(90,90)
title('Width/Distance for raw y data, no mask/aberrations');
xlabel 'Width (\mum)', ylabel 'Distance (m)', zlabel 'Measured Intensity'
set(gca,'Xlim',[-2.4*99.5 2.4*99.5],'YLim',[-12.5e-3 12.5e-3])

%% PART 3 %% 
% "Lock in" the value found for w0 into the x and y functions, so that when calculating 
% peak intensity, a different w0 isn't used for each value
fittedX_orig = funx(w0x,X);
fittedY_orig = funy(w0y,X);
fittedX_pixel = 1e-4*funx(w0x,pixel_size);
fittedY_pixel = 1e-4*funy(w0y,pixel_size);

% Find the radii of the beamsradiusX = (1:(p+1));
radiusX = (1:(p+1));
radiusY = (1:(p+1));
for i = 1:(p+1)
    radiusX(i) = abs(1./fittedX_orig(i) - a0x);
    radiusY(i) = abs(1./fittedY_orig(i) - a0y);
end

% Calculate the 1/e.^2 Gaussian distributions with appropriate peak intensities I0
e2xGaussian = zeros((p+1),200);
e2yGaussian = zeros((p+1),200);
for i = 1:(p+1)
    for j = 1:200
        e2xGaussian(i,j) = exp(-2*((fittedX_pixel(j)).^2)./fittedX_orig(i).^2)/fittedX_orig(i);
        e2yGaussian(i,j) = exp(-2*((fittedY_pixel(j)).^2)./fittedY_orig(i).^2)/fittedY_orig(i);
    end
end

% Plot the fitted x data
figure
surf(pixel_size,X,e2xGaussian,'LineStyle','None')
title('Width/Distance for fitted x data, no mask/aberrations');
xlabel 'Width (\mum)', ylabel 'Distance (m)', zlabel 'Normalized Intensity'
set(gca,'Xlim',[-2.4*99.5 2.4*99.5],'YLim',[-12.5e-3 12.5e-3]), view(90,90)
% Plot the fitted y data
figure
surf(pixel_size,X,e2yGaussian,'LineStyle','None')
title('Width/Distance for fitted y data, no mask/aberrations');
xlabel 'Width (\mum)', ylabel 'Distance (m)', zlabel 'Normalized Intensity'
set(gca,'Xlim',[-2.4*99.5 2.4*99.5],'YLim',[-12.5e-3 12.5e-3]), view(90,90)

%% Plot the bowties side-by-side
figure
subplot(1,2,1)
surf(pixel_size,X,Xdata,'LineStyle','None'), view(90,90)
title('Width/Distance for raw x data, no mask/aberrations');
xlabel 'Width (\mum)', ylabel 'z (m)', zlabel 'Measured Intensity'
set(gca,'XLim',[-2.4*99.5 2.4*99.5],'YLim',[-12.5e-3 12.5e-3])

subplot(1,2,2)
surf(pixel_size,X,e2xGaussian,'LineStyle','None')
title('Width/Distance for fitted x data, no mask/aberrations');
xlabel 'Width (\mum)', ylabel 'z (m)', zlabel 'Normalized Intensity'
set(gca,'Xlim',[-2.4*99.5 2.4*99.5],'YLim',[-12.5e-3 12.5e-3]), view(90,90)

%% Plot the simulated profile and lineouts
%%% ROW 1 %%%
% Plot the 2D beam profiles
figure
subplot(2,3,2)
surf(image_2D_xc,image_2D_yc,fresn_centre_2D,'LineStyle','None')
title('Simulated 2D plots')
set(gca,'Xlim',[-200 200],'Ylim',[-200 200]), view(2)

subplot(2,3,1)
surf(image_2D_xl,image_2D_yl,fresn_left_2D,'LineStyle','None')
set(gca,'Xlim',[-200 200],'Ylim',[-200 200]), view(2)

subplot(2,3,3)
surf(image_2D_xr,image_2D_yr,fresn_right_2D,'LineStyle','None')
set(gca,'Xlim',[-200 200],'Ylim',[-200 200]), view(2)

%%% ROW 2 %%%
% Plot the simulated lineouts
subplot(2,3,4);
plot(image_2D_xl,fresn_left_2D(end/2,:)./sqrt(2))
set(gca,'Xlim',[-2.4*99.5 2.4*99.5],'Ylim',[0 1])

subplot(2,3,5);
plot(image_2D_xc,fresn_centre_2D(end/2,:))
set(gca,'Xlim',[-2.4*99.5 2.4*99.5],'Ylim',[0 1]), title('Simulated lineouts')

subplot(2,3,6);
plot(image_2D_xr,fresn_right_2D(end/2,:)./sqrt(2))
set(gca,'Xlim',[-2.4*99.5 2.4*99.5],'Ylim',[0 1])

%% Plot the measured profiles and lineouts
load Centre_focal.txt
load Centre_minuszR.txt
load Centre_pluszR.txt

%%% ROW 1 %%%
figure
subplot(2,3,1)
surf(pixel_size,pixel_size,Centre_minuszR,'LineStyle','None'), view(2)
set(gca,'Xlim',[-99.5*2.4 99.5*2.4],'Ylim',[-99.5*2.4 99.5*2.4])
subplot(2,3,2)
surf(pixel_size,pixel_size,Centre_focal,'LineStyle','None'), view(2), title('Measured 2D plots')
set(gca,'Xlim',[-99.5*2.4 99.5*2.4],'Ylim',[-99.5*2.4 99.5*2.4])
subplot(2,3,3)
surf(pixel_size,pixel_size,Centre_pluszR,'LineStyle','None'), view(2)
set(gca,'Xlim',[-99.5*2.4 99.5*2.4],'Ylim',[-99.5*2.4 99.5*2.4])

%%% ROW 2 %%%
subplot(2,3,4)
plot(pixel_size,Xdata(24,:)./max(Xdata(24,:))./sqrt(2))
set(gca,'Xlim',[-99.5*2.4 99.5*2.4],'Ylim',[0 1])
subplot(2,3,5)
plot(pixel_size,Xdata(27,:)./max(Xdata(27,:))), title('Measured lineouts')
set(gca,'Xlim',[-99.5*2.4 99.5*2.4],'Ylim',[0 1])
subplot(2,3,6)
plot(pixel_size,Xdata(30,:)./max(Xdata(30,:))./sqrt(2))
set(gca,'Xlim',[-99.5*2.4 99.5*2.4],'Ylim',[0 1])

%% Plot the simulated and raw lineouts together
figure
subplot(1,3,1);
plot(pixel_size,Xdata(24,:)./max(Xdata(24,:))./sqrt(2));
hold on
plot(image_2D_xl,fresn_left_2D(end/2,:)./sqrt(2))
set(gca,'Xlim',[-2.4*99.5 2.4*99.5],'Ylim',[0 1])

subplot(1,3,2);
plot(pixel_size,Xdata(27,:)./max(Xdata(27,:)));
hold on
plot(image_2D_xc,fresn_centre_2D(end/2,:))
set(gca,'Xlim',[-2.4*99.5 2.4*99.5],'Ylim',[0 1]), title('Simulated & Measured lineouts');

subplot(1,3,3);
plot(pixel_size,Xdata(30,:)./max(Xdata(30,:))./sqrt(2));
hold on
plot(image_2D_xr,fresn_right_2D(end/2,:)./sqrt(2))
set(gca,'Xlim',[-2.4*99.5 2.4*99.5],'Ylim',[0 1])

%% Plot a 3x3 of the simulated and measured profiles and overlapped lineouts
figure
subplot(3,3,2)
surf(image_2D_xc,image_2D_yc,fresn_centre_2D,'LineStyle','None')
title('Simulated 2D plots')
set(gca,'Xlim',[-99.5*2.4 99.5*2.4],'Ylim',[-99.5*2.4 99.5*2.4]), view(2)
subplot(3,3,1)
surf(image_2D_xl,image_2D_yl,fresn_left_2D,'LineStyle','None')
set(gca,'Xlim',[-99.5*2.4 99.5*2.4],'Ylim',[-99.5*2.4 99.5*2.4]), view(2)
subplot(3,3,3)
surf(image_2D_xr,image_2D_yr,fresn_right_2D,'LineStyle','None')
set(gca,'Xlim',[-99.5*2.4 99.5*2.4],'Ylim',[-99.5*2.4 99.5*2.4]), view(2)

subplot(3,3,4)
surf(pixel_size,pixel_size,Centre_minuszR,'LineStyle','None'), view(2)
set(gca,'Xlim',[-99.5*2.4 99.5*2.4],'Ylim',[-99.5*2.4 99.5*2.4])
subplot(3,3,5)
surf(pixel_size,pixel_size,Centre_focal,'LineStyle','None'), view(2), title('Measured 2D plots')
set(gca,'Xlim',[-99.5*2.4 99.5*2.4],'Ylim',[-99.5*2.4 99.5*2.4])
subplot(3,3,6)
surf(pixel_size,pixel_size,Centre_pluszR,'LineStyle','None'), view(2)
set(gca,'Xlim',[-99.5*2.4 99.5*2.4],'Ylim',[-99.5*2.4 99.5*2.4])

subplot(3,3,7);
plot(pixel_size,Xdata(24,:)./max(Xdata(24,:))./sqrt(2));
hold on
plot(image_2D_xl,fresn_left_2D(end/2,:)./sqrt(2))
set(gca,'Xlim',[-2.4*99.5 2.4*99.5],'Ylim',[0 1])
subplot(3,3,8);
plot(pixel_size,Xdata(27,:)./max(Xdata(27,:)));
hold on
plot(image_2D_xc,fresn_centre_2D(end/2,:))
set(gca,'Xlim',[-2.4*99.5 2.4*99.5],'Ylim',[0 1]), title('Simulated & Measured lineouts');
subplot(3,3,9);
plot(pixel_size,Xdata(30,:)./max(Xdata(30,:))./sqrt(2));
hold on
plot(image_2D_xr,fresn_right_2D(end/2,:)./sqrt(2))
set(gca,'Xlim',[-2.4*99.5 2.4*99.5],'Ylim',[0 1])