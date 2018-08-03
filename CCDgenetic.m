%% Mirror and camera initialisation
close; clear; clc

% Add necessary directories to the MATLAB path
addpath(genpath('C:\Program Files\Alpao'));
addpath('C:\Users\Kevin\Downloads\Elliot\MRes_Dropbox\Experimental_Data_(Elliot)\exchangeData');

% Establish mirror connection and set mirror parameters
mirrorSN = 'BAX160';        % Set mirror serial name
dm = asdkDM(mirrorSN);      % Initialise new mirror object
nbAct = dm.Get('NbOfActuator');     % Get the number of actuators

% Initialise actuator data and send it to the mirror
dmdata = zeros(nbAct,1);
dm.Send(dmdata);

% Set up variables and algorithmic parameters
pixels = 1:200;
NumberOfUnknowns = length(dmdata);
StringLength = 12;
TotalStringLength = StringLength*NumberOfUnknowns;
Individuals = 25;
Generations = 150;
CrossOverProbability = 0.8;
InitialMutationRate = 0.005; % Each bit has a .5% chance to change
MutationRate = InitialMutationRate;
UpperRange = 0.025;
LowerRange = -0.025;
WidthGuess = 10;
Pressure = 1;       % Apply evolutionary pressure (increase mutation if algorithm stalls)
StallTolerance = ceil(Generations/10);    % Number of generations before applying pressure

% Apply elitism (automatically select a % of the best to pass through)
Elitism = 1;
if Elitism == 1
    EliteNumber = ceil(0.1*Individuals);
else
    EliteNumber = 0;
end

% Connect to camera
vid = imaq.VideoDevice('winvideo',1,'RGB32_3088x2076','ROI',[2210 688 200 200]);
ccddata = im2double(step(vid));

% Plot initial beam profile
zfig = figure(3);
surf(pixels,pixels,ccddata(:,:,2),'LineStyle','None'), view(2)
set(zfig,'Position',[850 550 400 400])

% Print a seeded message upon section completion
fprintf('Mirror Initialisation Successful (Checkrand: %f)\n',rand);

%% Genetic algorithm main program
% Initialise the fields of the structs and variables
Population(Individuals).Strings = [];
Population(Individuals).Individual = [];
Population(Individuals).Fitness = [];
NewStrings(Individuals).Strings = [];
MeanFitness = 1:Generations;
MaxFitness = 1:Generations;
FitnessStdev = 1:Generations;
BestIndividual = 1:TotalStringLength;
StallCounter = 0;

% First generation
for gen = 1
    % Create an initial population of individuals
    for indiv = 1:Individuals
        Population(indiv).Strings = binarystring(StringLength,NumberOfUnknowns);
        Population(indiv).Individual = fastbin2dec(Population(indiv).Strings,UpperRange,LowerRange,StringLength);
        
        % Send surface to mirror
        dmdata(:,1) = Population(indiv).Individual(:);
        dm.Send(dmdata);
        pause(0.25)
        
        % Grab image from the camera
        ccddata = im2double(step(vid));
        
        % Find the brightest pixel using a centre of gravity (cog) method
        cog = centreofgravity(ccddata);
        if isnan(cog)
            Population(indiv).Fitness = 0;
            continue
        end
        
        % Find the lineouts
        [xlineout,ylineout] = lineouts(pixels,ccddata);
        
        % Generate gaussian fits for the image lineouts
        [gaussx,gaussy,fitx,fity,resx,resy] = gaussfit(pixels,xlineout,ylineout);

        % Calculate the ratio between the widths of the fits
        if fitx(1) > fity(1)
            widthratio = fitx(1)/fity(1);
        else
            widthratio = fity(1)/fitx(1);
        end
        
        if fitx(1) <= WidthGuess/10
            Population(indiv).Fitness = 0;
            continue
        end
        
        if fity(1) <= WidthGuess/10
            Population(indiv).Fitness = 0;
            continue
        end
        
        % Create a pupil function
        beamdiameter = WidthGuess/2*(pi/2);
        [x,y] = meshgrid(pixels,pixels);
        pupil = ((x - max(pixels)/2).^2 + (y - max(pixels)/2).^2 <= beamdiameter.^2);
        selecteddata = ccddata(:,:,2);
        selecteddata(pupil) = 0;
        selecteddata(selecteddata<0.015) = 0;
        pupilsum = sum(sum(selecteddata))/50;
        
        % Create a fitness function
        Population(indiv).Fitness = 1.025/...   % Would be 1 but is accounting for noise
                                    (((fitx(1)/2 + fity(1)/2)/WidthGuess + pupilsum)...
                                     *(widthratio + (resx + resy)/5));
        
        % Display data for every individual
        clc
        fprintf('Individual = %i/%i, ',indiv,max(Individuals));
        fprintf('Generation = %i/%i\n',gen,max(Generations));
        fprintf('Xwidth = %f, Ywidth = %f\n',fitx(1),fity(1));
        fprintf('Lineout width ratios = %f\n',widthratio);
        fprintf('Xres = %.3f, Yres = %.3f\n',resx,resy);
        if max(max(ccddata(:,:,2))) >= 0.98
            fprintf('Warning! Camera saturated!')
            pause(0.5)
        end

        % Plot the x data
        xfig = figure(1);
        set(xfig,'Position',[50 550 400 400])
        plot(pixels,xlineout,pixels,gaussx)
        title('xlineout'), xlim([0 max(pixels)]), ylim([0 1])
        % Plot the y data
        yfig = figure(2);
        set(yfig,'Position',[450 550 400 400])
        plot(pixels,ylineout,pixels,gaussy)
        title('ylineout'), xlim([0 max(pixels)]), ylim([0 1])
        % Plot beam profile in 2D
        zfig = figure(3);
        set(zfig,'Position',[850 550 400 400])
        surf(pixels,pixels,ccddata(:,:,2),'LineStyle','None'), view(2)
    end
    
    % Find fitness statistics
    MeanFitness(gen) = mean([Population.Fitness]);
    SumFitness = sum([Population.Fitness]);
    [MaxFitness(gen),MaxFitnessIndex] = max([Population.Fitness]);
    FitnessStdev(gen) = std([Population.Fitness]);
    BestIndividual = Population(MaxFitnessIndex).Individual;
    
    % Display data for every generation (for loop for cleanness)
    clc
    fprintf('Individual = %i/%i, ',indiv,max(Individuals));
    fprintf('Generation = %i/%i\n',gen,max(Generations));
    fprintf('Xwidth = %f, Ywidth = %f\n',fitx(1),fity(1));
    fprintf('Lineout width ratios = %f\n',widthratio);
    fprintf('Xres = %.3f, Yres = %.3f\n',resx,resy);
    fprintf('MaxFitness = %f, ',MaxFitness(gen));
    fprintf('MeanFitness = %f\n',MeanFitness(gen));
    fprintf('StallCounter = %i, ',StallCounter);
    fprintf('MutationRate = %.3f\n',MutationRate);
end

% All other generations
for gen = 2:Generations
    if Elitism == 1
        % Apply Elitism, automatically sending the best individual through
        [~,index] = sortrows([Population.Fitness].');
        Population = Population(index(end:-1:1));
        clear index

        for indiv = 1:EliteNumber
            NewStrings(indiv).Strings = Population(indiv).Strings;
        end
    end
    
    % Generate offspring for the rest of the generation
    for indiv = (EliteNumber + 1):Individuals
        % Select individuals for crossover
        Parent1 = selection(Individuals,Population,SumFitness);
        Parent2 = selection(Individuals,Population,SumFitness);
        
        % Perform crossover and generate a temporary population
        if rand <= CrossOverProbability
            for j = 1:NumberOfUnknowns
                % Generate two random locations in the binary string for crossover
                CrossSite1 = round((StringLength - 1)*rand + 1);
                CrossSite2 = round((StringLength - 1)*rand + 1);
                % Grab the strings of the parents
                Selected1 = Population(Parent1).Strings(j,:);
                Selected2 = Population(Parent2).Strings(j,:);
                % Perform crossover
                if rand >= 0.5      % If rand >= 0.5, swap first half
                    NewStrings(indiv).Strings(j,:) = horzcat(Selected1(1:CrossSite1),Selected2(CrossSite1+1:StringLength));
                else                % Else swap second half
                    NewStrings(indiv).Strings(j,:) = horzcat(Selected1(1:CrossSite2),Selected2(CrossSite2+1:StringLength));
                end
                % If crossover is performed, randomly mutate the children
                for k = 1:StringLength
                    if rand <= MutationRate
                        if NewStrings(indiv).Strings(j,k) == 1
                            NewStrings(indiv).Strings(j,k) = 0;
                        else
                            NewStrings(indiv).Strings(j,k) = 1;
                        end
                    end
                end
            end
        % No crossover case
        else
            if rand >= 0.5
                NewStrings(indiv).Strings = Population(Parent1).Strings;
            else
                NewStrings(indiv).Strings = Population(Parent2).Strings;
            end
            % Randomly mutate the children
            for k = 1:StringLength
                if rand <= MutationRate
                    if NewStrings(indiv).Strings(j,k) == 1
                        NewStrings(indiv).Strings(j,k) = 0;
                    else
                        NewStrings(indiv).Strings(j,k) = 1;
                    end
                end
            end
        end
    end
    
    % Transfer elite population without recalculating fitness
    if Elitism == 1
        for indiv = 1:EliteNumber
            Population(indiv).Strings = NewStrings(indiv).Strings;
            Population(indiv).Individual = fastbin2dec(Population(indiv).Strings,UpperRange,LowerRange,StringLength);

            % Send surface to mirror
            dmdata(:,1) = Population(indiv).Individual(:);
            dm.Send(dmdata);
            pause(0.25)

            % Grab image from the camera
            ccddata = im2double(step(vid));

            % Display data every loop
            clc
            fprintf('Individual = %i/%i (Elite), ',indiv,max(Individuals));
            fprintf('Generation = %i/%i\n',gen,max(Generations));
            fprintf('Xwidth = %f, Ywidth = %f\n',fitx(1),fity(1));
            fprintf('Lineout width ratios = %f\n',widthratio);
            fprintf('Xres = %.3f, Yres = %.3f\n',resx,resy);
            fprintf('MaxFitness = %f, ',MaxFitness(gen-1));
            fprintf('MeanFitness = %f\n',MeanFitness(gen-1));
            fprintf('StallCounter = %i, ',StallCounter);
            fprintf('MutationRate = %.3f\n',MutationRate);
            if max(max(ccddata(:,:,2))) >= 0.98
                fprintf('Warning! Camera saturated!')
                pause(0.5)
            end

            % Plot the x data
            xfig = figure(1);
            set(xfig,'Position',[50 550 400 400])
            plot(pixels,xlineout,pixels,gaussx)
            title('xlineout'), xlim([0 max(pixels)]), ylim([0 1])
            % Plot the y data
            yfig = figure(2);
            set(yfig,'Position',[450 550 400 400])
            plot(pixels,ylineout,pixels,gaussy)
            title('ylineout'), xlim([0 max(pixels)]), ylim([0 1])
            % Plot beam profile in 2D
            zfig = figure(3);
            set(zfig,'Position',[850 550 400 400])
            surf(pixels,pixels,ccddata(:,:,2),'LineStyle','None'), view(2)
        end
    end
    
    % Transfer crossover children to next generation
    for indiv = (EliteNumber + 1):Individuals
        Population(indiv).Strings = NewStrings(indiv).Strings;
        Population(indiv).Individual = fastbin2dec(Population(indiv).Strings,UpperRange,LowerRange,StringLength);
        
        % Send surface to mirror
        dmdata(:,1) = Population(indiv).Individual(:);
        dm.Send(dmdata);
        pause(0.25)
        
        % Grab image from the camera
        ccddata = im2double(step(vid));
        
        % Find the brightest pixel using a centre of gravity (cog) method
        cog = centreofgravity(ccddata);
        if isnan(cog)
            Population(indiv).Fitness = 0;
            continue
        end
        
        % Find the lineouts
        [xlineout,ylineout] = lineouts(pixels,ccddata);
        
        % Generate gaussian fits for the image lineouts
        [gaussx,gaussy,fitx,fity,resx,resy] = gaussfit(pixels,xlineout,ylineout);
        
        % Calculate the ratio between the widths of the fits
        if fitx(1) > fity(1)
            widthratio = fitx(1)/fity(1);
        else
            widthratio = fity(1)/fitx(1);
        end
        
        if fitx(1) <= WidthGuess/10
            Population(indiv).Fitness = 0;
            continue
        end
        
        if fity(1) <= WidthGuess/10
            Population(indiv).Fitness = 0;
            continue
        end
        
        % Create a pupil function
        beamdiameter = WidthGuess/2*(pi/2);
        [x,y] = meshgrid(pixels,pixels);
        pupil = ((x - max(pixels)/2).^2 + (y - max(pixels)/2).^2 <= beamdiameter.^2);
        selecteddata = ccddata(:,:,2);
        selecteddata(pupil) = 0;
        selecteddata(selecteddata<0.015) = 0;
        pupilsum = sum(sum(selecteddata))/50;
        
        % Create a fitness function
        Population(indiv).Fitness = 1.025/...
                                    (((fitx(1)/2 + fity(1)/2)/WidthGuess + pupilsum)...
                                     *(widthratio + (resx + resy)/5));
        
        % Display data every loop
        clc
        fprintf('Individual = %i/%i, ',indiv,max(Individuals));
        fprintf('Generation = %i/%i\n',gen,max(Generations));
        fprintf('Xwidth = %f, Ywidth = %f\n',fitx(1),fity(1));
        fprintf('Lineout width ratios = %f\n',widthratio);
        fprintf('Xres = %.3f, Yres = %.3f, ',resx,resy);
        fprintf('PupilSum = %.3f\n',pupilsum);
        fprintf('MaxFitness = %f, ',MaxFitness(gen-1));
        fprintf('MeanFitness = %f\n',MeanFitness(gen-1));
        fprintf('StallCounter = %i, ',StallCounter);
        fprintf('MutationRate = %.3f\n',MutationRate);
        if max(max(ccddata(:,:,2))) >= 0.98
            fprintf('Warning! Camera saturated!')
            pause(0.5)
        end

        % Plot the x data
        xfig = figure(1);
        set(xfig,'Position',[50 550 400 400])
        plot(pixels,xlineout,pixels,gaussx)
        title('xlineout'), xlim([0 max(pixels)]), ylim([0 1])
        % Plot the y data
        yfig = figure(2);
        set(yfig,'Position',[450 550 400 400])
        plot(pixels,ylineout,pixels,gaussy)
        title('ylineout'), xlim([0 max(pixels)]), ylim([0 1])
        % Plot beam profile in 2D
        zfig = figure(3);
        set(zfig,'Position',[850 550 400 400])
        surf(pixels,pixels,ccddata(:,:,2),'LineStyle','None'), view(2)
    end
    
    % Find fitness statistics
    MeanFitness(gen) = mean([Population.Fitness]);
    SumFitness = sum([Population.Fitness]);
    [MaxFitness(gen),MaxFitnessIndex] = max([Population.Fitness]);
    FitnessStdev(gen) = std([Population.Fitness]);
    
    % Calculate whether a stall has occured
    if MaxFitness(gen) == MaxFitness(gen - 1)
        StallCounter = StallCounter + 1;
    else
        MutationRate = InitialMutationRate;
        StallCounter = 0;
    end
    
    % Display data for every generation
    clc
    fprintf('Individual = %i/%i, ',indiv,max(Individuals));
    fprintf('Generation = %i/%i\n',gen,max(Generations));
    fprintf('Xwidth = %f, Ywidth = %f\n',fitx(1),fity(1));
    fprintf('Lineout width ratios = %f\n',widthratio);
    fprintf('Xres = %.3f, Yres = %.3f, ',resx,resy);
    fprintf('PupilSum = %.3f\n',pupilsum);
    fprintf('MaxFitness = %f, ',MaxFitness(gen));
    fprintf('MeanFitness = %f\n',MeanFitness(gen));
    fprintf('StallCounter = %i, ',StallCounter);
    fprintf('MutationRate = %.3f\n',MutationRate);
    
    % Apply evolutionary pressure
    if (StallCounter == StallTolerance) && (Pressure == 1)
        if MutationRate <= 0.25     % Increasing beyond this would have a mutation rate
                                    % that is too high.
            MutationRate = MutationRate + InitialMutationRate;
            fprintf('Stress Applied!\n')
            pause(0.5)
            StallCounter = 0;
        end
    end
end

% Find best ever generation
[~,index] = sortrows([Population.Fitness].');
Population = Population(index(end:-1:1));
clear index
% Send best individual to mirror
dmdata(:,1) = Population(1).Individual;
dm.Send(dmdata);
pause(0.25)

% Grab image from the camera
ccddata = im2double(step(vid));
delete(vid);

% Find the brightest pixel using a centre of gravity (cog) method
cog = centreofgravity(ccddata);

% Find the lineouts
[xlineout,ylineout] = lineouts(pixels,ccddata);

% Generate gaussian fits for the image lineouts
[gaussx,gaussy,fitx,fity] = gaussfit(pixels,xlineout,ylineout);

% Save variables for later use
save('C:\Users\Kevin\Downloads\Elliot\MRes_Dropbox\Experimental_Data_(Elliot)\DMga\CorrectedSurface.mat','dmdata')

% Display data for the final generation
clc
fprintf('Individual = %i/%i, ',indiv,max(Individuals));
fprintf('Generation = %i/%i\n',gen,max(Generations));
fprintf('Xwidth = %f, Ywidth = %f\n',fitx(1),fity(1));
fprintf('Lineout width ratios = %f\n',widthratio);
fprintf('Xres = %.3f, Yres = %.3f\n',resx,resy);

% Plot the x data
xfig = figure(1);
set(xfig,'Position',[50 550 400 400])
plot(pixels,xlineout,pixels,gaussx)
title('xlineout'), xlim([0 max(pixels)]), ylim([0 1])
% Plot the y data
yfig = figure(2);
set(yfig,'Position',[450 550 400 400])
plot(pixels,ylineout,pixels,gaussy)
title('ylineout'), xlim([0 max(pixels)]), ylim([0 1])
% Plot beam profile in 2D
zfig = figure(3);
set(zfig,'Position',[850 550 400 400])
surf(pixels,pixels,ccddata(:,:,1),'LineStyle','None'), view(2)
% Plot evolution of mean fitness
figure(4);
plot(1:max(Generations),MeanFitness,1:max(Generations),MaxFitness)
xlabel('Generation'), ylabel('Fitness')
xlim([1 gen]), ylim([0 1])
legend('Average Fitness','Max Fitness')

% Print a seeded message upon section completion
fprintf('\nAlgorithm finished. Checkrand: %f\n',rand);

%% Functions and some brief descriptions
% Create a binary string representing an individual 
function string = binarystring(StringLength,NumberOfUnknowns)
    string = round(ones(NumberOfUnknowns,StringLength).*rand(NumberOfUnknowns,StringLength));
end

% Create a function which turns binary strings into decimal numbers.
% MATLAB's built-in 'bin2dec' is capable of this, but is too slow.
function decimalout = fastbin2dec(binaryin,UpperRange,LowerRange,StringLength)
    decimalout = binaryin*(2.^(size(binaryin, 2)-1:-1:0)).';
    
    % Ensure that the output is scaled between the specified problem range
    decimalout = decimalout*(UpperRange - LowerRange)/(2.^(StringLength)-1) + LowerRange;
    
    % Mirror is 16-bit and minimum actuator range is 0.0005, so ensure that
    % strings are rounded to this value
    accuracy = 0.0005;
    decimalout = round(decimalout/accuracy)*accuracy;
end

% Find the centre of gravity (cog) of the image
function cog = centreofgravity(ccddata)
    binaryImage = true(size(ccddata));
    measurements = regionprops(binaryImage,ccddata,'WeightedCentroid');
    cog = measurements(1).WeightedCentroid;
    cog = round(cog);
end

% Find the lineouts of the image
function [xlineout,ylineout] = lineouts(pixels,ccddata)
    xlineout = ccddata(max(pixels)/2,:,2);
    ylineout = ccddata(:,max(pixels)/2,2);
    % Normalise the lineouts
    xlineout = xlineout./max(xlineout);
    ylineout = transpose(ylineout./max(ylineout));
end

% Generate gaussian fits to the image lineouts
function [gaussx,gaussy,fitx,fity,resx,resy] = gaussfit(pixels,xlineout,ylineout)
    % Set options
    options = optimset('Display','off');
    
    % Generate gaussian fits
    gaussx = @(a,xlineout) exp(-(xlineout - a(2)).^2/(2*sqrt(2)*a(1)));
    gaussy = @(a,ylineout) exp(-(ylineout - a(2)).^2/(2*sqrt(2)*a(1)));
    [fitx,resx] = lsqcurvefit(gaussx,[60,max(pixels)/2],pixels,xlineout,[],[],options);
    [fity,resy] = lsqcurvefit(gaussy,[60,max(pixels)/2],pixels,ylineout,[],[],options);
    gaussx = gaussx(fitx,pixels);
    gaussy = gaussy(fity,pixels);
end

% Select individuals for crossover using fitness
% proportional selection (roulette wheel method)
function parent = selection(Individuals,Population,SumFitness)
    RouletteWheel = rand*SumFitness;
    Sum = 0;
    for i = 1:Individuals
        Sum = Sum + Population(i).Fitness;
        if Sum >= RouletteWheel
            parent = i;
            break
        end
    end
end
