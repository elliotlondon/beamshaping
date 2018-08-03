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
Individuals = 30;
Generations = 150;
CrossOverProbability = 0.8;
InitialMutationRate = 0.005;
MutationRate = InitialMutationRate;
UpperRange = 0.05;
LowerRange = -0.05;
GaussNumber = 4;
WidthGuess = 12;
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
ccdmax = max(max(ccddata(:,:,2)));

% Grab image from the camera
ccddata = im2double(step(vid));

% Create a 2D supergaussian function
[X,Y] = meshgrid(pixels,pixels);
X = X - max(pixels)/2;
Y = Y - max(pixels)/2;
sgauss = ccdmax*exp(-sqrt((X.^2 + Y.^2)/WidthGuess.^2).^(2*GaussNumber));

% Create a square 2D supergaussian function
[X,Y] = meshgrid(pixels,pixels);
X = X - max(pixels)/2;
Y = Y - max(pixels)/2;
ssgauss = ccdmax*exp(-(X.^(GaussNumber) + Y.^(GaussNumber))...
                     /WidthGuess.^(GaussNumber));

% Calculate euclidean distance between sgauss function and data
eudist = sqrt((sgauss - ccddata(:,:,2)).^2);

% Plot the euclidean distance
figure(2)
surf(eudist,'LineStyle','None'), view(2)

% Plot initial beam profile
figure(3)
surf(ccddata(:,:,2),'LineStyle','None'), view(2)
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
        ccdmax = max(max(ccddata(:,:,2)));
        
        % Create a supergaussian
        [X,Y] = meshgrid(pixels,pixels);
        X = X - max(pixels)/2;
        Y = Y - max(pixels)/2;
        sgauss = ccdmax*exp(-sqrt((X.^2 + Y.^2)/WidthGuess.^2).^(2*GaussNumber));
        
        % Calculate euclidean distance between sgauss function and data
        eudist = sqrt((sgauss - ccddata(:,:,2)).^2);
        eudist(eudist<0.01) = 0;
        eusum = sum(sum(eudist))/max(pixels);
        
        % Create a fitness function
        if eusum <= 1
            Population(indiv).Fitness = 1 - eusum;
        else
            % If < 0, individual is unsuccessful, but still has a chance
            % at crossover (to ensure inclusion of mutations).
            Population(indiv).Fitness = 0.05;
        end
        
        % Display data for every individual
        clc
        fprintf('Individual = %i/%i, ',indiv,max(Individuals));
        fprintf('Generation = %i/%i\n',gen,max(Generations));
        fprintf('eusum = %.3f\n',eusum);
        if max(max(ccddata(:,:,2))) >= 0.98
            fprintf('Warning! Camera saturated!')
        end
        
        % Plot the euclidean distance
        yfig = figure(2);
        set(yfig,'Position',[450 550 400 400])
        surf(eudist,'LineStyle','None'), view(2)
        % Plot beam profile in 2D
        zfig = figure(3);
        set(zfig,'Position',[850 550 400 400])
        surf(ccddata(:,:,2),'LineStyle','None'), view(2)
    end
    
    % Find fitness statistics
    MeanFitness(gen) = mean([Population.Fitness]);
    SumFitness = sum([Population.Fitness]);
    [MaxFitness(gen),MaxFitnessIndex] = max([Population.Fitness]);
    FitnessStdev(gen) = std([Population.Fitness]);
    BestIndividual = Population(MaxFitnessIndex).Individual;
    
    % Display data for every generation
    clc
    fprintf('Individual = %i/%i, ',indiv,max(Individuals));
    fprintf('Generation = %i/%i\n',gen,max(Generations));
    fprintf('MeanFitness = %f\n',MeanFitness(gen));
    fprintf('MaxFitness = %f\n',MaxFitness(gen));
end

% All other generations
for gen = 2:Generations
    if Elitism == 1
        % Apply Elitism, automatically sending the best individual through
        [~,index] = sortrows([Population.Fitness].');
        Population = Population(index(end:-1:1));

        for indiv = 1:EliteNumber
            NewStrings(indiv).Strings = Population(indiv).Strings;
            
            Population = Population(index);
        end
    end
    
    % Generate offspring for the rest of the generation
    % Go in steps of two, as each crossover generates two children.
    for indiv = (EliteNumber + 1):2:Individuals    
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
                if rand >= 0.5
                    NewStrings(indiv).Strings(j,:) = horzcat(Selected1(1:CrossSite1),Selected2(CrossSite1+1:StringLength));
                else
                    NewStrings(indiv).Strings(j,:) = horzcat(Selected2(1:CrossSite2),Selected1(CrossSite2+1:StringLength));
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
        % No crossover case
        else
            % Do a tournament to choose the new individual
            if Population(Parent1).Fitness > Population(Parent2).Fitness
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
    
    if Elitism == 1
        % Transfer elite population without recalculating fitness
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
            fprintf('MaxFitness = %f, ',MaxFitness(gen-1));
            fprintf('MeanFitness = %f\n',MeanFitness(gen-1));
            fprintf('StallCounter = %i, ',StallCounter);
            fprintf('MutationRate = %.3f\n',MutationRate);
            if max(max(ccddata(:,:,2))) >= 0.98
                fprintf('Warning! Camera saturated!')
                pause(0.5)
            end
            
            % Plot the euclidean distance
            yfig = figure(2);
            set(yfig,'Position',[450 550 400 400])
            surf(eudist,'LineStyle','None'), view(2)
            % Plot beam profile in 2D
            zfig = figure(3);
            set(zfig,'Position',[850 550 400 400])
            surf(ccddata(:,:,2),'LineStyle','None'), view(2)
        end
    end
    
    % Transfer temporary population to real population
    for indiv = (EliteNumber + 1):Individuals
        Population(indiv).Strings = NewStrings(indiv).Strings;
        Population(indiv).Individual = fastbin2dec(Population(indiv).Strings,UpperRange,LowerRange,StringLength);
        
        % Send surface to mirror
        dmdata(:,1) = Population(indiv).Individual(:);
        dm.Send(dmdata);
        pause(0.25)
        
        % Grab image from the camera
        ccddata = im2double(step(vid));
        ccdmax = max(max(ccddata(:,:,2)));
        
        % Create a supergaussian
        [X,Y] = meshgrid(pixels,pixels);
        X = X - max(pixels)/2;
        Y = Y - max(pixels)/2;
        sgauss = ccdmax*exp(-sqrt((X.^2 + Y.^2)/WidthGuess.^2).^(2*GaussNumber));
        
        % Calculate euclidean distance between sgauss function and data
        eudist = sqrt((sgauss - ccddata(:,:,2)).^2);
        eudist(eudist<0.01) = 0;
        eusum = sum(sum(eudist))/max(pixels);
        
        % Create a fitness function
        if eusum <= 1
            Population(indiv).Fitness = 1 - eusum;
        else
            % If < 0, individual is unsuccessful, but still has a chance
            % at crossover (to ensure inclusion of mutations).
            Population(indiv).Fitness = 0.05;
        end
        
        % Display data every loop
        clc
        fprintf('Individual = %i/%i, ',indiv,max(Individuals));
        fprintf('Generation = %i/%i\n',gen,max(Generations));
        fprintf('MaxFitness = %f, ',MaxFitness(gen-1));
        fprintf('MeanFitness = %f\n',MeanFitness(gen-1));
        fprintf('StallCounter = %i, ',StallCounter);
        fprintf('MutationRate = %.3f\n',MutationRate);
        if max(max(ccddata(:,:,2))) >= 0.98
            fprintf('Warning! Camera saturated!')
            pause(0.5)
        end

        % Plot the euclidean distance
        yfig = figure(2);
        set(yfig,'Position',[450 550 400 400])
        surf(eudist,'LineStyle','None'), view(2)
        % Plot beam profile in 2D
        zfig = figure(3);
        set(zfig,'Position',[850 550 400 400])
        surf(ccddata(:,:,2),'LineStyle','None'), view(2)
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
    fprintf('MaxFitness = %f, ',MaxFitness(gen));
    fprintf('MeanFitness = %f\n',MeanFitness(gen));
    fprintf('StallCounter = %i, ',StallCounter);
    fprintf('MutationRate = %.3f\n',MutationRate);
    
    % Apply evolutionary pressure
    if (StallCounter == StallTolerance) && (Pressure == 1)
        if MutationRate <= 0.25
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
ccdmax = max(max(ccddata(:,:,2)));

% Find the lineouts
[xlineout,ylineout] = lineouts(pixels,ccddata);

% Generate gaussian fits for the image lineouts
[gaussx,gaussy,resx,resy] = supergaussfit(pixels,xlineout,ylineout,WidthGuess,GaussNumber);
[gaussx2,gaussy2] = gaussfit(pixels,xlineout,ylineout,WidthGuess);

% Create a supergaussian
[X,Y] = meshgrid(pixels,pixels);
X = X - max(pixels)/2;
Y = Y - max(pixels)/2;
sgauss = ccdmax*exp(-sqrt((X.^2 + Y.^2)/WidthGuess.^2).^(2*GaussNumber));

% Calculate euclidean distance between sgauss function and data
ccdnoiseless = ccddata(:,:,2);
ccdnoiseless(ccdnoiseless<0.001) = 0;
eudist = sqrt((sgauss - ccdnoiseless).^2);
eusum = sum(sum(eudist))/max(pixels);

% Save variables for later use
save('C:\Users\Kevin\Downloads\Elliot\MRes_Dropbox\Experimental_Data_(Elliot)\DMgasuper\CorrectedSurface.mat','dmdata')

% Display data for the final generation
clc
fprintf('Individual = %i/%i, ',indiv,max(Individuals));
fprintf('Generation = %i/%i\n',gen,max(Generations));
fprintf('X resnorm = %.3f, Y resnorm = %.3f\n',resx,resy);

% Plot the x data
xfig = figure(1);
set(xfig,'Position',[50 550 400 400])
plot(pixels,xlineout,pixels,gaussx,pixels,gaussx2)
set(gca,'FontSize',20)
xlim([0 max(pixels)]), ylim([0 1])
xlabel('Pixels'), ylabel('Normalised integrated count rate')
legend('Raw x data','4th order supergaussian fit','Gaussian fit')
% Plot the y data
yfig = figure(2);
set(yfig,'Position',[450 550 400 400])
plot(pixels,ylineout,pixels,gaussy,pixels,gaussy2)
set(gca,'FontSize',20)
xlim([0 max(pixels)]), ylim([0 1])
xlabel('Pixels'), ylabel('Normalised integrated count rate')
legend('Raw y data','4th order supergaussian fit','Gaussian fit')
% Plot beam profile in 2D
zfig = figure(3);
set(zfig,'Position',[850 550 400 400])
surf(pixels,pixels,ccddata(:,:,2),'LineStyle','None'), view(2)
% Plot evolution of mean fitness
figure(4);
plot(1:max(Generations),MeanFitness,1:max(Generations),MaxFitness)
set(gca,'FontSize',20)
xlabel('Generation'), ylabel('Fitness')
xlim([1 gen]), ylim([0 1])
legend('Average Fitness','Maximum fitness')
% Plot the difference between the euclidean distance and the measured data
figure(5)
surf(eudist,'LineStyle','None'), view(2), colorbar
xlabel('x'), ylabel('y');

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

% Find the lineouts of the image
function [xlineout,ylineout] = lineouts(pixels,ccddata)
    xlineout = ccddata(max(pixels)/2,:,2);
    ylineout = ccddata(:,max(pixels)/2,2);
    % Normalise the lineouts
    xlineout = xlineout./max(xlineout);
    ylineout = transpose(ylineout./max(ylineout));
end

% Generate gaussian fits to the image lineouts
function [gaussx,gaussy] = gaussfit(pixels,xlineout,ylineout,WidthGuess)
    % Set options
    options = optimset('Display','off');
    
    % Generate gaussian fits
    gaussx = @(a,xlineout) exp(-2*((xlineout - a(1))/(a(2))).^2);
    gaussy = @(a,ylineout) exp(-2*((ylineout - a(1))/(a(2))).^2);
    fitx = lsqcurvefit(gaussx,[max(pixels)/2,WidthGuess],pixels,xlineout,[],[],options);
    fity = lsqcurvefit(gaussy,[max(pixels)/2,WidthGuess],pixels,ylineout,[],[],options);
    gaussx = gaussx(fitx,pixels);
    gaussy = gaussy(fity,pixels);
end

% Generate supergaussian fits to the image lineouts
function [gaussx,gaussy,resx,resy] = supergaussfit(pixels,xlineout,ylineout,WidthGuess,GaussNumber)
    % Set options
    options = optimset('Display','off');
    
    % Generate gaussian fits
    gaussx = @(a,xlineout) exp(-2*((xlineout - a(1))/a(2)).^GaussNumber);
    gaussy = @(a,ylineout) exp(-2*((ylineout - a(1))/a(2)).^GaussNumber);
    [fitx,resx] = lsqcurvefit(gaussx,[max(pixels)/2,WidthGuess],pixels,xlineout,[],[],options);
    [fity,resy] = lsqcurvefit(gaussy,[max(pixels)/2,WidthGuess],pixels,ylineout,[],[],options);
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
