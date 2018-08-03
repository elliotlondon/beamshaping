% ---------------- SET FIXED PARAMETERS ---------------- %
close; clear; clc

% Set up variables and algorithmic parameters
options = optimset('Display','off');
NumberOfUnknowns = 2;
StringLength = 8;
TotalStringLength = StringLength*NumberOfUnknowns;
Individuals = 25;
Generations = 50;
CrossOverProbability = 0.7; 
MutationProbability = 0.01;
Elitism = 1;
GridSize = 257;

if Elitism == 1
    EliteNumber = ceil(0.1*Individuals);
else
    EliteNumber = 0;
end

% Set up test function (waffle grid). This function has a maximum of 1,
% which the genetic algorithm aims to find via guessing a random point
% within the array and then iteratively improving until it finds the
% maximum.
[gridx,gridy] = meshgrid(linspace(-1,1,GridSize));
waffle = exp(-(gridx.^2 + gridy.^2)/1.5).*cos(pi*(gridx)*3.5);
waffle = waffle.*transpose(waffle);

h = figure(1);
surf(waffle,'LineStyle','None'),view(2)
xlim([1 GridSize]), ylim([1 GridSize])
xlabel('x'), ylabel('y'), zlabel('f(x,y)')

%% ---------------- GA MAIN PROGRAM ---------------- %
% Initialise the fields of the structs and variables
Population(Individuals).Strings = [];
Population(Individuals).Individual = [];
Population(NumberOfUnknowns).Fitness = [];
General(Generations).People(Individuals).Individual = [];
NewStrings(Individuals).Strings = [];
MeanFitness = 1:Generations;
FittestIndividual = 1:Generations;
FitnessStdev = 1:Generations;
successgen = 0;

% First generation
for gen = 1
    % Create an initial population of individuals
    for indiv = 1:Individuals
        Population(indiv).Strings = binarystring(StringLength,NumberOfUnknowns);
        Population(indiv).Individual = fastbin2dec(Population(indiv).Strings);
        
        % Find the fitness of each individual
        Population(indiv).Fitness = (waffle(Population(indiv).Individual(1),Population(indiv).Individual(2)) + 1)/2;
        
        % Show individuals on grid
        h = figure(2);
        titlesprint = sprintf('Generation = %i/%i',1,max(Generations));
        title(titlesprint)
        contour(waffle),view(2)
        hold on
        plot(Population(indiv).Individual(1),Population(indiv).Individual(2),'*','Color','b')
        
        General(gen).People(indiv).Individual = Population(indiv).Individual;
    end
    
    % Find fitness statistics
    MeanFitness(gen) = mean([Population.Fitness]);
    SumFitness = sum([Population.Fitness]);
    FittestIndividual(gen) = max([Population.Fitness]);
    FitnessStdev(gen) = std([Population.Fitness]);
    
    % Display statistics
    clc
    fprintf('Generation = %i/%i\n',1,max(Generations));
    fprintf('MeanFitness = %f\n',MeanFitness(gen));
    fprintf('MaxFitness = %f\n',FittestIndividual(gen));
    fprintf('SumFitness = %f\n',SumFitness);
end

clf

% All other generations
for gen = 2:Generations
    if Elitism == 1
        % Perform elitism (always pass best individual into next generation)
        for indiv = 1:EliteNumber
            % Sort array so that best individuals are at the top
            [~,index] = sortrows([Population.Fitness].');
            Population = Population(index(end:-1:1)); 
            
            % Pass elite individuals through
            NewStrings(indiv).Strings = Population(indiv).Strings;
            
            Population = Population(index);
        end
    end
        for indiv = EliteNumber+1:2:Individuals
            % Select individuals for crossover
            Parent1 = selection(Individuals,Population,SumFitness);
            Parent2 = selection(Individuals,Population,SumFitness);

            % Perform crossover and generate a temporary population
            if rand <= CrossOverProbability
                [NewStrings(indiv).Strings,NewStrings(indiv+1).Strings] = ...
                crossover(Population,Parent1,Parent2,StringLength,NumberOfUnknowns);
                for j = 1:NumberOfUnknowns
                    for k = 1:StringLength
                        % Mutate the temporary population
                        if rand <= MutationProbability
                            if NewStrings(indiv).Strings(j,k) == 1
                                NewStrings(indiv).Strings(j,k) = 0;
                                fprintf('Mutation occured...\n')
                            else
                                NewStrings(indiv).Strings(j,k) = 1;
                                fprintf('Mutation occured...\n')
                            end
                        end
                    end
                end
            else
                NewStrings(indiv).Strings = nocrossover(Population,Parent1);
                NewStrings(indiv+1).Strings = nocrossover(Population,Parent2);
                for j = 1:NumberOfUnknowns
                    for k = 1:StringLength
                        % Mutate the temporary population
                        if rand <= MutationProbability
                            if NewStrings(indiv).Strings(j,k) == 1
                                NewStrings(indiv).Strings(j,k) = 0;
                                fprintf('Mutation occured...\n')
                            else
                                NewStrings(indiv).Strings(j,k) = 1;
                                fprintf('Mutation occured...\n')
                            end
                        end
                    end
                end
            end
        end
    
    % Transfer the temporary population into the new population
    for indiv = 1:Individuals
        Population(indiv).Strings = NewStrings(indiv).Strings;
        Population(indiv).Individual = fastbin2dec(Population(indiv).Strings);
        
        % Find the fitness of each individual
        Population(indiv).Fitness = (waffle(Population(indiv).Individual(1),Population(indiv).Individual(2)) + 1)/2;
        
        % Show individuals on grid
        h = figure(2);
        titlesprint = sprintf('Generation = %i/%i',gen,max(Generations));
        title(titlesprint)
        contour(waffle),view(2)
        hold on
        plot(Population(indiv).Individual(1),Population(indiv).Individual(2),'*','Color','b')
        
        General(gen).People(indiv).Individual = Population(indiv).Individual;
    end
    
    % Find fitness statistics
    MeanFitness(gen) = mean([Population.Fitness]);
    SumFitness = sum([Population.Fitness]);
    [FittestIndividual(gen),MaxFitnessIndex] = max([Population.Fitness]);
    FitnessStdev(gen) = std([Population.Fitness]);
    
    % Display statistics
    clc
    fprintf('Generation = %i/%i\n',gen,max(Generations));
    fprintf('MeanFitness = %f\n',MeanFitness(gen));
    fprintf('MaxFitness = %f\n',FittestIndividual(gen));
    fprintf('SumFitness = %f\n',SumFitness);
    
    % Break if the algorithm is entirely successful
    if successgen == 0
        if FittestIndividual(gen) == 1
            successgen = gen;
        end
    end
    
    if gen < (Generations - 1)
        clf
    end
end

% Plot evolution of mean fitness
h = figure(3);
plot(1:Generations,MeanFitness)
hold on
plot(1:Generations, FittestIndividual)
xlabel('Generation'), ylabel('Fitness')
xlim([1 gen]), ylim([0 1])
set(gca,'FontSize',20)
legend('Average Fitness', 'Maximum Fitness')

% Save figure
saveas(h,sprintf('fitness.png'))

% Print a seeded message upon completion
fprintf('Algorithm finished. Checkrand = %f\n',rand);

%% Plotting

% Put contour maps in subplots
subplot(2,2,1)
for i = 1:Individuals
    title('Generation = 1')
    contour(waffle),view(2)
    hold on
    plot(General(1).People(i).Individual(1),General(1).People(i).Individual(2),'*','Color','b')
    xlabel('x'), ylabel('y')
end

subplot(2,2,2)
for i = 1:Individuals
    title('Generation = 10')
    contour(waffle),view(2)
    hold on
    plot(General(10).People(i).Individual(1),General(10).People(i).Individual(2),'*','Color','b')
    xlabel('x'), ylabel('y')
end

subplot(2,2,3)
for i = 1:Individuals
    title('Generation = 20')
    contour(waffle),view(2)
    hold on
    plot(General(20).People(i).Individual(1),General(20).People(i).Individual(2),'*','Color','b')
    xlabel('x'), ylabel('y')
end

subplot(2,2,4)
for i = 1:Individuals
    title('Generation = 30')
    contour(waffle),view(2)
    hold on
    plot(General(30).People(i).Individual(1),General(30).People(i).Individual(2),'*','Color','b')
    xlabel('x'), ylabel('y')
end

%% Functions and some brief descriptions
% Create a binary string representing an individual 
function string = binarystring(StringLength,NumberOfUnknowns)
    string = round(ones(NumberOfUnknowns,StringLength).*rand(NumberOfUnknowns,StringLength));
end

% Turns binary strings into decimal numbers. MATLAB's built-in 
% 'bin2dec' is capable of this, but is too slow.
function decimalout = fastbin2dec(binaryin)
    decimalout = binaryin*(2.^(size(binaryin, 2)-1:-1:0)).' + 1;
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

% Perform crossover between two parents to create a new individual
function [newindiv1,newindiv2] = crossover(Population,Parent1,Parent2,StringLength,NumberOfUnknowns)
    % Generate two random locations in the binary string for crossover
    CrossSite = round((StringLength*NumberOfUnknowns - 1)*rand + 1);
    % Grab the strings of the parents
    selected1 = horzcat(Population(Parent1).Strings(1,:),Population(Parent1).Strings(2,:));
    selected2 = horzcat(Population(Parent2).Strings(1,:),Population(Parent2).Strings(2,:));
    
    % Perform crossover
    c1 = horzcat(selected1(1:CrossSite),selected2(CrossSite+1:end));
    c2 = horzcat(selected2(1:CrossSite),selected1(CrossSite+1:end));
    
    % Change to correct array dimensions
    newindiv1(1,1:StringLength) = c1(1:StringLength);
    newindiv1(2,1:StringLength) = c1(StringLength+1:StringLength*NumberOfUnknowns);
    newindiv2(1,1:StringLength) = c2(1:StringLength);
    newindiv2(2,1:StringLength) = c2(StringLength+1:StringLength*NumberOfUnknowns);
end

% Transfer individuals in the case of no crossover
function newindiv1 = nocrossover(Population,Parent1)
    newindiv1 = Population(Parent1).Strings;
end