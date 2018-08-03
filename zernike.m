% This program creates circular plots of the Zernike polynomials.

%% Plot the results (for one)
% Initialise variables
x = linspace(1,-1,500);
[X,Y] = meshgrid(x,x);
[theta,r] = cart2pol(X,Y);
circ = r <= 1;
z = nan(size(X));

z(circ) = zernikepolynomials(40,0,r(circ),theta(circ));
figure('Units','normalized')
surf(x,x,z,'LineStyle','None'), view(2), axis off
set(gca,'XTick',[],'YTick',[])
title('Zernike function Z_5^1(r,\theta)')

%% Plot the results (if more than one)
% Initialise variables
x = linspace(1,-1,500);
[X,Y] = meshgrid(x,x);
[theta,r] = cart2pol(X,Y);
circ = r <= 1;
y = nan(size(X));

% Set the zernike functions to be found
n = [0  1  1  2  2  2  3  3  3  3  4  4  4  4  4];
m = [0 -1  1 -2  0  2 -3 -1  1  3 -4 -2  0  2  4];
Nplot = [5 13 15 21 23 25 29 31 33 35 37 39 41 43 45];

% Calculate the zernike functions
y = zernikepolynomials(n,m,r(circ),theta(circ));

figure('Units','normalized')
for i = 1:length(n)
    z(circ) = y(:,i);
    subplot(5,9,Nplot(i))
    surf(x,x,z,'LineStyle','None'), view(-60,45), axis square
    set(gca,'XTick',[],'YTick',[])
    title(['Z_{' num2str(n(i)) '}^{' num2str(m(i)) '}'])
    suptitle('Zernike polynomials (up to 4th order)')
end

%% Function for calculating the zernike polynomials on a unit sphere
function z = zernikepolynomials(n,m,r,theta)
    % Initialise variables
    n = n(:);
    m = m(:);
    r = r(:);
    theta = theta(:);
    length_r = length(r);
    m_abs = abs(m);                             % Find the absolute value of m
    z = zeros(length_r,length(n));             % Create an empty array for z
    rpowers = [];
    
    % Find powers of r
    for i = 1:length(n)
        rpowers = [rpowers m_abs(i):2:n(i)];          % sum of r.^(n - 2l)
    end
    rpowers = unique(rpowers);                        % Ignore repeated values of r
    
    % Calculate powers of r. Do this by raising each value in r_powers to the power of p.
    % If the 1st entry in r_powers is zero, skip this indes and add ones.
    if rpowers(1) == 0
        rpowern = arrayfun(@(p) r.^p,rpowers(2:end),'UniformOutput',false);
        rpowern = cat(2,rpowern{:});        % Concatenate into one array
        rpowern = [ones(length_r,1) rpowern];
    else
        rpowern = arrayfun(@(p) r.^p,rpowers,'UniformOutput',false);
        rpowern = cat(2,rpowern{:});                    % Concatenate into one array
    end
    
    % Calculate the polynomials
    for i = 1:length(n)
        l = 0:(n(i) - m_abs(i))/2;              % Calculate the values of l
        powers = n(i):-2:m_abs(i);              % Calculate powers in descending order
        for j = length(l):-1:1                  % Start with highest l (for factorials)
            % Calculate the polynomials and store them within R
            p = (1 - 2*mod(l(j),2))*...                        % (-1).^l
                prod(2:n(i) - l(j))/ ...                      % (n - l)!
                prod(2:l(j))/...                              % l!
                prod(2:(n(i) + m_abs(i))/2 - l(j))/...        % (n + m)/2 - 1
                prod(2:(n(i) - m_abs(i))/2 - l(j));           % (n - m)/2 - 1
            
            % circ = index(powers) if index(powers) = r_powers
            circ = (powers(j) == rpowers);
            % Calculate the values of the polynomials themselves, and add
            % each term in the series to z.
            z(:,i) = z(:,i) + p*rpowern(:,circ);
        end
    end
    
    % Choose whether to use cos or sin, dependent upon the sign of m
    circ_plus = m > 0;
    circ_minus = m < 0;
    
    if any(circ_plus)
        z(:,circ_plus) = z(:,circ_plus).*cos(theta*m_abs(circ_plus)');
    end
    if any(circ_minus)
        z(:,circ_minus) = z(:,circ_minus).*sin(theta*m_abs(circ_minus)');
    end
end