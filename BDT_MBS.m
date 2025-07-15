%% Assignement 2 FINE 452
% Matthieu Fouilloux, 261082501
% Giulia Gambaretto, 261103315

%% Part (i) a: Replicating the BDT excel example in Matlab
clear all

% Set parameters
T = 3; % End date, in years
N = 3; %Number of steps in the binomial tree
dt = T / N; %set time step size
m = (dt:dt:T)';

observedData = table(m, 'VariableNames',{'Maturity'});
% EXCEL REPLICATION
observedData.price = [0.9091, 0.8116, 0.7118]'; 
observedData.volatility = [NaN, 0.19, 0.18]';

%Create empty short rate tree. time i and node j
shortTree = NaN(N, N);

%Starting guesses for fminsearch
startR = 0.01; % Initial guess for the short rate
startParam = [0.01, 0.01]; % Initial guesses for rate and volatility


% Create first node in tree
i = 1;
thisPrice = observedData.price(i); % First price observation
thisMaturity = T - i; % First maturity observation

% No volatility matching for the first node
thisError = @(r) (thisPrice - (0.5*1+0.5*1)*exp(-r))^2; 

thisR = fminsearch(thisError, startR);

shortTree(i, 1) = thisR;


% Moving on in the tree

for i = 2:N
    i/N
    thisPrice = observedData.price(i);
    thisMaturity = T - i; 
    thisVolatility = observedData.volatility(i);

    thisPriceAndVol = [thisPrice, thisVolatility];

   
    thisError = @(param) sum((thisPriceAndVol - BondTree(param, i, ...
        shortTree)).^2)

    thisParam = fminsearch(thisError, startParam);

    % Build next step in tree with the right parameters
    thisBottomR = thisParam(1); % Short rate at bottom node
    thisSigma = thisParam(2); % Volatility at this step

    treeStep = exp(2 * thisSigma);
    shortTree(i, 1) = thisBottomR;
    for j = 2:(i)
        shortTree(i, j) = shortTree(i, j-1) * treeStep; 
    end
end
%% BondTree function

function PriceandVol = BondTree(param, maturity, shortTree)

    % Initialize output vector and bond price matrix
    PriceandVol = [NaN, NaN];
    BondPrice = NaN(maturity, maturity);

    % Construct the short rate tree
    shortTree(maturity, 1) = param(1);
    for i = 2:maturity
        shortTree(maturity, i) = shortTree(maturity, i-1) * ...
        exp(2 * param(2));
    end
    disp('Short Rate Tree:');
    disp(shortTree);

    % Initialize bond prices at maturity
    for t = 1:maturity
        % Discounted face value (assumes face value is 1)
        BondPrice(maturity, t) = exp(-shortTree(maturity, t)); 
    end

    % Calculate bond prices backward through the tree
    for i = maturity-1:-1:1
        for j = 1:i
            % Calculate bond price at each node using future cash
            % flows discounted by the short rate
            BondPrice(i, j) = (0.5 * BondPrice(i+1, j) + 0.5 ...
                * BondPrice(i+1, j+1)) * exp(-shortTree(i, j));
        end
    end

    % Display bond price matrix
    disp('Bond Price Matrix:');
    disp(BondPrice);

    % Calculate yields and volatility
    % Yield at the second period (node 1)
    yield_1 = -log(BondPrice(2, 1)) / (maturity-1); 
    % Yield at the second period (node 2)
    yield_2 = -log(BondPrice(2, 2)) / (maturity-1); 
    

    % Calculate the initial bond price at the root node
    Price = BondPrice(1, 1);

    % Calculate yield volatility
    vol = 0.5 * log(yield_2 / yield_1);

    % Return price and volatility
    PriceandVol = [Price, vol];

    % Optional: Display results for verification
    disp('Price and Volatility:');
    disp(PriceandVol);
    disp('BondPrice tree')
    disp(BondPrice);
    disp('shortTree')
    disp(shortTree)
end



%% Part (i) b: Estimating Yield and Volatility Term Structures

%import data for the treasury yield curve

yieldData = readtable('yieldAndVolatility2000.csv')

maturities = yieldData.m;  % Maturities from the data
observedYields = yieldData.y;  % Observed yields
observedVols = yieldData.vol; % Observed volatilities




%par vector contains, respectively, beta0, beta1, beta2, and tau
startPar=[1,1,1,1];


%Alias the error function with observed data and parameters.


nelsonSiegelYield=@(par, m) par(1) + (par(2) + par(3)) .* (1 - exp(-m ...
    ./ par(4))) ./ (m ./ par(4)) - par(3) .* exp(-m ./ par(4));
nelsonSiegelVol=@(par, m) par(1) + (par(2) + par(3)) .* (1 - exp(-m ...
    ./ par(4))) ./ (m ./ par(4)) - par(3) .* exp(-m ./ par(4));

% Error function for yield and volatiliry
errorNelsonSiegelYield = @(par) sum((observedYields - nelsonSiegelYield ...
    (par, maturities)).^2);
errorNelsonSiegelVol = @(par) sum((observedVols - nelsonSiegelVol ...
    (par, maturities)).^2);



%Run optimization
[parMinYield,errorMinYield] = fminsearch(errorNelsonSiegelYield, startPar);
[parMinVol,errorMinVol]   = fminsearch(errorNelsonSiegelVol, startPar);



%Alias new function fittedYield that takes as input some maturity m and 
%returns fitted yields at the estimated parameters
fittedFunctionYield= @(m) nelsonSiegelYield(parMinYield, m);
fittedFunctionVol  = @(m)  nelsonSiegelVol(parMinVol, m);


yieldData.nelsonSiegelYield= fittedFunctionYield(yieldData.m);
yieldData.nelsonSiegelVol= fittedFunctionVol(yieldData.m);

% Plot observed vs. fitted yields and volatilities

subplot(2,1,1);
plot(yieldData.m,yieldData.y,'*',yieldData.m,yieldData.nelsonSiegelYield);
legend('Observed Yields','Nelson Siegel Fitted Yields');
title('Observed vs. Fitted Yields');


subplot(2,1,2)
plot(yieldData.m,yieldData.vol,'*',yieldData.m,yieldData.nelsonSiegelVol)
legend('Observed Vols','Nelson Siegel Fitted Vols')
title('Observed vs. Fitted Volatilities')

%% Part (ii) a : Fitting to a Dynamic tree
% Now let's get 360 yield and volatility estimates

m360 = (1/12:1/12:30)'; % vector of maturities from 1 to 360

% calculate yields and volatilities for each maturities using fitted
% functions
estimatedYields360 = fittedFunctionYield(m360);
estimatedVols360 = fittedFunctionVol(m360);

% create table to visualize results
extendedYieldData = table(m360, estimatedYields360, estimatedVols360, ...
    'VariableNames', ...
    {'Maturity', 'Yield', 'Volatility'});
extendedYieldData.monthsteps = (1:1:360)';

N = 360
T = 30
dt = T/N

for i = 1:N
    extendedYieldData.prices(i) = 1/(1+(extendedYieldData.Yield(i)/2) ...
        )^(extendedYieldData.Maturity(i)*2);
end

shortTree=NaN(N,N);
startR = 0.01
startParam = [0.01,0.01]

% first node
i = 1;
thisPrice = extendedYieldData.prices(i);
thisMaturity = N - i
thisError = @(r) (extendedYieldData.prices(i) - ...
    (0.5*1+0.5*1)*exp(-r*(dt)))^2;
thisR=fminsearch(thisError,startR);
shortTree(i,1)=thisR;

% moving on in the rate tree

for i = 2:N
    i/N
    if i>= 3
        startParam = thisParam;
    end
    thisPrice = extendedYieldData.prices(i);
    thisMaturity = T - extendedYieldData.Maturity(i); 
    thisVolatility = extendedYieldData.Volatility(i);

    thisPriceAndVol = [thisPrice, thisVolatility];

   
    thisError = @(param) sum((thisPriceAndVol - ...
        BondTree2(param, i, shortTree, dt)).^2);

    thisParam = fminsearch(thisError, startParam);

    % Build next step in tree with the right parameters
    thisBottomR = thisParam(1); % Short rate at bottom node
    thisSigma = thisParam(2); % Volatility at this step


   treeStep = exp(2 * thisSigma * sqrt(dt));

    % Populate the short rate tree for this time step
    shortTree(i, 1) = thisBottomR;
    for j = 2:(i)
        shortTree(i, j) = shortTree(i, j-1) * treeStep;
        % Exponentially grow the rates upward in the tree
    end
end

%%
function PriceandVol = BondTree2(param, maturity, shortTree, dt)

    % Initialize output vector and bond price matrix
    PriceandVol = [NaN, NaN];
    BondPrice = NaN(maturity, maturity);

    % Construct the short rate tree
    shortTree(maturity, 1) = param(1);
    for i = 2:maturity
        shortTree(maturity, i) = shortTree(maturity, i-1) *...
            exp(2 * param(2)*sqrt(dt));
    end
    

    % Initialize bond prices at maturity
    for t = 1:maturity
        BondPrice(maturity, t) = exp(-shortTree(maturity, t)*dt); 
        % Discounted face value (assumes face value is 1)
    end

    % Calculate bond prices backward through the tree
    for i = maturity-1:-1:1
        for j = 1:i
            % Calculate bond price at each node using future cash
            % flows discounted by the short rate
            BondPrice(i, j) = (0.5 * BondPrice(i+1, j) + 0.5 * ...
                BondPrice(i+1, j+1)) * exp(-shortTree(i, j)*dt);
        end
    end

    % Calculate yields and volatility
    yield_1 = -log(BondPrice(2, 1)) / ((maturity-1)*(dt)); 
    yield_2 = -log(BondPrice(2, 2)) / ((maturity-1)*(dt)); 

    % Calculate the initial bond price at the root node
    Price = BondPrice(1, 1);


    % Calculate yield volatility
    vol = (0.5 * log(yield_2 / yield_1))/sqrt(dt);


    % Return price and volatility
    PriceandVol = [Price, vol];

end


%% part (ii) b : value a 30-year Mortgage-backed security

%% What is the fixed, monthly amount due?
% Fixed, monthly amount due
faceValue = 100000000; 
mortgageRate = 0.07; 
T = 30; 
N = 360;
dt = T/N
monthlyRate = exp(mortgageRate * dt) - 1; 

% Fixed monthly payment calculation
Pmt = (faceValue * monthlyRate) / (1 - (1 + monthlyRate)^(-N));

disp(['The fixed monthly payment is: $', num2str(Pmt, '%.2f')]);

%% After t monthly payments have been made, 
% what is the outstanding principal?
format bank;
maturity = 360
remainingBalance = faceValue;
remainingPrincipal = zeros(N, 1); 
remainingPrincipal(1) = faceValue;

for t = 2:maturity+1
    interest = remainingBalance * monthlyRate;
    principal = Pmt - interest;
    remainingBalance = remainingBalance - principal;
    remainingPrincipal(t) = remainingBalance;
end

disp(remainingPrincipal)
format short;



%% What is the value of the mortgage if borrowers refinance optimally? 
maturity = 360;
format bank;
MortgageTreeOptimal = NaN(maturity+1, maturity+1);
ContinuationValueOptimal = NaN(maturity, maturity);
OutstandingPrincipalOptimal = NaN(maturity+1, maturity+1);

MortgageTreeOptimal(maturity+1, :) = Pmt;

for i = maturity+1:-1:1
    for j = 1:i
        OutstandingPrincipalOptimal(i,j) = remainingPrincipal(i);
    end
end

for i = maturity:-1:1
    for j = 1:i
        if i == 1
            ContinuationValueOptimal(i,j) = (0.5 * MortgageTreeOptimal ...
                (i+1, j) + 0.5 * MortgageTreeOptimal(i+1, j+1)) * ...
            exp(-shortTree(i,j)*dt);
            MortgageTreeOptimal(i,j) = min(ContinuationValueOptimal(i,j),...
                OutstandingPrincipalOptimal(i,j));
        else
            ContinuationValueOptimal(i,j) = (0.5 * MortgageTreeOptimal ...
                (i+1, j) + 0.5 * MortgageTreeOptimal(i+1, j+1)) *...
            exp(-shortTree(i,j)*dt);
            MortgageTreeOptimal(i,j) = min(ContinuationValueOptimal ...
                (i,j), OutstandingPrincipalOptimal(i,j)) + Pmt;
        end
    end
end



mortgageValueOptimal = MortgageTreeOptimal(1,1);
disp(['The value of the mortgage with optimal refinancing is: $', ...
    num2str(mortgageValueOptimal, '%.2f')]);

format short;

%% What is the value of the mortgage if borrowers never refinance? 

maturity = 360;
format bank
MortgageTreeNonOptimal = NaN(maturity+1, maturity+1);
MortgageTreeNonOptimal(maturity+1, :) = Pmt;

for i = maturity:-1:1
    for j = 1:i
        if i == 1
            MortgageTreeNonOptimal(i,j) = (0.5 * MortgageTreeNonOptimal ...
                (i+1, j) + 0.5 * MortgageTreeNonOptimal(i+1, j+1)) * ...
            exp(-shortTree(i,j)*dt);
        else
            MortgageTreeNonOptimal(i,j) = (0.5 * MortgageTreeNonOptimal ...
                (i+1, j) + 0.5 * MortgageTreeNonOptimal(i+1, j+1))...
            * exp(-shortTree(i,j)*dt) + Pmt;
        end
    end
end

mortgageValueNonOptimal = MortgageTreeNonOptimal(1,1);
disp(['The value of the mortgage with non optimal refinancing is: $', ...
    num2str(mortgageValueNonOptimal, '%.2f')]);

format short;


%% Duration and convexity


factor = 0.4:0.02:2; % Range of factors

for k = 1:length(factor)
    new_short_tree = shortTree*factor(k);
    
    
    for i = maturity:-1:1 
        for j = 1:i
            if i == 1
                ContinuationValueOptimal(i,j) = (0.5 * ...
                    MortgageTreeOptimal (i+1, j) + 0.5 * ...
                    MortgageTreeOptimal(i+1, j+1)) * ...
                exp(-new_short_tree(i,j)*dt);
                MortgageTreeOptimal(i,j) = min(ContinuationValueOptimal ...
                    (i,j), OutstandingPrincipalOptimal(i,j));
            else
                ContinuationValueOptimal(i,j) = (0.5 * ...
                    MortgageTreeOptimal(i+1, j) + 0.5 * ...
                    MortgageTreeOptimal(i+1, j+1)) * ...
                exp(-new_short_tree(i,j)*dt);
                MortgageTreeOptimal(i,j) = min(ContinuationValueOptimal ...
                    (i,j), OutstandingPrincipalOptimal(i,j)) + Pmt;
            end
        end
    mortgageValueOptimal(k) = MortgageTreeOptimal(1,1);
        
    DurationOptimal(k) = (-1 / MortgageTreeOptimal(1,1)) * ...
    ((MortgageTreeOptimal(2,2) - MortgageTreeOptimal(2,1)) / ...
    ((new_short_tree(2,2) - new_short_tree(2,1))));
    DurationUpOptimal(k) = (-1 / MortgageTreeOptimal(2,2)) *...
    ((MortgageTreeOptimal(3,3) - MortgageTreeOptimal(3,2)) / ...
    ((new_short_tree(3,3) - new_short_tree(3,2))));
    DurationDownOptimal(k) = (-1 / MortgageTreeOptimal(2,1)) * ...
    ((MortgageTreeOptimal(3,2) - MortgageTreeOptimal(3,1)) / ...
    ((new_short_tree(3,2) - new_short_tree(3,1))));

    ConvexityOptimal(k) = DurationOptimal(k)^2 - ...
    (DurationUpOptimal(k) - DurationDownOptimal(k))/...
    ((new_short_tree(2, 2) - new_short_tree(2, 1)));

    end
    

    for i = maturity:-1:1 
        for j = 1:i
            if i == 1
                MortgageTreeNonOptimal(i,j) = (0.5 * ...
                    MortgageTreeNonOptimal(i+1, j) + 0.5 * ...
                    MortgageTreeNonOptimal(i+1, j+1)) *...
                exp(-new_short_tree(i,j)*dt);
            else
                MortgageTreeNonOptimal(i,j) = (0.5 * ...
                    MortgageTreeNonOptimal(i+1, j) + 0.5 * ...
                    MortgageTreeNonOptimal(i+1, j+1)) * ...
                exp(-new_short_tree(i,j)*dt) + Pmt;
            end
        end
    end

    mortgageValueNonOptimal(k) = MortgageTreeNonOptimal(1,1);
    DurationNonOptimal(k) = (-1 / MortgageTreeNonOptimal(1,1)) *...
        ((MortgageTreeNonOptimal(2,2) - MortgageTreeNonOptimal(2,1)) / ...
        ((new_short_tree(2,2) - new_short_tree(2,1))));
    DurationUpNonOptimal(k) = (-1 / MortgageTreeNonOptimal(2,2)) *...
    ((MortgageTreeNonOptimal(3,3) - MortgageTreeNonOptimal(3,2)) / ...
    ((new_short_tree(3,3) - new_short_tree(3,2))));
    DurationDownNonOptimal(k) = (-1 / MortgageTreeNonOptimal(2,1)) *...
    ((MortgageTreeNonOptimal(3,2) - MortgageTreeNonOptimal(3,1)) / ...
    ((new_short_tree(3,2) - new_short_tree(3,1))));
    ConvexityNonOptimal(k) = DurationNonOptimal(k)^2 - ...
    (DurationUpNonOptimal(k) - DurationDownNonOptimal(k))/...
    ((new_short_tree(2, 2) - new_short_tree(2, 1)));
    interest_rate(k) = -new_short_tree(1,1)

end



%% Report : Create a figure with 5 subplots

figure;

% Plot 1: Fitted Yield Curve and Observed Yields
subplot(2, 3, 1);
plot(yieldData.m,yieldData.y,'*', 'DisplayName', 'Observed Yields')
hold on;
plot(yieldData.m,yieldData.nelsonSiegelYield, '-', 'DisplayName', ...
    'Fitted Yield Curve');
hold off;
grid on;
title('Yields');
legend('Observed Yields','Nelson Siegel Fitted Yields', ...
    'FontName', 'Times New Roman');
xlabel('Years', 'FontName', 'Times New Roman');
ylabel('Yields', 'FontName', 'Times New Roman');
set(gca, 'FontName', 'Times New Roman');


% Plot 2: Fitted Volatility Curve and Observed Volatilities
subplot(2, 3, 2);
plot(yieldData.m, yieldData.vol, '*', 'DisplayName', ...
    'Observed Volatilities');
hold on;
plot(yieldData.m, yieldData.nelsonSiegelVol, '-', 'DisplayName', ...
    'Fitted Volatility Curve');
hold off;
title('Volatilities');
legend('Observed Vols','Nelson Siegel Fitted Vols', 'FontName', ...
    'Times New Roman');
xlabel('Years', 'FontName', 'Times New Roman');
ylabel('Volatility', 'FontName', 'Times New Roman');
set(gca, 'FontName', 'Times New Roman');


% Plot 3: Values of prepayable and non prepayable mortgage pools
subplot(2, 3, 3);
plot(mortgageValueOptimal, 'g', 'DisplayName', 'Prepayable');
hold on;
plot(mortgageValueNonOptimal, 'r', 'DisplayName', 'Non-Prepayable');
hold off;
title('Mortgage Values');
legend('Mortgage Value with Prepayment', ...
    'Mortgage Value with No Prepayment', 'FontName', 'Times New Roman');
xticks(0:9:length(factor));
xticklabels(string([round(-interest_rate(1:9:length(factor)),4), ...
    round(-interest_rate(length(factor)),4)]));
xlabel('Interest Rates', 'FontName', 'Times New Roman');
ylabel('Mortgage Value', 'FontName', 'Times New Roman');
set(gca, 'FontName', 'Times New Roman');


% Plot 4: Durations of prepayable and non prepayable mortgage pools
subplot(2, 3, 4);
plot(DurationOptimal, 'g', 'DisplayName', 'Prepayable');
hold on;
plot(DurationNonOptimal, 'r', 'DisplayName', 'Non-Prepayable');
hold off;
title('Durations as a Function of Interest Rates');
legend('Duration with Prepayment','Duration with No Prepayment', ...
    'FontName', 'Times New Roman');
xticks(0:9:length(factor));
xticklabels(string([round(-interest_rate(1:9:length(factor)),4), ...
    round(-interest_rate(length(factor)),4)]));
xlabel('Interest Rates', 'FontName', 'Times New Roman');
ylabel('Years', 'FontName', 'Times New Roman');
set(gca, 'FontName', 'Times New Roman');


% Plot 5: Convexities of prepayable and non prepayable mortgage pools
subplot(2, 3, 5);
plot(ConvexityOptimal, 'g', 'DisplayName', 'Prepayable');
hold on;
plot(ConvexityNonOptimal, 'r', 'DisplayName', 'Non-Prepayable');
hold off;
title('Convexities as a Function of Interest Rates');
legend('Convexity with Prepayment','Convexity with No Prepayment', ...
    'FontName', 'Times New Roman');
xticks(0:9:length(factor));
xticklabels(string([round(-interest_rate(1:9:length(factor)),4), ...
    round(-interest_rate(length(factor)),4)]));
xlabel('Interest Rates', 'FontName', 'Times New Roman');
set(gca, 'FontName', 'Times New Roman');
