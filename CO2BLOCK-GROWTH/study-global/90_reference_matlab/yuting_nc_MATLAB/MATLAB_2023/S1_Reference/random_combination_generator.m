% This script is used to find 1000 random combinations of growth rate (Rr)
% and target_s to calculate the corresponding peak year (Pr) and Storage
% resource required (Qr) 
% 

%% Testing data
% inputs/knowns for each Country 
start_year = 2030;% t2030
% current cumulative storage 
start_q = 0.12147; % change
 target_year = 2050; % t2050

% 4 variables: 
% range to search for equation to find peak year -> output 1 -->cannot be
% earlier than 2050 as it is outside the temporal limitation of this tool 
Pr = linspace(2050, 2350, 2000);
% range to search for equation storage resource --> output 2
Qr = [logspace(log10(1), log10(100000), 4000)];

% random growth rates - defined between 0.01 and 0.2 - input 1
Rr =  0.01 + (0.2-0.01)*rand(1,10000); 
%Rr = 0.1345; % single growth rate 


target_s = 0.001 + (1.28-0.001)*rand(1, 10000); % % change

% preallocate
min_peak = zeros(length(Qr),1); % minimum peak year found using the minimisation function between calculated and true data
qt_target = zeros(length(Qr),1);
% corresponding min storage
Qmin = zeros(length(Rr),1); % asscoiated minimum C to minimum peak year 
Prmin = zeros(length(Rr),1); 


for j=1:length(Rr)
    
    for i=1:length(Qr)
        
            for k=1:length(Pr)
            % calculate storage rate at 2030           
            p = ((Qr(i)-start_q)./(1+exp(Rr(j)*(Pr(k)- start_year))));
            % Evaluate fit
            fit_diff(k) = (p-start_q)^2; % minimum squared difference between calculated and Q(2030)  
            end
        [rowk] = find(fit_diff== min(fit_diff));
        min_peak(i) = Pr(rowk); % find the minium peak year
        
               %for h=1:length(target_s)
        % Find minium at give storage rate 
        qt_target(i) = (target_s(j) - ((Qr(i)-start_q).*Rr(j)*exp(Rr(j)*(Pr(rowk)-target_year))/(1+exp(Rr(j)*(Pr(rowk)-target_year)))^2))^2;
               %end
    end
    % Find minimum at given storage rate
    [ifit] = find(qt_target== min(qt_target));
    Qmin(j) = Qr(ifit);
    Prmin(j) =  min_peak(ifit);
    
end
% 

% recording the data in a table for the complete 1000 list - delete extra
% rows of data if necessary 
Growth = Rr.';
Storage = Qmin;
Peak = Prmin;
Target = target_s.';
T = table(Growth,Storage,Peak,Target);
todelete1 = T.Storage >= 404 ;% change
T(todelete1,:) = [];
size(T);
todelete2 = T.Storage == 1;% 
T(todelete2,:) = [];
size(T);
todelete3 = T.Peak >= 2100;
T(todelete3,:) = [];
todelete4 = T.Peak == 2050;
T(todelete4,:) = [];
size(T);
writetable(T, 'Canada2.txt') % change

% use this when you are running the code for a single target + growth rate
% combo 
%[ifit, jfit] = find(qt_target== min(min(qt_target)));
%peak_year = min_peak(ifit)
%storage_resource = Qr(ifit)
%storage_rate = Rr;