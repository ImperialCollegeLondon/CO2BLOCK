%This code concatenate the initial master table of more than 

load('US.txt')
years = US(:,1); % grwoth rate
Q = US(:,2); % MT - sotrage resource  
Pr = US(:,3); % Gt - peak year 
target = US(:,4);% target 

load('China.txt')
years2 = China(:,1); % years
Q2 = China(:,2); % MT - storage rate 
Pr2 = China(:,3); % Gt - cumulative storage 
target2 = China(:,4);

load('EU.txt')
years3 = EU(:,1); % years
Q3 = EU(:,2); % MT - storage rate 
Pr3 = EU(:,3); % Gt - cumulative storage 
target3 = EU(:,4);

load('UK.txt')
years4 = UK(:,1); % years
Q4 = UK(:,2); % MT - storage rate 
Pr4 = UK(:,3); % Gt - cumulative storage 
target4 = UK(:,4);

load('Thailand.txt')
years5 = Thailand(:,1); % years
Q5 = Thailand(:,2); % MT - storage rate 
Pr5 = Thailand(:,3); % Gt - cumulative storage 
target5 = Thailand(:,4);

load('Canada.txt')
years6 = Canada(:,1); % years
Q6 = Canada(:,2); % MT - storage rate 
Pr6 = Canada(:,3); % Gt - cumulative storage 
target6 = Canada(:,4);

load('Brazil.txt')
years7 = Brazil(:,1); % years
Q7 = Brazil(:,2); % MT - storage rate 
Pr7 = Brazil(:,3); % Gt - cumulative storage 
target7 = Brazil(:,4);

load('Australia.txt')
years8 = Australia(:,1); % years
Q8 = Australia(:,2); % MT - storage rate 
Pr8 =Australia(:,3); % Gt - cumulative storage 
target8 = Australia(:,4);

load('Korea.txt')
years9 = Korea(:,1); % years
Q9= Korea(:,2); % MT - storage rate 
Pr9 =Korea(:,3); % Gt - cumulative storage 
target9 = Korea(:,4);

load('Saudi.txt')
years10 = Saudi(:,1); % years
Q10 = Saudi(:,2); % MT - storage rate 
Pr10 =Saudi(:,3); % Gt - cumulative storage 
target10 = Saudi(:,4);

load('UAE.txt')
years11 = UAE(:,1); % years
Q11 = UAE(:,2); % MT - storage rate 
Pr11 =UAE(:,3); % Gt - cumulative storage 
target11 = UAE(:,4);

load('Indonesia.txt')
years12 = Indonesia(:,1); % years
Q12 = Indonesia(:,2); % MT - storage rate 
Pr12 =Indonesia(:,3); % Gt - cumulative storage 
target12 = Indonesia(:,4);

for i=1:length(target)
    
    total(i) = target(i)+target2(i)+target3(i)+target5(i)+target6(i)+target7(i)+target8(i)+target9(i)+target10(i)+...
        target12(i);
end
Sum = total.'

% Master database 
keySet = {'USgrowth','USstorage','USpeak','UStarget','CNgrowth','CNstorage','CNpeak','CNtarget',...
    'EUgrowth','EUstorage','EUpeak','EUtarget','UKgrowth','UKstorage','UKpeak','UKtarget',...
    'THgrowth','THstorage','THpeak','THtarget','CAgrowth','CAstorage','CApeak','CAtarget'...
    'BRgrowth','BRstorage','BRpeak','BRtarget','ASgrowth','ASstorage','ASpeak','AStarget',...
    'SKgrowth','SKstorage','SKpeak','SKtarget','SAgrowth','SAstorage','SApeak','SAtarget',...
    'UAEgrowth','UAEstorage','UAEpeak','UAEtarget','Indogrowth','Indostorage','Indopeak','Indotarget','Sum'};
T = table(years,Q,Pr,target,years2,Q2,Pr2,target2,years3,Q3,Pr3,target3,years4,Q4,Pr4,target4,...
    years5,Q5,Pr5,target5,years6,Q6,Pr6,target6,years7,Q7,Pr7,target7,years8,Q8,Pr8,target8,...
    years9,Q9,Pr9,target9,years10,Q10,Pr10,target10,years11,Q11,Pr11,target11,years12,Q12,Pr12,target12,Sum);
writetable(T, '12_country_database')



