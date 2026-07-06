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

load('Australia.txt')
years7 = Australia(:,1); % years
Q7 = Australia(:,2); % MT - storage rate 
Pr7 =Australia(:,3); % Gt - cumulative storage 
target7 = Australia(:,4);

load('Korea.txt')
years8 = Korea(:,1); % years
Q8= Korea(:,2); % MT - storage rate 
Pr8 =Korea(:,3); % Gt - cumulative storage 
target8 = Korea(:,4);

load('MiddleEast.txt')
years9 = MiddleEast(:,1); % years
Q9 = MiddleEast(:,2); % MT - storage rate 
Pr9 = MiddleEast(:,3); % Gt - cumulative storage 
target9 = MiddleEast(:,4);

load('Indonesia.txt')
years10 = Indonesia(:,1); % years
Q10 = Indonesia(:,2); % MT - storage rate 
Pr10 =Indonesia(:,3); % Gt - cumulative storage 
target10 = Indonesia(:,4);

total = zeros(length(1000),1);
growth = zeros(length(1000),1);
storage = zeros(length(1000),1);
peak = zeros(length(1000),1);


for i=1:length(target)
    
    total(i) = target(i)+target2(i)+target4(i)+target3(i)+target5(i)+target6(i)+target7(i)+target8(i)+target9(i)+target10(i);
end
Sum = total.';
save('S1','Sum')
% Master database 
keySet = {'USgrowth','USstorage','UStarget','CNgrowth','CNstorage','CNtarget',...
    'EUgrowth','EUstorage','EUtarget','UKgrowth','UKstorage','UKtarget',...
    'THgrowth','THstorage','THtarget','CAgrowth','CAstorage','CAtarget'...
    'ASgrowth','ASstorage','AStarget',...
    'SKgrowth','SKstorage','SKtarget','MEgrowth','MEstorage','MEpeak','MEtarget',...
  'Indogrowth','Indostorage','Indotarget','Sum'};
T = table(years,Q,target,years2,Q2,target2,years3,Q3,target3,years4,Q4,target4,...
    years5,Q5,target5,years6,Q6,target6,years7,Q7,target7,years8,Q8,target8,...
    years9,Q9,target9,years10,Q10,target10,Sum);
writetable(T, '10_country_database')

Histo = total.';
h = histogram(Histo,32,'FaceColor','#C875C4')
hold on 

xlabel('Total storage rate for 2050 [Gt/year]')
ylabel('Frequency')




