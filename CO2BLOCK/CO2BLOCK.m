% This tool provides estimate of the CO2 storage capacity of a geological
% reservoir under different scenarios of well number and distance.
% Wells are placed into a grid configuration with equal number of rows and 
% columns, or with number of row and columns that differ by 1.

% Please, read the file README and the User Guide for instructions.
% This software is free. Please cite CO2BLOCK as:
%   https://github.com/co2block/CO2BLOCK 
% De Simone and Krevor (2021).  A tool for first order estimates and optimisation of dynamic storage resource capacity in saline aquifers”. International Journal of Greenhouse Gas Control, 106, 103258.

clearvars; close all;

%%%%%%  INPUT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-- data file directory and name
fpath = '';                         % directory of the input data file
fname = 'example_data.xlsx';        % name of the input data file
%-

%-- setting parameters
correction = 'off' ;                % set on/off if you want to apply correction for superposition
dist_min = 2 ;                      % minimum inter-well distance [km]
dist_max = 'auto';                  % maximum inter-well distance [km]. Set a number or 'auto' if you prefer automatic calculation
nr_dist = 30 ;                      % number of inter-well distances to explore
nr_well_max = 'auto';               % maximum number of wells. Set a number or 'auto' if you prefer automatic calculation
rw = 0.2 ;                          % well radius [m]
time_yr = 30 ;                      % time of injection [years]
maxQ = 5 ;                          % maximum sustainable injection rate per well because of technological limitations [Mton/years]

%%%%%%%%%%% END OF INPUT DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%calculate
[d_list,well_list,d_max,Q_M_each,V_M,Table_Q,Table_V] = calculate(fpath,fname,correction,dist_min,...
    dist_max,nr_dist,nr_well_max,rw,time_yr,maxQ);


%%    
%%%%%%%%%   PLOTS  %%%%%%%%%%%%%%%%%%%%  
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');


figure();   % sustainable per well flow-rate Q_M_each
[C,h] = contour(d_list, well_list,real(Q_M_each),'color', [.6 .6 .6], 'linewidth', 1);  hold on; 
clabel(C,h, 'Fontsize', 12, 'FontWeight','bold', 'color',[.6 .6 .6]);  
plot(d_max, well_list,'linewidth',2,'color','r'); 
set(gca,'FontSize',12);
xlim([min(d_list),max(d_list)]);
title('Maximum sustainable per well flow-rate $Q_{M}$ (Mt/yr) ', 'Fontsize',16);
xlabel('inter-well distance $d$ (km)','Fontsize',18); 
ylabel('number of wells $n$','Fontsize',18); hold off;

figure();   % sustainable storage capacity V_M
levels_nr = 20; %logspace(log(0.1),log(max(V_M(:))),15);
[C,h] = contour(d_list, well_list,real(V_M),levels_nr,'color', [.6 .6 .6], 'linewidth', 1);  hold on;  
h.LevelList=round(h.LevelList,2) ;
clabel(C,h, 'Fontsize', 12, 'FontWeight','bold', 'color',[.6 .6 .6]);
plot(d_max, well_list,'linewidth',2,'color','r'); 
set(gca,'FontSize',12);
xlim([min(d_list),max(d_list)]);
title('Maximum sustainable storage $V_{M}$ (Gt) ', 'Fontsize',16);
xlabel('inter-well distance $d$ (km)','Fontsize',18); 
ylabel('number of wells $n$','Fontsize',18); hold off;

%% 
%%%%%%%%%%%%%%%%%%%%  WRITE TABLES  %%%%%%%%%%%%%%%%%%%% 
writetable(Table_Q,'Q_M_max_per_well_inj_rate.xls', 'WriteRowNames',true);
writetable(Table_V,'V_M_max_storage_capacity.xls', 'WriteRowNames',true);
