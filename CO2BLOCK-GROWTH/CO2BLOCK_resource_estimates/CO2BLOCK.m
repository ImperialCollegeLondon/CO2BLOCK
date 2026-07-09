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
fname = 'Global_complementary.xlsx';        % name of the input data file
%-

%-- setting parameters
nr_region = 12;                     % number of storage regions in the input file
correction = 'off' ;                % set on/off if you want to apply correction for superposition
dist_min = 2 ;                      % minimum inter-well distance [km]
dist_max = 'auto'; %'auto';          % maximum inter-well distance [km]. Set a number or 'auto' if you prefer automatic calculation
nr_dist = 200 ;                     % number of inter-well distances to explore
nr_well_max = 'auto';  %'auto';      % maximum number of wells. Set a number or 'auto' if you prefer automatic calculation
rw = 0.2 ;                          % well radius [m]
time_yr = 30 ;                     % time of injection [years]
maxQ = 10 ;                         % maximum sustainable injection rate per well because of technological limitations [Mton/years]
minQ = 1;                           % minimum injection rate per site [Mt/y]
                                    % assign 0 to exclude min rate constraints
%%%%%%%%%%% END OF INPUT DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%calculate
% 
% for region_no = 1:nr_region
%     [region_name,thick,area_res,perm,por,dens_c,~,visc_w,compr,p_lim,rc,~,~, ~] = read_data(fpath,fname,region_no);
% 
%     % the first column includes the constant factor in radius of influence
%     diffusivity_regions(region_no,1) = (2.25*perm/(visc_w*compr))^0.5;
%     diffusivity_regions(region_no,2) = rc;
% end


% for region_no = 1:nr_region
%     [region_name,thick,area_res,~,por,dens_c,~,~,compr,p_lim,~,~,~, ~] = read_data(fpath,fname,region_no);
% 
%     % static theoretical capacity (Gt)
%     Max_Vol_theory = area_res * 1e6 * thick * por * dens_c / 1e12;
% 
%     % capacity under closed tank model conditions (compressibility)
%     % note that the real cmpressibility that should be used here equals
%     % compr/phi. Thus, this phi cancels out by the other phi in the
%     % equation used for the convesion of the total volume to pore volume
%     Max_Vol_tank = area_res * 1e6 * thick * dens_c * compr * p_lim * 1e6 / 1e12;
% 
%     Tank_model_summary {region_no , 1} = region_name;
%     Tank_model_summary {region_no , 2} = Max_Vol_theory;
%     Tank_model_summary {region_no , 3} = Max_Vol_tank;
% 
%     % Convert to table and set column names
%     Tank_model_summary_table = cell2table(Tank_model_summary, ...
%         'VariableNames', {'Region', 'Max theoretical capacity [Gt]','Max tank capacity [Gt]'});
% 
%     writetable(Tank_model_summary_table,'Tank_model_summary_Anadarko.xls');
% end

for region_no = 1:nr_region
[region_name,d_list,well_list,d_max,Q_M_each,V_M,Table_Q,Table_V,p_sup_vec] = calculate(fpath,fname,region_no,correction,dist_min,...
    dist_max,nr_dist,nr_well_max,rw,time_yr,maxQ,minQ);

% Truncate to 31 characters if needed
region_name = region_name(1:min(end, 31));

% find maximum volume/rate for each of the injection scenarios
% Convert table to array (numerical part only)
dataArray_Q = table2array(Table_Q);
dataArray_Vol = table2array(Table_V);

% Find where Q < minQ
% rows and columns include different site arrangements and
% distances respectively
mask = dataArray_Q < minQ;

% makig single injection site exempt from the minQ criterion
mask(1,1) = false;

% Set masked entries to zero
dataArray_Q(mask) = 0;
dataArray_Vol(mask) = 0;

% Find the max value and its linear index
[max_Q, ~] = max(dataArray_Q(:));
[max_Vol, linearIdx] = max(dataArray_Vol(:));

% calculating the max Q anf vol in each row (for each scenario with a
% particular number of injection sites)
% this is for plotting the effect of the number of injection sites on the
% maximum storage capacity that can be achieved
Max_Q_vol_scenarios_Qmin{region_no}(:,1) = well_list';
Max_Q_vol_scenarios_Qmin{region_no}(:,2) = max(dataArray_Q, [], 2);
Max_Q_vol_scenarios_Qmin{region_no}(:,3) = max(dataArray_Vol, [], 2);

% Convert linear index to row and column
[row, col] = ind2sub(size(dataArray_Vol), linearIdx);

% calculation results for all regions
% Column 1: name of storage regions 
% Column 2: maximum achievable injection rate 
% Column 3: injection rate/site corresponding to maximum storage capacity
% Column 4: maximum storage capacity
% Column 5: corresponding number of sites
% Column 6: corresponding maximum distance
% column 7: injection rate per region 
Regional_storage_summary {region_no , 1} = region_name;
Regional_storage_summary {region_no , 2} = max_Q;
Regional_storage_summary {region_no , 3} = Table_Q{row,col};
Regional_storage_summary {region_no , 4} = max_Vol;
Regional_storage_summary {region_no , 5} = well_list(row);
Regional_storage_summary {region_no , 6} = round(d_list(col),1);
Regional_storage_summary {region_no , 7} = well_list(row)*Table_Q{row,col};

% Convert to table and set column names
Regional_storage_summary_table = cell2table(Regional_storage_summary, ...
    'VariableNames', {'Region', 'Max Q [Mt/y]','Optimum Q [Mt/y]','Max capacity [Gt]', 'Number of wells', 'Distance','Region Q [Mt/y]'});

%%%%%%%%%%%%%%%%%%%%  WRITE TABLES  %%%%%%%%%%%%%%%%%%%% 

% Generating tables of injection rate and volume with Q < Qmin set to zero
Table_Q_minQ = array2table(dataArray_Q, ...
    'VariableNames', Table_Q.Properties.VariableNames, ...
    'RowNames', Table_Q.Properties.RowNames);

Table_V_minQ = array2table(dataArray_Vol, ...
    'VariableNames', Table_V.Properties.VariableNames, ...
    'RowNames', Table_V.Properties.RowNames);

% writetable(Table_Q,['Q_M_max_per_well_inj_rate_',num2str(time_yr),'y','.xls'], 'Sheet', region_name, 'WriteRowNames',true);
% writetable(Table_V,['V_M_max_storage_capacity_',num2str(time_yr),'y','.xls'], 'Sheet', region_name, 'WriteRowNames',true);

% writetable(Table_Q_minQ,['Q_M_max_per_well_inj_rate_minQ=',num2str(minQ),'_',num2str(time_yr),'y','.xls'], 'Sheet', region_name, 'WriteRowNames',true);
% writetable(Table_V_minQ,['V_M_max_storage_capacity_minQ=',num2str(minQ),'_',num2str(time_yr),'y','.xls'], 'Sheet', region_name, 'WriteRowNames',true);

disp(['Calculations of region ',num2str(region_no),' completed'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% injectivity plot
if region_no == 1 || region_no == 22
    fig = figure('Name',['Injectivity_',region_name],'Color', 'w',...
        'Units','centimeters','Position', [2 2 9 9]);
    
    % --- Plot contour ---
    if region_no == 1
        levels= [0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7];
    end

    if region_no == 22
        levels= [0.5,1,4,8,12,16,19];
    end

    [C, h] = contour(d_list, well_list, real(Q_M_each),levels, ...
                     'Color', [.6 .6 .6], 'LineWidth', 0.5);

    hold on;
    
    % --- Apply contour labels ---
    ht = clabel(C,h, ...
            'FontSize', 8, ...
            'Interpreter', 'tex', ...
            'FontName', 'Helvetica',...
            'Color', [.6 .6 .6]);
        
    % --- d_max line ---
    plot(d_max, well_list, 'LineWidth', 1.5, 'Color', [0.00,0.45,0.74]);

    hold on

    scatter(d_list(col), well_list(row), 50, 'k', '+', 'LineWidth', 1.5); 

    labelText = sprintf('Capacity = %.2f Gt', max_Vol);

    text(d_list(col) * 1.05, ...   % small shift to the right
     well_list(row) *1.1, ...  % small shift upward
     labelText, ...
     'FontSize', 8, ...
     'FontName', 'Helvetica', ...
     'Interpreter', 'tex', ...
     'HorizontalAlignment', 'left', ...
     'VerticalAlignment', 'bottom');

    % --- Axes formatting ---
    set(gca, 'FontSize', 8, ...
             'FontName', 'Helvetica', ...
             'TickLabelInterpreter', 'tex', ...
             'Box', 'on');
    
    if region_no == 1
        xlim([2, 50]);
        ylim([1, 200]);

        text(-4, 207, '(a)', 'FontSize',10, 'Interpreter','tex')
    end

    if region_no == 22
        xlim([2, 80]);
        ylim([1, 1000]);

        text(-9, 1030, '(b)', 'FontSize',10, 'Interpreter','tex')
    end
   
    title(region_name, ...
          'FontSize', 8, 'Interpreter', 'tex', 'FontWeight','normal', ...
          'FontName', 'Helvetica');
    
    xlabel('Inter-site distance (km)', ...
           'FontSize', 8, 'Interpreter', 'tex', 'FontName', 'Helvetica');
    
    ylabel('Number of sites', ...
           'FontSize', 8, 'Interpreter', 'tex', 'FontName', 'Helvetica');
    
    inset = get(gca,'TightInset');
    set(gca,'LooseInset', inset + 0.005);

    exportgraphics(gcf,['Injectivity_',region_name,'.pdf'], 'ContentType','vector');
    exportgraphics(gcf, ['Injectivity_',region_name,'.png'], ...
    'Resolution', 600, 'BackgroundColor', 'white');
    saveas(gcf, ['Injectivity_',region_name,'.svg']);

    hold off;
end
end

writetable(Regional_storage_summary_table,['Regional_storage_summary_minQ=',num2str(minQ),'_',num2str(time_yr),'y','.xls']);

Total_capacity = sum(cellfun(@sum, Regional_storage_summary(:,4)));

disp(['Total storage capacity = ',num2str(Total_capacity),' [Gt]'])


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot storage resources
% Sort alphabetically by region name
[~, idx] = sort(Regional_storage_summary(:,1));
idx = flipud(idx);
Regional_storage_summary_sorted = Regional_storage_summary(idx, :);

% Time axis (in years)
t = linspace(0, time_yr, 500);  % 0 to 50 years, 100 points

% Set up figure
storage_capacity_plot = figure;
storage_capacity_plot.Color = [1,1,1];
storage_capacity_plot.Units = 'centimeters';
storage_capacity_plot.Position = [6 6 15.5 10.5];    % [x y w h]
hold on;

% Initialize cumulative storage
Bottom_curve = zeros(size(t));

cmap = parula(nr_region);

% Loop through each region
for i = 1:size(Regional_storage_summary_sorted, 1)
    Top_curve = Bottom_curve + Regional_storage_summary_sorted{i,3} * Regional_storage_summary_sorted{i,5}*t/1000;  
    
    % Fill area between curves using region color
    region_color = cmap(i, :);
    
    fill([t, fliplr(t)], [Bottom_curve, fliplr(Top_curve)], ...
         region_color, 'EdgeColor', 'none', 'FaceAlpha', 0.9, ...
         'DisplayName', Regional_storage_summary_sorted{i,1});
    
    % Plot solid top curve line, exclude from legend
    plot(t, Top_curve, '-', 'Color', [0 0 0], 'LineWidth', 1, 'HandleVisibility', 'off');

    % Update bottom_curve for next region
    Bottom_curve = Top_curve;

end

ylim([0 105]);

set(gca, 'LineWidth', 0.5,...
         'Layer', 'top', ...
         'FontName', 'Helvetica', ...
         'FontSize', 8, ...
         'FontWeight', 'normal', ...
         'FontAngle', 'normal', ...
         'Box', 'on', ...
         'XGrid', 'off', ...
         'YGrid', 'off');

% Final touches
xlabel('Time (y)', 'FontName', 'Helvetica', 'FontSize', 8, 'FontWeight', 'normal');
ylabel('Storage capacity [Gt]', 'FontName', 'Helvetica', 'FontSize', 8, 'FontWeight', 'normal');

h = findobj(gca, 'Type', 'Patch');
lgd = legend(h, flip(Regional_storage_summary_sorted(:,1)), ...
             'FontName', 'Helvetica', ...
             'FontSize', 8, ...
             'FontWeight', 'normal', ...
             'FontAngle', 'normal', ...
             'Location', 'eastoutside', 'Interpreter', 'tex');

set(lgd, 'Box', 'off');

% Adjust figure layout to make room for legend
set(gca, 'LooseInset', get(gca, 'TightInset'));

% Export as vector PDF for journal submission
exportgraphics(gcf, ['storage capacity_minQ=',num2str(minQ),'.pdf'], 'ContentType','vector');
exportgraphics(gcf, ['storage capacity_minQ=',num2str(minQ),'.png'], ...
    'Resolution', 600, 'BackgroundColor', 'white');
saveas(gcf, ['storage capacity_minQ=',num2str(minQ),'.svg']);






% 
% %% effect of injection duration on storage capacity
% count=1;
% nr_dist = 1000;
% for time_yr = 10:5:200   
%     % total storage capacity for different injection periods
%     storage_capacity_time(count,1) = time_yr;
% 
%     storage_capacity_time(count,2) = 0;
% 
%     % calculation for each of the regions
%     for region_no = 1:nr_region
%         [region_name,d_list,well_list,d_max,Q_M_each,V_M,Table_Q,Table_V,p_sup_vec] = calculate(fpath,fname,region_no,correction,dist_min,...
%             dist_max,nr_dist,nr_well_max,rw,time_yr,maxQ);
% 
%         % Truncate to 31 characters if needed
%         region_name = region_name(1:min(end, 31));
% 
%         % find maximum volume/rate for each of the injection scenarios
%         % Convert table to array (numerical part only)
%         dataArray_Vol = table2array(Table_V);
% 
%         % Find the max value and its linear index
%         [max_Vol, linearIdx] = max(dataArray_Vol(:));    
% 
%         % total storage capacity for each injection period
%         storage_capacity_time(count,2) = storage_capacity_time(count,2) + max_Vol;
%     end
% 
% 
%     count = count +1;
% end
% 
% 
% storage_capacity_vs_time_plot = figure;
% storage_capacity_vs_time_plot.Color = [1,1,1];
% storage_capacity_vs_time_plot.Units = 'centimeters';
% storage_capacity_vs_time_plot.Position = [6 6 10 9];    % [x y w h]
% 
% plot(storage_capacity_time(:,1),storage_capacity_time(:,2),...
%     '-', 'Color', [0 0 0], 'LineWidth', 1);



