% This code implements an integrated framework for assessing the scaleup of
% geological CO2 storage by coupling techno-economic growth models with
% physics-based reservoir pressure constraints. It identifies feasible
% deployment pathways and storage resources that satisfy both industrial
% growth and subsurface injectivity limitations.

% Please, read the file README for instructions. This software is free. 
% Please cite CO2BLOCK-GROWTH as:
% https://github.com/ImperialCollegeLondon/CO2BLOCK-GROWTH/
% and
% Kivi, I.R., Gao, X., Krevor, S. (2026). Coupled geophysical and 
% technoeconomic growth constraints on geological carbon storage scaleup 
% with an application to the UK. Under review 


clearvars; close all;

% -------------------------------------------------------------------------
%                              Input data
% -------------------------------------------------------------------------

%-- data file directory and name
fpath = '';                         % directory of the input data file
fname = 'UK_sites.xlsx';        % name of the input data file

outputDir = 'output';               % name of the folder to save the outputs

allocation_order = "manual";                    % ascend: allocating resources in the order of ascending storage capacity (smallest first)
                                                % descend: allocating resources in the order of descending storage capacity (largest first)
                                                % random: random allocation of resources
                                                % manual: assigning the order of the units manually 

% order of the units when assigning them to scaleup trajectories
% this is used only when allocation_order = "manual"
sort_idx = [6 7 8 9 10 11 12 13 1 18 20 22 23 24 26 27 3 4 5 15 2 14 15 16 17 19 21 25];                                                
               
% Loading resources from a saved files allows for rapidly testing a variety 
% of growth models
% To generate and save resources file, run the code to the beginning of the
% section "Allocation of storage resources to growth trajectories". The
% saved file can be called later for repeated screening of resource
% allocation to growth trajectories. 
storage_resource_calculation = 'savedfile';     % 'savedfile': allocation from already calculated storage resources
                                                % 'calculate': calculate storage resources and then do allocation

if storage_resource_calculation == 'savedfile'
    % resource allocation from saved calculations of storage resources
    % this is an example of resources calculated for UK offshore sites
    load(fullfile(outputDir,'Storage_resources_400y_initialparameters_freq5.mat'));

    % Logistic growth model with cumulative storage amount of 69 Gt and
    % 2050 target rate of 1600 Mt/y
    load(fullfile(outputDir,'Growthrate_model_Cum69Gt_peak1600Mty_175Mty_130y.mat'))
    
    inj_duration = 400;
    
    allocation_duration = 130;

    %-- setting parameters
    dist_min = 2 ;                      % minimum inter-well distance [km]
    dist_max = 'auto';                  % maximum inter-well distance [km]. Set a number or 'auto' if you prefer automatic calculation
    nr_dist = 100 ;                     % number of inter-well distances to explore
    nr_well_max = 'auto';               % maximum number of wells. Set a number or 'auto' if you prefer automatic calculation
    rw =   0.2;                         % well radius [m]
    maxQ = 20 ;                         % maximum sustainable injection rate per well because of technological limitations [Mton/years]
    minQ = 1;                           % minimum injection rate per site [Mt/y]
                                            % assign 0 to exclude min rate constraints
    
else

    % if storage resources are calculated for the first time, the required
    % data need to be defined here
    inj_duration = 400;
    
    allocation_duration = 130;

    % doing calculations of storage resources
    
    %-- setting parameters
    correction = 'off' ;                % set on/off if you want to apply correction for superposition
    dist_min = 2 ;                      % minimum inter-well distance [km]
    dist_max = 'auto';                  % maximum inter-well distance [km]. Set a number or 'auto' if you prefer automatic calculation
    nr_dist = 100 ;                     % number of inter-well distances to explore
    nr_well_max = 'auto';               % maximum number of wells. Set a number or 'auto' if you prefer automatic calculation
    rw =   0.2;                         % well radius [m]
    maxQ = 20 ;                         % maximum sustainable injection rate per well because of technological limitations [Mton/years]
    minQ = 1;                           % minimum injection rate per site [Mt/y]
                                        % assign 0 to exclude min rate constraints

    % Assigning maximum number of injection sites for each of the storage
    % units/regions
    % nr_well_max_list = [2;192;34;21;100;5;3;1;3;1;1;18;12;76;25;28;136;415;5;100;343;309;31;1;96;1;15];
       
    % number of storage regions in the input file
    nr_region = height(readtable(fname,'VariableNamingRule','preserve'));

    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    % time series
    storage_periods = 5:5:inj_duration;
    n_time = length(storage_periods);


% -------------------------------------------------------------------------
%                   Storage resource calculations
% -------------------------------------------------------------------------
    for region_no = 1:nr_region
        % to assign maximum number of sites for each unit/region, the line
        % below needs to be uncommented
        % nr_well_max = nr_well_max_list(region_no);

        [region_name,~,well_list,~,Q_all,V_all,~] = calculate_timeseries(fpath,fname,region_no,correction,dist_min,...
            dist_max,nr_dist,nr_well_max,rw,storage_periods,maxQ);
    
        % Truncate to 31 characters if needed
        region_name = region_name(1:min(end, 31));
 
        for it = 1:n_time

            dataArray_Q = Q_all(:,:,it);      % per-site rate [Mt/yr] table among all configurations at time it
            dataArray_Vol = V_all(:,:,it);    % total volume [Gt] table among all configurations at time it
        
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
            [max_Vol, linearIdx] = max(dataArray_Vol(:));
        
            % Convert linear index to row and column
            [row, col] = ind2sub(size(dataArray_Vol), linearIdx);
        
            % summary of optimum injection rates for different scenarios
            % Column 1: number of storage regions 
            % Column 2: name of storage regions 
            % Column >3: injection rate per region  
            % Total optimum injection rate for this time
            % dataArray_Q is per-well rate, so multiply by number of wells
            Rate_opt_all_scenarios {region_no , 1} = region_no;
            Rate_opt_all_scenarios {region_no , 2} = region_name;
            Rate_opt_all_scenarios {region_no,it+2} = well_list(row) * dataArray_Q(row,col);
        
            % summary of the number of the sites for different scenarios
            % Column 1: number of storage regions 
            % Column 2: name of storage regions 
            % Column >3: Optimum number of sites for this time  
            Site_no_all_scenarios {region_no , 1} = region_no;
            Site_no_all_scenarios {region_no , 2} = region_name;
            Site_no_all_scenarios {region_no , it+2} = well_list(row);
        end

        disp(['Calculations of region ',num2str(region_no),' completed'])
    end
      
% -------------------------------------------------------------------------
%                   writing summary of the calculations
% -------------------------------------------------------------------------
    Rate_opt_all_scenarios_varname{1} = 'Region_no';
    Rate_opt_all_scenarios_varname{2} = 'Region_name';
    
    for i = 1:n_time
        Rate_opt_all_scenarios_varname{i + 2} = ['Q [Mt/yr] for t= ' num2str(storage_periods(i)) ' y'];
    end
    
    % Convert to table and set column names
    Rate_opt_all_scenarios_table = cell2table(Rate_opt_all_scenarios, ...
            'VariableNames', Rate_opt_all_scenarios_varname);
    
    writetable(Rate_opt_all_scenarios_table,fullfile(outputDir,['Rate_opt_summary_minQ=',num2str(minQ),'_maxQ=',num2str(maxQ),'_sites=',num2str(nr_well_max),'_',num2str(inj_duration),'y','.xlsx']));

    % writing summary of site numbers per regions for all scenarios
    Site_no_all_scenarios_varname{1} = 'Region_no';
    Site_no_all_scenarios_varname{2} = 'Region_name';
    
    for i = 1:n_time
        Site_no_all_scenarios_varname{i + 2} = ['Q [Mt/yr] for t= ' num2str(storage_periods(i)) ' y'];
    end
    
    % Convert to table and set column names
    Site_no_all_scenarios_table = cell2table(Site_no_all_scenarios, ...
            'VariableNames', Site_no_all_scenarios_varname);
    
    writetable(Site_no_all_scenarios_table,fullfile(outputDir,['Site_no_summary_minQ=',num2str(minQ),'_maxQ=',num2str(maxQ),'_sites=',num2str(nr_well_max),'_',num2str(inj_duration),'y','.xlsx']));
end


% -------------------------------------------------------------------------
%           Allocation of storage resources to growth trajectories
% -------------------------------------------------------------------------
[Resource_assignment_table, storage_resource_total, step_curve_start, ...
 step_curve_end, region_name_array_sorted, site_no_step] = ...
    resource_allocation(Growthrate_model, sort_idx, Rate_opt_all_scenarios, ...
    Site_no_all_scenarios, storage_periods, allocation_duration, ...
    allocation_order, minQ, maxQ, nr_region, fpath, fname, correction, ...
    dist_min, dist_max, nr_dist, nr_well_max, rw, outputDir);

allocation_step = length(storage_resource_total);


% -------------------------------------------------------------------------
%                plotting storage resource assignment
% -------------------------------------------------------------------------
years = Growthrate_model(:,1);
total_rate = Growthrate_model(:,2);

Figure_storage_assignment = figure('Color','w');
set(Figure_storage_assignment, 'Units', 'centimeters', 'Position', [6 6 15.5 10.5]);  % 150 mm × 90 mm

hold on

% Plot growth curve (main black)
hGrowth = plot(years, total_rate, 'k-', 'LineWidth', 1.5, 'DisplayName','Growth rate model');

% Color map
cmap = flip(parula(allocation_step));

% Helper: function to get interpolated growth value at arbitrary time
growth_at = @(t) interp1(years, total_rate, t, 'linear', 'extrap');

% Loop regions
for i = 1:allocation_step 
    % getting parameters for resource allocation
    [constraining_curve_x, constraining_curve_y, x_LT, x_RT] = ...
    get_allocation_polygon(i, years, total_rate, storage_resource_total, ...
    step_curve_start, step_curve_end, allocation_duration);

    hFill(i) = fill(constraining_curve_x, constraining_curve_y, cmap(i,:), 'FaceAlpha', 0.8, 'EdgeColor','none', ...
         'DisplayName', region_name_array_sorted{i});
    hold on

    if i < allocation_step
        plot([x_LT x_RT], [storage_resource_total(i) storage_resource_total(i)], ...
         '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, 'HandleVisibility', 'off');

        hold on

        % plot([2134 2143], [storage_resource_total(i) storage_resource_total(i)], ...
        %  '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, 'HandleVisibility', 'off');
        % 
        % hold on
    end
    
    
    % text(2080 + allocation_duration/2-9, 0.5*(max(constraining_curve_y)+min(constraining_curve_y)), num2str(site_no_step(i)), 'FontSize',8, 'Interpreter','tex')
    
end

set(gca, 'LineWidth', 0.5,...
         'Layer', 'top', ...
         'FontName', 'Helvetica', ...
         'FontSize', 8, ...
         'FontWeight', 'normal', ...
         'FontAngle', 'normal', ...
         'Box', 'on', ...
         'XGrid', 'off', ...
         'YGrid', 'off');

xlabel('Time (yr)', 'FontName', 'Helvetica', 'FontSize', 8, 'FontWeight', 'normal');
ylabel('Storage rate (Mt yr^{-1})', 'FontName', 'Helvetica', 'FontSize', 8, 'FontWeight', 'normal');
% title('Projected CO_{2} Storage Growth and Regional Assignments', ...
%       'FontName', 'Helvetica', 'FontSize', 8, 'FontWeight', 'normal');

% subset figure lettering
text(2010, 2350, '(a)', 'FontSize',10, 'Interpreter','tex')

% text(2125, 1650, 'Number of sites', 'FontSize',8, 'Interpreter','tex')

text(2119.6,327.82,'S3: 175 Mt yr^{-1}', ...
    'FontSize',8,'FontName','Helvetica', ...
    'VerticalAlignment','bottom')

xlim([2025 2160]);
% choose y-limit to show some headroom
ylim([0 1700]);

xticks([2040 2060 2080 2100 2120 2140 2160])
yticks([0 200 400 600 800 1000 1200 1400 1600])

% legend
legendHandles = [flipud(hFill(:)); hGrowth];

legendLabels = [flip(region_name_array_sorted(:)); ...
                {'Growth rate model'}];

legend(legendHandles, legendLabels, ...
       'FontName', 'Helvetica', ...
       'FontSize', 8, ...
       'FontWeight', 'normal', ...
       'FontAngle', 'normal', ...
       'Box','off', ...
       'Location', 'eastoutside', ...
       'Interpreter', 'tex');

uistack(hGrowth, 'top');

set(gca, 'LooseInset', get(gca, 'TightInset'));

disp(['Number of injection siters required = ',num2str(sum(site_no_step))])

% Export as vector PDF for journal submission
exportgraphics(gcf, fullfile(outputDir,['storage_assignment_minQ=',num2str(minQ),'_maxQ=',num2str(maxQ),'_sites=',num2str(nr_well_max),'_',num2str(allocation_duration),'y','.pdf']), 'ContentType','vector');
exportgraphics(gcf, fullfile(outputDir,['storage_assignment_minQ=',num2str(minQ),'_maxQ=',num2str(maxQ),'_sites=',num2str(nr_well_max),'_',num2str(allocation_duration),'y','.png']), ...
    'Resolution', 600, 'BackgroundColor', 'white');
saveas(gcf, fullfile(outputDir,['storage_assignment_minQ=',num2str(minQ),'_maxQ=',num2str(maxQ),'_sites=',num2str(nr_well_max),'_',num2str(allocation_duration),'y','.svg']));



%% Assignment complete?
if storage_resource_total(end) < max(total_rate)
    disp('Storage resources is not enough to accomodate the projected storage rates in this scenario')
    disp(['peak rate achieved = ',num2str(storage_resource_total(end))])
end


% -------------------------------------------------------------------------
%         Different scenarios of random allocation of resources
% -------------------------------------------------------------------------
%% Different random allocation scenarios
% number of random allocation scenarios
random_no = 100;

failed_arrangements = NaN(random_no, 2);

years = Growthrate_model(:,1);
total_rate = Growthrate_model(:,2);
for k = 1:random_no
    storage_resource_step = [];
    site_no_step = [];
    storage_resource_step = [];
    storage_resource_total = [];
    step_curve_start = [];
    step_curve_end = [];
    period_step = [];
    cum_vol = [];
    step_start = [];
    step_end = [];

    % Helper: function to get interpolated growth value at arbitrary time
    growth_at = @(t) interp1(years, total_rate, t, 'linear', 'extrap');
    
    % Convert to numeric matrix of rates of resource use
    resource_rate_matrix = cell2mat(Rate_opt_all_scenarios(:,3:end)); % [nRegions x nPeriods]
    
    region_name_array = Rate_opt_all_scenarios(:,2);
    
    region_no_array = Rate_opt_all_scenarios(:,1);
    
    % Convert to numeric matrix of the number of injection sites
    Site_no_matrix = cell2mat(Site_no_all_scenarios(:,3:end)); % [nRegions x nPeriods]
        
    % Sort regions by total storage capacity in a rando way
    sort_idx = randperm(size(resource_rate_matrix, 1));

    sort_idx_random_list(k,:) = sort_idx;

    resource_rate_matrix_sorted = resource_rate_matrix(sort_idx,:);
    region_name_array_sorted = region_name_array(sort_idx);
    region_no_array_sorted = region_no_array(sort_idx);
    
    Site_no_matrix_sorted = Site_no_matrix(sort_idx,:);
    
    % remaining time to allocate storage resources
    t_remaining = allocation_duration; 
    
    % remaining rate from the growth rate model (-: allocated; +: need to be allocated in the next steps)
    % a vector of the same size as the growth rate model
    total_rate_remain = total_rate;
    
    for allocation_step = 1 : nr_region
    
        % index of first storage period > remaining time to allocate
        idx = find(storage_periods >= t_remaining, 1, 'first');   
    
        % possible cumulative injection volumes for the unit based on the
        % injection duration adn rate
        vol_possible = storage_periods(1:idx) .* resource_rate_matrix_sorted(allocation_step,1:idx);
    
        % storage resource allocated in the current step (Mt/y)
        storage_resource_step(allocation_step) = 0;
    
        breaking = 0;
        break_max_rate = 0;
    
        for i=idx:-1:1
            % remaining rate to be assigned to storage resources for this
            % provisional injection rate
            total_rate_remain_test = total_rate_remain - resource_rate_matrix_sorted(allocation_step,i);
        
            % Find leftmost index where total_rate >= storage_resource
            left_idx = find(total_rate_remain_test >= 0, 1, 'first');
        
            % Find rightmost index where total_rate >= storage_resource
            right_idx = find(total_rate_remain_test >= 0, 1, 'last');
    
            % cumulative storage rate (Mt/y)
            if allocation_step == 1
                storage_resource_total_provisional = resource_rate_matrix_sorted(allocation_step,i);
            else
                storage_resource_total_provisional = storage_resource_total(allocation_step-1) + resource_rate_matrix_sorted(allocation_step,i);
            end
    
            % finding times at which the growth curve is intersected by the 
            % considered injection rate
            if ~isempty(left_idx)    % this is to exclude the curve of the last assignment
                % start point (year) of the curve of this step
                if years(left_idx) > 2030
                    step_curve_start_provisional = years(left_idx-1) + (storage_resource_total_provisional - total_rate(left_idx-1)) * (years(left_idx) - years(left_idx-1)) / (total_rate(left_idx) - total_rate(left_idx-1));
                else
                    step_curve_start_provisional = 2030;
                end
            end
        
            if ~isempty(left_idx)    % this is to exclude the curve of the last assignment 
                % end point (year) of the curve of this step
                if years(right_idx) < years(end)
                    step_curve_end_provisional = years(right_idx) + (storage_resource_total_provisional - total_rate(right_idx)) * (years(right_idx+1) - years(right_idx)) / (total_rate(right_idx+1) - total_rate(right_idx));
                else
                    step_curve_end_provisional = years(end);
                end
            end
    
            % time and injection rates above the considered rate
            % used for the calculation of area beneath the growth curve
            if ~isempty(left_idx)
                area_time = [step_curve_start_provisional;years(left_idx:right_idx);step_curve_end_provisional];
                area_rate = [0;total_rate_remain_test(total_rate_remain_test>=0);0];
    
                % calculating the 
                % area beneath the growth curve for the considered
                % injection rate
                if allocation_step == 1
                    area_growthrate = trapz(years, total_rate) - trapz(area_time, area_rate);
                else
                    area_growthrate = trapz(years, total_rate) - cum_vol(allocation_step-1) - trapz(area_time, area_rate);
                end
            elseif allocation_step == 1
                area_growthrate = trapz(years, total_rate);
            else
                area_growthrate = trapz(years, total_rate) - cum_vol(allocation_step-1);
            end
    
            % storage resource allocated in the current step (Mt/y)
            if (vol_possible(i) > area_growthrate || abs(vol_possible(i)-area_growthrate)<0.01) && resource_rate_matrix_sorted(allocation_step,i) >= storage_resource_step(allocation_step) 
                % number of sites to achieve the desired rate
                site_no_step(allocation_step) = Site_no_matrix_sorted(allocation_step,i);
                
                storage_resource_step(allocation_step) = resource_rate_matrix_sorted(allocation_step,i);
            
                storage_resource_total(allocation_step) = storage_resource_total_provisional;
    
                % real period of the operation of the assigned storage resource
                period_step(allocation_step) = storage_periods(i);

                if ~isempty(left_idx)
                    step_curve_start(allocation_step) = step_curve_start_provisional;
                    step_curve_end(allocation_step) = step_curve_end_provisional;
                end

                % cumulative volume allocated
                if allocation_step == 1
                    cum_vol(allocation_step) = area_growthrate;
                else
                    cum_vol(allocation_step) = cum_vol(allocation_step-1) + area_growthrate;
                end

                % check if the assigned rate exceeded the maximum rate of the
                % model
                if storage_resource_total_provisional >= max(total_rate)
                    break_max_rate = 1;

                    % final rate remaining to be addressed (the last region)
                    remaining_rate_final = max(total_rate) - storage_resource_total(allocation_step-1);
    
    
                    % repeat storage resource calculations for the last
                    % allocated region to find the minimum required sites to
                    % address the remaining injection rate and volume
                    [~,~,well_list,~,~,~,Table_Q,Table_V,~] = calculate(fpath,fname,sort_idx(allocation_step),correction,dist_min,...
                         dist_max,nr_dist,nr_well_max,rw,storage_periods(i),maxQ);

                    % find maximum volume/rate for each of the injection scenarios
                    % Convert table to array (numerical part only)
                    dataArray_Q = table2array(Table_Q);
                    dataArray_Q_total = table2array(well_list'.*Table_Q);
    
                    dataArray_Vol = table2array(Table_V);
                
                    % Find where Q < minQ
                    % rows and columns include different site arrangements and
                    % distances respectively
                    mask = dataArray_Q < minQ;
                
                    % makig single injection site exempt from the minQ criterion
                    mask(1,1) = false;
            
                    % Set masked entries to zero
                    dataArray_Q_total(mask) = 0;
                    dataArray_Vol(mask) = 0;
    
                    % Logical mask for valid scenarios
                    valid_mask = (dataArray_Q_total > remaining_rate_final) & (dataArray_Vol*1000 > area_growthrate);
                    
                    % Find all valid indices
                    [row_idx, col_idx] = find(valid_mask);
    
                    % Find smallest row index (minimum number of sites)
                    min_row = min(row_idx);
                
                    % Among those rows, take the first valid distance
                    first_col = col_idx(find(row_idx == min_row, 1));
    
                    Q_last = dataArray_Q(min_row, first_col)*well_list(min_row);
                    V_last = dataArray_Vol(min_row, first_col);
    
                    site_no_step(allocation_step) = well_list(min_row);
                    
                    storage_resource_step(allocation_step) = Q_last;
                
                    storage_resource_total(allocation_step) = storage_resource_total(allocation_step-1)+Q_last;
                
                end
    
                if isempty(left_idx)
                    breaking = 1;
                end
            end
    
            if break_max_rate ==1
                break
            end
        end
    
        total_rate_remain = total_rate_remain - storage_resource_step(allocation_step);
     
        % time frame of the storage use for each step
        % time starting and ending injection into resources of each step
        if allocation_step == 1
            step_start(allocation_step) = 2030;
            step_end(allocation_step) = 2030+allocation_duration-1;
        else
            step_start(allocation_step) = step_curve_start(allocation_step-1);
            step_end(allocation_step) = step_curve_end(allocation_step-1);
        end
    
        % storage resource assignment
        % Column 1: storage resource assignment step
        % Column 2: number of storage regions 
        % Column 3: name of storage regions 
        % Column 4: start of the storage resource assignment step
        % Column 5: end of the storage resource assignment step
        % Column 6: real operation duration of the assigned storage resource
        % Column 7: injection rate increment in this assignment step
        % Column 8: cumulative injection rate in this assignment step
        % Column 9: number of sites to achieve the rate
        % Column 10: start of the corresponding curve
        % Column 11: end of the corresponding curve
    
        warning('off','MATLAB:table:ModifiedAndSavedVarnames');

        if breaking ~= 1  
            % time remained to be assigned to storage resources
            t_remaining = step_curve_end(allocation_step) - step_curve_start(allocation_step)+1;
    
        else
            break;
        end
    end

    if storage_resource_total(end) < max(total_rate)
        failed_arrangements(k,1) = k;
        failed_arrangements(k,2) = nr_region;
    end

    % total number of sites for each allocation scenario
    site_no_allocation_scenario(k,1) = sum(site_no_step);

    % number of regions used in the allocation scenario
    region_no_allocation_scenario(k,1) = allocation_step;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Random scenarios plot
figure('Color','w','Units','centimeters','Position',[5 5 14 9])

blue = [0 0.4470 0.7410];
red  = [0.8500 0.3250 0.0980];

x = 1:random_no;

% ---------------- LEFT AXIS ----------------
yyaxis left
h1 = plot(x, site_no_allocation_scenario, ...
    'Color', blue, 'LineWidth', 1);
ylabel('Number of injection sites','FontSize',8,'FontName','Helvetica')

ylim([600 2400])

yticks([600 800 1000 1200 1400 1600])

hold on

% Limited horizontal segments (LEFT axis)
plot([0 25],[1191 1191],'k--','LineWidth',0.8,'HandleVisibility','off')
plot([0 25],[777 777],'k--','LineWidth',0.8,'HandleVisibility','off')
plot([0 25],[988 988],'k--','LineWidth',0.8,'HandleVisibility','off')

% LEFT axis annotations
text(4.99,1208.9,'prioritising large units', ...
    'FontSize',7,'FontName','Helvetica', ...
    'VerticalAlignment','bottom')

text(0.24,748.35,'prioritising small units', ...
    'FontSize',7,'FontName','Helvetica', ...
    'VerticalAlignment','top')

text(23.8,1353.25,'planned projects path', ...
    'FontSize',7,'FontName','Helvetica', ...
    'VerticalAlignment','bottom')

text(34.5,2301,'failed allocations', ...
    'FontSize',7,'FontName','Helvetica', ...
    'VerticalAlignment','bottom')

annotation('arrow',[0.269919642857142 0.357610119047619],...
    [0.281291666666667 0.453976851851852],'HeadWidth',6,'HeadLength',6);

annotation('arrow',[0.485744047619047 0.484988095238095],...
    [0.887305555555556 0.927287037037037],'HeadWidth',6,'HeadLength',6);

hold off


% ---------------- RIGHT AXIS ----------------
yyaxis right
h2 = plot(x, region_no_allocation_scenario, ...
    'Color', red, 'LineWidth', 1);
ylabel('Number of geological units','FontSize',8,'FontName','Helvetica')

ylim([-20 33])

yticks([10 15 20 25 30])

hold on

% Limited horizontal segments (RIGHT axis)
plot([75 100],[13 13],'k--','LineWidth',0.8,'HandleVisibility','off')
plot([75 100],[27 27],'k--','LineWidth',0.8,'HandleVisibility','off')
% plot([75 100],[20 20],'k--','LineWidth',0.8,'HandleVisibility','off')

hold on

% showing failed arrangements
scatter(failed_arrangements(:,1),failed_arrangements(:,2), 20, ...
    'o', ...              % circle markers
    'MarkerEdgeColor', 'k', ...   % black edge
    'LineWidth', 1);            % marker edge thickness

% Text annotations (aligned with segments)
text(77,13.3,'prioritising large units', ...
    'FontSize',7,'FontName','Helvetica', ...
    'VerticalAlignment','bottom')

text(77,32,{'prioritising small units/','planned projects path'}, ...
    'FontSize',7,'FontName','Helvetica', ...
    'VerticalAlignment','top')

% text(82,21,'planned projects path', ...
%     'FontSize',7,'FontName','Helvetica', ...
%     'VerticalAlignment','bottom')

hold off


% ---------------- COMMON FORMATTING ----------------
xlabel('Random allocation scenario', ...
    'FontSize',8,'FontName','Helvetica')

legend([h1 h2], {'Injection sites','Geological units'}, ...
    'Location','northwest', ...
    'FontSize',8,'FontName','Helvetica', ...
    'Box','off')

set(gca, ...
    'FontSize',8, ...
    'FontName','Helvetica', ...
    'LineWidth',0.5, ...
    'TickDir','in')

set(gca,'LooseInset',max(get(gca,'TightInset'),0.02))


exportgraphics(gcf, fullfile(outputDir,['random allocation','.pdf']), 'ContentType','vector');
exportgraphics(gcf, fullfile(outputDir,['random allocation','.png']), ...
    'Resolution', 600, 'BackgroundColor', 'white');
saveas(gcf, fullfile(outputDir,['random allocation','.svg']));
