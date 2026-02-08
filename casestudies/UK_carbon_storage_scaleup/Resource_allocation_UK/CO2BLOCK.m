% This is adapted from the original CO2BLOCK version for storage resource
% allocation to geological CO2 storage scaleup trajectories in the UK.
% Details about the algorithm can be found in the readme file

clearvars; close all;

%%%%%%  INPUT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-- data file directory and name
fpath = '';                         % directory of the input data file
fname = 'UK_sites.xlsx';        % name of the input data file

% saving the code outputs into the 'results' folder
if ~exist('results', 'dir'); mkdir('results'); end

storage_resource_calculation = 'savedfile';     % use either 
                                                % 'savedfile': allocation from already calculated storage resources
                                                % the code in "Storage_capacity_UK" was run and the workspace was saved
                                                % or
                                                % 'calculate': calculate storage resources and then do allocation

if storage_resource_calculation == 'savedfile'
    % resource allocation from saved calculations of storage resources
    load('Storage_resources_400y_initialparameters.mat');

    % growth rate model for the target 175 Mt/y on 2050 and 100 years injection
    load(fullfile('Growth_Gompertz_models', 'Growthrate_model_Cum78Gt_175Mty_130y.mat'))
    
    % storage resources calculated for different durations up to inj_duration
    inj_duration = 400;
    
    % resource allocation is done for this duration
    allocation_duration = 130;

else
    % growth rate model for the target 175 Mt/y on 2050 and 100 years injection
    load(fullfile('Growth_Gompertz_models', 'Growthrate_model_Cum78Gt_175Mty_130y.mat'))
    
    inj_duration = 130;
    
    allocation_duration = 130;

    % doing calculations of storage resources
    period_count = 0;
    
    % loop on time of injection [years]
    for time_yr = 5:1:inj_duration
    
        %-- setting parameters
        nr_region = 27;                     % number of storage regions in the input file
        correction = 'off' ;                % set on/off if you want to apply correction for superposition
        dist_min = 2 ;                      % minimum inter-well distance [km]
        dist_max = 'auto';                  % maximum inter-well distance [km]. Set a number or 'auto' if you prefer automatic calculation
        nr_dist = 200 ;                     % number of inter-well distances to explore
        nr_well_max = 'auto';               % maximum number of wells. Set a number or 'auto' if you prefer automatic calculation
        rw = 0.2 ;                          % well radius [m]
        maxQ = 20 ;                         % maximum sustainable injection rate per well because of technological limitations [Mton/years]
        minQ = 1;                           % minimum injection rate per site [Mt/y]
                                            % assign 0 to exclude min rate constraints
    
    
    %%%%%%%%%%% END OF INPUT DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        period_count = period_count + 1;
        storage_periods(period_count) = time_yr;
    
        %%
        %calculate
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
            mask = dataArray_Q < minQ;
        
            % Set masked entries to zero
            dataArray_Q(mask) = 0;
            dataArray_Vol(mask) = 0;
        
            % Find the max value and its linear index
            [max_Q, ~] = max(dataArray_Q(:));
            [max_Vol, linearIdx] = max(dataArray_Vol(:));
        
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
            warning('off','MATLAB:table:ModifiedAndSavedVarnames');
            Regional_storage_summary {region_no , 1} = region_name;
            Regional_storage_summary {region_no , 2} = max_Q;
            Regional_storage_summary {region_no , 3} = dataArray_Q(linearIdx);
            Regional_storage_summary {region_no , 4} = max_Vol;
            Regional_storage_summary {region_no , 5} = well_list(row);
            Regional_storage_summary {region_no , 6} = round(d_list(col),1);
            Regional_storage_summary {region_no , 7} = well_list(row)*dataArray_Q(linearIdx);
        
            % Convert to table and set column names
            Regional_storage_summary_table = cell2table(Regional_storage_summary, ...
                'VariableNames', {'Region', 'Max Q [Mt/y]','Optimum Q [Mt/y]','Max capacity [Gt]', 'Number of wells', 'Distance','Region Q [Mt/y]'});
        
        
            % summary of optimum injection rates for different scenarios
            % Column 1: name of storage regions 
            % Column 2: name of storage regions 
            % Column >3: injection rate per region  
            Rate_opt_all_scenarios {region_no , 1} = region_no;
            Rate_opt_all_scenarios {region_no , 2} = region_name;
            Rate_opt_all_scenarios {region_no , period_count+2} = well_list(row)*dataArray_Q(linearIdx);
        
            % summary of the number of the sites for different scenarios
            % Column 1: name of storage regions 
            % Column 2: name of storage regions 
            % Column >3: site number  
            Site_no_all_scenarios {region_no , 1} = region_no;
            Site_no_all_scenarios {region_no , 2} = region_name;
            Site_no_all_scenarios {region_no , period_count+2} = well_list(row);
        end
    
        % writetable(Regional_storage_summary_table,['Regional_storage_summary_minQ=',num2str(minQ),'_',num2str(time_yr),'y','.xls']);
    
        Total_capacity = sum(cellfun(@sum, Regional_storage_summary(:,4)));
    
        disp(['Total storage capacity for ',num2str(time_yr),' y injection = ',num2str(Total_capacity),' [Gt]'])
    end
    
    % writing summary of optimum injection rates per regions for all scenarios
    Rate_opt_all_scenarios_varname{1} = 'Region_no';
    Rate_opt_all_scenarios_varname{2} = 'Region_name';
    
    for i = 1:period_count
        Rate_opt_all_scenarios_varname{i + 2} = ['Q [Mt/y] for t= ' num2str(storage_periods(i)) ' y'];
    end
    
    % Convert to table and set column names
    Rate_opt_all_scenarios_table = cell2table(Rate_opt_all_scenarios, ...
            'VariableNames', Rate_opt_all_scenarios_varname);
    
    writetable(Rate_opt_all_scenarios_table,['Rate_opt_summary_minQ=',num2str(minQ),'_maxQ=',num2str(maxQ),'_sites=',num2str(nr_well_max),'_',num2str(time_yr),'y','.xlsx']);


    % writing summary of site numbers per regions for all scenarios
    Site_no_all_scenarios_varname{1} = 'Region_no';
    Site_no_all_scenarios_varname{2} = 'Region_name';
    
    for i = 1:period_count
        Site_no_all_scenarios_varname{i + 2} = ['Q [Mt/y] for t= ' num2str(storage_periods(i)) ' y'];
    end
    
    % Convert to table and set column names
    Site_no_all_scenarios_table = cell2table(Site_no_all_scenarios, ...
            'VariableNames', Site_no_all_scenarios_varname);
    
    writetable(Site_no_all_scenarios_table,fullfile('results', ['Site_no_summary_minQ=',num2str(minQ),'_maxQ=',num2str(maxQ),'_sites=',num2str(nr_well_max),'_',num2str(time_yr),'y','.xlsx']));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Allocation of subsurface resources to growth trajectories
years = Growthrate_model(:,1);
total_rate = Growthrate_model(:,2);

% Helper: function to get interpolated growth value at arbitrary time
growth_at = @(t) interp1(years, total_rate, t, 'linear', 'extrap');

% Convert to numeric matrix of rates of resource use
resource_rate_matrix = cell2mat(Rate_opt_all_scenarios(:,3:end)); % [nRegions x nPeriods]

region_name_array = Rate_opt_all_scenarios(:,2);

region_no_array = Rate_opt_all_scenarios(:,1);

% Convert to numeric matrix of the number of injection sites
Site_no_matrix = cell2mat(Site_no_all_scenarios(:,3:end)); % [nRegions x nPeriods]

% Sort regions by total storage capacity (sum of sustainable rates over all periods)
[~, sort_idx] = sort(sum(resource_rate_matrix,2), 'descend');
% [~, sort_idx] = sort(sum(resource_rate_matrix,2));

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
            
            % check if the assigned rate exceeded the maximum rate of the
            % model
            if storage_resource_total_provisional >= max(total_rate)
                break_max_rate = 1;
            end
            
            % number of sites to achieve the desired rate
            site_no_step(allocation_step) = Site_no_matrix_sorted(allocation_step,i);
            
            storage_resource_step(allocation_step) = resource_rate_matrix_sorted(allocation_step,i);
        
            storage_resource_total(allocation_step) = storage_resource_total_provisional;

            if ~isempty(left_idx)
                step_curve_start(allocation_step) = step_curve_start_provisional;
                step_curve_end(allocation_step) = step_curve_end_provisional;
            end

            % real period of the operation of the assigned storage resource
            period_step(allocation_step) = storage_periods(i);

            % cumulative volume allocated
            if allocation_step == 1
                cum_vol(allocation_step) = area_growthrate;
            else
                cum_vol(allocation_step) = cum_vol(allocation_step-1) + area_growthrate;
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

    Resource_assignment {allocation_step , 1} = allocation_step;
    Resource_assignment {allocation_step , 2} = region_no_array_sorted(allocation_step);
    Resource_assignment {allocation_step , 3} = region_name_array_sorted(allocation_step);
    Resource_assignment {allocation_step , 4} = step_start(allocation_step);
    Resource_assignment {allocation_step , 5} = step_end(allocation_step);
    Resource_assignment {allocation_step , 6} = period_step(allocation_step);
    Resource_assignment {allocation_step , 7} = storage_resource_step(allocation_step);
    Resource_assignment {allocation_step , 8} = storage_resource_total(allocation_step);
    Resource_assignment {allocation_step , 9} = site_no_step(allocation_step);
    if breaking ~= 1
        Resource_assignment {allocation_step , 10} = step_curve_start(allocation_step);
        Resource_assignment {allocation_step , 11} = step_curve_end(allocation_step);

        % time remained to be assigned to storage resources
        t_remaining = step_curve_end(allocation_step) - step_curve_start(allocation_step)+1;

    else
        break;
    end
end

% Convert to table and set column names
Resource_assignment_table = cell2table(Resource_assignment, ...
     'VariableNames', {'Step', 'Region no','Region name','Start [y]', 'End [y]', 'Duration [y]','Step rate increment [Mt/y]','Step rate cumulative [Mt/y]','No_sites','Curve start [y]', 'Curve end [y]'});

writetable(Resource_assignment_table,fullfile('results', ['Resource_assignment_minQ=',num2str(minQ),'_maxQ=',num2str(maxQ),'_sites=',num2str(nr_well_max),'_',num2str(time_yr),'y','.xlsx']));



%% plotting storage resource assignment
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
    % define x boundaries for region i
    % LB: left bottom
    % LT: left top
    % RT: right top
    % RB: right bottom
    if i==1
        x_LB = 2030;
        x_RB = 2030+allocation_duration-1;
    else
        x_LB = step_curve_start(i-1);
        x_RB = step_curve_end(i-1);
    end

    idx1 = find(years >= x_LB, 1, 'first');
    idx4 = find(years >= x_RB, 1, 'first');

    % if i < allocation_step
    if storage_resource_total(i) < max(total_rate)        
        x_LT = step_curve_start(i);
        x_RT = step_curve_end(i);

        idx2 = find(years >= x_LT, 1, 'first');
        idx3 = find(years >= x_RT, 1, 'first');
    end

    % Build constraining curves for each region 
    
    if i==1
        if storage_resource_total(i)<total_rate(1) && storage_resource_total(i)<total_rate(end)
            constraining_curve_x = [x_LB;  x_LB;                       x_RB;                           x_RB];
            constraining_curve_y = [0;     storage_resource_total(i);  storage_resource_total(i);      0   ];  
        else if storage_resource_total(i)<total_rate(1) && storage_resource_total(i)>=total_rate(end)
            constraining_curve_x = [x_LB;  x_LB;                       x_RT;                         years(idx3:(idx4-1));        x_RB;                                               x_RB;    transpose(floor(x_RB):-1:ceil(x_LB))];
            constraining_curve_y = [0;     storage_resource_total(i);  storage_resource_total(i);    total_rate(idx3:(idx4-1));   min(growth_at(x_RB),storage_resource_total(i));     0;       transpose(zeros(size(floor(x_RB):-1:ceil(x_LB))))];  
        else
            constraining_curve_x = [x_LB;  x_LB;             years(idx1:(idx2-1));        x_LT;             transpose(ceil(x_LT):floor(x_RT));                               x_RT;                         years(idx3:(idx4-1));        x_RB;                                               x_RB;    transpose(floor(x_RB):-1:ceil(x_LB))];
            constraining_curve_y = [0;     growth_at(x_LB);  total_rate(idx1:(idx2-1));   growth_at(x_LT);  transpose(growth_at(x_LT)*ones(size(ceil(x_LT):floor(x_RT))));   storage_resource_total(i);    total_rate(idx3:(idx4-1));   min(growth_at(x_RB),storage_resource_total(i));     0;       transpose(zeros(size(floor(x_RB):-1:ceil(x_LB))))];  
        end
        end


            
        % constraining_curve_x = [x_LB;  x_LB;             years(idx1:(idx2-1));        x_LT;             transpose(ceil(x_LT):floor(x_RT));                               x_RT;                         years(idx3:(idx4-1));        x_RB;                                               x_RB;    transpose(floor(x_RB):-1:ceil(x_LB))];
        % constraining_curve_y = [0;     growth_at(x_LB);  total_rate(idx1:(idx2-1));   growth_at(x_LT);  transpose(growth_at(x_LT)*ones(size(ceil(x_LT):floor(x_RT))));   storage_resource_total(i);    total_rate(idx3:(idx4-1));   min(growth_at(x_RB),storage_resource_total(i));     0;       transpose(zeros(size(floor(x_RB):-1:ceil(x_LB))))];  
    else if storage_resource_total(i) < max(total_rate)
        if storage_resource_total(i)<total_rate(1) && storage_resource_total(i)<total_rate(end)
            constraining_curve_x = [x_LB;  x_LB;                     x_RB;                           x_RB];
            constraining_curve_y = [storage_resource_total(i-1);     storage_resource_total(i);  storage_resource_total(i);      storage_resource_total(i-1)];  
        else if storage_resource_total(i)<total_rate(1) && storage_resource_total(i)>=total_rate(end)
            constraining_curve_x = [x_LB;  x_LB;                     x_RT;                         years(idx3:(idx4-1));        x_RB;                                              x_RB;                           transpose(floor(x_RB):-1:ceil(x_LB))];
            constraining_curve_y = [storage_resource_total(i-1);     storage_resource_total(i);  storage_resource_total(i);    total_rate(idx3:(idx4-1));   min(growth_at(x_RB),storage_resource_total(i));    storage_resource_total(i-1);    transpose(storage_resource_total(i-1)*(ones(size(floor(x_RB):-1:ceil(x_LB)))))];
        else
            constraining_curve_x = [x_LB;                          x_LB;               years(idx1:(idx2-1));        x_LT;             transpose(ceil(x_LT):floor(x_RT));                               x_RT;                        years(idx3:(idx4-1));        x_RB;                                              x_RB;                           transpose(floor(x_RB):-1:ceil(x_LB))];
            constraining_curve_y = [storage_resource_total(i-1);   growth_at(x_LB);    total_rate(idx1:(idx2-1));   growth_at(x_LT);  transpose(growth_at(x_LT)*ones(size(ceil(x_LT):floor(x_RT))));   storage_resource_total(i);   total_rate(idx3:(idx4-1));   min(growth_at(x_RB),storage_resource_total(i));    storage_resource_total(i-1);    transpose(storage_resource_total(i-1)*(ones(size(floor(x_RB):-1:ceil(x_LB)))))];
        end
        end
            
        % constraining_curve_x = [x_LB;                          x_LB;               years(idx1:(idx2-1));        x_LT;             transpose(ceil(x_LT):floor(x_RT));                               x_RT;                        years(idx3:(idx4-1));        x_RB;                                              x_RB;                           transpose(floor(x_RB):-1:ceil(x_LB))];
        % constraining_curve_y = [storage_resource_total(i-1);   growth_at(x_LB);    total_rate(idx1:(idx2-1));   growth_at(x_LT);  transpose(growth_at(x_LT)*ones(size(ceil(x_LT):floor(x_RT))));   storage_resource_total(i);   total_rate(idx3:(idx4-1));   min(growth_at(x_RB),storage_resource_total(i));    storage_resource_total(i-1);    transpose(storage_resource_total(i-1)*(ones(size(floor(x_RB):-1:ceil(x_LB)))))];
        
    else
        constraining_curve_x = [x_LB;                          x_LB;               years(idx1:(idx4-1));        x_RB;                           transpose(floor(x_RB):-1:ceil(x_LB))];
        constraining_curve_y = [storage_resource_total(i-1);   growth_at(x_LB);    total_rate(idx1:(idx4-1));   storage_resource_total(i-1);    transpose(storage_resource_total(i-1)*ones(size(floor(x_RB):-1:ceil(x_LB))))];
    end
    end


    fill(constraining_curve_x, constraining_curve_y, cmap(i,:), 'FaceAlpha', 0.8, 'EdgeColor','none', ...
         'DisplayName', region_name_array_sorted{i});
    hold on

    if i < allocation_step
        plot([x_LT x_RT], [storage_resource_total(i) storage_resource_total(i)], ...
         '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, 'HandleVisibility', 'off');
    
        hold on
    end
    
    text(2030 + allocation_duration/2-9, 0.5*(max(constraining_curve_y)+min(constraining_curve_y)), num2str(site_no_step(i)), 'FontSize',8, 'Interpreter','tex')
    
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

xlabel('Time (y)', 'FontName', 'Helvetica', 'FontSize', 8, 'FontWeight', 'normal');
ylabel('Storage rate (Mt/y)', 'FontName', 'Helvetica', 'FontSize', 8, 'FontWeight', 'normal');
% title('Projected CO_{2} Storage Growth and Regional Assignments', ...
%       'FontName', 'Helvetica', 'FontSize', 8, 'FontWeight', 'normal');

% subset figure lettering
text(2013, 2000, '(a)', 'FontSize',10, 'Interpreter','tex')

xlim([2030 2030+allocation_duration]);
% choose y-limit to show some headroom
ylim([0 2000]);

% legend
[~, hObjects] = legend('show','FontName', 'Helvetica', ...
             'FontSize', 8, ...
             'FontWeight', 'normal', ...
             'FontAngle', 'normal', ...
             'Box','off',...
             'Location', 'eastoutside',...
             'Interpreter', 'tex');

uistack(hGrowth, 'top');

set(gca, 'LooseInset', get(gca, 'TightInset'));

% Export as vector PDF for journal submission
exportgraphics(gcf, fullfile('results', ['storage_assignment_minQ=',num2str(minQ),'_maxQ=',num2str(maxQ),'_sites=',num2str(nr_well_max),'_',num2str(time_yr),'y','.pdf']), 'ContentType','vector');
exportgraphics(gcf, fullfile('results', ['storage_assignment_minQ=',num2str(minQ),'_maxQ=',num2str(maxQ),'_sites=',num2str(nr_well_max),'_',num2str(time_yr),'y','.png']), ...
    'Resolution', 600, 'BackgroundColor', 'white');
saveas(gcf, fullfile('results', ['storage_assignment_minQ=',num2str(minQ),'_maxQ=',num2str(maxQ),'_sites=',num2str(nr_well_max),'_',num2str(time_yr),'y','.svg']));


%% Assignment complete?
if storage_resource_total(end) < max(total_rate)
    disp('Storage resources is not enough to accomodate the projected storage rates in this scenario')
    disp(['peak rate achieved = ',num2str(storage_resource_total(end))])
end
