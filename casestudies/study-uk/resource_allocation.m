function [Resource_assignment_table, ...
          storage_resource_total, ...
          step_curve_start, ...
          step_curve_end, ...
          region_name_array_sorted, ...
          site_no_step] = ...
          resource_allocation( ...
          Growthrate_model, ...
          sort_idx,...
          Rate_opt_all_scenarios, ...
          Site_no_all_scenarios, ...
          storage_periods, ...
          allocation_duration, ...
          allocation_order, ...
          minQ, maxQ, ...
          nr_region, ...
          fpath, fname, correction, ...
          dist_min, dist_max, nr_dist, ...
          nr_well_max, rw, ...
          outputDir)

%% Growth model

years = Growthrate_model(:,1);
total_rate = Growthrate_model(:,2);

%% Input matrices
% Convert to numeric matrix of rates of resource use
% [nRegions x nPeriods]
resource_rate_matrix = cell2mat(Rate_opt_all_scenarios(:,3:end));

region_name_array    = Rate_opt_all_scenarios(:,2);
region_no_array      = Rate_opt_all_scenarios(:,1);

% Convert to numeric matrix of the number of injection sites
% [nRegions x nPeriods]
Site_no_matrix = cell2mat(Site_no_all_scenarios(:,3:end));

%% Sort storage regions
% Sort regions by total storage capacity (sum of sustainable rates over all periods)
switch allocation_order
    case 'descend'
        [~, sort_idx] = sort(sum(resource_rate_matrix,2), 'descend');

    case 'ascend'
        [~, sort_idx] = sort(sum(resource_rate_matrix,2), 'ascend');

    case 'random'
        sort_idx = randperm(size(resource_rate_matrix,1));

    case 'manual'
        sort_idx = sort_idx;

    otherwise
        error('allocation_order must be ''descend'', ''ascend'', ''manual'', or ''random''.');
end

resource_rate_matrix_sorted = resource_rate_matrix(sort_idx,:);
region_name_array_sorted    = region_name_array(sort_idx);
region_no_array_sorted      = region_no_array(sort_idx);
Site_no_matrix_sorted       = Site_no_matrix(sort_idx,:);

%% Preallocation
% storage resource allocated in the current step (Mt/y)
storage_resource_step  = zeros(nr_region,1);
storage_resource_total = zeros(nr_region,1);
site_no_step           = zeros(nr_region,1);
period_step            = zeros(nr_region,1);
cum_vol                = zeros(nr_region,1);

step_start       = nan(nr_region,1);
step_end         = nan(nr_region,1);
step_curve_start = nan(nr_region,1);
step_curve_end   = nan(nr_region,1);

Resource_assignment = cell(nr_region,11);

%% Initial remaining demand
% remaining time to allocate storage resources
t_remaining = allocation_duration;

% remaining rate from the growth rate model (-: allocated; +: need to be allocated in the next steps)
% a vector of the same size as the growth rate model
total_rate_remain = total_rate;

%% Allocation loop

for allocation_step = 1:nr_region
    % index of first storage period > remaining time to allocate
    idx = find(storage_periods >= t_remaining, 1, 'first');

    if isempty(idx)
        idx = numel(storage_periods);
    end

    % possible cumulative injection volumes for the unit based on the
    % injection duration adn rate
    vol_possible = storage_periods(1:idx) .* ...
                   resource_rate_matrix_sorted(allocation_step,1:idx);

    breaking = 0;
    break_max_rate = 0;

    for i = idx:-1:1

        current_rate = resource_rate_matrix_sorted(allocation_step,i);

        % remaining rate to be assigned to storage resources for this
        % provisional injection rate
        total_rate_remain_test = total_rate_remain - current_rate;

        % Find leftmost index where total_rate >= storage_resource
        left_idx  = find(total_rate_remain_test >= 0, 1, 'first');

        % Find rightmost index where total_rate >= storage_resource
        right_idx = find(total_rate_remain_test >= 0, 1, 'last');

        % cumulative storage rate (Mt/y)
        if allocation_step == 1
            storage_resource_total_provisional = current_rate;
        else
            storage_resource_total_provisional = ...
                storage_resource_total(allocation_step-1) + current_rate;
        end

        %% Find intersection between cumulative rate and growth curve

        % finding times at which the growth curve is intersected by the 
        % considered injection rate
        if ~isempty(left_idx)   % this is to exclude the curve of the last assignment
            % start point (year) of the curve of this step
            if years(left_idx) > 2030
                step_curve_start_provisional = ...
                    years(left_idx-1) + ...
                    (storage_resource_total_provisional - total_rate(left_idx-1)) * ...
                    (years(left_idx) - years(left_idx-1)) / ...
                    (total_rate(left_idx) - total_rate(left_idx-1));
            else
                step_curve_start_provisional = 2030;
            end

            % end point (year) of the curve of this step
            if years(right_idx) < years(end) 
                step_curve_end_provisional = ...
                    years(right_idx) + ...
                    (storage_resource_total_provisional - total_rate(right_idx)) * ...
                    (years(right_idx+1) - years(right_idx)) / ...
                    (total_rate(right_idx+1) - total_rate(right_idx));
            else
                step_curve_end_provisional = years(end);
            end
        end

        %% Calculate required volume under growth model

        % time and injection rates above the considered rate
        % used for the calculation of area beneath the growth curve
        if ~isempty(left_idx)

            area_time = [step_curve_start_provisional; ...
                         years(left_idx:right_idx); ...
                         step_curve_end_provisional];

            area_rate = [0; ...
                         total_rate_remain_test(total_rate_remain_test >= 0); ...
                         0];

            % calculating the 
            % area beneath the growth curve for the considered
            % injection rate
            if allocation_step == 1
                area_growthrate = trapz(years,total_rate) - trapz(area_time,area_rate);
            else
                area_growthrate = trapz(years,total_rate) ...
                                - cum_vol(allocation_step-1) ...
                                - trapz(area_time,area_rate);
            end

        else

            if allocation_step == 1
                area_growthrate = trapz(years,total_rate);
            else
                area_growthrate = trapz(years,total_rate) ...
                                - cum_vol(allocation_step-1);
            end
        end

        %% Check whether this region can satisfy required volume

        % storage resource allocated in the current step (Mt/y)
        if (vol_possible(i) > area_growthrate || ...
            abs(vol_possible(i)-area_growthrate) < 0.01) && ...
            current_rate >= storage_resource_step(allocation_step)

            % number of sites to achieve the desired rate
            site_no_step(allocation_step) = ...
                Site_no_matrix_sorted(allocation_step,i);

            storage_resource_step(allocation_step) = current_rate;

            storage_resource_total(allocation_step) = ...
                storage_resource_total_provisional;

            % real period of the operation of the assigned storage resource
            period_step(allocation_step) = storage_periods(i);

            if ~isempty(left_idx)
                step_curve_start(allocation_step) = step_curve_start_provisional;
                step_curve_end(allocation_step)   = step_curve_end_provisional;
            end

            % cumulative volume allocated
            if allocation_step == 1
                cum_vol(allocation_step) = area_growthrate;
            else
                cum_vol(allocation_step) = ...
                    cum_vol(allocation_step-1) + area_growthrate;
            end

            %% Final-region correction if cumulative rate exceeds model maximum

            % check if the assigned rate exceeded the maximum rate of the
            % model
            if storage_resource_total_provisional >= max(total_rate)

                break_max_rate = 1;

                % final rate remaining to be addressed (the last region)
                if allocation_step == 1
                    remaining_rate_final = max(total_rate);
                else
                    remaining_rate_final = ...
                        max(total_rate) - storage_resource_total(allocation_step-1);
                end

                % repeat storage resource calculations for the last
                % allocated region to find the minimum required sites to
                % address the remaining injection rate and volume
                [~,~,well_list,~,~,~,Table_Q,Table_V,~] = calculate( ...
                    fpath, fname, sort_idx(allocation_step), correction, ...
                    dist_min, dist_max, nr_dist, nr_well_max, rw, ...
                    storage_periods(i), maxQ);

                % find maximum volume/rate for each of the injection scenarios
                % Convert table to array (numerical part only)
                dataArray_Q       = table2array(Table_Q);
                dataArray_Q_total = table2array(well_list'.*Table_Q);
                dataArray_Vol     = table2array(Table_V);

                % Find where Q < minQ
                % rows and columns include different site arrangements and
                % distances respectively
                mask = dataArray_Q < minQ;

                % makig single injection site exempt from the minQ criterion
                mask(1,1) = false;

                % Set masked entries to zero
                dataArray_Q_total(mask) = 0;
                dataArray_Vol(mask)     = 0;

                % Logical mask for valid scenarios
                valid_mask = ...
                    dataArray_Q_total > remaining_rate_final & ...
                    dataArray_Vol * 1000 > area_growthrate;

                % Find all valid indices
                [row_idx, col_idx] = find(valid_mask);

                if isempty(row_idx)

                    warning('No valid scenario found for the final allocated region.');

                else

                    % Find smallest row index (minimum number of sites)
                    min_row = min(row_idx);

                    % Among those rows, take the first valid distance
                    first_col = col_idx(find(row_idx == min_row,1));

                    Q_last = dataArray_Q(min_row,first_col) * well_list(min_row);

                    site_no_step(allocation_step) = well_list(min_row);
                    storage_resource_step(allocation_step) = Q_last;

                    if allocation_step == 1
                        storage_resource_total(allocation_step) = Q_last;
                    else
                        storage_resource_total(allocation_step) = ...
                            storage_resource_total(allocation_step-1) + Q_last;
                    end
                end
            end

            if isempty(left_idx)
                breaking = 1;
            end
        end

        if break_max_rate == 1
            break
        end
    end

    %% Update remaining growth-rate demand

    total_rate_remain = total_rate_remain - storage_resource_step(allocation_step);

    %% Assignment start/end time

    % time frame of the storage use for each step
    % time starting and ending injection into resources of each step
    if allocation_step == 1
        step_start(allocation_step) = 2030;
        step_end(allocation_step)   = 2030 + allocation_duration - 1;
    else
        step_start(allocation_step) = step_curve_start(allocation_step-1);
        step_end(allocation_step)   = step_curve_end(allocation_step-1);
    end

    %% Store assignment table row

    warning('off','MATLAB:table:ModifiedAndSavedVarnames');

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
    Resource_assignment{allocation_step,1} = allocation_step;
    Resource_assignment{allocation_step,2} = region_no_array_sorted(allocation_step);
    Resource_assignment{allocation_step,3} = region_name_array_sorted{allocation_step};
    Resource_assignment{allocation_step,4} = step_start(allocation_step);
    Resource_assignment{allocation_step,5} = step_end(allocation_step);
    Resource_assignment{allocation_step,6} = period_step(allocation_step);
    Resource_assignment{allocation_step,7} = storage_resource_step(allocation_step);
    Resource_assignment{allocation_step,8} = storage_resource_total(allocation_step);
    Resource_assignment{allocation_step,9} = site_no_step(allocation_step);

    if breaking ~= 1

        Resource_assignment{allocation_step,10} = step_curve_start(allocation_step);
        Resource_assignment{allocation_step,11} = step_curve_end(allocation_step);

        % time remained to be assigned to storage resources
        t_remaining = step_curve_end(allocation_step) ...
                    - step_curve_start(allocation_step) + 1;

    else
        break;
    end
end

%% Trim unused rows

last_step = allocation_step;

storage_resource_total = storage_resource_total(1:last_step);
site_no_step           = site_no_step(1:last_step);
step_curve_start       = step_curve_start(1:last_step);
step_curve_end         = step_curve_end(1:last_step);
region_name_array_sorted = region_name_array_sorted(1:last_step);

Resource_assignment = Resource_assignment(1:last_step,:);

%% Convert to table

% Convert to table and set column names
Resource_assignment_table = cell2table(Resource_assignment, ...
    'VariableNames', {'Step', ...
                      'Region no', ...
                      'Region name', ...
                      'Start [y]', ...
                      'End [y]', ...
                      'Duration [y]', ...
                      'Step rate increment [Mt/y]', ...
                      'Step rate cumulative [Mt/y]', ...
                      'No_sites', ...
                      'Curve start [y]', ...
                      'Curve end [y]'});

%% Save table

writetable(Resource_assignment_table, ...
    fullfile(outputDir, ...
    ['Resource_assignment_minQ=',num2str(minQ), ...
     '_maxQ=',num2str(maxQ), ...
     '_sites=',num2str(nr_well_max), ...
     '_',num2str(allocation_duration),'y.xlsx']));

end