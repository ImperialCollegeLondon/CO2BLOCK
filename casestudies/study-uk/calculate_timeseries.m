function [site_name,d_list,well_list,d_max,Q_all,V_all,p_sup_all]...
    = calculate_timeseries(fpath,fname,site_no,correction,dist_min,dist_max,nr_dist,nr_well_max,rw,time_yr_vec,maxQ)

    %read data
    [site_name,thick,area_res,perm,por,dens_c,visc_c,visc_w,compr,p_lim,rc,gamma,delta, omega] = read_data(fpath,fname,site_no);
    
    n_time = length(time_yr_vec);

    % initialize
    if strcmp(nr_well_max,'auto')                                               % calculate maximum well number if not set
        nr_well_max = floor(area_res/(dist_min^2));
    end
    
    if strcmp(dist_max,'auto')                                                  % calculate maximum interwell distance if not set
        dist_max = sqrt(2*area_res)/2;
    end
    
    d_list = linspace(dist_min,dist_max,nr_dist) ;                          % inter-well distance list
    M0 = perm/1e-13 ;                                                      % guess value for total injection rate [Mton/y] 

    if isfinite(rc)
        % quadrature points inside equivalent circular reservoir
        n_lookup_r = 40;       % radial positions for source well
        n_quad_r   = 35;       % radial integration points
        n_quad_th  = 60;      % angular integration points

        rr = rc * sqrt(((1:n_quad_r)-0.5)/n_quad_r);
        tt = linspace(0,2*pi,n_quad_th+1);
        tt(end) = [];

        [RR,TT] = meshgrid(rr,tt);
        sample_x = RR(:).*cos(TT(:));
        sample_y = RR(:).*sin(TT(:));
        sample_weight = ones(size(sample_x))/numel(sample_x);

        % lookup distance may exceed rc because well arrays can extend outside equivalent circle
        max_lookup_d = sqrt(2) * dist_max*1000 * sqrt(nr_well_max);
        lookup_d = linspace(0,max_lookup_d,n_lookup_r);
    end

% -------------------------------------------------------------------------
%            building all static well/distance geometry (only once)
% -------------------------------------------------------------------------

    well_list = [];   
    d_max = [];   
    w_id = 0 ;
    Q0_vec = [];
    well_geom = struct();                                                   % saving x-grid, y-grid, number of sites (w) and distance vector for each configuration

    for x_grid_num = 1:sqrt(nr_well_max)                                    % number of wells on a horizontal row
        if x_grid_num*(x_grid_num+1) < nr_well_max
            plus = 1;
        else
            plus = 0;
        end
        for y_grid_num = x_grid_num:x_grid_num+plus                         % number of wells on a  vertical row
            w_id = w_id +1 ;                                                % well number scenario ID
            w = x_grid_num*y_grid_num ;                                     % well number for each scenario

            well_list(w_id) = w;                                            % store in vector
            d_max(w_id) = sqrt(area_res/w);                                 % maximum interwell distance for each well number [km]
            
            Q0_vec(w_id) = M0*1e9/dens_c/365/86400/w;                                % injection rate per well [m3/s]

            well_geom(w_id).x_grid_num = x_grid_num;
            well_geom(w_id).y_grid_num = y_grid_num;
            well_geom(w_id).w = w;
            
            for d = 1:nr_dist
                distance = d_list(d)*1000;

                % calculate distances
                wells_coord_x = repmat((0:distance:distance*x_grid_num-1),[y_grid_num,1]);
                wells_coord_y = repmat(transpose(0:distance:distance*y_grid_num-1),[1,x_grid_num]);

                central_well_x = ceil(x_grid_num/2);     % posiiton of the central well in x-coord vector
                central_well_y = ceil(y_grid_num/2);     % posiiton of the central well in y-coord vector

                dist_vec_x = wells_coord_x - wells_coord_x(central_well_y,central_well_x);        % distance in x from central well [km]
                dist_vec_y = wells_coord_y - wells_coord_y(central_well_y,central_well_x) ;       % distance in y from central well [km]
                dist_vec   = sqrt(dist_vec_x.^2+dist_vec_y.^2) ;                                  % distance from central well [m]
                dist_vec(central_well_y,central_well_x) = rw ;                                    % assign wells radius to the central well

                well_geom(w_id).dist_vec{d} = dist_vec;

                if isfinite(rc)
                    % coordinates relative to geometric center of well array
                    array_center_x = mean(wells_coord_x(:));
                    array_center_y = mean(wells_coord_y(:));

                    well_x_domain = wells_coord_x - array_center_x;
                    well_y_domain = wells_coord_y - array_center_y;

                    well_radial_pos = sqrt(well_x_domain.^2 + well_y_domain.^2);

                    well_geom(w_id).well_radial_pos{d} = well_radial_pos;   
                end
            end
        end
    end
          
    % allocate outputs
    Q_all = nan(w_id,nr_dist,n_time);
    V_all = nan(w_id,nr_dist,n_time);
    p_sup_all = nan(w_id,nr_dist,n_time);

% -------------------------------------------------------------------------
%                            Time loop
% -------------------------------------------------------------------------
    for it = 1:n_time

        time_yr = time_yr_vec(it);       
        time = time_yr*86400*365;              %injection time [sec]

        R_influence = sqrt(2.246*perm*time/(visc_w*compr));   % pressure propagation radius for the time of injection

        b = zeros(1,w_id);
        p_sup_vec = zeros(w_id,nr_dist);

        % doing pressure calculations
        for nn = 1:w_id
            w = well_list(nn);
            Q0 = Q0_vec(nn);

            csi = sqrt(Q0*time/pi/por/thick);                              % average plume extension [m]
            psi = exp(omega)*csi;                                          % equivalent plume extension [m]
            p_c = (Q0*visc_w)/(2*pi*thick*perm)/1e6;                       % characteristic pressure  [MPa]

            if isfinite(rc)
                % lookup table for average open-boundary Nordbotten pressure
                % built once per well-number scenario, reused for all distances
                avg_Nord_lookup = zeros(size(lookup_d));
    
                for kk = 1:length(lookup_d)
    
                    xw = lookup_d(kk);
                    yw = 0;
    
                    r_samples = sqrt((sample_x-xw).^2 + (sample_y-yw).^2);
                    r_samples(r_samples < rw) = rw;
    
                    p_samples_unit = zeros(size(r_samples));
    
                    for jj = 1:length(r_samples)
                        % use open-boundary Nordbotten, so rc = inf
                        p_samples_unit(jj) = Nordbotten_solution( ...
                            r_samples(jj), R_influence, psi, inf, gamma);
                    end
    
                    avg_Nord_lookup(kk) = sum(p_samples_unit .* sample_weight);
                end
            end

            for d = 1:nr_dist
                distance = d_list(d)*1000;
                dist_vec = well_geom(nn).dist_vec{d};

                if isfinite(rc)
                    well_radial_pos = well_geom(nn).well_radial_pos{d};

                    % global closed-tank pressure
                    p_tank = (w*Q0*time)/(area_res*1e6*thick*compr)/1e6;     % [MPa]
    
                    p_sup = p_tank; 

                    for i = 1:w
                        r = dist_vec(i);
    
                        % open-boundary Nordbotten pressure at central well
                        p_Nord_center = Nordbotten_solution( ...
                            r, R_influence, psi, inf, gamma) * p_c;   % overpressure according to Nordbotten and Celia solution for overpressure [MPa]
    
                        % average open-boundary Nordbotten pressure over equivalent circular domain
                        d_source = well_radial_pos(i);
    
                        avg_Nord_unit = interp1( ...
                            lookup_d, avg_Nord_lookup, d_source, ...
                            'linear', 'extrap');
    
                        p_Nord_avg = avg_Nord_unit * p_c;
    
                        % add only zero-mean local contribution
                        p_sup = p_sup + (p_Nord_center - p_Nord_avg);
                    end

                    sup_error = 0 ;                                   % not considering superposition error for the case of closed boundary conditions
                    b(nn) = (visc_w-visc_c)/4/pi/perm/thick;  

                else
                    p_sup = 0;

                    for i = 1:w
                        r = dist_vec(i);
                        Delta_p = Nordbotten_solution(r,R_influence,psi,rc,gamma)*p_c;       % overpressure according to Nordbotten and Celia solution for overpressure [MPa]
                        p_sup = p_sup +  Delta_p ;                                 % superposed overpressure [MPa]
                    end
    
                    switch correction                                              % correction for superposition error (De Simone et al., GRL2019)
                        case 'off'
                            sup_error = 0 ;
                            b(nn) = (visc_w-visc_c)/4/pi/perm/thick  ;                   
                        case 'on'  
                            if w < 9 ||  R_influence*csi/(distance^2) < 1
                                sup_error = 0; 
                                b(nn) = (visc_w-visc_c)/4/pi/perm/thick  ;
                            else      
                                sup_error = w*delta/4 * log(R_influence*csi/(distance^2)) ;  
                                b(nn) = (visc_w-visc_c)/4/pi/perm/thick * (1+w/4) ; 
                            end
                    end
                end

                p_sup =  p_sup - sup_error*p_c ; 
                p_sup_vec(nn,d) =  p_sup  ;          
            end
        end

        % calculate injectable per well flow-rate Q_M_each and total...
        % storage capacity V_M for each scenario for this time
        b_mat = repmat(b',1,nr_dist);   
        q1 = repmat(Q0_vec',1,nr_dist); 
        well_mat = repmat(well_list',1,nr_dist); 
    
        p1 = p_sup_vec*1e6;
        p2 = p_lim*1e6;   

        q2 = - p2./b_mat./(lambertw(-1,-p2./q1./b_mat.*exp(-p1./q1./b_mat))) ;     % [m3/s] limit flow rate at each well with  non-linear multi-phase relationship
        Q_M_each = q2*86400*365*dens_c/1e9 ;                           % [ Mton/year]   limit flow rate at each well with non-linear multi-phase relationship 

        %rescale according to lower contstraints
        for dd = 1:nr_dist
            distance= d_list(dd)*1000;
            if Q_M_each(1,dd) > 0.9999*maxQ            % rescale for n_well = 1
               Q_M_each(1,dd) =  0.9999*maxQ ;
            end
            for nn = 2:w_id                             % rescale for n_well > 1
                if Q_M_each(nn,dd) > 0.9999*maxQ || q2(nn,dd) > (distance)^2*pi*por*thick/4.001/time 
                   Q_M_each(nn,dd) = min( 0.9999*maxQ, (distance)^2*pi*por*thick/4.001/time*86400*365*dens_c/1e9) ;
                end
            end
        end
    
        Q_M_tot = Q_M_each.*well_mat ;                     % [ Mton/year]  total sustainable flow rate   
        V_M = Q_M_tot.*time_yr/1000 ;                      % [Gton]        total sustainable injectable mass
    
        % upper constraint 
        d_mat = repmat(d_list,w_id,1); 
        d_max_mat = repmat(d_max,nr_dist,1); 
        d_max_check =  d_mat./d_max_mat' ;       
    
        % apply upper constraint
        possible = d_max_check < 1 ; 
        Q_poss_each = Q_M_each.* possible;
        V_poss = V_M .* possible;

        Q_all(:,:,it) = Q_poss_each;
        V_all(:,:,it) = V_poss;
        p_sup_all(:,:,it) = p_sup_vec;
    end
  
end
