%------------------------------------------------------------------------%
% Billy Orion Mazotti
% Department of Aerospace Engineering, University of Michigan, Ann Arbor
% Laboratory for Air Transportation, Infrastructure, and Connected 
% Environments (LATTICE)
%------------------------------------------------------------------------%

%{

Preparation Tasks:
1. Download data from OpenStreetMap as shape (.shp) files (recommend using
MyGeoData for ease of data conversion)
2. Use OSM_2_Formatted.m to convert raw data .shp files to .mat files used
for anlaysis in main.m

Running Steps:
1. Load Formatted Datasets 

%}

%% LOAD FORMATTED DATASETS
clear
clc
close all

%----------------------- SELECT Dataset Variables -----------------------%
city_choices = ["San Francisco","New York City","Los Angeles"];
city = city_choices(1);
altitude = 200;
%------------------------------------------------------------------------%



if city == "San Francisco"

    disp("San Francisco Datasets Loading ...")

    %------------------- SAN FRANCISCO DATASETS -------------------%
    %--- sensor locations ---%
    %building sensor location
    load_buildings_1 = load('FormattedDatasets/SF_formatted/2_3_23_SF_buidings_1_sensors_N160927.mat');
    S_building_sens_1 = load_buildings_1.S_out;
    xy_data_1 = [S_building_sens_1.XLocation,S_building_sens_1.YLocation];
    
    %remove duplicate coordiantes
    xy_data_all = [[xy_data_1(:,1)],[xy_data_1(:,2)]];
    xy_data_unique = unique(xy_data_all,'rows');
    
    S_building_sens.XLocation = xy_data_unique(:,1);
    S_building_sens.YLocation = xy_data_unique(:,2);
    
    
    %--- start/stop locations ---%
    %customer/stop locations
    load_buildings = load('FormattedDatasets/SF_formatted/2_3_23_SF_customers_N125446.mat');
    S_customer_stop = load_buildings.S_out;
    xy_data = cell2mat(struct2cell(S_customer_stop));
    xyCustomers = [xy_data(1:size(xy_data,1)/2),xy_data(1+size(xy_data,1)/2:end)];
    %vendor/start locations
    load_buildings = load('FormattedDatasets/SF_formatted/2_3_23_SF_vendors_N159.mat');
    S_vendor_start = load_buildings.S_out;
    xy_data = cell2mat(struct2cell(S_vendor_start));
    xyVendors = [xy_data(1:size(xy_data,1)/2),xy_data(1+size(xy_data,1)/2:end)];
    
    if altitude == 200
        %--- occupancy data ---%
        load_buildings_1 = load('FormattedDatasets/SF_formatted/2_3_23_SF_buidings_1_occupancy_H200ft_N154.mat');
        S_building_occ_1 = load_buildings_1.S_out;
        xy_data_1 = cell2mat(struct2cell(S_building_occ_1));
        
        xyBuildings_old = [xy_data_1(1:size(xy_data_1,1)/2),xy_data_1(1+size(xy_data_1,1)/2:end)];
        xyBuildings = unique(xyBuildings_old,'rows');
    end
    if altitude == 400
        %--- occupancy data ---%
        load_buildings_1 = load('FormattedDatasets/SF_formatted/2_3_23_SF_buidings_1_occupancy_H400ft_N48.mat');
        S_building_occ_1 = load_buildings_1.S_out;
        xy_data_1 = cell2mat(struct2cell(S_building_occ_1));
        
        xyBuildings_old = [xy_data_1(1:size(xy_data_1,1)/2),xy_data_1(1+size(xy_data_1,1)/2:end)];
        xyBuildings = unique(xyBuildings_old,'rows');
    end
    

    %--- city contours data ---%
    load_struct = load('FormattedDatasets/SF_formatted/SF_border_convex.mat');
    S_contour_convex = load_struct.S_out;
    load_struct = load('FormattedDatasets/SF_formatted/SF_border.mat');
    S_contours = load_struct.S_out;
    S_contours(2).contour = S_contours(1).contour;
    size(struct2table(S_contours),1);
    
    minX_map_lat = -130;
    minY_map_long = 0;
    lat2meters = 10000000/90;
    long2meters = 40075161.2/360;
    minX_map_m = minX_map_lat * lat2meters;
    minY_map_m = minY_map_long * long2meters;
    
    S_contour_convex.contour(:,1) = round(S_contour_convex.contour(:,1)*lat2meters - minX_map_m);
    S_contour_convex.contour(:,2) = round(S_contour_convex.contour(:,2)*long2meters - minY_map_m);
    
    for idx = 1:1:size(struct2table(S_contours),1)
        S_contours(idx).contour(:,1) = round(S_contours(idx).contour(:,1)*lat2meters - minX_map_m);
        S_contours(idx).contour(:,2) = round(S_contours(idx).contour(:,2)*long2meters - minY_map_m);
    end
    disp("San Francisco Datasets Loading Complete")


elseif city == "Los Angeles"

    disp("Los Angeles Datasets Loading ...")
    
    % ------------------- LOS ANGELES DATASETS ------------------- %
    %--- sensor locations ---%
    %building sensor location
    load_buildings_1 = load('FormattedDatasets/LA_formatted/2_2_23_LA_buidings_1_sensors_N538415.mat');
    S_building_sens_1 = load_buildings_1.S_out;
    xy_data_1 = [S_building_sens_1.XLocation,S_building_sens_1.YLocation];
    load_buildings_2 = load('FormattedDatasets/LA_formatted/2_2_23_LA_buidings_2_sensors_N542339.mat');
    S_building_sens_2 = load_buildings_2.S_out;
    xy_data_2 = [S_building_sens_2.XLocation,S_building_sens_2.YLocation];
    load_buildings_3 = load('FormattedDatasets/LA_formatted/2_2_23_LA_buidings_3_sensors_N695004.mat');
    S_building_sens_3 = load_buildings_3.S_out;
    xy_data_3 = [S_building_sens_3.XLocation,S_building_sens_3.YLocation];
    load_buildings_4 = load('FormattedDatasets/LA_formatted/2_2_23_LA_buidings_4_sensors_N227406.mat');
    S_building_sens_4 = load_buildings_4.S_out;
    xy_data_4 = [S_building_sens_4.XLocation,S_building_sens_4.YLocation];
    
    %remove duplicate coordiantes
    xy_data_all = [[xy_data_1(:,1);xy_data_2(:,1); ...
                    xy_data_3(:,1);xy_data_4(:,1)], ...
                    [xy_data_1(:,2);xy_data_2(:,2); ...
                    xy_data_3(:,2);xy_data_4(:,2)]];
    xy_data_unique = unique(xy_data_all,'rows');
    
    S_building_sens.XLocation = xy_data_unique(:,1);
    S_building_sens.YLocation = xy_data_unique(:,2);
    
    
    %--- start/stop locations ---%
    %customer/stop locations
    load_buildings = load('FormattedDatasets/LA_formatted/2_2_23_LA_customers_N360183.mat');
    S_customer_stop = load_buildings.S_out;
    xy_data = cell2mat(struct2cell(S_customer_stop));
    xyCustomers = [xy_data(1:size(xy_data,1)/2),xy_data(1+size(xy_data,1)/2:end)];
    %vendor/start locations
    load_buildings = load('FormattedDatasets/LA_formatted/2_3_23_LA_vendors_N453.mat');
    S_vendor_start = load_buildings.S_out;
    xy_data = cell2mat(struct2cell(S_vendor_start));
    xyVendors = [xy_data(1:size(xy_data,1)/2),xy_data(1+size(xy_data,1)/2:end)];
    
    
    %--- occupancy data ---%
    load_buildings_1 = load('FormattedDatasets/LA_formatted/2_2_23_LA_buidings_1_occupancy_H200ft_N20.mat');
    S_building_occ_1 = load_buildings_1.S_out;
    xy_data_1 = cell2mat(struct2cell(S_building_occ_1));
    load_buildings_2 = load('FormattedDatasets/LA_formatted/2_2_23_LA_buidings_2_occupancy_H200ft_N224.mat');
    S_building_occ_2 = load_buildings_2.S_out;
    xy_data_2 = cell2mat(struct2cell(S_building_occ_2));
    load_buildings_3 = load('FormattedDatasets/LA_formatted/2_2_23_LA_buidings_3_occupancy_H200ft_N144.mat');
    S_building_occ_3 = load_buildings_3.S_out;
    xy_data_3 = cell2mat(struct2cell(S_building_occ_3));
    load_buildings_4 = load('FormattedDatasets/LA_formatted/2_2_23_LA_buidings_4_occupancy_H200ft_N21.mat');
    S_building_occ_4 = load_buildings_4.S_out;
    xy_data_4 = cell2mat(struct2cell(S_building_occ_4));
    
    xyBuildings_old = [xy_data_1(1:size(xy_data_1,1)/2),xy_data_1(1+size(xy_data_1,1)/2:end);
                    xy_data_2(1:size(xy_data_2,1)/2),xy_data_2(1+size(xy_data_2,1)/2:end);
                    xy_data_3(1:size(xy_data_3,1)/2),xy_data_3(1+size(xy_data_3,1)/2:end);
                    xy_data_4(1:size(xy_data_4,1)/2),xy_data_4(1+size(xy_data_4,1)/2:end)];
    xyBuildings = unique(xyBuildings_old,'rows');
    
    
    %--- city contours data ---%
    load_struct = load('FormattedDatasets/LA_formatted/LA_border_convex.mat');
    S_contour_convex = load_struct.S_out;
    load_struct = load('FormattedDatasets/LA_formatted/LA_border.mat');
    S_contours(1).contour = load_struct.S_out(1).contour;
    S_contours(2).contour = load_struct.S_out(1).contour;
    
    minX_map_lat = -130;
    minY_map_long = 0;
    lat2meters = 10000000/90;
    long2meters = 40075161.2/360;
    minX_map_m = minX_map_lat * lat2meters;
    minY_map_m = minY_map_long * long2meters;
    
    S_contour_convex.contour(:,1) = round(S_contour_convex.contour(:,1)*lat2meters - minX_map_m);
    S_contour_convex.contour(:,2) = round(S_contour_convex.contour(:,2)*long2meters - minY_map_m);
    for idx = 1:1:size(struct2table(S_contours),1)
        S_contours(idx).contour(:,1) = round(S_contours(idx).contour(:,1)*lat2meters - minX_map_m);
        S_contours(idx).contour(:,2) = round(S_contours(idx).contour(:,2)*long2meters - minY_map_m);
    end
    disp("Los Angeles Datasets Loading Complete")

    
elseif city == "New York City"

    disp("New York City Datasets Loading ...")

    % ------------------- NEY YORK CITY DATASETS ------------------- %
    %--- sensor locations ---%
    %building sensor location
    load_buildings_1 = load('FormattedDatasets/NYC_formatted/2_3_23_NYC_buidings_1_sensors_N685639.mat');
    S_building_sens_1 = load_buildings_1.S_out;
    xy_data_1 = [S_building_sens_1.XLocation,S_building_sens_1.YLocation];
    load_buildings_2 = load('FormattedDatasets/NYC_formatted/2_3_23_NYC_buidings_2_sensors_N688024.mat');
    S_building_sens_2 = load_buildings_2.S_out;
    xy_data_2 = [S_building_sens_2.XLocation,S_building_sens_2.YLocation];
    
    %remove duplicate coordiantes
    xy_data_all = [[xy_data_1(:,1);xy_data_2(:,1)], ...
                    [xy_data_1(:,2);xy_data_2(:,2)]];
    
    xy_data_unique = unique(xy_data_all,'rows');
    S_building_sens.XLocation = xy_data_unique(:,1);
    S_building_sens.YLocation = xy_data_unique(:,2);
    
    
    %--- start/stop locations ---%
    %customer/stop locations
    load_buildings = load('FormattedDatasets/NYC_formatted/2_3_23_NYC_customers_N958747.mat');
    S_customer_stop = load_buildings.S_out;
    xy_data = cell2mat(struct2cell(S_customer_stop));
    xyCustomers = [xy_data(1:size(xy_data,1)/2),xy_data(1+size(xy_data,1)/2:end)];
    %vendor/start locations
    load_buildings = load('FormattedDatasets/NYC_formatted/2_3_23_NYC_vendors_N838.mat');
    S_vendor_start = load_buildings.S_out;
    xy_data = cell2mat(struct2cell(S_vendor_start));
    xyVendors = [xy_data(1:size(xy_data,1)/2),xy_data(1+size(xy_data,1)/2:end)];
    
    
    %--- occupancy data ---%
    load_buildings_1 = load('FormattedDatasets/NYC_formatted/2_3_23_NYC_buidings_1_occupancy_H200ft_N1388.mat');
    S_building_occ_1 = load_buildings_1.S_out;
    xy_data_1 = cell2mat(struct2cell(S_building_occ_1));
    load_buildings_2 = load('FormattedDatasets/NYC_formatted/2_3_23_NYC_buidings_2_occupancy_H200ft_N55.mat');
    S_building_occ_2 = load_buildings_2.S_out;
    xy_data_2 = cell2mat(struct2cell(S_building_occ_2));
    
    xyBuildings_old = [xy_data_1(1:size(xy_data_1,1)/2),xy_data_1(1+size(xy_data_1,1)/2:end);
                    xy_data_2(1:size(xy_data_2,1)/2),xy_data_2(1+size(xy_data_2,1)/2:end)];
    xyBuildings = unique(xyBuildings_old,'rows');
    
    
    %--- city contours data ---%
    load_struct = load('FormattedDatasets/NYC_formatted/NYC_border_convex.mat');
    S_contour_convex = load_struct.S_out;
    load_struct = load('FormattedDatasets/NYC_formatted/NYC_border.mat');
    S_contours = load_struct.S_out;
    
    minX_map_lat = -100;
    minY_map_long = 15;
    lat2meters = 10000000/90;
    long2meters = 40075161.2/360;
    minX_map_m = minX_map_lat * lat2meters;
    minY_map_m = minY_map_long * long2meters;
    
    S_contour_convex.contour(:,1) = round(S_contour_convex.contour(:,1)*lat2meters - minX_map_m);
    S_contour_convex.contour(:,2) = round(S_contour_convex.contour(:,2)*long2meters - minY_map_m);
    for idx = 1:1:size(struct2table(S_contours),1)
        S_contours(idx).contour(:,1) = round(S_contours(idx).contour(:,1)*lat2meters - minX_map_m);
        S_contours(idx).contour(:,2) = round(S_contours(idx).contour(:,2)*long2meters - minY_map_m);
    end
    disp("New York City Datasets Loading Complete")
end


%---------------------- VIEW SIMULATED ENVIORNMENT ----------------------%

disp("Normalizing Datasets ...")
[entireMapEdges_local,xyBuildings_local,xyCustomers_local_unfiltered, ...
    xyVendors_local_unfiltered, S_building_sens_local, ...
    S_contour_convex_local,S_contours_local] ...
                    = NormalizeData(S_building_sens, xyCustomers, ...
                                        xyVendors,xyBuildings,S_contour_convex,S_contours);
disp("Normalizing Complete!")


disp("Generating Occupancy Map ...")
% max display map HxW ~20,000x20,000 
mapHeight = entireMapEdges_local(3);
mapWidth = entireMapEdges_local(4);
mapResolution = 0.7;
omap = binaryOccupancyMap(mapHeight+100,mapWidth+100,mapResolution);
% load in building obstacles
setOccupancy(omap,xyBuildings_local,1);
disp("Generating Occupancy Map Complete!")

% confirm customer and vendor locations are not on occupied spaces
xyCustomers_local = CheckLocationOccupancy(omap,xyCustomers_local_unfiltered);
xyVendors_local = CheckLocationOccupancy(omap,xyVendors_local_unfiltered);
customer_cng_per = sum(sum(xyCustomers_local-xyCustomers_local_unfiltered,2) ~=0)/length(xyCustomers_local);
xyCustomer_change = max(sqrt(sum((xyCustomers_local - xyCustomers_local_unfiltered).^2,2)));
vendor_cng_per = sum(sum(xyVendors_local-xyVendors_local_unfiltered,2) ~=0)/length(xyVendors_local);
xyVendor_change = max(sqrt(sum((xyVendors_local - xyVendors_local_unfiltered).^2,2)));



% View Occupancy Map
figure()

hold on
% show occupancy map
show(omap,"local")
% show pseudo occupancy map
scatter(xyBuildings_local(:,1),xyBuildings_local(:,2),1,'k')

% show city contour
h = plot(S_contour_convex_local.contour(:,1),S_contour_convex_local.contour(:,2),'k--','LineWidth',1);
for idx = 1:1:size(struct2table(S_contours),1)
    plot(S_contours_local(idx).contour(:,1),S_contours_local(idx).contour(:,2),'k','LineWidth',1)
end
hold off

xlabel('X [m]','Interpreter','Latex','FontSize',15)
ylabel('Y [m]','Interpreter','Latex','FontSize',15)
ax = ancestor(h, 'axes');
ax.XAxis.Exponent = 4;
ax.YAxis.Exponent = 4;
title('')
sgtitle(city + " - " +altitude+ "ft",'interpreter','latex','FontSize',15)
axis equal

customer_vendor_title = city + " Customers and Vendors";



% View Vendors and Customers
figure()

hold on
% show city contour
h = plot(S_contour_convex_local.contour(:,1),S_contour_convex_local.contour(:,2),'k--','LineWidth',1);
for idx = 1:1:size(struct2table(S_contours),1)
% for idx = 1:1:1
    plot(S_contours_local(idx).contour(:,1),S_contours_local(idx).contour(:,2),'k','LineWidth',1)
end

scatter(xyCustomers_local(:,1),xyCustomers_local(:,2),5,'filled', ...
    'DisplayName','customers','MarkerFaceColor',[0 1 1],'MarkerEdgeColor',[0 1 1])
scatter(xyVendors_local(:,1),xyVendors_local(:,2),50,'^', ...
    'DisplayName','vendors','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0 0 0])
hold off

xlabel('X [m]','Interpreter','Latex','FontSize',15)
ylabel('Y [m]','Interpreter','Latex','FontSize',15)
ax = ancestor(h, 'axes');
ax.XAxis.Exponent = 4;
ax.YAxis.Exponent = 4;
sgtitle(city + " Customers and Vendors",'interpreter','latex','FontSize',15)
axis equal

num_customers = size(xyCustomers_local,1);
num_vendors = size(xyVendors_local,1);
disp("Number of Vendor/Start Locations: " + num2str(size(xyVendors_local,1)))
disp("Number of Customers/Stop Locations: " + num2str(size(xyCustomers_local,1)))


%% Experiment

generate_od = true;
generate_rrt = false;


%%%%%%%%%%%%%%%%%%%%%%%% --- SENSOR PARAMETERS --- %%%%%%%%%%%%%%%%%%%%%%%
% D_rad_BLA = 250;            % [m] bluetooth legacy advertising
% D_rad_BLR = 1000;           % [m] bluetooth long range
% D_rad_WFN = 2000;           % [m] Wi-Fi NAN
% D_rad_WFB = 2000;           % [m] Wi-Fi Beacon

D_rad = 1000;          % sensor radius for simulation %%%
% select number of receivers
num_receivers = 30;                                  %%%
% select number of trials
num_receiver_distribs = 1;                             
num_paths = 1;


MAX_RANGE_RULE = false;
% max_range = 10000;        % mean([4,10,14.5,20])*1e3 = 12125
max_range = 6500;           % average distance found for SF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure start and stop for visual
sf_start = [13270.4 9920 0];
sf_goal = [9973.11 408 0];
la_start = [19400 38600 0];
la_goal = [51827 34511 0];
nyc_start = [30960.1 30050.6 0];
nyc_goal = [39302.1 18650.6 0];


% od records
avg_cvg_record_od = zeros([num_paths,num_receiver_distribs]);
dist_record_od = zeros([num_paths,num_receiver_distribs]);

% rrt records
avg_cvg_record_rrt = zeros([num_paths,num_receiver_distribs]);
dist_record_rrt = zeros([num_paths,num_receiver_distribs]);
iter_time_record_rrt = zeros([num_paths,num_receiver_distribs]);
path_efficiency_record_rrt = zeros([num_paths,num_receiver_distribs]);

% other constants
num_vendors = size(xyVendors_local,1);
num_customers = size(xyCustomers_local,1);

tStart = tic;
for distrib_idx = 1:1:num_receiver_distribs   % vary distributions

    % generate new distribution
    [xyBuildingSensorLocations] = ...
            GenerateSensorLocations(S_building_sens_local,num_receivers);

    ttl_cvg_for_distrib_od = 0;
    ttl_cvg_for_distrib_rrt = 0;
    for path_idx = 1:1:num_paths    % vary start and goal points

        disp("Path "+num2str(path_idx + num_paths*(distrib_idx-1)) ...
                    + " out of " +num2str(num_receiver_distribs*num_paths))
        
        % select start and goal points
        start_idx_rand = round(rand(1)*num_vendors);
        % edge case when start idx is 0
        if start_idx_rand == 0
            start_idx_rand = 1;
        end
        start = xyVendors_local(start_idx_rand,:);
        search_for_goal = true;
        
        if MAX_RANGE_RULE
            while search_for_goal
                goal_idx_rand = round(rand(1)*num_customers);
                % edge case when goal idx is 0
                if goal_idx_rand == 0
                    goal_idx_rand = 1;
                end
                goal = xyCustomers_local(goal_idx_rand,:);
                min_dist = sqrt(sum((start-goal).^2,2));
                if min_dist <= max_range
                    search_for_goal = false;
                end
            end
        else
            goal_idx_rand = round(rand(1)*num_customers);
            % edge case when goal idx is 0
            if goal_idx_rand == 0
                goal_idx_rand = 1;
            end
            goal = xyCustomers_local(goal_idx_rand,:);
        end
        
        start = [start,0];
        goal = [goal, 0];
        
        % --- Origin Destination ---%
        if generate_od
            % generate path
            num_trajecotry_pnts = 1000;
            path_coordinates_od = OriginDesitnation(start,goal,num_trajecotry_pnts);
            path_coordinates_detected_od = DetectPoints(path_coordinates_od,xyBuildingSensorLocations,D_rad);
            distance_od = sqrt((start(1)-goal(1))^2+(start(2)-goal(2))^2);
            % analyze path
            cvg_percent_od = size(path_coordinates_detected_od,1)/size(path_coordinates_od,1);
            ttl_cvg_for_distrib_od = ttl_cvg_for_distrib_od + cvg_percent_od;
            avg_cvg_for_dist_od = ttl_cvg_for_distrib_od/path_idx;
            dist_record_od(path_idx,distrib_idx) = distance_od;
            avg_cvg_record_od(path_idx,distrib_idx) = avg_cvg_for_dist_od;
        end
        
        % --- RRT* --- %
        if generate_rrt
            smoothPath = false;
            tStart_rrt = tic;
            [path_coordinates,explore_coordinates] = RRTStar(start,goal,omap, smoothPath);
            iter_time_record_rrt(path_idx,distrib_idx) = toc(tStart_rrt);
            path_coordinates_detected = DetectPoints(path_coordinates,xyBuildingSensorLocations,D_rad);
            distance_rrt = ComputeTotalDistance(path_coordinates);
            % analyze path
            cvg_percent = size(path_coordinates_detected,1)/size(path_coordinates,1);
            ttl_cvg_for_distrib_rrt = ttl_cvg_for_distrib_rrt + cvg_percent;
            avg_cvg_for_dist_rrt = ttl_cvg_for_distrib_rrt/path_idx;
            dist_record_rrt(path_idx,distrib_idx) = distance_rrt;
            avg_cvg_record_rrt(path_idx,distrib_idx) = avg_cvg_for_dist_rrt;
            % compare path
            distance_od = sqrt((start(1)-goal(1))^2+(start(2)-goal(2))^2);
            path_efficiency_rrt = distance_rrt/distance_od;
            path_efficiency_record_rrt(path_idx,distrib_idx) = path_efficiency_rrt;
        end
        
    end
end
tEnd = toc(tStart);
disp("Run Complete")

%--- ANALYZE EXPERIMENT ---%
figure()
hold on
convergence_line_od = zeros([num_paths,1]);
convergence_record = zeros([num_paths-1,1]);
for distrib_idx = 1:1:num_receiver_distribs
    if generate_od
        plot([1:1:num_paths],avg_cvg_record_od(:,distrib_idx)*100,'DisplayName','OD')
    end
    convergence_line_od = convergence_line_od + avg_cvg_record_od(:,distrib_idx)*100;
    if generate_rrt
        plot([1:1:num_paths],avg_cvg_record_rrt(:,distrib_idx)*100,'DisplayName','RRT*')
    end
end
% % convergence_line_od = convergence_line_od/num_receiver_distribs)
% for idx = 1:1:length(convergence_record)
%     convergence_record(idx) = abs(convergence_line_od(idx+1)-convergence_line_od(idx));
% end
% plot([1:1:num_paths-1],convergence_record)
hold off
xlabel("Number of Paths For "+num2str(num_receiver_distribs) ...
                +" Sensor Distribution","Interpreter","latex",'FontSize',15)
ylabel('Average Coverage \%',"Interpreter","latex",'FontSize',15)
sgtitle('Convergence of Average Trajectory Coverage Over 10 Distributions', ...
    "Interpreter","latex",'FontSize',15)

disp(' ')
disp(' ')
disp("%%%%%%% Experiment Parameters %%%%%%%")
disp("Range of Reciever: " + num2str(D_rad) + " meters")
disp("Number of Receiver Distributions: " +  num_receiver_distribs)
disp("Number of Paths Per Distribution: " + num_paths)
disp("Number of Recievers: " +  num_receivers)
disp(' ')
disp("%%%%%%% OD Coverage %%%%%%%")
disp("Average Coverage %: "+ mean(avg_cvg_record_od(end,:))*100)
% disp("Std of Average Coverages %: " + num2str(std(avg_cvg_record_od(end,:))*100))
% disp("Max Diff %: " + num2str(((max(avg_cvg_record_od(end,:))*100) - (mean(avg_cvg_record_od(end,:))*100))) )
% disp("Min Diff %: " + num2str(((min(avg_cvg_record_od(end,:))*100) - (mean(avg_cvg_record_od(end,:))*100))) )
disp("Avg Distance [m]: " + num2str(mean(dist_record_od,'all')))
disp("Std of Avg Distance [m]: " + num2str(std(dist_record_od,0,'all')))
disp("Max Distance [m]: " + num2str(max(dist_record_od,[],'all')))
disp("Min Distance [m]: " + num2str(min(dist_record_od,[],'all')))

disp(' ')
if generate_rrt
    disp("%%%%%%% RRT* Coverage %%%%%%%")
    disp("Average Coverage %: "+ num2str(mean(avg_cvg_record_rrt(end,:))))
%     disp("Std of Average Coverages %: " + num2str(std(avg_cvg_record_rrt(end,:))*100))
%     disp("Max Diff %: " + num2str((max(avg_cvg_record_rrt(end,:))*100 - mean(avg_cvg_record_rrt(end,:))*100)) )
%     disp("Max Diff %: " + num2str((min(avg_cvg_record_rrt(end,:))*100 - mean(avg_cvg_record_rrt(end,:))*100)) )
    disp("Avg Distance [m]: " + num2str(mean(dist_record_rrt,'all')))
    disp("Std of Avg Distance [m]: " + num2str(std(dist_record_rrt,0,'all')))
    disp("Avg Efficiency: " + num2str(mean(path_efficiency_record_rrt,'all')))
    disp("Std Efficiency: " + num2str(std(path_efficiency_record_rrt,0,'all')))
    disp("Max Distance [m]: " + num2str(max(dist_record_rrt,[],'all')))
    disp("Min Distance [m]: " + num2str(min(dist_record_rrt,[],'all')))

end
disp("Run Time (seconds): " + tEnd)
disp(' ')
disp(' ')

%% VISUALIZE PATH

if generate_od
    figure()
    fig_2_title = city+" - OD";
    
    hold on
    %--- show occupancy map ---%
    % omap.show
    
    %--- show pseudo occupancy map ---%
    scatter(xyBuildings_local(:,1),xyBuildings_local(:,2),2,'k')
    
    %--- show vendors and customers ---%
    % scatter(xyVendors_local(:,1),xyVendors_local(:,2),50,'^', ...
    %     'DisplayName','vendors','MarkerFaceColor',[1 0 0])
    % scatter(xyCustomers_local(:,1),xyCustomers_local(:,2),5,'filled', ...
    %     'DisplayName','customers','MarkerFaceColor',[0 1 1])
    
    %---% show city contour ---%
    plot(S_contour_convex_local.contour(:,1),S_contour_convex_local.contour(:,2),'k--','LineWidth',1)
    for idx = 1:1:size(struct2table(S_contours),1)
    % for idx = 1:1:1
        plot(S_contours_local(idx).contour(:,1),S_contours_local(idx).contour(:,2),'k','LineWidth',1)
    end
    
    %--- show receiver circles ---%
    viscircles(xyBuildingSensorLocations, ...
                D_rad*ones([size(xyBuildingSensorLocations,1),1]), ...
                'linewidth',1,'LineStyle','-','Color','red');
    
    %--- plot path ---%
    h = scatter(path_coordinates_od(:,1),path_coordinates_od(:,2),10,'b','filled','DisplayName','path undetected');
    if ~isempty(path_coordinates_detected_od)
        scatter(path_coordinates_detected_od(:,1),path_coordinates_detected_od(:,2),10,'r','filled','DisplayName','path detected')
    end

    
    %--- show origin and destination ---%
    scatter(start(1),start(2),75,'d','LineWidth',2,'MarkerEdgeColor',[0 0 0], ...
        'MarkerFaceColor',[0 1 0],'DisplayName','origin')
    scatter(goal(1),goal(2),75,'d','LineWidth',2,'MarkerEdgeColor',[0 0 0], ...
        'MarkerFaceColor',[1 1 0],'DisplayName','destination')
    
    hold off
    title('')
    sgtitle(fig_2_title,'interpreter','latex','FontSize',15)
    % legend({'','','','','','','','path undetected','path detected','orign','destination'},'location','northwest')
    ax = ancestor(h, 'axes');
    ax.XAxis.Exponent = 4;
    ax.YAxis.Exponent = 4;
    axis equal
end

if generate_rrt
    figure()
    fig_3_title = city+" - RRT*";
    
    hold on
    %--- show occupancy map ---%
    % omap.show
    
    %--- show pseudo occupancy map ---%
    scatter(xyBuildings_local(:,1),xyBuildings_local(:,2),2,'k')
    
    %--- show vendors and customers ---%
    % scatter(xyVendors_local(:,1),xyVendors_local(:,2),50,'^', ...
    %     'DisplayName','vendors','MarkerFaceColor',[1 0 0])
    % scatter(xyCustomers_local(:,1),xyCustomers_local(:,2),5,'filled', ...
    %     'DisplayName','customers','MarkerFaceColor',[0 1 1])
    
    %---% show city contour ---%
    plot(S_contour_convex_local.contour(:,1),S_contour_convex_local.contour(:,2),'k--','LineWidth',1)
    for idx = 1:1:size(struct2table(S_contours),1)
    % for idx = 1:1:1
        plot(S_contours_local(idx).contour(:,1),S_contours_local(idx).contour(:,2),'k','LineWidth',1)
    end
    
    %--- show receiver circles ---%
    viscircles(xyBuildingSensorLocations, ...
                D_rad*ones([size(xyBuildingSensorLocations,1),1]), ...
                'linewidth',1,'LineStyle','-','Color','red');
    
    %--- plote path ---%
    h = scatter(path_coordinates(:,1),path_coordinates(:,2),10,'b','filled','DisplayName','RRT* Undetected')
    if ~isempty(path_coordinates_detected)
        scatter(path_coordinates_detected(:,1),path_coordinates_detected(:,2),10,'r','filled','DisplayName','RRT* Detected')
    end
%     plot(explore_coordinates(:,1),explore_coordinates(:,2),'DisplayName','RRT* Explored')
    
    %--- show origin and destination ---%
    scatter(start(1),start(2),75,'d','LineWidth',2,'MarkerEdgeColor',[0 0 0], ...
        'MarkerFaceColor',[0 1 0],'DisplayName','origin')
    scatter(goal(1),goal(2),75,'d','LineWidth',2,'MarkerEdgeColor',[0 0 0], ...
        'MarkerFaceColor',[1 1 0],'DisplayName','destination')
    
    hold off
    title('')
    sgtitle(fig_3_title,'interpreter','latex','FontSize',15)
    % legend({'','','','','','','','path undetected','path detected','orign','destination'},'location','northwest')
    ax = ancestor(h, 'axes');
    ax.XAxis.Exponent = 4;
    ax.YAxis.Exponent = 4;
    axis equal
end



%% FUNCTIONS

function [xySensorLocations] = GenerateSensorLocations(S_sens,num_receivers)
%     numberOfAllSensors = size(struct2table(S_sens),1);
%     probability = rand([numberOfAllSensors,1]);
%     sensorMask = probability <= pecentageOfSensors;
%         
%     % use binary array to select points
%     xySensorLocations_idx = [1:1:size(struct2table(S_sens),1)];
%     xSensorLocations = [S_sens(xySensorLocations_idx((sensorMask))).LocationX];
%     ySensorLocations = [S_sens(xySensorLocations_idx((sensorMask))).LocationY];
%     xySensorLocations = horzcat(xSensorLocations',ySensorLocations');

    num_all_possible = size(S_sens.XLocation,1);
    xy_all_possible = zeros([num_all_possible,2]);
    for idx = 1:1:num_all_possible
        xy_all_possible(idx,:) = [S_sens.XLocation(idx),S_sens.YLocation(idx)];
    end
    
    xySensorLocations = zeros([num_receivers,2]);
    for idx = 1:1:num_receivers
        
        select_idx = round(rand(1)*size(xy_all_possible,1));
        if select_idx == 0
            select_idx = 1;
        end
        xySensorLocations(idx,:) = [xy_all_possible(select_idx,1),xy_all_possible(select_idx,2)];
        xy_all_possible(select_idx,:) = [];

    end
end


function path_coordinates = OriginDesitnation(start,goal,num_trajectroy_pnts)
    
    x_coords = linspace(start(1), goal(1), num_trajectroy_pnts);
    y_coords = linspace(start(2), goal(2), num_trajectroy_pnts);
    path_coordinates = [x_coords(:), y_coords(:)];
    
end


function path_coordinates_detected = DetectPoints(path_coordinates,xyBuildingSensorLocations,D_rad)
    
    num_trajecotry_pnts = size(path_coordinates,1);

    % comptue trajectory coverage
    x1 = path_coordinates(:,1);
    y1 = path_coordinates(:,2);
    x2 = xyBuildingSensorLocations(:,1);
    y2 = xyBuildingSensorLocations(:,2);
    
    Z = zeros(num_trajecotry_pnts,1);
    for traj_pnt = 1:1:num_trajecotry_pnts
        Z(traj_pnt) = min((x1(traj_pnt)-x2).^2 + (y1(traj_pnt)-y2).^2);
    end
%     % vecotrized distance equation (that runs slower?)
%     Z = min(x1.^2 - 2*x1*x2' + x2'.^2 + y1.^2 - 2*y1*y2' + y2'.^2,[],2);
    
    path_coordinates_detected = [x1(Z<=D_rad.^2),y1(Z<=D_rad.^2)];
end

function total_distance = ComputeTotalDistance(path_coordinates)
    total_distance = 0;
    for idx = 1:1:size(path_coordinates,1)-1
        line_distance = sqrt((path_coordinates(idx,1)-path_coordinates(idx+1,1))^2 ...
        + (path_coordinates(idx,2)-path_coordinates(idx+1,2))^2);
        total_distance = total_distance + line_distance;
    end
end

function xy_out = CheckLocationOccupancy(omap,xy_locations)
    occ_check_cust = checkOccupancy(omap,xy_locations);
    % seperate positions on unoppupied and occupied space
    xy_unoccupied = xy_locations(occ_check_cust==0,:);
    xy_out = zeros([size(xy_locations)]);
    num_unoccupied_initially = size(xy_unoccupied,1);
    xy_out(1:num_unoccupied_initially,:) = xy_unoccupied;
    xy_occupied = xy_locations(occ_check_cust==1,:);
    % randomly move positions on occupied spaces until they are no longer on 
    % occupied spaces
    if size(xy_occupied,1) > 0
        for idx = 1:1:size(xy_occupied,1)
            occupied = true;
            pos = xy_occupied(idx,:);
            while occupied
                % pick random direction and allow max movement of 2 steps
                r = 2;
                theta_rand = rand(1)*2*pi - pi;
                xy_rand = [r*cos(theta_rand),r*sin(theta_rand)];
                pos = pos + xy_rand;
                occupied = checkOccupancy(omap,pos);
            end
            xy_out(num_unoccupied_initially+idx,:) = pos;
        end
    end
end






