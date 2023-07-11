%% Building Receiver/Sensor Locations (HNoneft files)
% Building Data Conversion from Shape Read to Mat File
% asumes geometry of entires are polygons
clear
clc
close all

% --- SF CONSTANTS & CORRECTIONS --- % (37.77N, -122.41W)
% minX_map_lat = -130;
% minY_map_long = 0;

% --- LA CONSTANTS & CORRECTIONS --- % (34.05N, -118.24W)
% minX_map_lat = -130;
% minY_map_long = 0;

% --- NYC CONSTANTS & CORRECTIONS --- % (40.71N, -74.0060W)
minX_map_lat = -100;
minY_map_long = 15;

city_section = 2;
city = "NYC";

filePath = "OSM_datasets/"+city+"/"+city+"_buildings_"+ ...
    num2str(city_section) + "/CutDown/building_polygon_HNoneft.mat";

convexContourPath = "FormattedDatasets/"+city+"_formatted/"+city+"_border_convex.mat"


oneContour = true;  % necessary condition to consider all convex options
lat2meters = 10000000/90;
long2meters = 40075161.2/360;
minX_map_m = minX_map_lat * lat2meters;
minY_map_m = minY_map_long * long2meters;
sensor_percentage = 1;


load_building = load(convexContourPath);
S_contour = load_building.S_out;
[S_out,num] = PolygonXYLocation(filePath,minX_map_m,minY_map_m,lat2meters,long2meters,S_contour,oneContour);

outputFilePath = "FormattedDatasets/"+city+"_formatted/"+getDateLabel()+"_"+city+"_buidings_"+ ...
    num2str(city_section) + "_sensors_N"+num2str(num)+".mat"

save(outputFilePath, 'S_out','-v7.3')

disp("Receivers Saved!")


%% Building Occupancy Map Data
clear
clc
close all


% --- SF CONSTANTS & CORRECTIONS --- % (37.77N, -122.41W)
% minX_map_lat = -130;
% minY_map_long = 0;

% --- LA CONSTANTS & CORRECTIONS --- % (34.05N, -118.24W)
% minX_map_lat = -130;
% minY_map_long = 0;

% --- NYC CONSTANTS & CORRECTIONS --- % (40.71N, -74.0060W)
minX_map_lat = -100;
minY_map_long = 15;

city = "NYC";
city_section = 2;
altitude_ft = -1;

if altitude_ft < 0
    shapeFilePath = "OSM_datasets/"+city+"_buildings_"+ ...
        num2str(city_section) + "/CutDown/building_polygon_HNoneft.mat"
else
    shapeFilePath = "OSM_datasets/"+city+"_buildings_"+ ...
        num2str(city_section) + "/CutDown/building_polygon_H"+num2str(altitude_ft)+"ft.mat"
end
convexContourPath = "FormattedDatasets/"+city+"_formatted/"+city+"_border_convex.mat"


numTableRows_raw = 0;
lat2meters = 10000000/90;
long2meters = 40075161.2/360;
minX_map_m = minX_map_lat * lat2meters;
minY_map_m = minY_map_long * long2meters;
map_minXY_meters = [minX_map_m,minY_map_m;minX_map_m,minY_map_m];
altitude_m = altitude_ft/3.281;

load_building = load(convexContourPath);
S_contour_convex = load_building.S_out;
load(shapeFilePath)
[S_out,num_buildings] = OccupancyMapMesh(shapeFilePath,numTableRows_raw,minX_map_m, ...
                            minY_map_m,map_minXY_meters,lat2meters,long2meters,altitude_m, ...
                            S_contour_convex)

if altitude_ft < 0
    outputFilePath = "FormattedDatasets/"+city+"_formatted/" + ...
        "2_3_23_"+city+"_buidings_"+ num2str(city_section) + ...
        "_occupancy_HNoneft_N"+ ...
        num2str(num_buildings)+".mat"
else    
    outputFilePath = "FormattedDatasets/"+city+"_formatted/"+getDateLabel()+"_"+city+"_buidings_"+ ...
        num2str(city_section) + "_occupancy_H" + ...
        num2str(altitude_ft)+"ft_N"+num2str(num_buildings)+".mat"
end

save(outputFilePath, 'S_out','-v7.3')

disp("Occupancy Map Data Saved!")

%% Building Polygon Map Data
clear
clc
close all


% --- SF CONSTANTS & CORRECTIONS --- % (37.77N, -122.41W)
% shapeFilePath = "OSM_datasets/SF_buildings/building_polygon.shp";
% minX_map_lat = -130;
% minY_map_long = 0;
% outputFilePath = "FormattedDatasets/BuildingPolygon_SF_" + ...
%                 "BuildingMygeodata_H"+num2str(altitude_ft)+"ft_1_31_23_.mat"

% --- LA CONSTANTS & CORRECTIONS --- % (34.05N, -118.24W)
minX_map_lat = -130;
minY_map_long = 0;

city = "NYC"
num_case = "4"
altitude_ft = 400;



decimate = false;
makeMeters = false;
numTableRows_raw = 0;

filePath = "FormattedDatasets/"+city+"_formatted/GeoFence/2_3_23_"+city+"_geofence_"+num_case+"_H400ft.mat";
convexContourPath = "FormattedDatasets/"+city+"_formatted/"+city+"_border_convex.mat";
load_building = load(convexContourPath);
S_contour_convex = load_building.S_out;
lat2meters = 10000000/90;
long2meters = 40075161.2/360;
minX_map_m = minX_map_lat * lat2meters;
minY_map_m = minY_map_long * long2meters;
altitude_m = altitude_ft/3.281;

[S_out,num_buildings] = PolygonMapCoords(filePath,numTableRows_raw,minX_map_m, ...
                            minY_map_m,lat2meters,long2meters,altitude_m,decimate,makeMeters,S_contour_convex);

outputFilePath = "FormattedDatasets/"+city+"_formatted/PolygonMap/"+getDateLabel()+"_"+city+"_buidings_"+num_case+"_polygon_H" + ...
                num2str(altitude_ft)+"ft_N"+num2str(num_buildings)+".mat"

% save file
save(outputFilePath, 'S_out','-v7.3')

disp("Polygon Map Data Saved!")

%% ANALYSIS
clear
clc

shapeFilePath = "OSM_datasets/SF_buildings/building_polygon.shp";
S_raw = shaperead(shapeFilePath);

close all
% number of rows in struct
numTableRows_raw = size(struct2table(S_raw),1)

histogram_heights = zeros([numTableRows_raw,1]);

for idx = 1:1:numTableRows_raw
    height = str2double(S_raw(idx).height);
    if ~isnan(height)
        histogram_heights(idx) = height;
    end
end

% histogram_heights = rmmissing(histogram_heights);
histogram_heights = histogram_heights(histogram_heights>0);
num_nonzero_heights = size(histogram_heights,1)
perentage_heights = num_nonzero_heights/numTableRows_raw*100;

figure(99)
sgtitle({"Heights of Buildingsin SF","Buildings With Heights: " + ...
    num2str(perentage_heights) + "% | AVG: " + num2str(mean(histogram_heights)) ... 
    + "m | STD: " + num2str(std(histogram_heights))+ "m"},'FontSize',12)
xline(121.92,'r','linewidth',2)
hold on
histogram(histogram_heights,500)
hold off
set(gca,'XScale','log')
set(gca,'YScale','log')
legend({'400ft Amazon Cruise Alt'},'location','north')
xlabel("Building Height [m]")
ylabel("Number of Buildings per Bin")

%% Vendor/Start Locations (NOT CONVEX)
clear
clc

% --- SF CONSTANTS & CORRECTIONS --- % (37.77N, -122.41W)
% minX_map_lat = -130;
% minY_map_long = 0;

% --- LA CONSTANTS & CORRECTIONS --- % (34.05N, -118.24W)
% minX_map_lat = -130;
% minY_map_long = 0;

% --- NYC CONSTANTS & CORRECTIONS --- % (40.71N, -74.0060W)
minX_map_lat = -100;
minY_map_long = 15;

city = "NYC";


shapeFilePath_polygon = "OSM_datasets/"+city+"_shopgeneral/shop_general_polygon.shp";
shapeFilePath_point = "OSM_datasets/"+city+"_shopgeneral/shop_general_point.shp";
load_struct = load("FormattedDatasets/"+city+"_formatted/"+city+"_border.mat");
S_contour = load_struct.S_out;

oneContour = false;
lat2meters = 10000000/90;
long2meters = 40075161.2/360;
minX_map_m = minX_map_lat * lat2meters;
minY_map_m = minY_map_long * long2meters;

S_out_point = PointXYLocation(shapeFilePath_point,minX_map_m,minY_map_m,lat2meters,long2meters,S_contour);
[S_out_polygon,num] = PolygonXYLocation(shapeFilePath_polygon,minX_map_m,minY_map_m,lat2meters,long2meters,S_contour,oneContour);
% 
S_out.XLocation = [S_out_point.XLocation;S_out_polygon.XLocation];
S_out.YLocation = [S_out_point.YLocation;S_out_polygon.YLocation];

%--- Save Data ---%
num_vendors = size(S_out.XLocation,1);
outputFilePath = "FormattedDatasets/"+city+"_formatted/"+getDateLabel()+"_vendors_N"+num2str(num_vendors)+".mat"
save(outputFilePath, 'S_out')

disp("Vendor Data Saved!")


%% Customer/Stop Locations
clear
clc

% --- SF CONSTANTS & CORRECTIONS --- % (37.77N, -122.41W)
% minX_map_lat = -130;
% minY_map_long = 0;
% pop_density_sqr_miles = 18629.1;
% city = "SF";

% % --- LA CONSTANTS & CORRECTIONS --- % (34.05N, -118.24W)
% minX_map_lat = -130;
% minY_map_long = 0;
% pop_density_sqr_miles = 8304.2;
% city = "LA";

% --- NYC CONSTANTS & CORRECTIONS --- % (40.71N, -74.0060W)
minX_map_lat = -100;
minY_map_long = 15;
pop_density_sqr_miles = 29303.2;
city = "NYC";

shapeFilePath = "OSM_datasets/"+city+"/"+city+"_residential/residential_polygon.shp";
load_struct = load("FormattedDatasets/"+city+"_formatted/"+city+"_border.mat");
S_contour = load_struct.S_out;
pop_density_sqr_meters = pop_density_sqr_miles/(5280^2)*3.28084^2;
numTableRows_raw = 0;
lat2meters = 10000000/90;
long2meters = 40075161.2/360;
minX_map_m = minX_map_lat * lat2meters;
minY_map_m = minY_map_long * long2meters;
map_minXY_meters = [minX_map_m,minY_map_m;minX_map_m,minY_map_m];


[S_out,num_customers] = CustomerStopLocations(shapeFilePath, ...
                                numTableRows_raw, ...
                                minX_map_m,minY_map_m,map_minXY_meters,lat2meters,long2meters, ...
                                pop_density_sqr_meters,S_contour)


outputFilePath = "FormattedDatasets/"+city+"_formatted/"+getDateLabel()+"_"+city+"_customers_N"+num2str(num_customers)+".mat"
save(outputFilePath, 'S_out')

disp("Customer Data Saved!")



%% FUNCTIONS
clear
clc
close all

%--- Test for GenerateLatticePoints ---%
% r=10; % radius
% C=[1 1];
% theta=0:2*pi/9:2*pi; % the angle
% m=r*[cos(theta')+C(1) sin(theta')+C(2)]; % the points you asked
% xContourPoints = m(:,1)'
% yContourPoints = m(:,2)'
% BBoxPoints = [min(xContourPoints),min(yContourPoints);
%         max(xContourPoints),max(yContourPoints)]
% [xLattice,yLattice] = GenerateLatticePoints(BBoxPoints,xContourPoints,yContourPoints);


% Polygon Map Data
function [S_out,num_buildings] = PolygonMapCoords(filePath,numTableRows_raw,minX_map_m, ...
                                    minY_map_m,lat2meters,long2meters,altitude,decimate,makeMeters,S_contour_convex)
   

    disp("Loading OSM File ...")
    [fp,nm,ext] = fileparts(filePath);
    if ext == ".mat"
        load_buildings = load(filePath);
        S_raw = load_buildings.S_out;
        disp("Loading OSM Shape File Complete!")
    elseif ext == ".shp"
        S_raw = shaperead(filePath);
        disp("Loading OSM Mat File Complete!")
    end
        
    % number of rows in struct
    if numTableRows_raw == 0
        numTableRows_raw = size(struct2table(S_raw),1);
    end
     
    
    % convert coordiantes to meters and normalize map dat to global 
    % reference point map_min_meters
    disp("Normalizing Map Data Coordinates ...")
    S_temp = struct;
    if makeMeters
        temp_idx = 1;
        for raw_idx = 1:1:numTableRows_raw
            height = str2double(S_raw(raw_idx).height);
            if ~isnan(height) && height > altitude
                S_temp(temp_idx).X = ...
                    round(rmmissing(S_raw(raw_idx).X) * ...
                    lat2meters - minX_map_m);
                S_temp(temp_idx).Y = ...
                    round(rmmissing(S_raw(raw_idx).Y) * ...
                    long2meters - minY_map_m);
                temp_idx = temp_idx + 1;
            end
        end
    else
        temp_idx = 1;
        for raw_idx = 1:1:numTableRows_raw
            height = str2double(S_raw(raw_idx).height);
            if ~isnan(height) && height > altitude
                disp(raw_idx)
                S_temp(temp_idx).X = rmmissing(S_raw(raw_idx).X);
                S_temp(temp_idx).Y = rmmissing(S_raw(raw_idx).Y);
                temp_idx = temp_idx + 1;
            else
                disp(raw_idx)
            end
    
        end
    end
    
    S_contour_convex.contour(:,1) = round(S_contour_convex.contour(:,1)*lat2meters - minX_map_m);
    S_contour_convex.contour(:,2) = round(S_contour_convex.contour(:,2)*long2meters - minY_map_m);
    disp("Normalizing Map Complete ...")
    numTableRows_temp = size(struct2table(S_temp),1);
    
    
    % confirm building locations are within convex contour of city
    S_temp2 = struct;
    temp2_idx = 1;
    for temp_idx = 1:1:numTableRows_temp
        xLocations = S_temp(temp_idx).X;
        yLocations = S_temp(temp_idx).Y;
        in_city_mask = inpolygon(xLocations,yLocations, ...
                            S_contour_convex.contour(:,1), ...
                            S_contour_convex.contour(:,2));
        if length(xLocations) - length(xLocations(in_city_mask==1)) == 0
            S_temp2(temp2_idx).X = xLocations;
            S_temp2(temp2_idx).Y = yLocations;
            temp2_idx = temp2_idx + 1;
        end
    end
    
    S_temp = S_temp2;
    numTableRows_temp = size(struct2table(S_temp),1);
    num_buildings = numTableRows_temp;
    disp("Number of Buildings Above " + num2str(altitude*3.281) ...
        + "ft and in Contour: " + num_buildings)

    disp("Collecting Polygon Coordinates ...")
    % count number of polygon coordinates
    num_poly_coords = 0;
    for idx = 1:numTableRows_temp
        num_poly_coords = num_poly_coords + size(S_temp(idx).X,2);   
    end
    if decimate
        vertices_RPS = [];
    else
        vertices_RPS = zeros([num_poly_coords,3]);
    end
    
    % save and index polygon points for visibility graph
    row_idx = 1;
    count = 0;
    for polygon_idx = 1:numTableRows_temp
        if decimate
            if size(S_temp(polygon_idx).X,2) > 10 % decimate large polygons
                polygon = DecimatePoly([[S_temp(polygon_idx).X';S_temp(polygon_idx).X(1)], ...
                                [S_temp(polygon_idx).Y';S_temp(polygon_idx).Y(1)]],[0.5 2], 0);
                count = count + 1;
            elseif size(S_temp(polygon_idx).X,2) > 20 % decimate large polygons
                polygon = DecimatePoly([[S_temp(polygon_idx).X';S_temp(polygon_idx).X(1)], ...
                                    [S_temp(polygon_idx).Y';S_temp(polygon_idx).Y(1)]],[0.25 2], 0);
                count = count + 1;
            else
                polygon = [[S_temp(polygon_idx).X';S_temp(polygon_idx).X(1)], ...
                                    [S_temp(polygon_idx).Y';S_temp(polygon_idx).Y(1)]];
            end
        else
            polygon = [[S_temp(polygon_idx).X';S_temp(polygon_idx).X(1)], ...
                                    [S_temp(polygon_idx).Y';S_temp(polygon_idx).Y(1)]];
        end
        num_pnts = size(polygon,1);
        vertices_RPS(row_idx:row_idx+num_pnts-1,:) = ...
                                   [polygon, polygon_idx*ones(num_pnts,1)];
        row_idx = row_idx + num_pnts;
    end
    disp("Collecting Polygon Coordinates Complete!")
    
    
    S_out(1).vertices_RPS = vertices_RPS;
    disp("Number of Decimations: " + num2str(count))

end


function [S_out,num_customers] = CustomerStopLocations(shapeFilePath,numTableRows_raw, ...
                                        minX_map_m,minY_map_m,map_minXY_meters,lat2meters,long2meters, ...
                                        density,S_contour)
    
    S_raw = shaperead(shapeFilePath);
    
    if numTableRows_raw == 0
        numTableRows_raw = size(struct2table(S_raw),1);
    end
    
    % convert coordiantes to meters and normalize map dat to global 
    % reference point map_min_meters
    disp("Normalizing Map Data Coordinates ...")
    temp_idx = 1;
    S_temp.X = [];
    S_temp.Y = [];
    S_temp.BoundingBox = [];
    for raw_idx = 1:1:numTableRows_raw
        S_temp(temp_idx).BoundingBox = ...
            round(S_raw(raw_idx).BoundingBox.* ...
            [lat2meters,long2meters]) - map_minXY_meters;
        S_temp(temp_idx).X = ...
            round(rmmissing(S_raw(raw_idx).X) * ...
            lat2meters - minX_map_m);
        S_temp(temp_idx).Y = ...
            round(rmmissing(S_raw(raw_idx).Y) * ...
            long2meters - minY_map_m);
        temp_idx = temp_idx + 1;

    end
    % convert contour coordinates to meters and normalize map to minXY
    for contour_idx = 1:1:size(struct2table(S_contour),1)
        S_contour(contour_idx).contour(:,1) = ...
            round(S_contour(contour_idx).contour(:,1)*lat2meters - minX_map_m);
        S_contour(contour_idx).contour(:,2) = ...
            round(S_contour(contour_idx).contour(:,2)*long2meters - minY_map_m);
    end
    disp("Normalizing Map Data Coordinates Complete!")
    

    % compute mesh coords
    disp("Computing Mesh Coordinates ...")
    for idx = 1:numTableRows_raw
        disp(num2str(idx) + " of " + num2str(numTableRows_raw))
        [S_temp(idx).XMesh,S_temp(idx).YMesh] = ...
            GenerateLatticePoints(S_temp(idx).BoundingBox, ...
                                    S_temp(idx).X, ...
                                    S_temp(idx).Y);
    end
    disp("Computing Mesh Coordinates Complete!")
    
    numXYMesh_coords = 0;
    for idx = 1:numTableRows_raw
        numXYMesh_coords = numXYMesh_coords + size(S_temp(idx).XMesh,1);
        
    end
    allLatticeX = zeros([numXYMesh_coords,1]);
    allLatticeY = zeros([numXYMesh_coords,1]);
    
    disp("Collecting Mesh Coordinates ...")
    idx_pos = 1;
    for idx = 1:numTableRows_raw
        space = size(S_temp(idx).XMesh,1);
        allLatticeX(idx_pos:idx_pos+space-1) = S_temp(idx).XMesh;
        allLatticeY(idx_pos:idx_pos+space-1) = S_temp(idx).YMesh;
        idx_pos = idx_pos + size(S_temp(idx).XMesh,1);
    end
    disp("Collecting Mesh Coordinates Complete!")

    % select drone stops
    disp("Choosing Customer Stop Locations ...")
    [xLocations,yLocations] = ...
        CustomerVendorDistribution(density,allLatticeX,allLatticeY);
    

    % confirm customer locations are within city boundaries
    in_city_mask = zeros([size(xLocations),1]);
    for contour_idx = 1:1:size(struct2table(S_contour),1)
        in_contour_mask = inpolygon(xLocations,yLocations, ...
                            S_contour(contour_idx).contour(:,1), ...
                            S_contour(contour_idx).contour(:,2));
        in_city_mask = in_city_mask + in_contour_mask;
    end
    xLocations = xLocations(in_city_mask>0);
    yLocations = yLocations(in_city_mask>0);
    disp("Choosing Customer Stop Locations Complete!")
    
    num_customers = size(xLocations,1);
    disp("Numer of Customers/Vendors Chosen: " + num2str(num_customers))
    S_out.XLocation = xLocations;
    S_out.YLocation = yLocations;
end


% Occupancy Map Data
function [S_out,num_buildings] = OccupancyMapMesh(filePath,numTableRows_raw,minX_map_m, ...
                                    minY_map_m,map_minXY_meters,lat2meters,long2meters, ...
                                    altitude,S_contour_convex)
    %{
        Args:
            shapeFilePath = "SF_residential/residential_polygon.shp";
            numTableRows:   first X number of regions from raw database to 
                            extract data from; set to 0 to extract from 
                            all rows          
            minX_map_lat:   minimum lattitue you'd expect (in degrees) to
                            ensure that all lattitude values are postive
            minY_map_long:  minimum longitude you'd expect (in degrees) to
                            ensure that all lattitude values are postive
            outputFilePath: file destination to save output struct "S_out"
            LocationType:   customer, vendor, or building
        Returns:
            S_out: for idx ranging 1 to numTables*percentage
                - S_out(idx).XMesh
                - S_out(idx).YMesh
    %}

    % save shape file to struct
    disp("Loading OSM File ...")
    [fp,nm,ext] = fileparts(filePath);
    if ext == ".mat"
        load_buildings = load(filePath);
        S_raw = load_buildings.S_out;
        disp("Loading OSM Shape File Complete!")
    elseif ext == ".shp"
        S_raw = shaperead(filePath);
        disp("Loading OSM Mat File Complete!")
    end
        
    % number of rows in struct
    if numTableRows_raw == 0
        numTableRows_raw = size(struct2table(S_raw),1);
    end
     
    
    % convert coordiantes to meters and normalize map dat to global 
    % reference point map_min_meters
    disp("Normalizing Map Data Coordinates ...")
    temp_idx = 1;
    S_temp.X = [];
    S_temp.Y = [];
    S_temp.BoundingBox = [];
    for raw_idx = 1:1:numTableRows_raw
        if altitude < 0
            S_temp(temp_idx).BoundingBox = ...
                round(S_raw(raw_idx).BoundingBox.* ...
                [lat2meters,long2meters]) - map_minXY_meters;
            S_temp(temp_idx).X = ...
                round(rmmissing(S_raw(raw_idx).X) * ...
                lat2meters - minX_map_m);
            S_temp(temp_idx).Y = ...
                round(rmmissing(S_raw(raw_idx).Y) * ...
                long2meters - minY_map_m);
            temp_idx = temp_idx + 1;
        else
            height = str2double(S_raw(raw_idx).height);
            if ~isnan(height) && height >= altitude
                S_temp(temp_idx).BoundingBox = ...
                    round(S_raw(raw_idx).BoundingBox.* ...
                    [lat2meters,long2meters]) - map_minXY_meters;
                S_temp(temp_idx).X = ...
                    round(rmmissing(S_raw(raw_idx).X) * ...
                    lat2meters - minX_map_m);
                S_temp(temp_idx).Y = ...
                    round(rmmissing(S_raw(raw_idx).Y) * ...
                    long2meters - minY_map_m);
                temp_idx = temp_idx + 1;
            else
                disp("!!! DELETED HEIGHT FOR SOME REASEON ???")
                disp(height)
                disp(S_raw(raw_idx).height)
            end
        end
    end
    S_contour_convex.contour(:,1) = round(S_contour_convex.contour(:,1)*lat2meters - minX_map_m);
    S_contour_convex.contour(:,2) = round(S_contour_convex.contour(:,2)*long2meters - minY_map_m);
    disp("Normalizing Map Data Coordinates Complete!")
    numTableRows_temp = size(struct2table(S_temp),1);
    num_buildings = numTableRows_temp;
    disp("Number of Buildings Above " + num2str(altitude*3.281) ...
        + "ft: " + num_buildings)
   
    % compute mesh coords
    disp("Computing Mesh Coordinates ...")
    for idx = 1:numTableRows_temp
        [S_temp(idx).XMesh,S_temp(idx).YMesh] = ...
            GenerateLatticePoints(S_temp(idx).BoundingBox, ...
                                    S_temp(idx).X, ...
                                    S_temp(idx).Y);
    end

    % count number of x and y coordinates
    disp("Initilziing Mesh Arrays ...")
    numXYMesh_coords = 0;
    for idx = 1:numTableRows_temp
        numXYMesh_coords = numXYMesh_coords + size(S_temp(idx).XMesh,1);   
    end
    allLatticeX = zeros([numXYMesh_coords,1]);
    allLatticeY = zeros([numXYMesh_coords,1]);
    
    disp("Adding Mesh Coordinates to Arrays ...")
    idx_pos = 1;
    for idx = 1:numTableRows_temp
        space = size(S_temp(idx).XMesh,1);
        try
            allLatticeX(idx_pos:idx_pos+space-1) = S_temp(idx).XMesh;
            allLatticeY(idx_pos:idx_pos+space-1) = S_temp(idx).YMesh;
        catch ME
            disp(ME)
            disp("Size of LeftX")
            disp(size(allLatticeX(idx_pos:idx_pos+space-1)))
            disp("Size of RightX")
            disp(size(S_temp(idx).XMesh))
            disp(ME)
            disp("Size of LeftY")
            disp(size(allLatticeY(idx_pos:idx_pos+space-1)))
            disp("Size of RightY")
            disp(size(S_temp(idx).YMesh))
        end
        idx_pos = idx_pos + size(S_temp(idx).XMesh,1);
    end

    allLatticeX = allLatticeX(allLatticeX ~= 0);
    allLatticeY = allLatticeY(allLatticeY ~= 0);

    % ensure points are within the convex contour
    xLocations = allLatticeX;
    yLocations = allLatticeY;

    % confirm building locations are within convex contour of city    
    in_city_mask = inpolygon(xLocations,yLocations, ...
                        S_contour_convex.contour(:,1), ...
                        S_contour_convex.contour(:,2));

    xLocations = xLocations(in_city_mask>0);
    yLocations = yLocations(in_city_mask>0);

    allLatticeX = xLocations;
    allLatticeY = yLocations;

    % store mesh arrays in output struct
    S_out(1).ObstacleLocationX = allLatticeX;
    S_out(1).ObstacleLocationY = allLatticeY;

end



function [xLocations,yLocations] = CustomerVendorDistribution(density,xLattice,yLattice)
    %{
        Randomly select points in customer regions (polygons) to serve as 
        customer sensor and delivery locations

        Args:
            density: points (locations)/polygon_area (area [m^2])
            xLattice: integer array of x coords of the occupancy points i.e.
            the lattice/polygon points
            ,: integer array of y coords of the occupancy points i.e.
            the lattice/polygon points
        Returns:
            xLocations: integer array of x coords of customer locations
            yLocations: integer array of y coords of customer locations

    %}
    
    totalNumberOfPoints = size(xLattice,1);
    
    % create binary array for random index selection
    D_uniform_rand = rand([totalNumberOfPoints,1]);
    D_uniform_hasD = D_uniform_rand <= density;

    % use binary array to select points
    xLocations = xLattice(D_uniform_hasD);
    yLocations = yLattice(D_uniform_hasD);
    
    % view for testing
    close all
    figure(98)
    hold on
    scatter(xLocations,yLocations,40,'b','filled')
    hold off
    xlabel('Longitude [meters]')
    ylabel('Lattitude [meters]')

    axis('equal')
end


function [xLattice,yLattice] = GenerateLatticePoints(bbox,xContourPoints,yContourPoints)
    %{
        generate lattice points for each building to surve as occupancy
        points, accounting for buildings with irregular footprints
        Args:
            bbox: array of building footprint [xmin,ymin;xmax,ymax]
            xPoints: array of postive x coords of building contour points
            yPoints: array of postive y coords of building contour points
        Returns:
            xLattice: integer array of x coords of the occupancy points
            yLattice: integer array of y coords of the occupancy points

    %}
    % scale down large coordinate values (found was necessary for values 
    % on the order of 1e12)
    minXY = [min(xContourPoints),min(yContourPoints)];
    bbox = bbox - minXY;
    xContourPoints  = xContourPoints - minXY(1);
    yContourPoints = yContourPoints  - minXY(2);

    % find all points within and on the edge of the bounding box
    [xAllPoints,yAllPoints] = meshgrid(bbox(1):bbox(2),bbox(3):bbox(4));

    % find all points within and on the edge of irregular building 
    % footprint + scale-up the coordiantes to their original magnitudes
    inPoints = inpolygon(xAllPoints,yAllPoints,xContourPoints,yContourPoints);
    xLattice = xAllPoints(inPoints) + minXY(1);
    yLattice = yAllPoints(inPoints) + minXY(2);

%     % view for testing
%     figure(99)
%     hold on
%     scatter(xAllPoints,yAllPoints,40,'r','filled')
%     scatter(xAllPoints(inPoints),yAllPoints(inPoints),40,'g','filled')
%     plot(xContourPoints,yContourPoints,'b','linewidth',3)
%     scatter(xContourPoints,yContourPoints,40,'b','filled')
%     hold off

end

function [S_out,num] = PointXYLocation(shapeFilePath,minX_map_m,minY_map_m,lat2meters,long2meters,S_contour)
    
    disp("Loading OSM Shape File ...")
    S_raw = shaperead(shapeFilePath);
    disp("Loading OSM Shape File Complete")
    numTableRows_raw = size(struct2table(S_raw),1);

    % convert coordiantes to meters and normalize map dat to to minXY 
    % reference point map_min_meters
    disp("Normalizing Map Data Coordinates ...")
    temp_idx = 1;
    for raw_idx = 1:1:numTableRows_raw        
        S_temp(temp_idx).X = ...
            round(rmmissing(S_raw(raw_idx).X) * ...
            lat2meters - minX_map_m);
        S_temp(temp_idx).Y = ...
            round(rmmissing(S_raw(raw_idx).Y) * ...
            long2meters - minY_map_m);
        temp_idx = temp_idx + 1;
    end
    % convert contour coordinates to meters and normalize map to minXY
    for contour_idx = 1:1:size(struct2table(S_contour),1)
        S_contour(contour_idx).contour(:,1) = ...
            round(S_contour(contour_idx).contour(:,1)*lat2meters - minX_map_m);
        S_contour(contour_idx).contour(:,2) = ...
            round(S_contour(contour_idx).contour(:,2)*long2meters - minY_map_m);
    end
    disp("Normalizing Map Complete!")
    
    xy_point_locations = zeros(size(struct2table(S_temp),1),2);
    for idx = 1:1:size(struct2table(S_temp),1)
        xy_point_locations(idx,1) = S_temp(idx).X;
        xy_point_locations(idx,2) = S_temp(idx).Y;   
    end

    
    xLocations = xy_point_locations(:,1);
    yLocations = xy_point_locations(:,2);
    disp("Number of Points In and Out of Contour(s): " ...
        + num2str(length(xLocations)))
    % confirm locations are within city boundaries
    in_city_mask = zeros([size(xLocations),1]);
    for contour_idx = 1:1:size(struct2table(S_contour),1)

        in_contour_mask = inpolygon(xLocations,yLocations, ...
                            S_contour(contour_idx).contour(:,1), ...
                            S_contour(contour_idx).contour(:,2));
        in_city_mask = in_city_mask + in_contour_mask;
    end
    xLocations = xLocations(in_city_mask>0);
    yLocations = yLocations(in_city_mask>0);
    num = length(xLocations);
    disp("Number of Points In Contour(s): " + num2str(num))

    S_out.XLocation = xLocations;
    S_out.YLocation = yLocations;
end

function [S_out,num] = PolygonXYLocation(filePath,minX_map_m,minY_map_m,lat2meters,long2meters,S_contour,oneContour)
    
    % save shape file to struct
    disp("Loading OSM File ...")
    [fp,nm,ext] = fileparts(filePath);
    if ext == ".mat"
        load_buildings = load(filePath);
        S_raw = load_buildings.S_out;
        disp("Loading OSM Shape File Complete!")
    elseif ext == ".shp"
        S_raw = shaperead(filePath);
        disp("Loading OSM Mat File Complete!")
    end
    numTableRows_raw = size(struct2table(S_raw),1)

    % convert coordiantes to meters and normalize map dat to global 
    % reference point map_min_meters
    disp("Normalizing Map Data Coordinates ...")
    temp_idx = 1;
    for raw_idx = 1:1:numTableRows_raw        
        S_temp(temp_idx).X = ...
            round(rmmissing(S_raw(raw_idx).X) * ...
            lat2meters - minX_map_m);
        S_temp(temp_idx).Y = ...
            round(rmmissing(S_raw(raw_idx).Y) * ...
            long2meters - minY_map_m);
        temp_idx = temp_idx + 1;
    end
    % convert contour coordinates to meters and normalize map to minXY
    if oneContour
        S_contour.contour(:,1) = ...
            round(S_contour.contour(:,1)*lat2meters - minX_map_m);
        S_contour.contour(:,2) = ...
            round(S_contour.contour(:,2)*long2meters - minY_map_m);
    else
        for contour_idx = 1:1:size(struct2table(S_contour),1)
            S_contour(contour_idx).contour(:,1) = ...
                round(S_contour(contour_idx).contour(:,1)*lat2meters - minX_map_m);
            S_contour(contour_idx).contour(:,2) = ...
                round(S_contour(contour_idx).contour(:,2)*long2meters - minY_map_m);
        end
    end
    disp("Normalizing Map Complete!")

    xy_polygon_locations = zeros(size(struct2table(S_temp),1),2);
    for idx = 1:1:size(struct2table(S_temp),1)
        xy_polygon_locations(idx,1) = mean(S_temp(idx).X);
        xy_polygon_locations(idx,2) = mean(S_temp(idx).Y);   
    end

    xLocations = xy_polygon_locations(:,1);
    yLocations = xy_polygon_locations(:,2);
    disp("Number of Points In and Out of Contour(s): " ...
        + num2str(length(xLocations)))
    % confirm locations are within city boundaries
    in_city_mask = zeros([size(xLocations),1]);
    if oneContour
        in_contour_mask = inpolygon(xLocations,yLocations, ...
                                S_contour.contour(:,1), ...
                                S_contour.contour(:,2));
        in_city_mask = in_city_mask + in_contour_mask;
    else
        for contour_idx = 1:1:size(struct2table(S_contour),1)
            in_contour_mask = inpolygon(xLocations,yLocations, ...
                                S_contour(contour_idx).contour(:,1), ...
                                S_contour(contour_idx).contour(:,2));
            in_city_mask = in_city_mask + in_contour_mask;
        end
    end
    xLocations = xLocations(in_city_mask>0);
    yLocations = yLocations(in_city_mask>0);
    num = length(xLocations);

    disp("Number of Valid Points found: " + num2str(num))

    S_out.XLocation = xLocations;
    S_out.YLocation = yLocations;

end

function date_label = getDateLabel()

    datetime_1 = strsplit(string(datetime),"-")
    datetime_2 = strsplit(datetime_1(3)," ")
    datetime_3 = strsplit(datetime_2(2),":");
    date_label = ["Y"+datetime_2(1)+"_M"+datetime_1(2)+"_D"+datetime_1(1)+...
        "_h"+datetime_3(1)+"_m"+datetime_3(2)+"_s"+datetime_3(3)]

end

