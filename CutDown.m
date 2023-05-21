%% cut down information of large building shapfe files
clear
clc

city_section = "1";   % select city section 
city = "SF"
disp("Loading Original Dataset ...")
shapeFilePath = "OSM_datasets/"+city+"_buildings_"+ ...
    city_section + "/building_polygon.shp";
S_raw = shaperead(shapeFilePath);
disp("Loading Original Dataset Complete!")

disp("Copying Desired Files ...")
numTableRows_raw = size(struct2table(S_raw),1);

for in_idx = 1:1:numTableRows_raw
    S_out(in_idx).X = S_raw(in_idx).X;
    S_out(in_idx).Y = S_raw(in_idx).Y;
    S_out(in_idx).height = S_raw(in_idx).height;
    S_out(in_idx).BoundingBox = S_raw(in_idx).BoundingBox;
end
disp("Copying Desired Files Complete!")

disp("Saving CutDown Struct ...")
outputFilePath = "OSM_datasets/"+city+"_buildings_"+ ...
    city_section + "/CutDown/building_polygon.mat"

save(outputFilePath, 'S_out')
disp("Saving CutDown Struct Complete!")


%% filter hieghts
clear
clc

city = "SF";
city_section = 1;   % select city section 
altitude_ft = -1;

load_buildings = load("OSM_datasets/"+city+"_buildings_"+ ...
    num2str(city_section) + "/CutDown/building_polygon.mat");
S_raw = load_buildings.S_out;

% choose -1 for all footprints
altitude_m = altitude_ft/3.281;

numTableRows = size(struct2table(S_raw),1);
% filter heights
temp_idx = 1;
S_out.X = [];
S_out.Y = [];
S_out.height = '';
S_out.BoundingBox = [];
for raw_idx = 1:1:numTableRows
    height = str2double(S_raw(raw_idx).height);
    if altitude_m<0
        S_out(temp_idx).X = S_raw(raw_idx).X;
        S_out(temp_idx).Y = S_raw(raw_idx).Y;
        S_out(temp_idx).BoundingBox = S_raw(raw_idx).BoundingBox;
        S_out(temp_idx).height = S_raw(raw_idx).height;
        temp_idx = temp_idx + 1;
    
    elseif ~isnan(height) && height >= altitude_m
        S_out(temp_idx).X = S_raw(raw_idx).X;
        S_out(temp_idx).Y = S_raw(raw_idx).Y;
        S_out(temp_idx).BoundingBox = S_raw(raw_idx).BoundingBox;
        S_out(temp_idx).height = S_raw(raw_idx).height;
        temp_idx = temp_idx + 1;
    end
end

if altitude_ft < 0
    outputFilePath = "OSM_datasets/"+city+"_buildings_"+ ...
    num2str(city_section) + "/CutDown/" + ...
                    "building_polygon_HNoneft.mat"
else
    outputFilePath = "OSM_datasets/"+city+"_buildings_"+ ...
    num2str(city_section) + "/CutDown/" + ...
                "building_polygon_H"+num2str(altitude_ft)+"ft.mat"
end
save(outputFilePath, 'S_out')
%% check the magic
clear
clc

load("OSM_datasets/NYC_buildings_2/CutDown/building_polygon_H400ft.mat")

%% save all city borders/contours to mat
clear
clc
close all

load_building = shaperead("OSM_datasets/LA_border/border_level8_polygon.shp");
S_raw = load_building;

xyLocations = zeros([size(S_raw.X,2),2]);
xyLocations(:,1) = S_raw.X;
xyLocations(:,2) = S_raw.Y;

S_in.contour = xyLocations;

% find all contours
temp_idx = 1;
contour_counter_idx = 1;
S_out(1).contour = [];
for in_idx = 1:1:size(S_in.contour,1)
    if ~isnan(S_in.contour(in_idx))
        S_out(temp_idx,:).contour(contour_counter_idx,:) = S_in.contour(in_idx,:);
        contour_counter_idx = contour_counter_idx + 1;
    else
        temp_idx = temp_idx + 1;
        contour_counter_idx = 1;
    end
end

% select the correct contour to save
figure(1)
hold on
for idx = 1:1:size(struct2table(S_out),1)
    plot(S_out(idx).contour(:,1),S_out(idx).contour(:,2))
end
hold off
axis equal



outputFilePath ="FormattedDatasets/LA_formatted/LA_border.mat"
save(outputFilePath, 'S_out')

%% convex city border
clear
clc
close all

load_building = shaperead("OSM_datasets/LA_border/border_level8_polygon.shp");
S_raw = load_building;


P = rmmissing([S_raw.X',S_raw.Y']);
k = convhull(P);

figure(1)
hold on
plot(P(k,1),P(k,2))

hold off
axis equal

S_out.contour = [P(k,1),P(k,2)];

outputFilePath ="FormattedDatasets/LA_formatted/LA_border_convex.mat";
save(outputFilePath, 'S_out');
