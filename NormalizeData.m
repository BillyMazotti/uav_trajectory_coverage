
% NORMALIZE DATA FUNCTIOJN

function [entireMapEdges_local,xyBuildings,xyCustomers, ...
    xyVendors,S_building_sens,S_contour_convex,S_contours] ...
    = NormalizeData(S_building_sens,xyCustomers,xyVendors, ...
                            xyBuildings,S_contour_convex,S_contours)

    % find edges of customer, vendor, and building sensor maps 
    % (in global coords)
    %building
    [minX_map_V,minY_map_V,maxX_map_V,maxY_map_V] = FindMapEdges(S_building_sens);
    buildingSensMapEdges = [minX_map_V,minY_map_V,maxX_map_V,maxY_map_V];
    
    
    % find edges of customer, vendor, and building location/occupancy maps 
    % (in global coords)
    %customer
    minXY_map_C = [min(xyCustomers(:,1)),min(xyCustomers(:,2))];
    maxXY_map_C = [max(xyCustomers(:,1)),max(xyCustomers(:,2))];
    customerLocsMapEdges = [minXY_map_C,maxXY_map_C];
    
    %vendor
    minXY_map_V = [min(xyVendors(:,1)),min(xyVendors(:,2))];
    maxXY_map_V = [max(xyVendors(:,1)),max(xyVendors(:,2))];
    vendorLocsMapEdges = [minXY_map_V,maxXY_map_V];
    
    %building
    minXY_map_B = [min(xyBuildings(:,1)),min(xyBuildings(:,2))];
    maxXY_map_B = [max(xyBuildings(:,1)),max(xyBuildings(:,2))];
    buildingLocsMapEdges = [minXY_map_B,maxXY_map_B];
    
    % find edges of city contour 
    % (in global coords)
    minXY_map_contour = [min(S_contour_convex.contour(:,1)),min(S_contour_convex.contour(:,2))];
    maxXY_map_contour = [max(S_contour_convex.contour(:,1)),max(S_contour_convex.contour(:,2))];
    contourMapEdges = [minXY_map_contour,maxXY_map_contour];
    
    
    %total/encompassing map edges
    edges = vertcat(buildingSensMapEdges,customerLocsMapEdges, ...
                    vendorLocsMapEdges,buildingLocsMapEdges,contourMapEdges);
    entireMapEdges_global = [min(edges(:,1:2),[],1),max(edges(:,3:4),[],1)];
    
    
    % normalize building sensor struct
    for idx = 1:size(struct2table(S_building_sens),1)        
        S_building_sens.XLocation(idx) = ...
                S_building_sens.XLocation(idx) - entireMapEdges_global(1);
        S_building_sens.YLocation(idx) = ...
                S_building_sens.YLocation(idx) - entireMapEdges_global(2);
    end
    
    % normalize building occupancy and vendor/customer location arrays
    xyBuildings = xyBuildings - [entireMapEdges_global(1),entireMapEdges_global(2)];
    xyCustomers = xyCustomers - [entireMapEdges_global(1),entireMapEdges_global(2)];
    xyVendors = xyVendors - [entireMapEdges_global(1),entireMapEdges_global(2)];
    
    % normalize city contour
    S_contour_convex.contour = S_contour_convex.contour - [entireMapEdges_global(1),entireMapEdges_global(2)];
    for idx = 1:1:size(S_contours,2)
        S_contours(idx).contour = S_contours(idx).contour - [entireMapEdges_global(1),entireMapEdges_global(2)];
    end
    % find eges of normalized Customer, Vendor, Building sensor map
    [minX_map_B,minY_map_B,maxX_map_B,maxY_map_B] = FindMapEdges(S_building_sens);
    buildingSensMapEdges = [minX_map_B,minY_map_B,maxX_map_B,maxY_map_B];
    %occupancy/locations
    customerLocsMapEdges = [min(xyCustomers(:,1)),min(xyCustomers(:,2)), ...
                                max(xyCustomers(:,1)),max(xyCustomers(:,2))];
    vendorLocsMapEdges = [min(xyVendors(:,1)),min(xyVendors(:,2)), ...
                                max(xyVendors(:,1)),max(xyVendors(:,2))];
    buildingLocsMapEdges = [min(xyBuildings(:,1)),min(xyBuildings(:,2)), ...
                                max(xyBuildings(:,1)),max(xyBuildings(:,2))];


    minXY_map_contour = max(S_contour_convex.contour,[],1);
    maxXY_map_contour = min(S_contour_convex.contour,[],1);
    contourMapEdges = [minXY_map_contour,maxXY_map_contour];
    
    edges = vertcat(buildingSensMapEdges,customerLocsMapEdges, ...
                    vendorLocsMapEdges,buildingLocsMapEdges,contourMapEdges);
    entireMapEdges_local = [min(edges(:,1:2),[],1),max(edges(:,3:4),[],1)];

end