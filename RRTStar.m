% RRT* FUNCTION
function [path_coordinates,explore_coordinates] = RRTStar(start,goal,omap,smoothPath)

    disp("Generating Path with RRT* ...")
    tStart_path = tic;
    
    ss = stateSpaceSE2;
    sv = validatorOccupancyMap(ss);
    sv.Map = omap;
    
    sv.ValidationDistance = 10;
    ss.StateBounds = [omap.XWorldLimits; omap.YWorldLimits; [-pi pi]];
    planner = plannerRRTStar(ss,sv);
    
    % Maximum length of a motion allowed in the tree, specified as a scalar.
    planner.MaxConnectionDistance = 5;
    planner.BallRadiusConstant = 100;
    planner.ContinueAfterGoalReached = false;
    planner.MaxNumTreeNodes = 1e10;
    planner.MaxIterations = 1e5;
    planner.GoalBias = 0.05;
    
    [refPath,solnInfo] = plan(planner,start,goal);
    explore_coordinates = [solnInfo.TreeData(:,1),solnInfo.TreeData(:,2)];

    if ~isempty(refPath.States)
        disp("Generating Path with RRT* - SUCCESS")
    else
        disp("Generating Path with RRT* - FAIL")
    end
    
    toc(tStart_path)
    % interpolate path
    if ~smoothPath
        path_coordinates = horzcat(refPath.States(:,1),refPath.States(:,2));
    else
        tStart2 = tic;
        disp('Smoothing Path ...')
        
        % attempt 1
        waypoints = refPath.States;
        nWayPoints = refPath.NumStates;
        
        % Calculate the distance between waypoints
        distance = zeros(1,nWayPoints);
        for i = 2:nWayPoints
            distance(i) = norm(waypoints(i,1:3) - waypoints(i-1,1:3));
        end
        
        % Assume a UAV speed of 3 m/s and calculate time taken to reach 
        % each waypoint
        UAVspeed = 10;
        timepoints = cumsum(distance/UAVspeed);
        nSamples = nWayPoints*2;
        
        
        % Compute states along the trajectory
        path_coordinates = minsnappolytraj(waypoints',timepoints,nSamples, ...
                                        MinSegmentTime=2,MaxSegmentTime=20)';

        disp('Smoothing Path Complete!')
        toc(tStart2)        
    end

end