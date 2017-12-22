function [ projectedPoints ] = projectionOnSphere( pointsOfModel, center, radius)
% PROJECTIONONSPHERE
% The function is used to project model on the sphere with concrete center and radius.
%
% INPUTS:
%   pointsOfModel - matrix N x 3, where N i number of points of model;
%   represents cartesian coordinates
%   center - ideal center of the sphere
%   radius - ideal radius of the sphere
%
% OUTPUTS:
%   projectedPoints - points projected on sphere with defined radius and center.
%
% EXAMPLE
%
% BRIEF EXPLANATION:
%
% SEE ALSO:
%
% Author: Vaclava Piorecka (vaclava.piorecka@fbmi.cvut.cz, vaclava.piorecka@nudz.cz)
% Date:   
% 2017-10-01    creation of function

    numOfPoints = size(pointsOfModel,1);
    
    % Shift all points in a coordinate system centered at the center of the sphere 
    X = pointsOfModel(:,1) - ones(numOfPoints,1)*center(1,1);   
    Y = pointsOfModel(:,2) - ones(numOfPoints,1)*center(1,2);   
    Z = pointsOfModel(:,3) - ones(numOfPoints,1)*center(1,3);   
    
    % Length of shift vector 
    norm = sqrt(X.^2 + Y.^2 + Z.^2);
    
    % Scale the vector to the radius of the sphere 
    scaleData(:,1) = (radius.*X(:,1))./norm;
    scaleData(:,2) = (radius.*Y(:,1))./norm;
    scaleData(:,3) = (radius.*Z(:,1))./norm;

    % Transform to the original coordinate system
    projectedPoints(:,1) = scaleData(:,1) - ones(numOfPoints,1)*center(1,1);   
    projectedPoints(:,2) = scaleData(:,2) - ones(numOfPoints,1)*center(1,2);   
    projectedPoints(:,3) = scaleData(:,3) - ones(numOfPoints,1)*center(1,3);  

end

