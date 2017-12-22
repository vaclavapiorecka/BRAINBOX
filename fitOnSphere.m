function [ center, radius ] = fitOnSphere( pointsOfModel )
% FITONSPHERE 
% The function is used to calculate ideal projection of the model on the
% sphere. The output is ideal centre and ideal radius of the sphere.
%
% INPUTS:
%   pointsOfModel - matrix N x 3, where N i number of points of model;
%   represents cartesian coordinates
%
% OUTPUTS:
%   center - ideal center of the sphere
%   radius - ideal radius of the sphere
%
% EXAMPLE
%
% BRIEF EXPLANATION:
%   Minimizes function of the sphere: r^2 = (x-x0)^2+(y-y0)^2+(z-z0)^2,
%   where r is the radius, [x0,y0,z0] is the center of the sphere and
%   [x,y,z] are cartesian coordinates of model points.
%
% SEE ALSO:
%
% Author: Vaclava Piorecka (vaclava.piorecka@fbmi.cvut.cz, vaclava.piorecka@nudz.cz)
% Date:   
% 2017-10-01    creation of function


initialC = [mean(pointsOfModel(:,1)), mean(pointsOfModel(:,2)),mean(pointsOfModel(:,3))];
initialR = 1;

funSphere = @(outputRaC)sum((pointsOfModel(:,1)-(outputRaC(:,1))).^2+(pointsOfModel(:,2)-(outputRaC(:,2))).^2+(pointsOfModel(:,3)-(outputRaC(:,3))).^2-(initialR)^2)^2;
options = optimset('Display','none','MaxIter',1000);
[outputC] = fminsearch(funSphere,initialC,options);

funSphere = @(outputR)sum((pointsOfModel(:,1)-(outputC(:,1))).^2+(pointsOfModel(:,2)-(outputC(:,2))).^2+(pointsOfModel(:,3)-(outputC(:,3))).^2-(outputR)^2)^2;
options = optimset('Display','none','MaxIter',1000);
[outputR] = fminsearch(funSphere,initialR,options);

center = outputC;
radius = outputR;
end

