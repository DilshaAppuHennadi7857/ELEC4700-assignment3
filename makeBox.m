function makeBox(left,right,top,bottom,areaL,areaW)
%MAKEBOX Create and plot a box given sides and region
%   areaL and areaW define the total area of the region (workspace). Given
%   values for left, right, top, and bottom will be treated as fractional
%   values of the workspace - i.e. if right = 0.5, the right side of the
%   box will be located at x = 0.5*areaL. Thus, all values for left, right,
%   top and bottom must be between [0,1]. Locations of the sides are used
%   to determine the locations of the four verticies:
%
%       1--------2
%       |        |
%       |        |
%       |        |
%       3--------4
% 
%	Function does not turn hold on or off, that must be done in script.
%	This is just to avoid any confusion in plotting.

% ver1 = [left*areaL, top*areaW];
% ver2 = [right*areaL, top*areaW];
% ver3 = [left*areaL, bottom*areaW];
% ver4 = [right*areaL, bottom*areaW];

plot([left*areaL right*areaL], [top*areaW top*areaW], 'b') % 1 1nd 2
plot([right*areaL right*areaL], [top*areaW bottom*areaW], 'b') % 2 and 4
plot([left*areaL left*areaL], [top*areaW bottom*areaW], 'b') % 1 and 3
plot([left*areaL right*areaL], [bottom*areaW bottom*areaW], 'b') % 3 and 4

end

