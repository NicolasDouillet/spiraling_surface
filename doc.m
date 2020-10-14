%% spiraling_suface
%
% Function to compute, display, and save a meshed spiraling surface.
%
% Author & support : nicolas.douillet (at) free.fr, 2008-2020.
%
%% Syntax
% spiraling_suface(spiralid);
%
% spiraling_suface(spiralid, param);
%
% spiraling_suface(spiralid, param, grid);
%
% spiraling_suface(spiralid, param, grid, option_display);
%
% [V, T] = spiraling_suface(spiralid, param, grid, option_display);
%
%% Description
%
% spiraling_suface(spiralid) computes and displays the logarithmic
% spiraling surface with default parameter set on the default grid.
%
% spiraling_suface(spiralid, param) uses param for the spiral parameter vector.   
%
% spiraling_suface(spiralid, param, grid) uses grid for the spiral support grid vector.
%
% spiraling_suface(spiralid, param, grid, option_display) displays it when
% option_display is set to logical *true/1 (default), and doesn't
% when it is set to  logical false/0.
%
% [V, T] = spiraling_suface(spiralid, param, grid, option_display) stores the resulting
% vertices coordinates in the array V, and the corresponding triplet indices list in the array T. 
%
%% See also
%
% <https://fr.mathworks.com/matlabcentral/fileexchange/75109-galaxy-model galaxy_model>
%
%% Input arguments
%
% - *spiralid* : integer scalar double ine the range |[1; 9]|, the spiral identification number.
%
% # : logarithmic spiral
% # : Archimede spiral
% # : Exponential / Fermat spiral
% # : power function spiral
% # : gamma function spiral
% # : Inverse Archimede spiral
% # : Inverse logarithmic spiral
% # : Inverse exponential spiral
% # : Inverse gamma spiral
%
% - *grid* = [xmin, xmax, ymin, ymax, resolution_X, resolution_Y], real vector double, the spiral grid support parameters. numel(grid) = 6.
%
% - *param* = [amplitude, way, curv_coeff, shift, pow_coeff, nb_arms, phi_origin], real vector double, the spiral parameters. numel(param) = 7.
%
% - *option_display* : either logical *true/false or numeric *1/0.
%
%% Output arguments
%
%        [ |  |  |]
% - V = [Vx Vy Vz], real matrix double, the vertex coordinates. Size(V) = [nb_vertices,3].
%        [ |  |  |]
%
%        [ |  |  |]
% - T = [T1 T2 T3], positive integer matrix double, the triangulation. Size(T) = [nb_triangles,3].
%        [ |  |  |]
%
%% Example #1
% Computes and displays the log spiraling
% surface with default parameters on the default grid.
spiraling_suface;

%% Example #2
% Computes and display the Archimede spiraling
% surface with customized parameters, on the default grid.
param = [1, 1, 0.5, 0, 0, 1, 0.5]; % = [amplitude, way, curv_coeff, shift, power_coeff, nb_arms, phi_origin]
spiraling_suface(2,param);

%% Example #3
% Computes and display the square root spiraling
% surface with customized spiral and grid parameters.
param = [1, 1, 16, 0, 0.5, 2 0];
grid = [-45 45 -80 80 0.5 0.5]; % = [xmin, xmax, ymin, ymax, resolution_X, resolution_Y]
spiraling_suface(4,param,grid);