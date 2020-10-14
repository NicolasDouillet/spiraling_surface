function [V, T] = spiraling_suface(spiralid, param, grid, option_display)
%% spiraling_suface : function to compute, display, and save a meshed spiraling surface.
%
% Author & support : nicolas.douillet (at) free.fr, 2008-2020.
%
%
% Syntax
%
% spiraling_suface(spiralid);
% spiraling_suface(spiralid, param);
% spiraling_suface(spiralid, param, grid);
% spiraling_suface(spiralid, param, grid, option_display);
% [V, T] = spiraling_suface(spiralid, param, grid, option_display);
%
%
% Description
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
%
% Input arguments
%
% - spiralid : integer scalar double ine the range |[1; 9]|, the spiral identification number.
%
% 1 : logarithmic spiral
% 2 : Archimede spiral
% 3 : Exponential / Fermat spiral
% 4 : power function spiral
% 5 : gamma function spiral
% 6 : Inverse Archimede spiral
% 7 : Inverse logarithmic spiral
% 8 : Inverse exponential spiral
% 9 : Inverse gamma spiral
%
% - grid = [xmin, xmax, ymin, ymax, resolution_X, resolution_Y], real vector double, the spiral grid support parameters. numel(grid) = 6.
%
% - param = [amplitude, way, curv_coeff, shift, pow_coeff, nb_arms, phi_origin], real vector double, the spiral parameters. numel(param) = 7.
%
% - option_display : either logical *true/false or numeric *1/0.
%
%
% Output arguments
%
%       [ |  |  |]
% - V = [Vx Vy Vz], real matrix double, the vertex coordinates. Size(V) = [nb_vertices,3].
%       [ |  |  |]
%
%       [ |  |  |]
% - T = [T1 T2 T3], positive integer matrix double, the triangulation. Size(T) = [nb_triangles,3].
%       [ |  |  |]
%
%
% Example #1 : computes and displays the log spiraling
% surface with default parameters on the default grid.
% spiraling_suface;
%
% Example #2 : computes and display the Archimede spiraling
% surface with customized parameters, on the default grid.
% param = [1, 1, 0.5, 0, 0, 1, 0.5]; % = [amplitude, way, curv_coeff, shift, power_coeff, nb_arms, phi_origin]
% spiraling_suface(2,param);
%
% Example #3 : computes and display the square root spiraling
% surface with customized spiral and grid parameters.
% param = [1, 1, 16, 0, 0.5, 2 0];
% grid = [-45 45 -80 80 0.5 0.5]; % = [xmin, xmax, ymin, ymax, resolution_X, resolution_Y]
% spiraling_suface(4,param,grid);


%% Input parsing
assert(nargin < 5,'Too many input arguments.');    
    
if nargin < 4
    
    option_display = true;
    
    if nargin < 3 % default grid
        
        % Real numbers
        xmin = -30;
        xmax =  30; % xmax > xmin
        ymin = -30;
        ymax =  30; % ymax > ymin
        
        % Positive real numbers
        resolution_Y = 1; % < 0.5*(ymax-ymin)
        resolution_X = 1; % < 0.5*(xmax-xmin)
        
        if nargin < 2 % default param
            
            if ~nargin % nargin < 1
                
                spiralid = 1;
                
            else
                
                assert(isnumeric(spiralid) && spiralid == floor(spiralid) && spiralid > 0 && spiralid < 10,'spiralid must be an integer in the range |[1; 9]|.');
                
            end
            
            default_param = [1,    1,    8,     2,    0,   2, 0;... % log
                             1,   -1,    1,     0,    0,   2, 0;... % archi
                             1,    1,    0.04,  2,    0,   2, 0;... % exp
                             1,    1,    16,    0,    0.5, 2  0;... % pow
                             1,    1     0.1,   0.3,  0,   2, 0;... % gamma                                                           
                             1,    1,    0.002, 0.02, 0    2, 0;... % inv_archi                                                                                      
                             1,    90,   0.1,   3     0,   2, 0;... % inv_log
                             1,    60,   0.08,  0.72, 0,   2, 0;... % inv_exp                             
                             1,    30    0.05,  0.3,  0,   2, 0];   % inv_gamma
                                         
            param = default_param(spiralid,:);
            
        else % customized param
            
            assert(isreal(param) && numel(param) == 7,'Input argument param must be a real six elements vector.');
            
        end
        
    else % if nargin > 2 -> customized grid
        
        assert(isreal(grid) && numel(grid) == 6,'Input argument grid must be a vector containing six elements : [xmin, xmax, ymin, ymax, resolution_X, resolution_Y].');
        assert(grid(2) > grid(1) && grid(4) > grid(3),'Input argument grid vector must be such that xmax > xmin and ymax > ymin.');
        
        xmin = grid(1);
        xmax = grid(2);
        ymin = grid(3);
        ymax = grid(4);
        resolution_Y = grid(5);
        resolution_X = grid(6);
                
    end
    
    [y,x] = meshgrid(ymin:resolution_X:ymax,xmin:resolution_Y:xmax);
    
else
    
    assert(islogical(option_display) || isnumeric(option_display),'option_display parameter type must be either logical or numeric.');
    
end


%% Body

z = compute_spiral(spiralid,param,y,x);

Sx = size(x,1);
Sy = size(y,2);

V = cat(2,x(:),y(:),z(:));
T = build_triangulation(Sx,Sy);

if option_display
    
    clr_map = 'winter';
    dispiral(V,T,clr_map);
    
end


end % spiraling_suface


%% compute_spiral subfunction
function [z] = compute_spiral(spiral_type, param, y, x)


amplitude  = param(1); % sinus amplitude
way        = param(2); % spiral rotation way
curv_coeff = param(3); % curvature coefficient
shift      = param(4); % radius shift
pow_coeff  = param(5); % power coefficient
nb_arms    = param(6); % number of arms
phi_origin = param(7); % phase at the origin

r = sqrt(x.^2+y.^2);
Phi = angle(x+y*1i);

switch spiral_type
    
    case 1
        
        z = amplitude*sin(way*curv_coeff*log(r+shift)+nb_arms*Phi+phi_origin*pi); % 'log'
        
    case 2
        
        z = amplitude*sin(way*(curv_coeff*r+shift)+nb_arms*Phi+phi_origin*pi); % 'archi'
        
    case 3
        
        z = amplitude*sin(way*exp(curv_coeff*r+shift)+nb_arms*Phi+phi_origin*pi); % 'Fermat' / 'exp'                      
        
    case 4
        
        z = amplitude*sin(way*(curv_coeff*r+shift).^pow_coeff+nb_arms*Phi+phi_origin*pi); % 'pow'
        
    case 5
        
        z = amplitude*sin(way*gamma(curv_coeff*r+shift)+nb_arms*Phi+phi_origin*pi); % 'gamma'        
        
    case 6
        
        z = amplitude*sin(way./(curv_coeff*r+shift)+nb_arms*Phi+phi_origin*pi); % 'inv_archi'
        
    case 7
        
        z = amplitude*sin(way./log(curv_coeff*r+shift)+nb_arms*Phi+phi_origin*pi); % 'inv_log'
        
    case 8
        
        z = amplitude*sin(way./exp(curv_coeff*r+shift)+nb_arms*Phi+phi_origin*pi); % 'inv_exp'
        
    case 9
        
        z = amplitude*sin(way./gamma(curv_coeff*r+shift)+nb_arms*Phi+phi_origin*pi); % 'inv_gamma'                                                                
                
end


end % compute_spiral


%% build_triangulation subfunction
function [T] = build_triangulation(Sx, Sy)


r1 = cat(2,1,repelem(2:Sx-1,2),Sx);
r1 = reshape(r1,[2,Sx-1])';
R1 = cat(2,r1,(1:Sx-1)'+Sx); % 1st triangle row indices
% size(R1,1) = Sx-1

r2 = cat(2,1+Sx,repelem(2+Sx:2*Sx-1,2),2*Sx);
r2 = reshape(r2,[2,Sx-1])';
R2 = cat(2,(2:Sx)',fliplr(r2)); % 2nd triangle row indices

T = repmat(cat(1,R1,R2),[Sy-1,1]);
% size(T) = 2*(Sx-1)*(Sy-1)

T = T + Sx*repelem((0:Sy-2)',2*(Sx-1),3);


end % build_triangulation


%% dispiral subfunction
function [] = dispiral(V, T, clr_map)


f0 = [0 0 0];
f1 = [1 1 1];

h = figure(1);
set(h,'Position',get(0,'ScreenSize'));

trimesh(T,V(:,1),V(:,2),V(:,3));

ax = gca;
set(ax,'Color',f0,'XColor',f1,'YColor',f1,'ZColor',f1);
ax.Clipping = 'off';
set(gcf,'Color',f0);
view(-45,85);

shading interp;
colormap(clr_map);
axis square, axis tight;


end % dispiral