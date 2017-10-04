%Lab1 ANM FEM block

clear all; close all;
format long;
ge = unitsquare();
hmax = 0.05;
beta = [10; 10];
% for i = 1:12

[uu, AA, ~, ~] = My2DPoissonSolver(ge, hmax, beta);
%     norm_u(i) = uu'*AA*uu;
%     hmax = 0.125 -((0.125-0.03)/10)*i;
% end
norm_u = uu'*AA*uu;
err_norm = abs(1/45 - norm_u);
[p2, e2, t2]  = initmesh(ge, 'hmax', hmax);
x2 = p2(1,:);
y2 = p2(2,:);
u = sin(pi*x2).*sin(pi*y2);
%figure,
%pdesurf(p2, t2, u');


% help functions
% unit square geometry
function geom = unitsquare()
geom = [2 0 1 0 0 1 0;...
    2 1 1 0 1 1 0;...
    2 1 0 1 1 1 0;...
    2 0 0 1 0 1 0]';
end

% forcing-function
function f = funct(x, y)
f = transpose(2.*(x - x.^2) + 2.*(y - y.^2));
% f = transpose(2*(pi^2)*sin(pi*x).*sin(pi*y));
% f = ones(size(x))';
end

% penalty term for Robin BC
function ka = kappa(x,y)
if (x == 1 || x == 0)
    ka = 10^6;
elseif (y == 1 || y == 0)
    ka = 10^6;
else
    ka = 0;
end
end

% Poisson solver 2D
function [U, A, R, T] = My2DPoissonSolver(geom, hmax, beta)
[p, e, t] = initmesh(geom, 'hmax', hmax);
[A, R, b, r, T] = assemble2D(p, e, t, beta);
U = (A+R)\(b+r);
%U = (A+R+T)\(b+r);
figure, pdesurf(p, t, U);
% figure, pdemesh(p, e, t);
%figure, spy(A);
norm_A = norm(A-A',inf);
eig_A = eig(A);
eig_AR = eig(A+R);
%figure, plot(real(eig_A), imag(eig_A), 'b*', real(eig_AR), imag(eig_AR), 'g^')
if norm_A == 0
    disp(['A is symmetric, with smallest eigenvalue = ' num2str(min(eig_A))...
        '. The smallest eigenvalue of A+R is ' num2str(min(eig_AR))]);
end

end

% boundary term for Robin BC
function g = diric(x,y)

if (x == 1 || y == 1)
    g = 0;
elseif (x == 0 || y == 0 )
    g = 0;
else
    g = 1;
end
%g = 0;
end

function phi = hatfunction(x,y, area, beta)
b_ = [y(2)-y(3); y(3)-y(1); y(1)-y(2)]/2/area;
c_ = [x(3)-x(2); x(1)-x(3); x(2)-x(1)]/2/area;
a_ = [x(2)*y(3) - x(3)*y(2); x(3)*y(1)-x(1)*y(3);...
    x(1)*y(2) - x(2)*y(1)]/2/area;
phi = beta(1)*(a_'*b_ + b_'*b_*x' + c_'*b_*y')...
    + beta(2)*(a_'*c_ + b_'*c_*x' + c_'*c_*y');
end

% vector & matrix assembler
function [A, R, b, r, T] = assemble2D(p, e, t, beta)

% number of nodes
N_nodes = size(p,2);

% allocating space for matrices & vectors
A = sparse(N_nodes, N_nodes);
R = sparse(N_nodes, N_nodes);
T = sparse(N_nodes, N_nodes);
b = zeros(N_nodes, 1);
r = zeros(N_nodes, 1);

% looping over triangles
for K = 1:size(t,2)
    nodesi = t(1:3, K);
    x = p(1,nodesi);
    y = p(2,nodesi);
    
    % Calculating area of triangle K
    area = polyarea(x,y);
    
    % Hat function gradients
    b_= [y(2)-y(3); y(3)-y(1); y(1)-y(2)]/2/area;
    c_= [x(3)-x(2); x(1)-x(3); x(2)-x(1)]/2/area;
    
    % Hat function
    ph = hatfunction(x, y, area, beta);
    
    
    % Triangle element contribution to stiffness matrix A
    Ak = (b_*b_' + c_*c_')*area;
    Tk = ph*area;
    % forcing-function's contribution to load vector b
    bk = funct(mean(x),mean(y)).*area/3;
    
    A(nodesi, nodesi) = A(nodesi, nodesi) + Ak;
    T(nodesi, nodesi) = T(nodesi, nodesi) + Tk;
    b(nodesi) = b(nodesi) + bk;
    
end



for E = 1:size(e,2)
    nodes = e(1:2,E);
    x = p(1,nodes);
    y = p(2,nodes);
    
    ds = sqrt((x(1)-x(2))^2 + (y(1)-y(2))^2);
    k = kappa(mean(x), mean(y));
    gD = diric(mean(x), mean(y));
    R(nodes, nodes) = R(nodes, nodes) + k*[2 1; 1 2]*ds/6;
    r(nodes) = r(nodes) + k*gD*[1; 1]*ds/2;
    
    
end


end
