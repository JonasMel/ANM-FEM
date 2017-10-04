clear all; close all;
%P1


hmax = 0.1;
g = unitsquare();
[p, e, t] = initmesh(g, 'hmax', hmax);
figure(1);
pdegplot(g);
figure(2);
pdemesh(p, e, t);


x = p(1,:);
y = p(2,:);
u = x.*y;
figure(3);
pdesurf(p, t, u');

phi = zeros(length(p),1);
phi(floor(length(p)/2),1) = 1;
figure(4);
pdesurf(p, t, phi);

x2 = [0; 0; 1];
y2 = [0; 1; 0];


% Calculating area of triangle K
area = polyarea(x2,y2);

% Hat function gradients
b_=[y2(2)-y2(3); y2(3)-y2(1); y2(1)-y2(2)]/2/area;
c_=[x2(3)-x2(2); x2(1)-x2(3); x2(2)-x2(1)]/2/area;

% Triangle element contribution to stiffness matrix A
Ak = (b_*b_' + c_*c_')*area

% forcing-function's contribution to load vector b
bk = eff(x2,y2).*area/3

function f = eff(x2,y2)
    f =3 + 2*x2;
end

function geom = unitsquare()
geom = [2 0 1 0 0 1 0;...
    2 1 1 0 1 1 0;...
    2 1 0 1 1 1 0;...
    2 0 0 1 0 1 0]';
end

