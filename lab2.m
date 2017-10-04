clear all; close all;
geom  = unitsquare();
epsilon = [1 0.1 0.075 0.05];
% for i = 1:length(epsilon)
    
MyNewtonSolver(geom, epsilon(2));
% end

function geom = unitsquare()
geom = [2 0 1 0 0 1 0;...
    2 1 1 0 1 1 0;...
    2 1 0 1 1 1 0;...
    2 0 0 1 0 1 0]';
end


function MyNewtonSolver(geom, eps)
[p,e,t] = initmesh(geom,'hmax',0.1);
u = zeros(size(p,2),1);
for k = 1:5
    [J,r] = jacres(p,e,t,u, eps);
    d = J\r;
    u = u+d;
    norm_d(k) = 
end
figure, pdesurf(p,t,u);
end


function [J, r] = jacres(p, e, t, u, eps)

% triangle corner nodes
i = t(1,:);
j = t(2,:);
k = t(3,:);

% triangle midpoints
xc = (p(1,i) + p(1,j) + p(1,k))/3;
yc = (p(2,i) + p(2,j) + p(2,k))/3;

% evaluate u, a, a' and f
uu = (u(i)+u(j)+u(k))/3;
aa = feval('a',uu, eps);
% (numerical differentiation:)
da = (a(uu+1.e-8, eps)-a(uu, eps))/1.e-8;
ff = feval('f',xc,yc);
de = eu(uu);%feval('eu', u);
size_ff = size(ff);
size_aa = size(aa);
size_da = size(da);
% assemble lumped Jacobian and residual
%[Aa, Mm, b] = assema(p,t,aa',de',ff);
%[Ada, ~, ~] = assema(p,t,da',0,0);
[Aa, Mm, b] = assema(p,t,aa',-de',de');
size_Aa = size(Aa);
% size_Ada = size(Ada);
size_u = size(u);
% J = diag(Ada*u)+Aa;
J = Aa-Mm;
r = b-Aa*u;


% enforce BCs (homogeneous Dirichlet)
for i = 1:size(e,2)
    n = e(1,i);
    J(n,:) = 0;
    J(n,n) = 1;
    r(n) = 0;
end
end

function val = a(u, eps)
val = 1;%eps + eps*u.^2;

end

function val = eu(u)
val = exp(-u);
end

function val = f(x,y)
val  = 1;

end

function val = f2(u)
val = exp(-u);
end
