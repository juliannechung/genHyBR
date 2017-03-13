clear all
clc

%Tests the Toeplitz product Qx in 1D,2D and 3D against full matrix computations


%%%%%%%%%%%%%%%%%%%%%%1d test%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xmin = 0;
xmax = 1;
n = 10;
k = @(r) exp(-r./5.0);
theta = 1.0;

%Build row
Qr = createrow(xmin,xmax,n,k,theta);

%By the brute force approach
x = linspace(xmin,xmax,n);
[X,Y] = meshgrid(x,x);
R = abs(X-Y)./theta;
Q = k(R);
Q1 = Q(1,:);

%Test the rows
err1 = norm(Qr - Q1)./norm(Q1);
fprintf('The 2-norm error in the first row in 1D case is %g\n', err1);

for i = 1:3;
    x = rand(n,1);
    Qx = Q*x;
    Qxtoeplitz = toeplitzproduct(x,Qr,n);
    err = norm(Qx-Qxtoeplitz)./norm(Qxtoeplitz);
    fprintf('The 2-norm error in the mat-vec product in 1D case is %g\n', err);
end

%%%%%%%%%%%%%%%%%%%%%%2d test%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xmin = [0 0];
xmax = [1 1];
n = [10, 4];
k = @(r) exp(-r./5.0);
theta = [1.0 1.0];

%Build row
Qr = createrow(xmin,xmax,n,k,theta);

%By the brute force approach
x = linspace(xmin(1),xmax(1),n(1));    
y = linspace(xmin(2),xmax(2),n(2)); 
[xx,yy] = meshgrid(x,y);
[x1,y1] = meshgrid(xx(:), xx(:));
[x2,y2] = meshgrid(yy(:), yy(:));

R = sqrt((x1-y1).^2./theta(1) + (x2-y2).^2./theta(2));
Q = k(R);
Q1 = Q(1,:);

%Test the rows
err1 = norm(Qr - Q1)./norm(Q1);
fprintf('The 2-norm error in the first row in 2D case is %g\n', err1);

for i = 1:3;
    x = rand(prod(n(1:2)),1);
    Qx = Q*x;
    Qxtoeplitz = toeplitzproduct(x,Qr,n);
    err = norm(Qx-Qxtoeplitz)./norm(Qxtoeplitz);
    fprintf('The 2-norm error in the mat-vec product in 2D case is %g\n', err);
end

%%%%%%%%%%%%%%%%%%%%%%3d test%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xmin = [0 0 0];
xmax = [1 1 1];
n = [10, 20, 2];
k = @(r) exp(-r./5.0);
theta = [1.0 1.0 1.0];

%Build row
Qr = createrow(xmin,xmax,n,k,theta);

%By the brute force approach
x = linspace(xmin(1),xmax(1),n(1));    
y = linspace(xmin(2),xmax(2),n(2)); 
z = linspace(xmin(3),xmax(3),n(3)); 

[xx,yy,zz] = ndgrid(x,y,z);
[x1,y1] = meshgrid(xx(:), xx(:));
[x2,y2] = meshgrid(yy(:), yy(:));
[x3,y3] = meshgrid(zz(:), zz(:));

R = sqrt((x1-y1).^2./theta(1) + (x2-y2).^2./theta(2) + (x3-y3).^2./theta(3));
Q = k(R);
Q1 = Q(1,:);

%Test the rows
err1 = norm(Qr - Q1)./norm(Q1);
fprintf('The 2-norm error in the first row in 3D case is %g\n', err1);

for i = 1:3;
    x = rand(prod(n(1:3)),1);
    Qx = Q*x;
    Qxtoeplitz = toeplitzproduct(x,Qr,n);
    err = norm(Qx-Qxtoeplitz)./norm(Qxtoeplitz);
    fprintf('The 2-norm error in the mat-vec product in 3D case is %g\n', err);
end

