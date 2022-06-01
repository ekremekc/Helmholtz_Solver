clear all
clc
%%
N = 4

x1 = linspace(+1,0,N+1)'
dx = 1/N
%% Generate the positions of the gridpoints for P0 functions, x0
x0 = conv2(x1,[0.5;0.5],'valid')
%% Generate 1st order difference matrix for a P1 function (makes a P0 function)
D01 = eye(N+1) - diag(ones(1,N),+1) 
%%
D01 = D01(1:N,:) 
%%
D01 = D01/dx
% %% Generate the mass matrix for two P1 functions, M11
% M11 = 4*eye(N+1) + diag(ones(1,N),+1) + diag(ones(1,N),-1) 
% %%
% M11(1,1) = 2
% %%
% M11(N+1,N+1) = 2
% %%
% M11 = M11*(dx/6)
%%
M11 = zeros(N+1)
for i=1:N
    M11(i,i)=M11(i,i)+dx/3
    M11(i+1,i)=M11(i+1,i)+dx/6
    M11(i,i+1)=M11(i,i+1)+dx/6
    M11(i+1,i+1)=M11(i+1,i+1)+dx/3
end
%% Generate the mass matrix for two P0 functions, M00
M00 = eye(N)*dx
%% Generate the mass matrix for a P0 function * a P1 function, M01
M01 = eye(N+1) + diag(ones(1,N),+1); M01 = M01(1:N,:)/2*dx; 
%%
k0 = 10.  % b.c. at x = 0:  du/dx = k0 * u(0)
k1 = -10 % b.c. at x = 1:  du/dx = k1 * u(1)
%% Load in the gridpoints and matrices (see mat_FE for comments)
% Set f (P1 function)
f = (x1-0.5).^2
%% Create matrix A in the bulk
A = D01' * M00 * D01
%% Apply boundary condition at 0
A(N+1,N+1) = A(N+1,N+1) + k0 % Change sign because x is ordered backwards 
%% Apply boundary condition at 1
A(1,1)     = A(1,1)     - k1 % Change sign because x is ordered backwards          
%% Find u
u = A\M11 * f % u = inv(A) * mat.M11 * f;
%% plot
plot(x1,u,'k.-')
du = D01*u
slope = (u(2)-u(1))/(x1(2)-x1(1))