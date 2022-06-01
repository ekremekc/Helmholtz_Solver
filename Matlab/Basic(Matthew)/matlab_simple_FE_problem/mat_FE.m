function [pos,mat] = mat_FE(N)
% mat_FE
%
% Generate matrices and vectors for the finite element scheme
%
% INPUT: 
% N          N+1 is the number of gridpoints
%
% OUTPUTS:
% pos.x0     positions of the gridpoints x0_i
% pos.x1     positions of the gridpoints x1_i
% mat.D01    first derivative matrix: (dp/dx)0 = D01*p1
% mat.M00    mass matrix: <f0,g0> = f0^H * M00 * g0
% mat.M01    mass matrix: <f0,g1> = f0^H * M01 * g1
% mat.M11    mass matrix: <f1,g1> = f1^H * M11 * g1
%
% NOTES
% *0 are col vectors representing P0 functions. They contain N elements
% *1 are col vectors representing P1 functions. They contain N+1 elements
% *00 are matrices that map P0 to P0. They contain (N,N) elements
% *01 are matrices that map P1 to P0. They contain (N+1,N) elements
% *11 are matrices that map P1 to P1. They contain (N+1,N+1) elements

%% Generate the discretization points and the building block matrices
% Generate the positions of the gridpoints for P1 functions, x1
x1 = linspace(+1,0,N+1)'; dx = 1/N;
% Generate the positions of the gridpoints for P0 functions, x0
x0 = conv2(x1,[0.5;0.5],'valid');
% Generate 1st order difference matrix for a P1 function (makes a P0 function)
D01 = eye(N+1) - diag(ones(1,N),+1); D01 = D01(1:N,:); D01 = D01/dx;
% Generate the mass matrix for two P1 functions, M11
M11 = 4*eye(N+1) + diag(ones(1,N),+1) + diag(ones(1,N),-1); M11(1,1) = 2; M11(N+1,N+1) = 2; M11 = M11*(dx/6);
% Generate the mass matrix for two P0 functions, M00
M00 = eye(N)*dx;
% Generate the mass matrix for a P0 function * a P1 function, M01
M01 = eye(N+1) + diag(ones(1,N),+1); M01 = M01(1:N,:)/2*dx; 

%% Wrap into mat structure
pos.x0     = x0;
pos.x1     = x1;
mat.D01    = D01;
mat.M11    = M11;
mat.M00    = M00;
mat.M01    = M01;

end