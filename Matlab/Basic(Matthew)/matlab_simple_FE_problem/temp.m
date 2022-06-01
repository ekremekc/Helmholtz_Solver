function [] = temp()
%
% Use the weak form to solve the 1D Poisson equation:
%
%  d2u
%  ---- = f(x)              (1)
%  dx^2
%                                         du                         du
% with Robin boundary conditions k0 * u = -- at x = 0   and k1 * u = -- at x = 1  
%                                         dx                         dx
% 
% Define an inner product:
%
%                /1
%  < f , g > =   |   f * g * dx
%                /0
%
% With this inner product, multiply (1) by test function v(x):
% 
%        d2u
%  < v , ---- > = < v , f >
%        dx^2
%
% Perform integration by parts on the LHS:
%
%  [     du ]1     /1  dv   du          /1
%  [ v * -- ]   -  |   -- * -- * dx  =  |  v * f * dx            (2)
%  [     dx ]0     /0  dx   dx          /0
%
% Note that the boundary term can be written:
%
%  [     du ]1     
%  [ v * -- ]   =  v(1) * k1 * u(1) - v(0) * k0 * u(0)
%  [     dx ]0     
% 
% Now set u and v to be piecewise linear between gridpoints (i.e. P1 functions). 
% With the matrices defined in mat_FE.m the integrals are given (exactly) by:
%
%  /1  dv   du     
%  |   -- * -- * dx = (D01 * v)^T * M00 * (D01 * u) = v^T * (D01^T * M00 * D01) * u  
%  /0  dx   dx    
%
%  /1
%  |  v * f * dx    = v^T * M11 * f
%  /0
%
% We can therefore write (2) as:
%
%  v(1) * k1 * u(1) - v(0) * k0 * u(0) - v^T * (D01^T * M00 * D01) * u = v^T * M11 * f  
%
% The LHS of this can be created in matrix form by the following steps:
%  1. Calculate A := (D01^T * M00 * D01)
%  2. Add k1 to the element of A that corresponds to [v(1),u(1)]
%  3. subtract k0 from the element of A that corresponds to [v(0),u(0)]
%
% We can therefore write (2) as:
%
%  v^T * A * u = v^T * M11 * f              (3)
%
% The solution u is found by satisfying (3) for arbitrary v:
%
%  u = inv(A) * M11 * f                     (4)
% 
% Aha! This is NOT like an eigenvalue problem so v does not have the
% physical meaning I was expecting from eigenvalue problems. 
% I may need to perturb (1) in order to find a solution to v and the
% meaning to v. But (1) is linear so the problem seems trivial. 
% LEAVE FOR THE MOMENT. 

% Set the number of finite elements
N = 40;
% Set the Robin boundary conditions
k0 = 10.;  % b.c. at x = 0:  du/dx = k0 * u(0)
k1 = -10; % b.c. at x = 1:  du/dx = k1 * u(1)
% Load in the gridpoints and matrices (see mat_FE for comments)
[pos,mat] = mat_FE(N);
% Set f (P1 function)
f = fun_f(pos.x1);
% Create matrix A in the bulk
A = mat.D01' * mat.M00 * mat.D01;
% Apply boundary condition at 0
A(N+1,N+1) = A(N+1,N+1) + k0; % Change sign because x is ordered backwards 
% Apply boundary condition at 1
A(1,1)     = A(1,1)     - k1; % Change sign because x is ordered backwards          
% Find u
u = A\mat.M11 * f; % u = inv(A) * mat.M11 * f;

% Plot u and set bottom axis to zero
plot(pos.x1,u,'k.-'); Ylim = ylim; ylim([0,Ylim(2)]); xlabel('x'); ylabel('u'); grid on
title(['Weak form solution to u with du(0)/dx = ',num2str(k0),'*u(0) and du(1)/dx = ',num2str(k1),'*u(1)'])

end

function f = fun_f(x)

% forcing = 'uniform';
forcing = 'parabola';

switch forcing
    case 'uniform'
        f = 1 + 0*x;
    case 'parabola'
        f = (x-0.5).^2;
end

end
