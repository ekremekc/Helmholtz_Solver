clc 
clear all

tol = 10e-10;
ds = 2*tol;

s0 = 0.6;
k=0;
min_iter = 100;
while abs(ds)>tol & k<min_iter
    
    A = fun_A(s0);
    
    lambda  = eig(A);
    
    s = lambda(dsearchn(lambda,s0));
    
    ds = abs(s-s0);
    
    s0 = s;
    k=k+1;

end

disp(["Eigenvalue is: ",num2str(s)])
ac =fun_AC(s)

function [A] = fun_A(s)

A = [ 0 , 1                 ; ...
      0 , -0.5*s + s + sin(s)];
      
end  

function[set] = fun_AC(s)

A = [0 , 1;...
     0 , -0.5*s+s+sin(s)];

C = [1 , 0;...
       0 , 1];

set = A-s*C;
disp(['Determinant of (A-sC) = ',num2str(det(set))])
end


function [L] = fun_L(s)

L = [ -s , 1                 ; ...
      0 , -0.5*s + sin(s)];
      
end 
