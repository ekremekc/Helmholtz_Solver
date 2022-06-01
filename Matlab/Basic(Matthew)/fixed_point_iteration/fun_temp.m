function [] = fun_temp()

% Guess s
s = 0.6;
for nn = 1:20
    % Calculate A(s)
    A = fun_A(s);
    % Find the eigenvalues of A
    lam = eig(A)
    % Choose the eigenvalue closest to s
    s = lam(dsearchn(lam,s));
    % Display
    disp(num2str(s,16))
end
L = fun_LdL(s);
disp(['Determinant of L(s) = ',num2str(det(L))])
A=fun_A(2.61799387799104)
det(A)
end




function [L,dLds] = fun_LdL(s)

L = [ -s , 1           ; ...
       0 , (-0.5 +sin(s))];

dLds = [ -1 , 0 ; ...
          0 , cos(s)];
      
end      
      

function [A] = fun_A(s)

A = [ 0 , 1                 ; ...
      0 , -0.5 + s + sin(s)];
      
end      
      

