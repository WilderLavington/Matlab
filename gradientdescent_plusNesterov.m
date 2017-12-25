function [x,iterations] = gradientdescent_plusNesterov(f,gradf,tol,iter,xinit)
%initialization
t_old = rand;
%lambda 
lambda_0 = 1;
%y parameter
y_old = xinit;
%x parameter
x_old = xinit;
%primary loop
best_guess = xinit;
iterations = iter;
for i = 1:iter
    %backwards linesearch 
    t = linesearch(f,gradf,x_old,t_old);
    t_old = t;
    %initial jump by gradient decent
    y_new = x_old - (t)*gradf(x_old);
    %variables for Nesterov 
    lambda_1 = (1+(1+4*lambda_0^2)^.5)/2;
    gamma = (1-lambda_0)/lambda_1;
    lambda_0 = lambda_1;
    %next step 
    x_new = (1-gamma)*y_new+gamma*y_old;
    %stopping conditions
    if norm(gradf(x_new)) < tol
        iterations = i;
        break
    end
    % store info, update variables
    x_old = x_new;
    y_old = y_new;
end
x = y_new;
fprintf('Iterations (Nesterov): ')
disp(i)
end

