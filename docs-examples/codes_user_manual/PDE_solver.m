function [u,I] = PDE_solver(x,y,mu,sigma,PDE_rhs)
    % INPUTS: 
    % x: equispaced discretization nodes
    % y: random variables (vector of size 1 x N_rv)
    % mu: mean of the random field to be constructed
    % sigma=[sigma_1 ... sigma_N_rv]: variance of each part of the random field 
    % PDE_rhs: function handle of the rhs of the problem, PDE_rhs = @(x) ... where x is a row vector
    % OUTPUTS:
    % u: solution at x
    % I: spatial integral of u 

    N = length(x)-1; 
    h = (x(end)-x(1))/N;

    x_mid = x(1:end-1)+h/2; 
    
    % construct the random field
    N_rv = length(y); 
    what_reg = @(x) max(1,ceil(x*N_rv));
    
    RF = zeros(N_rv,1);
    for i = 1:N_rv
        RF(i) = mu + sigma(i)*y(i);
    end 
    a = RF(what_reg(x_mid));
     
    % approx integrals by mid-point rule
    f_mid = PDE_rhs(x_mid);
    f_int = f_mid(2:end-1);  
    rhs = h/2*[f_int(1);f_int(1:end-1)'+f_int(2:end)';f_int(end)];

    A = spdiags([-1/h*[a(2:end-1);0],1/h*(a(1:end-1)+a(2:end)),-1/h*[0;a(2:end-1)]],[-1 0 1],N-1,N-1);
    u = [0;A\rhs;0]; % adding boundary terms

    % trapezoidal rule
    I = [h/2 h*ones(1,N-1) h/2]*u; 
end
