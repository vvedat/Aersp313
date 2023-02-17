function [t, y] = rk2(f, y0, tf, h) %function call for rk2
%{
        Implementation of RK2 algorithm
    Input:
    f:  The function to be integrated, f(t, y)
    y0: Initial condition, 2-element vector, [displacement, velocity]
    tf: Final time, a scalar
    h:  Step size, a scalar
    Output:
    t: Time steps where solutions are computed, N-element vector
        N is number of steps including initial condition, same below
    y: Numerical solutions, Nx2 array
        Columns 1 and 2 contain the displacement and the velocity,
        respectively.
%}
    n = int32(ceil(tf/h))+1;  % determine the number of steps
    t = linspace(0, tf, n) ;  % generate the time step vector
    y = zeros(n,2);           % allocate the array for numerical solutions
    y(1,:) = y0';             % The first row of y is the initial condition
    y_n1 = [0,0];             % initiate values for new y
    
    for ii = 1:(n-1)          % interate loop until 1 less than number of steps
        %----------------------------------------------------------------------
        % TODO: Implement the loop
        %----------------------------------------------------------------------
        k1 = f(t(ii),y(ii,:));  %define k1 for each iteration
        k2 = f(t(ii)+(h/2) , y(ii,:)+(h/2).*k1'); %solve for k2 for each iteration
        y_n1 = y(ii,:)+h.*k2';  %solve for new y value for each iteration
        
        
        y(ii+1,:) = y_n1(1,:); %redefine y for next iteration
    end
end