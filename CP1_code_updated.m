clc;        % Clear command line
clear;      % Clear workspace variables
close all   % Closes open figures 
%% Part 2

z1 = .05;  % Zeta for part 2
o1 = 2;     % Omega for part 2
T0 = 1;     % T0 for part 2
t=0:0.1:8;   % time interval 0.1 step size
%{
    Define the f(t, y) function for Part 2.
    Input:
    t: Time step, a scalar
    y: Numerical solution of the current step, 2-element vector, [displacement, velocity]
    Output:
    A 2-element vector as derivative of y
%}

%--------------------------------------------------------------------------
% TODO: Implement the function f(t,y)
%--------------------------------------------------------------------------
f1 = @(t,y) [y(2);(-2*z1*o1*y(2))-(o1^2*y(1))+T0]; %define function f1 from analytical soln.

part1= 1-exp(-z1*o1.*t).*cos(o1.*t.*sqrt(1-z1^2)); %break down f_real for readability
part2= 3.*T0.*z1./(o1^2*sqrt(1-z1^2)); %break down f_real for readability
part3= exp(-z1*o1.*t).*sin(o1.*t.*sqrt(1-z1^2)); %break down f_real for readability
f_real = @(t)(T0/o1^2).*(part1)-(part2).*(part3); %particular soln
y_real= f_real(t)'; %take transpose of f_real

%---------------------------------------------------------------------
% Set up for Part 2
%---------------------------------------------------------------------

tf = 8;                    % Final time for part 2
h=0.1;                     % Step size
t=0:0.1:tf;                % Define array with initial time, step, and final time
y0 = [0.,0.];              % Initial condition
%--------------------------------------------------------------------------
% TODO: Compute the analytical solution
% HINT: 1. Use the vector calculation wisely to avoid for-loops
%       2. Define variables for coefficients to simplify the expressions
%--------------------------------------------------------------------------
ya1 = 0;    % The displacement
ya2 = 0;    % The velocity

%--------------------------------------------------------------------------
% Part 2(1) & (2)
%--------------------------------------------------------------------------
% This generates the numerical solution using RK2
[t2,y2] = rk2_skeleton(f1,y0,tf,h); %call rk2 to find particular soln.

% Generate the numerical solution using ode45 in MATLAB.
t3 = 0; % define initial time
y3 = 0; % define initial position
[tsol, ysol] = ode45(f1, t,y0); %call ode45 to find exact soln.
%--------------------------------------------------------------------------
% TODO: Process and plot the results
%--------------------------------------------------------------------------
plot(tsol,ysol,'.', t2,y2,'*',t,y_real,'-' ) %graph analytical soln, particular soln, and exact soln.

%--------------------------------------------------------------------------
% Part 2(3)
%--------------------------------------------------------------------------                                  
Nh = 8;                            % The number of step sizes to check
yr = [ ya1(end); ya2(end) ];       % Exact solution at last time step

hs = zeros(Nh,1);                  % Allocate the array for step sizes
es = zeros(Nh,1);                  % Allocate the array for errors.
for ii = 1:Nh                     % iterate loop 8 times
   plt = 0.1/2^(ii-1);            % Halving the step size each time 

    %----------------------------------------------------------------------
    % TODO: Implement the rest of the loop
    [t2,y2] = rk2_skeleton(f1,y0,tf,h); %call rk2 function to solve particular soln.
    hs(ii)= plt; %Half step size with each iteration
    
%     part4 = sqrt(ysol(ii,1).^2+ysol(ii,2).^2)
%     part5 = sqrt(y2(ii,1).^2+y2(ii,2).^2)
%     err = abs(part4-part5)
    
    err = sum(abs(y2(:,1)-y_real))/size(y_real,1); %calculate error
    es(ii)= err; %redefine error for next iteration
    %----------------------------------------------------------------------
    % You need to call rk2 using the step size h
    % And then record the step size (in hs) and the 
    % numerical error (in es)
end

%--------------------------------------------------------------------------
% TODO: Process and plot the results
%--------------------------------------------------------------------------
% Here is an example code for the log-log plot revealing the order of
% accuracy
er = es(end)*(hs/hs(end)).^2;
fig2 = figure(2);
plt = loglog(hs, es, 'b-',...
       hs, er, 'r--');
legend([plt(1), plt(2)], {'RK2', '2nd-order reference'},...
        'Location', 'SouthEast')
xlabel('Step Size')
ylabel('Error')
exportgraphics(fig2, './pics/p23.png', 'Resolution', 300,...
                                       'BackgroundColor', 'white');
%% Part 3 (1)

z3 = 0.05;  % Initial zeta for part 3
o3 = 100;   % Omega for part 3
h = 0.001;  % Time step for part 3   

f3 = @(t,y) [0.0; 0.0];   % Function to be integrated for part 3


%--------------------------------------------------------------------------
% Setup for Part 3
%--------------------------------------------------------------------------

tf3 = 0.2;      % Final time for Part 3
y0 = [0.,0.];   % Time step for Part 3
T0 = 1;         % T0 for part 3
zmax = 2.0;     % Max zeta value
dz = 0.1;       % Time steps for zeta values

%--------------------------------------------------------------------------
% TODO: Analyze, design and plot the results.
%--------------------------------------------------------------------------
