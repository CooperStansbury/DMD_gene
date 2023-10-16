%% Test SINDy

n = 3;
t = 1000;
X = zeros(n,t);

% Define the 3rd order polynomial function
% poly3rdOrder = @(t, y) [y(1)^3+; y(2)^3; y(3)^3; y(4)^3; y(5)^3];

%% GTP generated test

% Define the 3rd order polynomial function with additional terms
poly3rdOrder = @(t, y) [
    y(1)^3 + y(1)^2*y(3);
    y(2)^3 + y(1)*y(2)*y(3);
    y(3)^3 + y(2)^2*y(4);
    y(4)^3 + y(3)*y(4)*y(5);
    y(5)^3 + y(1)*y(4)*y(5)
];

% Initial conditions
y0 = rand(5,1);
% y0 = [1; 2; 3; 4; 5]; % Replace with your initial conditions

% Time span
tspan = [0:0.01:100]; % Replace with your desired time span

% Solve the ODE using ode45
[t, y] = ode45(poly3rdOrder, tspan, y0);

eps = SINDy(y',3,1,true);
