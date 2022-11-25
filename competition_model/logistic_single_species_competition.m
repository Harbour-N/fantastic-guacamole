%%% single species logistic competition model



% n inteval of intrest
n = linspace(0,3,100);

% plot n agaginst f(n) = dn/dt
figure()
t = 0; % need to set a t even though not used for this part
plot(n,logistic(t,n))

% t interval on which to numerically solve the ode
t_interval = [0, 5];

% Initial condition
IC = 3;

% this syntax also you to then use deval to eval at specific times
sol = ode45(@(t,y) logistic(t,y), t_interval, IC);

% this syntax also works
% [t,u] = ode45(@(t,y) logistic(t,y), t_interval, IC); 

t = linspace(0,5,100);
y = deval(sol,t);

% t = sol.x;
% y =  sol.y;

figure()
plot(t,y)

% plot solution with multiple different IC on same plot

IC = [0.1,0.8,1.3,2];
figure()
hold on

for i = 1:length(IC)

    sol = ode45(@(t,y) logistic(t,y), t_interval, IC(i));
    % evaluate all solutions at specific times to plot
    y = deval(sol,t);
    plot(t,y)


end

% plot a dotted line at n=1 the steady state
plot(t,ones(length(t)), 'r--')
xlabel("t")
ylabel("n(t)")
title("Logistic growth with different IC")









%%% functions

% function to calulate rhs of logistic growth
function rhs = logistic(t,n) % need this 'junk' t variable for ode45

    rhs  = n.*(1 - n);

end







