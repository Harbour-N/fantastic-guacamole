%%%% competition model

% proliefration rate
rho = 1;

% competition rates. must both be less than 1 or both greater than 1 in
% order for a coexisting steady state to exist
a12 = 1.2;
a21 = 0.8;

% function to plot phase plane 
two_speicies_competition_model_phase_plane_plot(a12,a21,rho);



t_interval = [0,100];

%set IC 
u0 = 0.101;
v0 = 0.1;
IC = [u0,v0];

sol = ode45(@(t,y) competition_model(t,y,a12,a21,rho), t_interval, IC);

%test = deval(sol,60);

figure()
plot(sol.y(1,:),sol.y(2,:))
xlim([0,1])
ylim([0,1])
title("sample trajectory")
xlabel("u")
ylabel("v")
hold off

figure()
plot(sol.x,sol.y(1,:))
hold on 
plot(sol.x,sol.y(2,:))
xlabel("t")
ylabel("population density")
legend("u", "v")
title("populations of u and v as functions of time")




%% minimum spanning tree stuff

% run a numeber of test, same IC -> assume all tumours originate from same
% location, same a12, a21 rho -> assume this is the true underlying
% mechanism of state transitions (i.e. this model simulates tumour growth accurately). run for a random ammount of time, to
% replicate biopsy samples take at unknown times since tumour originated

samples = 100;
sol_at_random_times = zeros(2,length(samples));

for i = 1:samples

    t_eval = randi([0,100]);
    sol_at_random_times(:,i)  = deval(sol,t_eval);


end


figure()
plot(sol_at_random_times(1,:), sol_at_random_times(2,:), 'o')
xlim([0,1])
ylim([0,1])
title("sample trajectory at random times")
xlabel("u")
ylabel("v")
hold off

% change IC, i.e. most tumours will have unique composition level of each
% of the two tissue states at the original time of biopsy. due to natural
% variation or microenviromnetal effects. (also we dont know what the
% initialastion in terms of tissue states would have been when the tumour
% 1st appeared)

samples = 100;
for i = 1:samples

    t_eval = randi([0,100]);
    u0 = 0.1*rand; % IC varies close to origin, i.e. tumour starts small
    v0 = 0.1*rand; 
    IC = [u0,v0];
    sol = ode45(@(t,y) competition_model(t,y,a12,a21,rho), t_interval, IC);
    sol_at_random_times_and_IC(:,i)  = deval(sol,t_eval);


end


figure()
plot(sol_at_random_times_and_IC(1,:), sol_at_random_times_and_IC(2,:), 'o')
xlim([0,1])
ylim([0,1])
xlabel("u")
ylabel("v")
hold off




figure()
% x and y coordinates, at which nodes occur
x = sol_at_random_times_and_IC(1,:);
y = sol_at_random_times_and_IC(2,:);

% make a fully connected graph
g  = graph(ones(length(x),length(x)), 'omitselfloops');

% set edge weights equal to distance. Code from: https://uk.mathworks.com/matlabcentral/answers/896827-how-to-make-complete-graph-from-co-ordinates
g.Edges.Weight = zeros(size(g.Edges.EndNodes, 1),1);
for i=1:size(g.Edges, 1)
    g.Edges.Weight(i) = norm([x(g.Edges.EndNodes(i,2))-x(g.Edges.EndNodes(i,1)) ...
        y(g.Edges.EndNodes(i,2))-y(g.Edges.EndNodes(i,1))]);
end

% plot abstract graph i.e. showing all connections but not in phase space
plot(g)

% plot fully connected graph with nodes in the correct position in phase
% space, and weights corresponding to actual distance in phase space
figure()
%h = plot(g,'XData', x,'YData', y, 'EdgeLabel',g.Edges.Weight);
h1 = plot(g,'XData', x,'YData', y);

% calulate minimum spanning tree
T = minspantree(g);

% plot abstract min spanning tree
%figure()
%plot(T)

% plot minimum  spanning tree in phase space 
figure()
%plot(T,'XData', x,'YData', y, 'EdgeLabel',T.Edges.Weight)
h2 = plot(T,'XData', x,'YData', y);
h2.NodeLabel = {}; % use this to remove node labels



%%
%%% functions 

% function to plot phase plane
function two_speicies_competition_model_phase_plane_plot(a12,a21,rho)

    u = 0:0.01:1;
    v = 0:0.01:1;

    % u-nullclines
    null_u = @(u) (1-u)/a12; % non-trivial


    % v-nullclines
    null_v = @(u) 1 - a21*u; % non-trivial

    trivial_nullcline = zeros(1,length(u));

    [U,V] = meshgrid(u,v);

    dUdt = U.*(1 - U - a12.*V) ;
    dVdt = rho.*V.*(1 - V - a21.*U);

    figure()
    h1 = streamslice(U,V,dUdt,dVdt);
    set(h1, 'Color', 'k')
    hold on
    set(groot, 'DefaultTextinterpreter', 'Latex')

    h2 = plot(u,null_u(u),'r');
    plot(trivial_nullcline,v,'r')

    h3 = plot(u,null_v(u), 'b');
    plot(u,trivial_nullcline,'b')

    xlim([0,1])
    ylim([0,1])
    xlabel("u")
    ylabel("v")
    title("u-v phase space trajectories")


    % plot star at steady states
    h4 = plot(0,0,'g*');
    plot(1,0,'g*')
    plot(0,1,'g*')
    plot((1 - a12) / (1 - a12*a21),(1 - a21) / (1 - a12*a21), 'g*')

    legend([h2,h3,h4], "u-nullcline", "v-nulllcine", "steady states")

    hold off

end


%%% functions

% function for the competion model
function dydt = competition_model(t,y,a21,a12,rho)


    u = y(1);
    v = y(2);
    %dydt = [y(1)*(1- y(1) - a12*y(2)); rho*y(2)*(1 - y(2) - a21*y(1))]; % must return a column vector
    dydt = [u*(1- u - a12*v); rho*v*(1 - v - a21*u)];

end








