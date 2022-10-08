% AE4803 Homework #2 DDP Cart Pendulum
close all;
clear;
syms x1 x2 theta1 theta2 u;
gamma = 0.5;
f = [x2;
    1/(1+0.05 * sin(theta1)^2) * (u + 0.05*sin(theta1) * (theta2^2 + 9.81 * cos(theta1)));
    theta2;
    1/(1+0.05 * sin(theta1)^2) * (-u*cos(theta1)-0.05*theta2^2*cos(theta1)*sin(theta1) - 1.05*9.81*sin(theta1))];
state = [x1,x2,theta1,theta2];
n = length(f);
m = length(u);

dt = 0.01;
A = jacobian(f, state);
B = jacobian(f, u);

Q = 0.1 * eye(n);

Q_f = 100 * eye(4);
Q_f(1,1) = 10000;
Q_f(2,2) = 1000;
Q_f(3,3) = 15000;
Q_f(4,4) = 1000;
R = .1;

% Set horizon and target
horizon = 50 * 10;
target = [0; 0; pi; 0];
    
% Initialize Initial Trajectory
x_nom = zeros(n, horizon);
u_nom = zeros(m, horizon);

% Set number of iterations
num_iter = 100;

% Create variables for Vxx, Vx, V
Vxx = zeros(n,n,horizon);
Vx = zeros(n,horizon);
V = zeros(1, horizon);

% Create variables for A (Phi) and B at each time step
A_f = zeros(n,n,horizon);
B_f = zeros(n, horizon);

% Constants
% Constants
l_xx = Q;
l_uu = R;
l_ux = zeros(1,4);
Q_k = dt * l_xx;
R_k = dt * l_uu;
P_k = dt * l_ux;

for k = 1:num_iter
    for i = 1:horizon

        l0 = u_nom(i)' *R *u_nom(i)  + (x_nom(:,i)-target)'*Q*(x_nom(:,i)-target);

        l_x = Q * (x_nom(:,i)-target);
        l_u = R * u_nom(:,i);

        q0(i) = dt * l0;
        q_k(:,i) = dt * l_x;
        r_k(:,i) = dt * l_u;
            
        [A_k,B_k] = state_control_q2(x_nom(:,i), u_nom(i));
        A_f(:,:,i) = eye(4,4) + dt * A_k;
        B_f(:,i) = dt * B_k;
    end

    % Set the Vxx, Vx, V
    Vxx(:,:,horizon)= Q_f;
    Vx(:,horizon) = Q_f * (x_nom(:,horizon) - target);
    V(horizon) = 0.5 * (x_nom(:,horizon) - target)' * Q_f * (x_nom(:,horizon) - target);
    % Back Propogate
    for j = (horizon-1):-1:1

        H = R_k + B_f(:,j)' * Vxx(:,:,j+1) * B_f(:,j);
        G = P_k + B_f(:,j)' * Vxx(:,:,j+1) * A_f(:,:,j);
        g = r_k(:,j) +  B_f(:,j)' * Vx(:,j+1);


        inv_H = inv(H);
        %feedback
        L_k(:,:,j)= - inv_H * G;
        %feedforward
        l_k (:,j) = - inv_H *g;

        % Old ones
        Vxx(:,:,j) = Q_k+ A_f(:,:,j)' * Vxx(:,:,j+1) * A_f(:,:,j) + L_k(:,:,j)' * H * L_k(:,:,j) + L_k(:,:,j)' * G + G' * L_k(:,:,j);
        Vx(:,j)= q_k(:,j) +  A_f(:,:,j)' *  Vx(:,j+1) + L_k(:,:,j)' * g + G' * l_k(:,j) + L_k(:,:,j)'*H * l_k(:,j);
        V(:,j) = q0(j) + V(j+1)   +   0.5 *  l_k (:,j)' * H * l_k (:,j) + l_k (:,j)' * g;
    end
    % Calculate x_bar
    dx = zeros(n,1);
    u_new = zeros(1,horizon);
    for i=1:(horizon-1)
        du = l_k(:,i) + L_k(:,:,i) * dx;
        dx = A_f(:,:,i) * dx + B_f(:,i) * du;
        u_new(i) = min(max(u_nom(i) + gamma * du,-10),10);
    end
    u_nom = u_new;

    for i=1:(horizon-1)

    F = [x_nom(2,i);
    1/(1+0.05 * sin(x_nom(3,i))^2) * (u_nom(i) + 0.05*sin(x_nom(3,i)) * (x_nom(4,i)^2 + 9.81 * cos(x_nom(3,i))));
    x_nom(4,i);
    1/(1+0.05 * sin(x_nom(3,i))^2) * (-u_nom(i)*cos(x_nom(3,i))-0.05*x_nom(4,i)^2*cos(x_nom(3,i))*sin(x_nom(3,i)) - 1.05*9.81*sin(x_nom(3,i)))];
       
    x_nom(:,i+1) = x_nom(:,i) + F * dt; 
        
    end
    
    hold on
    [a] = calcCost(x_nom,u_nom,target,Q,Q_f,R,dt);
    disp(a);
end

%%
figure(2)
plot(1:horizon, x_nom)
legend("x","\dot{x}","\theta","\dot{\theta}")
figure(3)
plot(1:horizon, u_nom)

%% Section 2
figure(2)
x = x_nom';
for k = 2:length(x)

    plot(x(k,1) + cos(x(k,3)-pi/2), sin(x(k,3)-pi/2),'b.')
    hold on
    rectangle('Position',[x(k,1)-0.2 -0.1 0.4 0.2], 'FaceColor','r')
    plot([x(k,1),x(k,1) + cos(x(k,3)-pi/2)], [0,sin(x(k,3)-pi/2)],'k')
    axis([-3 3 -3 3])
    hold off
    drawnow
end
hold on;
plot(x(:,1) + cos(x(:,3)-pi/2),sin(x(:,3)-pi/2))



