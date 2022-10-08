% AE4803 Homework #2 DDP Dubin's vehicle
close all;
clear;
syms x y theta u;

f = [cos(theta);sin(theta); u];
state = [x,y,theta];
dt = 0.02;
A = eye(3) + dt *jacobian(f, state);
B = dt * jacobian(f, u);

gamma = 0.5;

Q_f = 500 * eye(3);
Q = 1 * eye(3);
R = .1;

% Set horizon and target
horizon = 50 * 3;
target = [2;2; pi/2];
    
% Initialize Initial Trajectory
x_nom = zeros(3, horizon);
u_nom = zeros(1, horizon);

% Set number of iterations
num_iter = 10;

% Create variables for Vxx, Vx, V
Vxx = zeros(3,3,horizon);
Vx = zeros(3,horizon);
V = zeros(1, horizon);

% Create variables for A (Phi) and B at each time step
A_f = zeros(3,3,horizon);
B_f = zeros(3, horizon);

% Constants
l_xx = Q;
l_uu = R;
l_ux = zeros(1,3);
Q_k = dt * l_xx;
R_k = dt * l_uu;
P_k = dt * l_ux;
for k = 1:num_iter
    for i = 1:horizon - 1
        x_temp = x_nom(:,i) - target;
        l0 = u_nom(i)' *R *u_nom(i) + x_temp' * Q * x_temp;
        l_x = Q * x_temp;
        l_u = R * u_nom(i);
        q0(i) = dt * l0;
        q_k(:,i) = dt * l_x;
        r_k(:,i) = dt * l_u;

        
        A_f(:,:,i) = subs(A,{x,y,theta, u},{x_nom(1,i),x_nom(2,i),x_nom(3,i), u_nom(i)});
        B_f(:,i) = subs(B,{x,y,theta,u},{x_nom(1,i),x_nom(2,i),x_nom(3,i), u_nom(i)});

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
    dx = zeros(3,1);
    u_new = zeros(1,horizon);
    for i=1:(horizon-1)
        du = l_k(:,i) + L_k(:,:,i) * dx;
        dx = A_f(:,:,i) * dx + B_f(:,i) * du;
        u_new(i) = u_nom(i) + du;
    end
    
    for i=1:(horizon-1)
        theta1 = x_nom(3,i);
        F = [cos(theta1);sin(theta1);u_new(i)];

        x_nom(:,i+1) = x_nom(:,i) + F * dt;
        
    end
    

    u_nom = u_new;


    hold on
    [a] = calcCost(x_nom,u_nom,target,Q,Q_f,R,dt);
    disp(a);
end

%%
close all;
plot(x_nom(1,:),x_nom(2,:),'k-');
grid on
hold on
xlabel("X position [m]")
ylabel("Y Position [m]")
title("Dubin's Vehicle Trajectory")
plot(2,2,'ro')

