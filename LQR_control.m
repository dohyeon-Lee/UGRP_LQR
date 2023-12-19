clc;
clear;
%% system matrix
l = 1/21; %1/20;
m = 1.0;
c = 0.58;
g = 9.8;

A = [0 1 0 0;
     0 0 0 0;
     0 0 0 1;
     0 0 -g/l -c/m];
B = [0 1 0 1/l]';
C = [1 0 0 0;
     0 0 1 0];

%% check controllability
Uc = ctrb(A,B);
left_rank = length(A) - rank(Uc);
if left_rank > 0
    disp("system2 is uncontrollable");
    disp("controllability matrix : ");
    disp(Uc);
else 
    disp("system2 is controllable");
end
% controllable 하면 stabilizable 하다.
%% check observability
Uo = obsv(A,C);
left_rank = length(A) - rank(Uo);
if left_rank > 0
    disp("system2 is unobservable");
    disp("observability matrix : ");
    disp(Uo);
else 
    disp("system2 is observable");
end
%% Design LQR controller

% inital condition 
X_init = [0 0 30*(pi/180) 0]';
% LQR controller parameter
Q = [0.1 0;0 0.01];
R = 0.005;
% solve Ricatti equation
[P,L,K] = care(A,B,C'*Q*C,R)
%% simulation
activate_time = 10;
Hz = 50;
t = linspace(0,activate_time,activate_time*Hz);

dt = 1/Hz;
X_packet = zeros(4,activate_time*Hz);
for i = 1:activate_time*Hz
    X_packet(:,i) = expm((A-B*K)*t(i))*X_init;
end
U_packet = -K*X_packet;

% real input(velocity)
V_packet = cumtrapz(t, U_packet);

% desired twist
Vd_zero = zeros(1,activate_time*Hz);
Vd = [Vd_zero; Vd_zero; Vd_zero; Vd_zero; V_packet; Vd_zero]; % wx wy wz vx vy vz
%% calculate cost direct
% method 1
Y_packet = C*X_packet;
cost_1 = 0;
for i = 1:activate_time*Hz  
    cost_1 = cost_1 + (Y_packet(i)*Q*Y_packet(i) + U_packet(i)*R*U_packet(i))*dt;
end
disp("LQR cost calculated by integrate : ");
disp(cost_1);
% method 2
cost_2 = X_init'*P*X_init;
disp("LQR ideal cost:");
disp(cost_2);
%% graph
figure(1);
subplot(511)
grid on;
hold on;
plot(t, X_packet(2,:), 'r','LineWidth',2);
plot(t, X_packet(4,:), 'm','LineWidth',2);
legend('dx2','dtheta');
hold off;
subplot(512);
plot(t, X_packet(1,:), 'b','LineWidth',2);  
legend('x2');
subplot(513)
plot(t, X_packet(3,:), 'c','LineWidth',2);
legend('theta');
subplot(514)
plot(t, U_packet(1,:), 'b','LineWidth',2); 
legend('ddx2 (m/s^2)');
subplot(515)
plot(t, V_packet(1,:), 'r','LineWidth',2); 
legend('dx2 (m/s)');

% make csv
writeU = U_packet(1,:)';
writetheta = X_packet(3,:)'-pi/2;
writetheta_dot = X_packet(4,:)';
write = [writeU writetheta writetheta_dot];
label = ["u(t)" "theta" "theta_dot"];
write = [label;write];
filename = 'test0.csv';
writematrix(write, filename);

% make robot test csv
writeV = V_packet(1,:)';
filename = 'velocity_test.csv';
writematrix(writeV, filename);
 
filename = 'acceleration_test.csv';
writematrix(writeU, filename);
