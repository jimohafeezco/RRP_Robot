clear; clc;
syms t

% joint velocities when variables are changing 
a2=1;
d1=1;
px = sin(t);
py= cos(2*t);
pz= sin(3*t);

r= sqrt(px^2 +py^2);
s= pz-d1;

q1(t) = atan2(px, py);
q2(t)= atan2(r,s)+pi/2;
d3(t)= sqrt(r^2+s^2);

q1dot(t) = diff(q1, t);
q2dot(t) = diff(q2,t);
d3dot(t) = diff(d3, t);

figure
subplot(3,1,1);
hold on
fplot(@(t) q1dot(t));

fplot(@(t) q1(t));
legend("q_1", "Vel(q1)",'Interpreter','latex')
title("Joint Velocities with changing q")
hold off
subplot(3,1,2)
hold on
fplot(@(t) q2dot(t));
fplot(@(t) q2(t));
legend("q_2", "Vel(q2)",'Interpreter','latex')

hold off
subplot(3,1,3);
hold on
fplot(@(t) d3dot(t));

fplot(@(t) d3(t));
legend("d_3", "Vel(d3)",'Interpreter','latex')

hold off

% Solving jacobian starts here

syms q1 q2 q3 d1 a2 d3 
% 
% a2=1;
% d1=1;
% L(1) = Link([qtest(1), d(1), a(1), alpha(1), 0]);

% T = trotz(q)*transl(0,0,d)*transl(a,0,0)*trotx(alpha);
T1 = trotz(q1) * transl(0, 0, d1) * transl(0, 0, 0) * trotx(pi/2);
R1 = T1(1:3, 1:3);
O1 = T1(1:3, 4);

T2 = trotz(q2) * transl(0, 0, 0) * transl(0, 0, 0) * trotx(pi/2);
R2 = T2(1:3, 1:3);
O2 = T2(1:3, 4);

T3 = trotz(0) * transl(0, 0, a2+d3) * transl(0, 0, 0) * trotx(0);
R3 = T3(1:3, 1:3);
O3 = T3(1:3, 4);

H = T1*T2*T3;
R = H(1:3, 1:3);
% O = end effector
O = H(1:3, 4);

%% Analytical Jacobian

q = [q1; q2; d3];

Jv = [diff(O, q1), diff(O, q2), diff(O, d3)];
%Jv = jacobian(O, q);

Jw = [0, sin(q1), 0;
      0, -cos(q1), 0;
      1, 0, 0];
  
J = [Jv; Jw]

%% Geometrical Jacobian

O0 = [0 0 0]';
z0 = [0 0 1]';

z1 = R1*z0;
z2 = R1*R2*z0;

% The first joint (revolute)
Jv1 = cross(z0, O - O0);
Jw1 = z0;

% The second joint (revolute)
Jv2 = cross(z1, O - O1);
Jw2 = z1;

% The third joint (prismatic)
Jv3 = z2;
Jw3 = 0;

Jv = simplify([Jv1, Jv2, Jv3]);
Jw = simplify([Jw1, Jw2, Jw3]);
J = [Jv; Jw];


% clearvars -except d1 a2 q1 q2 d3 J Jv 
%% Analize determinant

det_Jv = simplify(det(Jv));




syms q1(t) q2(t) d1 a2 d3(t)

a2=1;
d1 = 1;

q1=sin(t);
q2=cos(2*t);
d3=sin(3*t);
% 
X_c = cos(q1)*(a2+d3)/cos(q2);
Y_c = sin(q1)*(a2+d3)/cos(q2);
Z_c = d1+(a2+d3)*tan(q2);
% X_c = O(1);
% Y_c = O(2);
% Z_c= O(3);

X_dot(t) = simplify(diff(X_c,t))
Y_dot(t) = simplify(diff(Y_c,t));
Z_dot(t) = simplify(diff(Z_c,t));


figure;
hold on
fplot(@(t) X_dot(t));
fplot(@(t) Y_dot(t));
fplot(@(t) Z_dot(t));
% title("q_")
legend("q_1", "q_2","d3",'Interpreter','latex')
grid on
hold off



% clf; 
syms px py pz real
time = 0:0.01:15;

px = 2*a2*sin(t);
py = 2*a2*cos(2*t);
pz = d1*sin(3*t);

O = [px py pz]';
O_time = subs(O, {d1, a2, t}, {1, 1, time});
O_time = double(O_time);

syms q_first q_second

% First inverse kinematic solution
q_first = [atan2(py, px), pi/2 + atan2(pz - d1, sqrt(px^2 + py^2)), sqrt(px^2 + py^2 + (pz - d1)^2) - a2]';
q_IK_1 = subs(q_first, {d1, a2, t}, {1, 1, time});
q_IK_1 = double(q_IK_1);



fig = figure(3);
set(gcf,'color','w');
plot(time, q_IK_1, 'LineWidth', 2)
grid on
xlabel('Time')
ylabel('Joint variables')
%title('First variant of inverse kinematics')
legend('q_1', 'q_2', 'd_3', 'Interpreter','latex')
frame = getframe(fig);
im = frame2im(frame);
[img,map] = rgb2ind(im,256);

q_IDK = q_IK_1(:, 1);
for i = 1:(length(time) - 1)
    timenow = time(1, i);
    q_IDK(:, i+1) = q_IDK(:, i) + qdot(timenow, q_IDK(:, i)) * timenow;
end


% 
% syms t d1  a2 q1 q2 d3_t theta1 theta2 d3 real
% 
% % 
% inv_Jv = simplify(inv(Jv));
% pdot = diff(O);
% qdot = simplify(inv_Jv * pdot);
% qdot = subs(qdot, {d1, a2}, {1, 1})
% 
% q_IDK = q_IK_1(:, 1);
% for i = 1:(length(time) - 1)
%     timenow = time(1, i);
%     q_IDK(:, i+1) = q_IDK(:, i) + qdot(timenow, q_IDK(:, i)) * timenow;
% end

% q1=sin(t);
% q2=cos(2*t);
% d3=sin(3*t);
% d1=1;
% a2=1;
% 
% qdot = matlabFunction(qdot, 'Vars', {t,[q1 q2 d3]'});
% 
% 
% q_IDK = q_IK_1(:, 1);
% for i = 1:(length(time) - 1)
%     timenow = time(1, i);
%     q_IDK(:, i+1) = q_IDK(:, i) + qdot(timenow, q_IDK(:, i)) * timenow;
% end
% 
% 
% q0 = q_IK_1(:, 1);
% [time, q_IDK] = ode45(@(t,q_IDK)qdot(t, q_IDK), time, q0);
% q_IDK = q_IDK';
% % qdot = matlabFunction(qdot, 'Vars', {t,[q1, q2 d3]'});
