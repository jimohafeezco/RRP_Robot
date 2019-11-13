clc;clear;clf;
syms q1 q2 d3 t

a2=1;
d1 = 1;
R0_0 = [0;0;1];
R0_1 = [-sin(q1);cos(q1);0];
R0_2 = [cos(q1)*sin(q2);sin(q1)*sin(q2);cos(q2)];
d0_0 = [0;0;0];
d0_1 = [0;0;d1];
d0_3 = [d3*cos(q1)*sin(q2) - a2*sin(q1);a2*cos(q1) + d3*sin(q1)*sin(q2);d1 + d3*cos(q2)];


% % finding Jacobian through partial differentiation approach
x = d3*cos(q1)*sin(q2) - sin(q1);
y = cos(q1) + d3*sin(q1)*sin(q2);
z =d1 + d3*cos(q2);
J = jacobian([x; y; z], [q1 q2 d3]);



r1 = cross(R0_0,(d0_3-d0_0));
r2 = simplify(cross(R0_1,(d0_3-d0_1)));
r3 = simplify(R0_2);

r4 = R0_0;
r5 = R0_1;
r6 = d0_0;

jacobi = [r1,r2,r3;r4,r5,r6];


JacobiLin = [r1, r2, r3];


px = 2*a2*sin(t);
py= 2*a2*cos(2*t);
pz= d1*sin(3*t);

r= sqrt(px^2 +py^2);
s= pz-d1;

q1(t) = atan2(px, py);
q2(t)= atan2(r,s)+pi/2;
d3(t)= sqrt(r^2+s^2);

q1dot(t) = diff(q1, t);
q2dot(t) = diff(q2,t);
d3dot(t) = diff(d3, t);

subplot(4,1,1);
hold on
fplot(@(t) q1dot(t));

fplot(@(t) q1(t));
legend("q_1", "Vel(q1)",'Interpreter','latex')

hold off
subplot(4,1,2)
hold on
fplot(@(t) q2dot(t));
fplot(@(t) q2(t));
legend("q_2", "Vel(q2)",'Interpreter','latex')

hold off
subplot(4,1,3);
hold on
fplot(@(t) d3dot(t));

fplot(@(t) d3(t));
legend("d_3", "Vel(d3)",'Interpreter','latex')

hold off

kkk=JacobiLin;
% 
invJac = inv(JacobiLin);
% 

px_dot(t) = diff(px,t);
py_dot(t) = diff(py,t);
pz_dot(t) = diff(pz,t);

v= [px_dot(t);py_dot(t);pz_dot(t)];

qq=invJac*v;


