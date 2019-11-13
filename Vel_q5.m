% clear all; clc;
syms q1(t) q2(t) d1 a2 d3(t)

a2=1;
d1 = 1;

q1=sin(t);
q2=cos(2*t);
d3=sin(3*t);
% 
% X_c = cos(q1)*(a2+d3)/cos(q2);
% Y_c = sin(q1)*(a2+d3)/cos(q2);
% Z_c = d1+(a2+d3)*tan(q2);
X_c = d3*cos(q1)*sin(q2) - a2*sin(q1)
Y_c = a2*cos(q1) + d3*sin(q1)*sin(q2)
Z_c= d1 + d3*cos(q2)


X_dot(t) = simplify(diff(X_c,t))
Y_dot(t) = simplify(diff(Y_c,t))
Z_dot(t) = simplify(diff(Z_c,t))


subplot(4,1,1);
fplot(@(t) X_dot(t));
title("Velocity of Joint1")
subplot(4,1,2);
fplot(@(t) Y_dot(t));
title("Velocity of Joint1")
subplot(4,1,3);
fplot(@(t) Z_dot(t));
title("Velocity of Joint1")
subplot(4,1,4);
hold on
fplot(@(t) X_dot(t));
fplot(@(t) Y_dot(t));
fplot(@(t) Z_dot(t));
title("Joints Velocity")

hold off
% fplot(@(t) X_dot(t),'b',@(t) Y_dot(t),'g');
