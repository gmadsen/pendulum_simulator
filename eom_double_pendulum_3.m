clc
clear all
syms tau1 m1 l1 theta1 theta1dot  g
syms tau2 m2 l2 theta2 theta2dot
syms theta1ddot theta2ddot

eqn1 = (m1 + m2)*l1^2*theta1ddot + m2*l2^2*(theta1ddot + theta2ddot) ...
    + m2*l1*l2*(2*theta1ddot + theta2ddot)*cos(theta2) ...
    - m2*l1*l2*(2*theta1dot + theta2dot)*theta2dot*sin(theta2) ...
    + (m1 + m2)*l1*g*sin(theta1) + m2*g*l2*sin(theta1 + theta2) == tau1



eqn2 = m2*l2^2*(theta1ddot + theta2ddot) + m2*l1*l2*theta1ddot*cos(theta2) ...
    + m2*l1*l2*theta1dot^2*sin(theta2) + m2*g*l2*sin(theta1 + theta2) == tau2

solx = solve(eqn1,theta1ddot);
soly = solve(eqn2,theta2ddot);


 f2 = subs(soly,theta1ddot,solx);

 f2_1 = solve(f2,theta2ddot);

 f1 = subs(solx,theta2ddot,f2_1);
 
 f1 = simplify(f1)
 f2_1 = simplify(f2_1,'Steps',100)
 
 eqs = [eqn1,eqn2];
 vars = [theta1ddot,theta2ddot];
 [a,b] = solve(eqs,vars)
 
 simplify(a,'Steps',50)
 simplify(b,'Steps',50)
