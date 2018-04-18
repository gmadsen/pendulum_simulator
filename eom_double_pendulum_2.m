clc
clear all
syms tau1 m1 l1 theta1 theta1dot  g
syms tau2 m2 l2 theta2 theta2dot 
syms theta1ddot theta2ddot

eqn1 = -tau1/l1 + (m1 + m2)*l1^2*theta1ddot ...
    + m2*l1*l2*theta2ddot*cos(theta1 - theta2) ...
    + m2*l1*l2*theta2dot^2 * sin(theta1 - theta2) ...
    + (m1 + m2)*g*l1*sin(theta1) == 0;



eqn2 = -tau2/(m2*l2) + l2*theta2ddot ...
    + l1*theta1ddot*cos(theta1 - theta2) ...
    - l1*theta1dot^2*sin(theta1 - theta2) ...
    + g*sin(theta2) == 0;

solx = solve(eqn1,theta1ddot);
soly = solve(eqn2,theta2ddot);


 f2 = subs(soly,theta1ddot,solx);
 
 f2_1 = solve(f2,theta2ddot)
 
 f1 = subs(solx,theta2ddot,f2_1)

 



    