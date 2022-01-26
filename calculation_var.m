clear all;clc;close;

syms v m g a b u c S E I R N

 N = 26634116; %Population size
 %~~~~~~~~~~~~~~~parameter~~~~~~~~~~~~~~~~~~~~
 u = 0.144; %rate of natural birth
 v = 0.051; %rate of natural death
 a = 0.010; %rate of disease-related death
 m = 0.974; %rate of recovery
 g = 1.428; %prob of changing from E to I
 b = 6.47; %transmission rate
 %b = 2.91;
% ~~~~~~~~~~~~~~~%R0~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% F = [0 b; 0 0]; %jacobian at DFE
% V = [m+v 0;-m g+a+v]; 
% Vinverse = inv(V); 
% P = F*Vinverse 
 
%~~~~~~~~~~~~~~%DEE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
eqn1 = u*N-v*S-b*I*S/N == 0;
eqn2 = b*I*S/N-(m+v)*E == 0;
eqn3 = m*E-(g+a+v)*I == 0;
eqn4 = g*I-v*R == 0;

%sol = solve([eqn1, eqn2, eqn3],[s, e, i]);
sol = solve([eqn1, eqn2, eqn3, eqn4],[S, E, I, R]);
SSol = double(sol.S)
ESol = double(sol.E)
ISol = double(sol.I)
RSol = double(sol.R)

%~~~~~~~~~~~~~%Jacobian~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% J = jacobian([u*N-v*S-b*I*S/N, b*I*S/N-(m+v)*E, m*E-(g+a+v)*I, g*I-v*R],[S,E,I,R]);
% disp (J)

  