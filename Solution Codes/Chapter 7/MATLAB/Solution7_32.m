function Solution7_32
%
%       Simple FEA program to predict the small deflection of a tensioned cable
%       pinned vertically at both ends and subected to a transverse force
%
%       This code solves problem 7.32 from the text
%       A.F. Bower, 'Solved Problems in Applied Mechanics of Solids'  
%       CRC press, Baton Rouge, 2026
%
%   It was downloaded from
%   https://github.com/albower/Applied_Mechanics_of_Solids
%
%% Parameters
close all
n_elements = 30;  % Number of elements
n_nodes = n_elements+1; % Number of nodes
T = 1;  % Cable tension
L = 1;  % Cable length
Le = L/n_elements; % Length of one element
q0 = 1; % Transverse load
%% Set up the finite element stiffness matrix and force vector;
%  Also calculate the exact deflection at the nodes
r = zeros(1,n_nodes);
x = zeros(n_nodes,1);
u_exact = zeros(n_nodes,1);
K = zeros(n_nodes);
for i=2:n_nodes-1
    K(i,i) = 2;
    K(i,i-1) = -1;
    K(i,i+1) = -1;
    r(i) = 1;
    x(i) = (i-1)*L/(n_nodes-1);
    u_exact(i) = q0*x(i)*(L-x(i))/(2*T);
end
K = K*T/Le;
r = r*q0*Le;
u_exact(1) = 0;
u_exact(n_nodes) = 0;
x(1) = 0;
x(n_nodes) = L;

% Modify the equation system to prescribe zero displacement at the ends
r(1) = 0;
r(n_nodes) = 0;
K(1,1) = 1;
K(n_nodes,n_nodes) = 1;
    
% Solve the equations
u = K\r';

% Plot the FEA solution for deflection along with the exact solution

figure1 = figure('NumberTitle','off','Name','Figure','Color',[1 1 1]); 
axes3 = axes('Parent',figure1,...
    'Position',[0.175 0.207142857142857 0.73 0.717857142857145]);
    hold(axes3,'on');
    set(axes3,'FontName','Times','FontSize',16,...
        'LineWidth',2,'XGrid','on','YGrid','on',...
        'Box','On')
plot(x,u,'DisplayName','FEA','Marker','o','LineStyle','none','Color','k')
plot(x,u_exact,'DisplayName','Exact','LineStyle','-','Color','k','LineWidth',2)

title('Transverse deflection of a tensioned cable');
 ylabel({'Deflection $wT/(q_0L^2)$'},'Interpreter','latex','FontSize',16);
 xlabel({'Position $x/L$'},'Interpreter','latex','FontSize',16);
 legend1 = legend(axes3,'show');
 set(legend1,...
    'Position',[0.369345250628175 0.340079370897922 0.278869035086111 0.221825391006842],...
    'Interpreter','latex',...
    'FontSize',14,...
    'FontName','Times New Roman');



end
