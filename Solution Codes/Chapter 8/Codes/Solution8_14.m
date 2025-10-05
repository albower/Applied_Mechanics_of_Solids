function Solution8_14

%
%       This code solves problem 8.14 the text
%       A.F. Bower 'Solved Problems in Applied Mechanics of Solids'  
%       CRC press, Baton Rouge, 2026
%
%       It was downloaded from
%       https://github.com/albower/Applied_Mechanics_of_Solids
%

    close all
    dt = 0.025;
    beta1 = 0.5;
    beta2 = 0.;
    nsteps = 1000;

    [time,uplot,vplot] = solve_oscillator(dt,beta1,beta2,nsteps);
    
    figure1 = figure('NumberTitle','off','Name','Figure','Color',[1 1 1]);
    axes3 = axes('Parent',figure1,...
        'Position',[0.175 0.207142857142857 0.73 0.717857142857145]);
    hold(axes3,'on');
    set(axes3,'FontName','Times','FontSize',16,...
        'LineWidth',2,'XGrid','on','YGrid','on',...
        'Box','On')
    plot(time,uplot,'LineStyle','-','Color','k','LineWidth',2,'Displayname','FEA')
    xlabel({'Time $t$'},'Interpreter','latex','FontSize',16);
    ylabel({'Displacement $u$'},'Interpreter','latex','FontSize',16);
    
    
    figure1 = figure('NumberTitle','off','Name','Figure','Color',[1 1 1]);
    axes3 = axes('Parent',figure1,...
        'Position',[0.175 0.207142857142857 0.73 0.717857142857145]);
    hold(axes3,'on');
    set(axes3,'FontName','Times','FontSize',16,...
        'LineWidth',2,'XGrid','on','YGrid','on',...
        'Box','On')
    
    nsteps = 300;
    [time,uplot,vplot] = solve_oscillator(dt,beta1,beta2,nsteps);
    eplot = (uplot.^2 + vplot.^2)/2;
    plot(time,eplot,'LineStyle','--','Color','k','LineWidth',2,'Displayname','\beta_1 = 1/2, \beta_2=0')
    xlabel({'Time $t$'},'Interpreter','latex','FontSize',16);
    ylabel({'Energy $E$'},'Interpreter','latex','FontSize',16);
    
    
    beta1 = 0.5;
    beta2 = 0.5;
    [time,uplot,vplot] = solve_oscillator(dt,beta1,beta2,nsteps);
    eplot = (uplot.^2 + vplot.^2)/2;
    plot(time,eplot,'LineStyle','-.','Color','k','LineWidth',2,'Displayname','\beta_1 = 1/2, \beta_2=1/2')
    
    ylim([0 0.55])
    
    beta1 = 1;
    beta2 = 0.5;
    [time,uplot,vplot] = solve_oscillator(dt,beta1,beta2,nsteps);
    eplot = (uplot.^2 + vplot.^2)/2;
    plot(time,eplot,'LineStyle','-','Color','k','LineWidth',2,'Displayname','\beta_1 = 1, \beta_2=1')
    
    legend1 = legend(axes3,'show');
    set(legend1,...
        'Position',[0.40877977343985 0.290079370897926 0.444642846605608 0.219444438625892],...
        'Interpreter','tex',...
        'FontSize',14,...
        'FontName','Times New Roman');

      
end
function [time,uplot,vplot] = solve_oscillator(dt,beta1,beta2,nsteps)

% Function to predict motion of a harmonic oscillator
% using the Newmark time integration algorithm
%
%   dt   Time step
%   beta1, beta2   Newmark parameters
%   nsteps    No. time steps.

    u = 1;   % Initial displacement
    v = 0;   % Initial velocity

    a = -u;  % Initial acceleration


    MK = 1+dt^2*beta2/2.;

    uplot = zeros(nsteps,1);  % Displacement history
    vplot = uplot;            % Velocity history
    time = uplot;             % Time values
    uplot(1) = u;
    time(1) = 0.;
    for i = 1:nsteps-1
       a1 = -1/MK*(u+dt*v+(1-beta2)*dt^2*a/2);  % Updated accel
       v1 = v + dt*((1-beta1)*a + beta1*a1);    % Updated velocity
       u1 = u + dt*v + dt^2*((1-beta2)*a + beta2*a1)/2.; % Updated displacement
       time(i+1) = time(i) + dt;
       uplot(i+1) = u;
       vplot(i+1) = v;
       u = u1;
       v = v1;
       a = a1;
    end

end
