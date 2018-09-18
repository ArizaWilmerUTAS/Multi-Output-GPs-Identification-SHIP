function [i] = plotgpnship(i, t, sys, y,nsys1,nsys2, std,name)
% Function plots results of simulation or one-step-ahead
% prediction. 
%
%% Syntax
% [i] = plotgp(i, t, sys, y, std);
%
%% Description
% In the upper 2/3 of the figure the output together with 95%
% confidence band and target is plotted, in the lower third of the figure
% the error with corresponding 95% confidence band is plotted. Can be used
% to plot in new or currently active figure. 
%
% Input: 
% * i   ... the figure handle, if i==0 function plots in current figure 
% * t   ... the time vector (x-axis) 
% * sys ... the target output 
% * y   ... the predicted output 
% * std ... the predicted standard deviation 
% 
% Output: 
% * i ... the figure handle 
%
% See Also:
% plotgpy, plotgpe
%
% Examples:
% demo_example_gp_simulation.m
%
%%
% * Written by Dejan Petelin

% check sizes of vectors 
sz = size(t);
if((size(t,1)~=sz(1) | size(sys,1)~=sz(1) |size(y,1)~=sz(1) | size(std,1)~=sz(1))...
        | (size(t,2)~=sz(2) | size(sys,2)~=sz(2) |size(y,2)~=sz(2) | size(std,2)~=sz(2)))
    warning(['figure ', num2str(i), ': vectors: t, tt, y, std must be same size']);
    disp(strcat(['t: ', num2str(size(t))])); 
    disp(strcat(['sys: ', num2str(size(sys))])); 
    disp(strcat(['y: ', num2str(size(y))])); 
    disp(strcat(['std: ', num2str(size(std))])); 
    out = -1; 
    return; 
end

ix_plot = 1:length(t);   
% reduce vector if borders are wanted in figure, 
% e.g. ixplot = 2:length(t)-1;

xfill = [t(ix_plot); flipdim(t(ix_plot),1)]; 
yfill = [y(ix_plot)+2*std(ix_plot);flipdim(y(ix_plot)-2*std(ix_plot),1)]; 


% if i==0 use current axis (figure) 
if (i~=0)
figure(i);
end 

% upper part of figure: y 
subplot(3,1,1:2)
fill(xfill, yfill, [7 7 7]/8, 'EdgeColor', [0 0 0]/8);
hold on
plot(t,nsys1, 'm-', 'LineWidth',2); 
plot(t,nsys2, 'g-', 'LineWidth',2);
plot(t,y, 'r', 'LineWidth',2); 
plot(t,sys, 'b:', 'LineWidth',2); 


%% introduce new data here

%%
set(gcf,'Color',[1 1 1]);
set(gca,'FontSize',32);
%ylim([-1.5 1.5]);

hold off 
grid on 
xlabel('t','FontSize', 32'); 
ylabel(name,'FontSize', 32'); 
set(gca,'FontSize',32)
%title(name, 'FontSize', 9)
legend({'\mu \pm 2\sigma', 'RNN1', 'NarxNN','\mu','system'},'FontSize',20); %
AX=axis; 
AX(1:2)=[t(1) t(end)];  
axis(AX); 

% lower part of figure: e  
subplot(3,1,3)
xfill = [t(ix_plot); flipdim(t(ix_plot),1)]; 
yfill = [2*std(ix_plot);zeros(size(std(ix_plot)))]; 
fill(xfill, yfill, [7 7 7]/8, 'EdgeColor', [0 0 0]/8);

hold on 
plot(t,abs(y-sys), 'k-', 'LineWidth',1); 
set(gca,'FontSize',24);
hold off 
legend('2\sigma', '|e|', 'Location','NorthEast'); 
grid 
xlabel('t'); 
ylabel('e'); 
AX=axis; 
AX(1:2)=[t(1) t(end)]; 
axis(AX); 



return 






