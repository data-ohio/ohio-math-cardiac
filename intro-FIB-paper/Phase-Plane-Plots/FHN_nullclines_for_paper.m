% FHN_Nullclines
% include the vector field and an impulse trajectory
% Set up for the first figure of the paper

% % healthy cell parameters
% sigma = 1.0;
% epsilon = .004;
% a = 0.02;
% b = 2.9;
% c = 1.0;

% % Functional Long QT cell 
% sigma = 0.6;
% epsilon = .004;
% a = 0.01;
% b = 4.5;
% c = 1.0;
 
% % Fast Epsilon parameters
% sigma = 1.6;
% epsilon = .006;
% a = 0.02;
% b = 2.9;
% c = 1.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unused settings %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% aged cell parameters
% a = 0.2;
% b = 1.5;
% c = .9;

% Dysfunctional Long QT cell
% a = 0.1;
% b = 5.5;
% c = 1.0;

% % Epsilon Long QT parameters
% a = 0.1;
% b = 1.5;
% c = 1.0;
% epsilon = .005

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

% healthy cell, Long QT, Fast Epsilon Parameters
sigma = [ 1.0 0.6 1.6];
epsilon = [ .004 .004 .006 ];
a = [ 0.02 0.01 0.02 ]; 
b = [ 2.9 4.5 2.9 ];
c = 1.0;

% time span segmented to allow for arrows
tspans = [ 0 1.5 ; 1.5 2.3; 2.3 31 ; 31 79.7; 79.7 80.4; 80.4 119; 119 300 ;
           0 1.4 ; 1.4 2.5; 2.5 35 ; 35 190; 190 192; 192 230; 230 350 ;
           0 1.0 ; 1.0 1.5; 1.5 21 ; 21 52.5; 52.5 53; 53 79; 79 200 ];

casename = ["healthy", "fast", "slow" ]

% program parameters
linewidthvalue = 2;
linestyle = '-'; 

% define the differential equations for Fitz-Hugh Nagumo
% with an initial impulse
function dydt = odefcn(t,y,a,b,c,epsilon,sigma)
  dydt = zeros(2,1);
   dydt(1) = sigma*(y(1).*(y(1)-a).*(c-y(1)) - y(2)) + 4*(t>0) - 4*(t>.1);
    dydt(2) = epsilon*(y(1) - b*y(2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin the calculations %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% three cases
for k = 1:3

% calculate values of nullclines
v_1 = -0.4:.001:1.2;

% plot line
wline = v_1/b(k);
fig = figure();
plot(v_1,wline,'LineWidth',linewidthvalue)
hold on

% plot cubic
wcubic = v_1.*(v_1-a(k)).*(c-v_1);
plot(v_1,wcubic,'LineWidth',linewidthvalue)
xlim([-0.4,1.2])
ylim([-.05,0.2])
axh = gca; % use current axes
line(get(axh,'XLim'), [0 0], 'Color','k', 'LineStyle', linestyle,'LineWidth',linewidthvalue);
line([0 0], get(axh,'YLim'), 'Color','k', 'LineStyle', linestyle,'LineWidth',linewidthvalue);
fontsize(fig, 16, "points")

% add the vector field
[X Y] = meshgrid(-0.35:.1:1.2,-.1:.02:1.2);
V_1 = sigma(k)*(X.*(X-a(k)).*(c-X) - Y);
V_2 = epsilon(k)*(X - b(k)*Y);
quiver(X,Y,V_1,V_2,8,'LineWidth',1);
%hold off
set(axh,'xtick',[])
set(axh,'ytick',[])

% Solve & add a trajectory to the plot
y0 = [0;0];
for j = 1:6
    [T Y] = ode45(@(t,y) odefcn(t,y,a(k),b(k),c,epsilon(k),sigma(k)),tspans(7*(k-1)+j,:),y0);
    h = plot(Y(:,1),Y(:,2),'k','LineWidth',linewidthvalue);
    hold on
    line2arrow(h,'headwidth',25,'headlength',15)
    y0 = Y(end,:)';
end
[T Y] = ode45(@(t,y) odefcn(t,y,a(k),b(k),c,epsilon(k),sigma(k)),tspans(7*k,:),y0);
plot(Y(:,1),Y(:,2),'k','LineWidth',1)
filename = strcat('phase_plane_', casename(k), '_nolabel.png');
exportgraphics(gcf,filename)
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the time series of the trajectory
tmax = 400;
tmin = -50;
y0 = [ 0; 0];
[T Y] = ode45(@(t,y) odefcn(t,y,a(k),b(k),c,epsilon(k),sigma(k)),[tmin tmax],y0,odeset('MaxStep',.05));
figure
plot(T,Y(:,1),'LineWidth',linewidthvalue)
hold on
plot(T,Y(:,2),'LineWidth',linewidthvalue)
axh = gca; % use current axes
set(axh,'xtick',[])
set(axh,'ytick',[])
% xlabel('t')
% ylabel('v_1',"Rotation",0)
axh = gca; % use current axes
line(get(axh,'XLim'), [0 0], 'Color','k', 'LineStyle', linestyle,'LineWidth',linewidthvalue);
line([0 0], get(axh,'YLim'), 'Color','k', 'LineStyle', linestyle,'LineWidth',linewidthvalue);
fontsize(fig, 16, "points")
filename = strcat('volt_time_', casename(k), '_nolabel.png');
exportgraphics(gcf,filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % end of cases