%% Plot bifurcation diagram

% Get the x values (parameter) and y values (state) of the bifurcation
% curve.
par_stst=getpar(stst_branch_wbifs,contpar);
y_stst=getx(stst_branch_wbifs,no_y);

% Get the x values (parameter) and y values (state) of the bifurcation
% points.
par_bf_hopf = bgetpar(stst_branch_wbifs,contpar,'hopf');
y_bf_hopf = bgetx(stst_branch_wbifs,no_y,'hopf');

% Create a figure.
fig_bf = figure('Name',append('Bifurcation diagram for ',namepar),'NumberTitle','off');
clf;
ax1=gca;
hold on;

% Plot the stable branch.
plot(ax1,par_stst(nunst_stst==0),y_stst(nunst_stst==0),'g.','DisplayName','stable');
% Plot the unstable branch.
plot(ax1,par_stst(nunst_stst>0),y_stst(nunst_stst>0),'b.','DisplayName','unstable');
% Plot the bifurcation point.
plot(ax1,par_bf_hopf,y_bf_hopf,'ks','DisplayName','hopf');

legend(ax1,'location','north');
xlabel(ax1,namepar);
ylabel(ax1,label_y);
title(ax1,sprintf(append('Bifurcation diagram for ',namepar)));

