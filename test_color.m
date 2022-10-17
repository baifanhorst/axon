x = linspace(0,1,101);
y = linspace(0,1,101);

[X,Y] = meshgrid(x,y);

J = X.*Y;

figure(1);
ax = gca;
fig = pcolor(X,Y,J);
set(fig, 'EdgeColor', 'none');
colorbar(ax);

x_cir=zeros(1,101);
y_cir=zeros(1,101);
for i=1:101
    x_cir(i) = 0.5 + 0.2*cos(2*pi/100*(i-1));
    y_cir(i) = 0.5 + 0.2*sin(2*pi/100*(i-1));
end


plot(ax, x_cir, y_cir, 'r');
