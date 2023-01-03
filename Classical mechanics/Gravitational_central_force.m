m_a = 3e-30;
M = 1; 
G = 39.4;

[x,y] = meshgrid(-1:.01:1,-1:.01:1) ; %
 VectorX = (-G.*M.*m_a./(((x.^2+y.^2)))).*(x./((x.^2+y.^2).^(1/2)));
 VectorY = (-G.*M.*m_a./(((x.^2+y.^2)))).*(y./((x.^2+y.^2).^(1/2))); 
 quiver(x,y,VectorX,VectorY)
 
hold on
 %% Sistema fuerza central gravitatoria

 m_a = 3e-30;
 M = 1; 
 G = 39.4;
 a_x=[];
 a_y=[];
 a_x(1) = (((-G*M*m_a)/((x(1)^2+y(1)^2)))*(x(1)/((x(1)^2+y(1)^2)^(1/2))))/m_a;
 a_y(1) = (((-G*M*m_a)/((x(1)^2+y(1)^2)))*(y(1)/((x(1)^2+y(1)^2)^(1/2))))/m_a;
 dt=.000001;


%% Position Verlet
m_a = 3e-30;
M = 1; 
G = 39.4;

[x,y] = meshgrid(-1:.01:1,-1:.01:1) ; %
 VectorX = (-G.*M.*m_a./(((x.^2+y.^2)))).*(x./((x.^2+y.^2).^(1/2)));
 VectorY = (-G.*M.*m_a./(((x.^2+y.^2)))).*(y./((x.^2+y.^2).^(1/2))); 
 quiver(x,y,VectorX,VectorY)
 
hold on
xpv=[]; 
ypv=[];

xpv(1)=2; %Posición inicial x 
ypv(1)=1; 
 v_x(1)=1.6255; %Velocidad inicial ex = 0.5
 v_y(1)=-2.6255; 
x_mid(1) = xpv(1) + (1/2)*dt*v_x(1);
y_mid(1) = ypv(1) + (1/2)*dt*v_y(1);
for i = 1:5000000 %%%%Iteraciones del método Verlet
    a_x(i+1) = (((-G*M*m_a)/((xpv(i)^2+ypv(i)^2)))*(xpv(i)/((xpv(i)^2+ypv(i)^2)^(1/2))))/m_a;
    a_y(i+1) = (((-G*M*m_a)/((xpv(i)^2+ypv(i)^2)))*(ypv(i)/((xpv(i)^2+ypv(i)^2)^(1/2))))/m_a;

    x_mid = xpv(i) + dt*v_x(i)/2;
    y_mid = ypv(i) + dt*v_y(i)/2;
    v_x(i+1) = v_x(i) + dt*a_x(i); 
    v_y(i+1) = v_y(i) + dt*a_y(i); 
    xpv(i+1) = x_mid + (dt*v_x(i+1))/2;
    ypv(i+1) = y_mid + dt*v_y(i+1)/2; 
end
plot(xpv,ypv)
title("Trayectoria")
xlabel("x")
ylabel("y")

hold off
% Excentricidad 
rt = sqrt(xpv.^2 + ypv.^2); 
rt_max = abs(max(rt));
rt_min = abs(min(rt));
a = (rt_max+rt_min)/2;
b = sqrt(rt_max*rt_min);
ex = sqrt(1- (b^2/a^2))

% Segunda ley de kepler
figure(3)
px = v_x*m_a;
py = v_y*m_a;
pz = zeros(1,length(px));

T_px= transpose(px);
T_py= transpose(py);
T_pz= transpose(pz);
T_x= transpose(xpv);
T_y= transpose(ypv);
p = [T_px T_py T_pz];
rk = [T_x T_y T_pz];

L= cross(rk,p);
dA=(cross(rk,p)/2*m_a);
plot(dA,'LineWidth',1)
title("Segunda ley de Kepler")
xlabel("Tiempo")
ylabel("Velocidad Areal")
% Tercera ley de kepler
P = 2*pi*sqrt(a^3/(G*M));
C = a^3 / P^2