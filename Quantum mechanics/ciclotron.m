tic
R = 0.05;   %Radio D's
A = 0.005;  %Distancia entre D
B = [0,0,1];      %Campo magnético
Vol = 5000; %Voltaje
q = 1.6e-19;%Carga de un protón
m = 1.667e-27; %Masa de un proton
w = q*norm(B)/(m);
E = Vol/A;
dt = 1e-13;

r = [0,0,0];  %Posiciones iniciales
x = [];
y = [];
z = [];
ED = 0;
E0 = [];
x(1) = 0;
y(1) = 0;
z(1) = 0;
v = [0,0,0];  %Velocidadades iniciales
F = [0,0,0];
t = 0;        %Tiempo inicial
i = 2;

while norm(r) < R
    if abs(r(1)) < A/2
        F(1) = q*E*sin(w*t); %Campo eléctrico
        k = 1;
    else
        F = q* cross(v,B);  %Campo magnético
        k = 2;
    end
    
    t = t+dt;

    v = v + (F/m).*dt;
    r = r + v.*dt;

    x(i) = r(1);
    y(i) = r(2);
    z(i) = r(3);
    vf(i) = norm(v);
    E0(i) = (0.5)*m*vf(i)^2;
    if k == 2
        ED = ED + (E0(i) - E0(i-1)) ;
    end
    i = i+1;
end

E0 = E0 /q ;
ED = ED /q;
E0(end)

figure(1)
th = linspace(pi/2, -pi/2, 100);
xR = R*cos(th)+dt/2;
yR = R*sin(th);
plot(xR,yR,'r','LineWidth',2);
axis equal;
hold on
plot([dt/2 dt/2], [-R R],'r','LineWidth',2)
xL = -(R*cos(th)+dt/2) ;
yL = -(R*sin(th)) ;
hold on
plot(xL,yL,'r','LineWidth',2)
hold on
plot([-dt/2 -dt/2], [-R R],'r','LineWidth',2)
xlim([-R-2*dt R+2*dt])
ylim([-R-2*dt R+2*dt])
hold on
plot3(x,y,z)
hold off

figure(2)
plot(E0)
title('Energía cinética de la partícula')
xlabel('Tiempo')
ylabel('Energía cinética')

figure(3)
plot(vf)
title('Velocidad de la partícula')
xlabel('Tiempo')
ylabel('m/s')