clc
clear all;
% Constante de gravitación universal
G = 4*pi^2;

% Masa 1, 2 y 3
m1 = 1;
m2 = 9.5e-4;
m3 = 2.8e-8;

% Constantes
k1 = G*m1;
k2 = G*m2;
k3 = G*m3;

n = 500;
t = linspace(0,20,n);
dt = t(2)-t(1);

% Vectores

% Masa 1

x1 = zeros(length(n),1);
y1 = zeros(length(n),1);
z1 = zeros(length(n),1);
vx1 = zeros(length(n),1);
vy1 = zeros(length(n),1);
vz1 = zeros(length(n),1);
ax1 = zeros(length(n),1);
ay1 = zeros(length(n),1);
az1 = zeros(length(n),1);

% Masa 2

x2 = zeros(length(t),1);
y2 = zeros(length(t),1);
z2 = zeros(length(t),1);
vx2 = zeros(length(t),1);
vy2 = zeros(length(t),1);
vz2 = zeros(length(t),1);
ax2 = zeros(length(t),1);
ay2 = zeros(length(t),1);
az2 = zeros(length(t),1);

% Masa 3

x3 = zeros(length(t),1);
y3 = zeros(length(t),1);
z3 = zeros(length(t),1);
vx3 = zeros(length(t),1);
vy3 = zeros(length(t),1);
vz3 = zeros(length(t),1);
ax3 = zeros(length(t),1);
ay3 = zeros(length(t),1);
az3 = zeros(length(t),1);

% Distancia entre los cuerpos

r12 = zeros(length(t),1);
r13 = zeros(length(t),1);
r23 = zeros(length(t),1);

% Condiciones iniciales

% Distancias

x1(1) = 0;
y1(1) = 0;
z1(1) = 0;

x2(1) = 5.23;
y2(1) = 0;
z2(1) = 0;

x3(1) = 5.303;
y3(1) = 0.02;
z3(1) = 0.02;

% Velocidades
vx1(1) = 0;
vy1(1) = 0;
vz1(1) = 0;

vx2(1) = 0;
vy2(1) = 2.7366;
vz2(1) = 0;

vx3(1) = -0.3;
vy3(1) = 2.1;
vz3(1) = 0;

% Distancias entre cuerpos 

r12(1) = ((x2(1)-x1(1))^2 + (y2(1)-y1(1))^2 + (z2(1)-z1(1))^2)^(1/2);
r13(1) = ((x3(1)-x1(1))^2 + (y3(1)-y1(1))^2 + (z3(1)-z1(1))^2)^(1/2);
r23(1) = ((x3(1)-x2(1))^2 + (y3(1)-y2(1))^2 + (z3(1)-z2(1))^2)^(1/2);

% Aceleraciones

ax1(1) = k2/r12(1)^3 * (x2(1)-x1(1)) + k3/r13(1)^3 * (x3(1)-x1(1));
ay1(1) = k2/r12(1)^3 * (y2(1)-y1(1)) + k3/r13(1)^3 * (y3(1)-y1(1));
az1(1) = k2/r12(1)^3 * (z2(1)-z1(1)) + k3/r13(1)^3 * (z3(1)-z1(1));

ax2(1) = -((k1/r12(1)^3 * (x2(1)-x1(1)) + k3/r23(1)^3 * (x3(1)-x2(1))));
ay2(1) = -((k1/r12(1)^3 * (y2(1)-y1(1)) + k3/r23(1)^3 * (y3(1)-y2(1))));
az2(1) = -((k1/r12(1)^3 * (z2(1)-z1(1)) + k3/r23(1)^3 * (z3(1)-z2(1))));

ax3(1) = -(k1/r13(1)^3 * (x3(1)-x1(1))) + k2/r23(1)^3 * (x3(1)-x2(1));
ay3(1) = -(k1/r13(1)^3 * (y3(1)-y1(1))) + k2/r23(1)^3 * (y3(1)-y2(1));
az3(1) = -(k1/r13(1)^3 * (z3(1)-z1(1))) + k2/r23(1)^3 * (z3(1)-z2(1));


for i = 1:length(t)-1
    
    x1(i+1) = x1(i) + vx1(i)*dt + (0.5).*(dt^2)*ax1(i);
    y1(i+1) = y1(i) + vy1(i)*dt + (0.5).*(dt^2)*ay1(i);
    z1(i+1) = z1(i) + vz1(i)*dt + (0.5).*(dt^2)*az1(i);
    
    x2(i+1) = x2(i) + dt*vx2(i) + (0.5).*(dt^2)*ax2(i);
    y2(i+1) = y2(i) + dt*vy2(i) + (0.5).*(dt^2)*ay2(i);
    z2(i+1) = z2(i) + dt*vz2(i) + (0.5).*(dt^2)*az2(i);
    
    x3(i+1) = x3(i) + dt*vx3(i) + (0.5).*(dt^2)*ax3(i);
    y3(i+1) = y3(i) + dt*vy3(i) + (0.5).*(dt^2)*ay3(i);
    z3(i+1) = z3(i) + dt*vz3(i) + (0.5).*(dt^2)*az3(i);
    
    r12(i+1) = ((x2(i+1)-x1(i+1))^2 + (y2(i+1)-y1(i+1))^2 + (z2(i+1)-z1(i+1))^2)^(1/2);
    r13(i+1) = ((x3(i+1)-x1(i+1))^2 + (y3(i+1)-y1(i+1))^2 + (z3(i+1)-z1(i+1))^2)^(1/2);
    r23(i+1) = ((x3(i+1)-x2(i+1))^2 + (y3(i+1)-y2(i+1))^2 + (z3(i+1)-z2(i+1))^2)^(1/2);
   
    ax1(i+1) = (k2/r12(i+1)^3) * (x2(i+1)-x1(i+1)) + (k3/r13(i+1)^3) * (x3(i+1)-x1(i+1));
    ay1(i+1) = (k2/r12(i+1)^3) * (y2(i+1)-y1(i+1)) + (k3/r13(i+1)^3) * (y3(i+1)-y1(i+1));
    az1(i+1) = (k2/r12(i+1)^3) * (z2(i+1)-z1(i+1)) + (k3/r13(i+1)^3) * (z3(i+1)-z1(i+1));
    
    ax2(i+1) = -((k1/r12(i+1)^3) * (x2(i+1)-x1(i+1)) + (k3/r23(i+1)^3) * (x3(i+1)-x2(i+1)));
    ay2(i+1) = -((k1/r12(i+1)^3) * (y2(i+1)-y1(i+1)) + (k3/r23(i+1)^3) * (y3(i+1)-y2(i+1)));
    az2(i+1) = -((k1/r12(i+1)^3) * (z2(i+1)-z1(i+1)) + (k3/r23(i+1)^3) * (z3(i+1)-z2(i+1)));
    
    ax3(i+1) = -(k1/r13(i+1)^3) * (x3(i+1)-x1(i+1)) + (k2/r23(i+1)^3) * (x3(i+1)-x2(i+1));
    ay3(i+1) = -(k1/r13(i+1)^3) * (y3(i+1)-y1(i+1)) + (k2/r23(i+1)^3) * (y3(i+1)-y2(i+1));
    az3(i+1) = -(k1/r13(i+1)^3) * (z3(i+1)-z1(i+1)) + (k2/r23(i+1)^3) * (z3(i+1)-z2(i+1));
    
    vx1(i+1) = vx1(i) + dt*(ax1(i) + ax1(i+1));
    vy1(i+1) = vy1(i) + dt*(ay1(i) + ay1(i+1));
    vz1(i+1) = vz1(i) + dt*(az1(i) + az1(i+1));
    
    vx2(i+1) = vx2(i) + dt*(ax2(i) + ax2(i+1));
    vy2(i+1) = vy2(i) + dt*(ay2(i) + ay2(i+1));
    vz2(i+1) = vz2(i) + dt*(az2(i) + az2(i+1));
    
    vx3(i+1) = vx3(i) + dt*(ax3(i) + ax3(i+1));
    vy3(i+1) = vy3(i) + dt*(ay3(i) + ay3(i+1));
    vz3(i+1) = vz3(i) + dt*(az3(i) + az3(i+1));
end

plot3(x1,y1,z1,'r')
hold on
plot3(x2,y2,z2,'b')
hold on
plot3(x3,y3,z3,'g')
xlabel('x')
ylabel('y')
zlabel('z')
