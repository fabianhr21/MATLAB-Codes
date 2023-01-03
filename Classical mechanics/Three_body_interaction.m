clear
m_1 = 1;
m_2 = 9.5e-4;
m_3 = 2.8e-8;
G = 4*pi^2;


% Coordenada del origen
r_0 = [0,0,0];

% Condiciones iniciales 
% ASTRO 1 (SOL)
x1 = 0;
y1 = 0;
z1 = 0;
vx1 = 0;
vy1 = 0;
vz1 = 0;

% ASTRO 2 (JUPITER)
x2 = 5.23;
y2 = 0;
z2 = 0;
vx2 = 0;
vy2 = 2.7366;
vz2 = 0;

% ASTRO 3 (COMETA)
x3 = 5.303;
y3 = 0.02;
z3 = 0.02;
vx3 = -0.3;
vy3 = 2.1;
vz3 = 0;clear
m_1 = 1;
m_2 = 9.5e-4;
m_3 = 2.8e-8;
G = 4*pi^2;


% Coordenada del origen
r_0 = [0,0,0];

% Condiciones iniciales 
% ASTRO 1 (SOL)
x1 = 0;
y1 = 0;
z1 = 0;
vx1 = 0;
vy1 = 0;
vz1 = 0;

% ASTRO 2 (JUPITER)
x2 = 5.23;
y2 = 0;
z2 = 0;
vx2 = 0;
vy2 = 2.7366;
vz2 = 0;

% ASTRO 3 (COMETA)
x3 = 5.303;
y3 = 0.02;
z3 = 0.02;
vx3 = -0.3;
vy3 = 2.1;
vz3 = 0;


r_1(:,:,1) = [x1,y1,z1];
r_2(:,:,1) = [x2,y2,z2];
r_3(:,:,1) = [x3,y3,z3];
v_1(1,:,1) = [vx1,vy1,vz1];
v_2(1,:,1) = [vx2,vy2,vz2];
v_3(1,:,1) = [vx3,vy3,vz3]; 

r12 = [];
r23 = [];
r32 = [];
r12(1) = norm(r_2(1,:,1)-r_1(1,:,1));
r13(1) = norm(r_3(1,:,1)-r_1(1,:,1));
r32(1) = norm(r_2(1,:,1)-r_3(1,:,1));

% Modulos de velocidad
v1t = [];
v2t = [];
v3t = [];
v1t(1) = norm(v_1(1,:,1));
v2t(1) = norm(v_2(1,:,1));
v3t(1) = norm(v_3(1,:,1));
dt = 0.0001;

%% Velocity verlet 
for i = 1:120000
    
    r_21 = (r_2(1,:,i) - r_1(1,:,i)) ./ norm(r_2(1,:,i)-r_1(1,:,i));
    r_31 = (r_3(1,:,i) - r_1(1,:,i)) ./ norm(r_3(1,:,i)-r_1(1,:,i));
    f_21 = G*m_1*m_2 / (norm(r_2(1,:,i)-r_1(1,:,i)))^2;
    f_31 = G*m_1*m_3 / (norm(r_3(1,:,i)-r_1(1,:,i)))^2;
    a_1(1,:,i+1) = (f_21.*r_21 + f_31.*r_31) ./ m_1;

    r_12 = (r_1(1,:,i) - r_2(1,:,i)) ./ norm(r_1(1,:,i)-r_2(1,:,i));
    r_32 = (r_3(1,:,i) - r_2(1,:,i)) ./ norm(r_3(1,:,i)-r_2(1,:,i));
    f_12 = G*m_2*m_1 / (norm(r_1(1,:,i)-r_2(1,:,i)))^2;
    f_32 = G*m_2*m_3 / (norm(r_3(1,:,i)-r_2(1,:,i)))^2;
    a_2(1,:,i+1) = (f_12.*r_12 + f_32.*r_32) ./ m_2;

    r_23 = (r_2(1,:,i) - r_3(1,:,i)) ./ norm(r_2(1,:,i)-r_3(1,:,i));
    r_13 = (r_1(1,:,i) - r_3(1,:,i)) ./ norm(r_1(1,:,i)-r_3(1,:,i));
    f_23 = G*m_3*m_2 / (norm(r_2(1,:,i)-r_3(1,:,i)))^2;
    f_13 = G*m_3*m_1 / (norm(r_1(1,:,i)-r_3(1,:,i)))^2;
    a_3(1,:,i+1) = (f_23.*r_23 + f_13.*r_13) ./ m_3;
    
    
    v_1m = v_1(1,:,i) + (1/2).*dt.*a_1(1,:,i);
    v_2m = v_2(1,:,i) + (1/2).*dt.*(a_2(1,:,i));
    v_3m = v_3(1,:,i) + (1/2).*dt.*(a_3(1,:,i));
    r_1(1,:,i+1) = r_1(1,:,i) + dt .* v_1m;
    r_2(1,:,i+1) = r_2(1,:,i) + dt .* v_2m;
    r_3(1,:,i+1) = r_3(1,:,i) + dt .* v_3m;
    v_1(1,:,i+1) = v_1m + (1/2).*dt.*a_1(1,:,i+1);
    v_2(1,:,i+1) = v_2m + (1/2).*dt.*a_2(1,:,i+1);
    v_3(1,:,i+1) = v_3m + (1/2).*dt.*a_3(1,:,i+1);
    
    v1t(i+1) = norm(v_1(1,:,i+1));
    v2t(i+1) = norm(v_2(1,:,i+1));
    v3t(i+1) = norm(v_3(1,:,i+1));
    M(1,:,i) = (m_1*r_1(1,:,i)+m_2*r_2(1,:,i)+m_3*r_3(1,:,i))/(m_1+m_2+m_3);
    r12(i+1) = norm(r_2(1,:,i+1)-r_1(1,:,i+1));
    r13(i+1) = norm(r_3(1,:,i+1)-r_1(1,:,i+1));
    r32(i+1) = norm(r_2(1,:,i+1)-r_3(1,:,i+1));
    MagM(i) = norm(M(1,:,i));
end

r_1x = squeeze(r_1(:,1,:));
r_1y = squeeze(r_1(:,2,:));
r_1z = squeeze(r_1(:,3,:));
r_2x = squeeze(r_2(:,1,:));
r_2y = squeeze(r_2(:,2,:));
r_2z = squeeze(r_2(:,3,:));
r_3x = squeeze(r_3(:,1,:));
r_3y = squeeze(r_3(:,2,:));
r_3z = squeeze(r_3(:,3,:));

%Vectores de masa
Mx = squeeze(M(:,1,:));
My = squeeze(M(:,2,:));
Mz = squeeze(M(:,3,:));

% Criterio Choque
choque = 0.0004; 

if min(r12) < choque
    disp("Choque 1 y 2")
end
if min(r13) < choque
    disp("Choque 1 y 3")
end
if min(r32) < choque
    disp("Choque 3 y 2")
end
minr12 = min(r12)
minr13 = min(r13)
minr32 = min(r32)

% Energía total mecánica
T = (1/2)*m_1*v1t.^2 + (1/2)*m_2*v2t.^2 + (1/2)*m_3*v3t.^2;
%U = (-G*m_1*m_2) ./ sqrt((r_2x-r_1x).^2+(r_2y-r_1y).^2+(r_2z-r_1z).^2 ) + (-G*m_1*m_3) ./ sqrt((r_3x-r_1x).^2+(r_3y-r_1y).^2+(r_3z-r_1z).^2 )+ (-G*m_2*m_3) ./ sqrt((r_3x-r_2x).^2+(r_3y-r_2y).^2+(r_3z-r_2z).^2 );
U = ((-G*m_1*m_2) ./ r12) + ((-G*m_1*m_3) ./ r13) + ((-G*m_2*m_3) ./ r32);
%U = transpose(U);
E0 = T+U;
figure(1)
plot(E0)
title('Energía total mecánica del sistema')

figure (2)
plot(T)

figure(4)
plot(U)

% Grafica trayectoria
figure(3)
view([45 45 45 ]) 
hold on
plot3(0,0,0,'r*')
plot3(r_2x,r_2y,r_2z,r_3x,r_3y,r_3z)
%title('Trayectoria de un sistema de 3 cuerpos')
xlabel('x')
ylabel('y')
zlabel('z')
legend('M1','M2','M3')
grid on
