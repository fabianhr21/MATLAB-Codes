% Mecánica cuántica barrera de potencial finito
% Última fecha de modificación: 20 de noviembre del 2022

clear all; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                      VALORES A EDITAR                   %%%%%%

V_0 = 15; % Potencial de la barrera % 15
a = 1; % Separación de las paredes de potencial a partir del origen.
t_0 = 0.01; % Poner valores de tiempo que se quiere (Valor mínimo 0.01)
infot = '|\Psi(x,t)|^{2} para tiempo t = 0'; % Actualizar con valor de t_0

a_p = 0.05; % Constante de la exponencial de la función en tiempo cero. % 0.05
% A mayor valor de a_p, phik se alarga y |PSI|^2 se hace mas puntiaguida. 
lcte = 3; % Constante de velocidad % 6
Elims = 14; % Rango máximo de la energía para calcular

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dE = 0.1; % Paso para energía
nE = [0:dE:Elims];

hbar = 1; %1.054571817.*10.^(-34);
m = 1; %9.1093837.*10.^(-31); % Masa del electrón
lims = 60; % Limites de x para la simulación [-lims, lims]
dx = 0.1; % Paso de discretización.
x = [-lims:dx:lims]; 
A = (2*a/pi).^(1/4); % Constante de normalización A para un perfil gaussiano

% Separamos x en tres secciones y a cada una le asignaremos un psi_n
x_1 = [-lims:dx:-a-dx];
x_2 = [-a:dx:a-dx];
x_3 = [a:dx:lims];

% Definimos el intervalo de tiempo
tlim = 12;
t = [0:0.01:tlim];

psi = zeros(length(x), length(nE));
T = zeros(1,Elims);
R = zeros(1,Elims);
tic
cont = 1;
for E = 0:dE:Elims
   
    l = sqrt(2.*m.*(E-V_0))./hbar;
    k = sqrt(2.*m.*E)./hbar;
    kapp = sqrt(2.*m.*(V_0-E))./hbar; % Utilizado cuando E<V
    
    % El potencial V_0 abarca el intervalo [-a, a]
    if E > V_0
        % Antes de la barrera: psi_1 = Aexp(ikx) + Bexp(-ikx)
        % En la barrera: psi_2 = Csin(lx) + Dcos(lx)
        % Después de la barrera: psi_3 = Fexp(ikx)
        % De esta manera, psi = [psi_1 psi_2 psi_3] abarca todo el espacio
        
        F = A*exp(-2*1i.*k.*a)./(cos(2*l.*a)-1i.*(k.^2+l.^2).*sin(2.*l.*a)./(2.*k.*l));
        B = 1i*sin(2.*l.*a)./(2.*k.*l).*(l.^2-k.^2).*F;
        C = F.*exp(1i.*k.*a).*(sin(l.*a)+1i.*k.*cos(l.*a)./l);
        D = F.*exp(1i.*k.*a).*(cos(l.*a)-1i.*k.*sin(l.*a)./l);
    
        psi_1 = A.*exp(1i.*k.*x_1) + B.*exp(-1i.*k.*x_1);
        psi_2 = C.*sin(l.*x_2) + D.*cos(l.*x_2);
        psi_3 = F.*exp(1i.*k.*x_3);
    elseif E == V_0
        % Antes de la barrera: psi_1 = Aexp(ikx) + Bexp(-ikx)
        % En la barrera: psi_2 = C + Dx
        % Después de la barrera: psi_3 = Fexp(ikx)
        % De esta manera, psi = [psi_1 psi_2 psi_3] psi abarca todo el espacio
    
        F = A.*exp(-2.*1i.*k.*a)./(1-1i.*k.*a);
        D = 1i.*k.*F.*exp(1i.*k.*a);
        C = F.*exp(1i.*k.*a)-D.*a;
        B = A.*exp(-2.*1i.*k.*a) - F;
    
        psi_1 = A.*exp(1i.*k.*x_1) + B.*exp(-1i.*k.*x_1);
        psi_2 = C + D.*x_2;
        psi_3 = F.*exp(1i.*k.*x_3);
    else 
        % Antes de la barrera: psi_1 = Aexp(ikx) + Bexp(-ikx)
        % En la barrera: psi_2 = Cexp(kapp*.x) + Dexp(-kapp.*x)
        % Después de la barrera: psi_3 = Fexp(ikx)
        % De esta manera, psi = [psi_1 psi_2 psi_3] psi abarca todo el espacio
    
        F = A.*exp(-2.*1i.*k.*a)./(cosh(2.*kapp.*a) + 1i.*(kapp.^2-k.^2).*sinh(2.*kapp.*a)./(2.*kapp.*k));
        C = F.*exp(1i.*k.*a-kapp.*a).*(1+1i.*k./kapp)./2;
        D = F.*exp(1i.*k.*a+kapp.*a).*(1-1i.*k./kapp)./2;
        B = (C.*exp(-kapp.*a) + D.*exp(kapp.*a) - A.*exp(-1i.*k.*a)).*exp(-1i.*k.*a);
        
        psi_1 = A.*exp(1i.*k.*x_1) + B.*exp(-1i.*k.*x_1);
        psi_2 = C.*exp(kapp.*x_2) + D.*exp(-kapp.*x_2);
        psi_3 = F.*exp(1i.*k.*x_3);
    
    end
    T(cont) = (conj(F).*F)/A^2;
    R(cont) = (conj(B).*B)/A^2;

    psi(:,cont) = [psi_1 psi_2 psi_3];
    cont = cont + 1;
end

%%% Coeficientes de transmisión y reflexión



r = 20; % Constante de desplazamiento
k = sqrt(2.*m.*nE)./hbar; % Redifino k pero esta vez como vector

desp = exp(1i*k*r); % Factor de desplazamiento de la onda de phik
vel_wav = exp(-(k-lcte).^2./(4*a_p)); % Factor de velocidad de la onda de phik

phik = 1/(2*pi*a_p)^(1/4)*desp.*vel_wav; % Vector de coeficientes.

% Creamos el vector PSI(x,t)
PSI = zeros(length(x),length(t));

% Definimos el intervalo de tiempo
tlim = 12;
t = [0:0.01:tlim];

% Armamos PSI(x,t) integrando phik, eigenestado y componente temporal en
% todo el rango de energía (0 a inf para nuestro caso).
for tt = 1:length(t)
    for kk = 1:length(x)
        PSI(kk,tt) = trapz(phik.*psi(kk,:).*exp((-1i*hbar*k.^2.*tt.*(0.01))./(2*m)));
    end
end


% Definimos una función de densidad de probabilidad de una vez y la
% normalizamos junto con PSI normal:
PSI2 = abs(PSI).^2;

for ni = 1:length(t)
    PSI(:,ni) = normalize(PSI(:,ni), "norm", 2);
    PSI2(:,ni) = normalize(PSI2(:,ni), "norm", 1).*10;
end

toc

%%%%%%%%%%%%%%%% Gráficas %%%%%%%%%%%%%%%%

% Gráfica de el eigenestado(x,k) y de (PSI(x,t))^2
figure(1)
subplot(1, 2, 1)
title("Graficas del eigenestado(x,k) y la funcion de onda(x,t)", 'fontsize', 17)
surf(x,k,real(psi)');
[caz, cel] = view;
v_sp1 = [5 -10 5];
[caz, cel] = view(v_sp1);
xlabel("x", 'FontSize', 17)
ylabel("k(Energía)", 'FontSize', 17)
zlabel("\psi(x,k)", 'FontSize', 17)
shading interp

subplot(1, 2, 2)
title("Graficas del eigenestado y la función de onda", 'fontsize', 17)
surf(x, t, PSI2')
clb = colorbar;
clb.Label.String = '|\Psi(x,t)|^{2}';
clb.Label.FontSize = 17;
view(2)
xlim([-lims lims])
xlabel("x", 'FontSize', 17)
ylabel("t", 'FontSize', 17)
shading interp

sgtitle("Gráficas del eigenestado(x,k) y el módulo cuadrado de la función de onda(x,t)", 'fontsize', 17)

% Gráficas de evolución temporal de PSI(x)^2 en tres y dos dimensiones
y_limit = 1.3*max(PSI2(:,1));
p_a = [-a y_limit y_limit];
p_b = [-a -y_limit y_limit];
p_c = [-a -y_limit -y_limit/5];
p_d = [-a y_limit -y_limit/5];

p_a2 = [a y_limit y_limit];
p_b2 = [a -y_limit y_limit];
p_c2 = [a -y_limit -y_limit/5];
p_d2 = [a y_limit -y_limit/5];

pts = [p_a; p_b; p_c; p_d; p_a];
pts_2 = [p_a2; p_b2; p_c2; p_d2; p_a2];

% Límites para graficacion en 2D
xlims = lims - 40;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Para ver evolucion temporal descomentar y pegar esto en la barra de comando

% figure(2)
% for y = 1:length(t)-1
%     subplot(1,2,1)
%     axis ([-xlims xlims -y_limit y_limit -y_limit y_limit])
%     plot3(x, imag(PSI(:,y).^2), real(PSI(:,y).^2))
%     hold on
%     plot3(x, imag(PSI2(:,y)), real(PSI2(:,y)))
%     plot3(pts(:,1), pts(:,2), pts(:,3))
%     plot3(pts_2(:,1), pts_2(:,2), pts_2(:,3))
%     hold off
% 
%     subplot(1,2,2)
%     axis ([-xlims xlims -y_limit/5 y_limit])
%     hold on
%     plot(x, real(PSI(:,y).^2))
%     plot(x, PSI2(:,y))
%     xline(-a, 'r:')
%     xline(a, 'r:')
% 
%     drawnow limitrate
%     cla
% end
% hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Gráficas de fotos en tiempos específicos

figure(3)


% Tile 4


plot(x, real(PSI(:,1100).^2))
hold on
plot(x, PSI2(:,1100))
xline(-a, 'r:')
xline(a, 'r:')
xlabel("x", 'Fontsize', 15)
ylabel('|\Psi(x,t)|^{2}', 'Fontsize', 15)
legend("Real(\Psi(x,t)^{2})", "|\Psi(x,t)|^{2}")
title('\Psi(x,t)^{2} transmitida y reflejada')
axis ([-xlims xlims -y_limit/5 y_limit])
% Figura individual para cierto tiempo
figure(4)
tlo_2 = tiledlayout(1,2);
title(tlo_2, infot, 'FontSize', 17)

% Tile 1
nexttile
hold on
plot(x, real(PSI(:,t_0*100).^2))
plot(x, PSI2(:,t_0*100))
xline(-a, 'r:')
xline(a, 'r:')
xlabel("x", 'Fontsize', 15)
ylabel('|\Psi(x,t)|^{2}', 'Fontsize', 15)
legend("Real(\Psi(x,t)^{2})", "|\Psi(x,t)|^{2}")
hold off

% Tile 2
nexttile
axis ([-lims lims -y_limit y_limit -y_limit y_limit])
plot3(x, imag(PSI(:,t_0*100).^2), real(PSI(:,t_0*100).^2))
hold on
plot3(x, imag(PSI2(:,t_0*100)), real(PSI2(:,t_0*100)))
plot3(pts(:,1), pts(:,2), pts(:,3))
plot3(pts_2(:,1), pts_2(:,2), pts_2(:,3))
xlabel("x", 'Fontsize', 15)
ylabel("i", 'Fontsize', 15)
zlabel('|\Psi(x,t)|^{2}', 'Fontsize', 15)
legend("\Psi(x,t)^{2}", "|\Psi(x,t)|^{2}", 'color,''north')
hold off

figure(5)
plot(nE,T,nE,R)
title("Coeficientes de transmisión y reflexión")
xlabel('E')
ylabel("Coeficiente")
legend("Transmisión", "Reflexión")


