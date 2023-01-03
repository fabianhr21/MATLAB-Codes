clc;clear
close all

%%%%%%% Definición de parámetros %%%%%%%%
pasos = 100; % Si se cambian pasos aqui, cambiar en función 'SED'
ti = -10;
tf = 10;
t = linspace(ti,tf,pasos); %Ventana en la que observamos el pulso
len = length(t);

u_buena = exp(-t.^2); %Primer pulso, u inicial [exp(-t.^2)]
u = u_buena;
h = (tf-ti)/pasos; %la delta x

%%%%%% inicio de cálculos para la propagación %%%%%%%%%%55
dz = h^2/2;
% z = 0:dz:3;
% lz = length(z);

U_RK = RG04(@SED,u,0,5,dz); %Función de RK4

u_r = U_RK(end,2:len+1); % Último vector de resultados

U_calculada = abs(u_r); %Valores reales

figure(1)
plot(t,u_buena, 'g', t, (U_calculada)', 'r');
figure(2)
imagesc(abs(U_RK(:,2:len+1)'));

%imagesc para graficar
% du/dz = iN*u

%%%%%%%%%%%%%%%%%%%%%%%%%% Función del sistema de ED %%%%%%%%%%%%%%%%%%%%%%%
function du = SED(u_func)
%%%% Definición de parámetros
pasos = 100;
ti = -10;
tf = 10;
t = linspace(ti,tf,pasos); %Ventana en la que observamos el pulso
len = length(t);
h = (tf-ti)/pasos;

u = u_func;  

d2u = zeros(1,len);

d2u(1,1) = (u(3) - 2*u(2) + u(1))/h^2; %segunda derv. hacia adelante [1]
for l = 2:len-1 %Segunda derivada centrada [2:len-1]
d2u(1,l) = (u(l+1) - 2*u(l) + u(l-1))/h^2;
end
d2u(1,len) = (u(len) - 2*u(len-1) + u(len-2))/h^2; %segundaderv. hacia atras [len]

abu = (abs(u).^2).*u;
N = ((1i*d2u./2)+(1i*abu));
du = N;
end


%%%%%%%%%%%%%%% RK4 %%%%%%%%%%%%
function U = RG04(func,u0,ti,tf,h)
stps = int16((tf-ti)/h); % Número de pasos
nd = length(u0); % Se utiliza para definir el número de ED en el sistema 
                 % y por lo tanto las columnas de resultados

U = zeros(stps,nd+1); % arreglo para almacenar los resultados
                       % se añade un uno para agregar la columna del tiempo 

ut = u0; 

tt = ti;
for j =1:stps 
    U(j,1) = tt; % U es una matriz, con columnas soluciones 
%                del sistema de ecuaciones
    
    for l = 2:nd+1 %%%%%  Define las columnas de la matriz. Cada columna es
                   %%%%%  una solución del sistema de ecuaciones

    U(j,l) = ut(l-1); % Como primeras soluciones van las condiciones iniciales
%                       despues se agregan las soluciones de RK4
    end
    
    k1 = func(ut);
    k2 = func(ut + h.*k1./2);
    k3 = func(ut + h.*k2./2);
    k4 = func(ut + h.*k3);
    
    ut = ut + h*(k1 + 2.*k2 +2.*k3 + k4)/6;
    
    tt = tt+h;
end

end

