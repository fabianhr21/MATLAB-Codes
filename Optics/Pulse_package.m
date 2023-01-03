clear
clc
close all
N=2^8;L=50;dx=L/N;dk=2*pi/L;
x=[(-N/2):1:(N/2-1)]*dx;
k=[(-N/2):1:(N/2-1)]*dk;
kshift=fftshift(k);kshift2=kshift.^2;

A_n = 5;        % Amplitud
W_n = 1/A_n;    % Anchura de pulsos, relación An = 1/Wn para minimizar distorsión
q = 2;       % Distancia entre puntos, mientras mayor se pueden ver todos los pulsos en la misma grafica
beta = -0.25;   % Bit-rate Gb/s
T_o = 10;       % Reescalamiento de anchura ps
L_d = (T_o)^2 ./ abs(beta); % Longitud de dispersión de fibra óptica
B = 20;                     % Parametro de dispersión ps2/km
z_real = 10800;             %Distancia en km de propagación 

u0 = zeros(1,length(x));
for n = -2:1:2
    uo=A_n.*sech((x + n.*q)./W_n);
    u0(1,:) = u0 + uo;
end

%figure(1);
%plot(x,abs(u0))

dz=dx.^2/4;zfinal=z_real/L_d;
pasos=ceil(zfinal/dz);

utotal=zeros(N,pasos+1);
utotal(:,1)=uo;
un=u0;

for cuenta=1:1:pasos
    F_NL=fft(exp(1i*dz*abs(un).^2).*un);
    F_D=exp(-1i*kshift2*dz/2).*F_NL;
    un=ifft(F_D);
    utotal(:,cuenta)=un;
end

figure(2);
imagesc(abs(utotal));

figure(1);
plot(x,abs(un),'r');