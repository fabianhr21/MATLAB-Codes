l_0= 1; %
k=1;
[x,y] = meshgrid(-5:.5:5,-5:.5:5) ; %to create the points I want to take a look at
 VectorX = -k.*(((x.^2+y.^2).^(1/2))-l_0).*(x./((x.^2+y.^2).^(1/2))); %is the Vector in the "x" direction
 VectorY = -k.*(((x.^2+y.^2).^(1/2))-l_0).*(y./((x.^2+y.^2).^(1/2))); %is the Vector in the "y" direction
 quiver(x,y,VectorX,VectorY)
 hold on
% Sistema Masa-Resorte sin Fuerza Disipativa
 m=1; %masa
 l_0= 1; %longitud del resorte
 k=1; %constante del resorte
 x=[]; 
 y=[];
 x(1)=1; %Posición inicial x
 y(1)=0; %&  y
 v_x=[];
 v_y=[];
 v_x(1)=.15; %Velocidad inicial
 v_y(1)=0.219; %
 a_x=[];
 a_y=[];
 r=(x(1)^2+y(1)^2)^(1/2);
 a_x(1)=-k*((r)-l_0)*(x(1)/(r))/m; %Aceleración inicial
 a_y(1)=-k*((r)-l_0)*(y(1)/(r))/m; %"" 
 dt=.01;
 
 for i=1:6000
     a_x(i+1)= (-k*(((x(i)^2+y(i)^2)^(1/2))-l_0)*(x(i)/((x(i)^2+y(i)^2)^(1/2)))/m);
     a_y(i+1)= (-k*(((x(i)^2+y(i)^2)^(1/2))-l_0)*(y(i)/((x(i)^2+y(i)^2)^(1/2)))/m);
     v_x(i+1)= v_x(i)+a_x(i)*dt;
     v_y(i+1)= v_y(i)+a_y(i)*dt;
     x(i+1)= x(i)+v_x(i)*dt;
     y(i+1)= y(i)+v_y(i)*dt;
 end
comet(x,y)
hold off

%% Conservación de la energía total mecánica
v  = (v_x.^2 + v_y.^2).^(1/2);
p  = ((x).^2 + (y).^2).^(1/2);
E0 = [];
U  = (1/2).*k.*(p-l_0).^2;
K  = (1/2).*m.* v .^2;
E0 = K + U;


plot(E0)
title("Energía Total")

%% Sistema Masa-Resorte con Fuerza Disipativa
l_0= 1; %
k=1;
[x,y] = meshgrid(-5:.5:5,-5:.5:5) ; %to create the points I want to take a look at
 VectorX = -k.*(((x.^2+y.^2).^(1/2))-l_0).*(x./((x.^2+y.^2).^(1/2))); %is the Vector in the "x" direction
 VectorY = -k.*(((x.^2+y.^2).^(1/2))-l_0).*(y./((x.^2+y.^2).^(1/2))); %is the Vector in the "y" direction
 quiver(x,y,VectorX,VectorY)
 hold on

 m=1; %masa
 l_0= 1; %longitud del resorte
 k=1; %constante del resorte
 x=[]; 
 y=[];
 x(1)=1; %Posición inicial x
 y(1)=0; %&  y
 v_x=[];
 v_y=[];
 v_x(1)=.15; %Velocidad inicial
 v_y(1)=.27; %
 a_x=[];
 a_y=[];
 r=(x(1)^2+y(1)^2)^(1/2);
 a_x(1)=-k*((r)-l_0)*(x(1)/(r))/m; %Aceleración inicial
 a_y(1)=-k*((r)-l_0)*(y(1)/(r))/m; %"" 
 dt=.001;
 b = 0.02;                         %Constante de amortiguación
a_dx = [];
a_dy = [];
a_dx(1) = -((b)*v_x(1))/m; %Fuerza Disipativa
a_dy(1) = -((b)*v_y(1))/m;

 for i=1:600000
     a_x(i+1)= (-k*(((x(i)^2+y(i)^2)^(1/2))-l_0)*(x(i)/((x(i)^2+y(i)^2)^(1/2)))/m) + a_dx(i);
     a_y(i+1)= (-k*(((x(i)^2+y(i)^2)^(1/2))-l_0)*(y(i)/((x(i)^2+y(i)^2)^(1/2)))/m) + a_dy(i);
     v_x(i+1)= v_x(i)+a_x(i)*dt;
     v_y(i+1)= v_y(i)+a_y(i)*dt;
     x(i+1)= x(i)+v_x(i)*dt;
     y(i+1)= y(i)+v_y(i)*dt;
     a_dx(i+1) = -((b)*v_x(i))/m; %Fuerza Disipativa
     a_dy(i+1) = -((b)*v_y(i))/m;
 end
plot(x,y)
hold off

%% Conservación de la energía total mecánica
v  = (v_x.^2 + v_y.^2).^(1/2);
p  = ((x).^2 + (y).^2).^(1/2);
E0 = [];
U  = (1/2).*k.*(p-l_0).^2;
K  = (1/2).*m.* v .^2;
E0 = K + U;

plot(E0)
title("Energía Total")
