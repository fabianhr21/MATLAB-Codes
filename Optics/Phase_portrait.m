clc; clear
close all
%%% Retrato del movimiento armonico

n = 20;
x = linspace(-1,1,n);
z = linspace(-10,10,n);
dx = linspace(-1,1,n);
%%% Parámetro de la ecuación 
k = 1/2;
%%% Malla que convina los vectores x, dx para hacer matriz
[X, Dx] = meshgrid(x,dx);
%%% Sistema de ecuaciones 
x1 = Dx;
%dx1 = -w^2.*X;
dx1 = 2.*X.*(k-X.^2);
%%% Función del gráfico
quiver(X,Dx,x1,dx1)