clc
clear
close all

%% Parámetros del sistema
a = 0.2; b = 0.2; c = 5.7;

% Condiciones iniciales
xi_0 = [0 0 0];
xi_1 = [1 1 1];
xi_2 = [-1 -1 -1];

% Ventana temporal
t_span = linspace(0, 200, 5000);

% Sistema de Rössler

rossler = @(t,x) [-x(2) - x(3);
           x(1) + a*x(2);
           b + x(3)*(x(1) - c)];

% Resolución
[t, x] = ode45(rossler, t_span, xi_0);
[~, y] = ode45(rossler, t_span, xi_1);
[~, z] = ode45(rossler, t_span, xi_2);

%% Diagrama de fase
figure(1);
plot3(x(:, 1), x(:, 2), x(:, 3), 'g-', 'LineWidth', 1.5); hold on;
plot3(y(:, 1), y(:, 2), y(:, 3), 'b--', 'LineWidth', 1.5);
plot3(z(:, 1), z(:, 2), z(:, 3), 'r-', 'LineWidth', 1.5);
plot3(xi_0(1), xi_0(2), xi_0(3), 'r*', 'MarkerSize',10,'MarkerFaceColor', 'g'); % Inicio X01
plot3(xi_1(1), xi_1(2), xi_1(3), 'go', 'MarkerFaceColor', 'b'); % Inicio x1
plot3(xi_2(1), xi_2(2), xi_2(3), 'bo', 'MarkerFaceColor', 'r'); % Inicio x2
xlabel('x'); ylabel('y'); zlabel('z');
title('Diagrama de Fase del Sistema No Lineal');
legend('Trayectoria 1', 'Trayectoria 2', 'Trayectoria 3', 'Inicio 1', 'Inicio 2', 'Inicio 3');
grid on; hold off;

% Gráfico de trayectorias temporales sistema no linealizado
figure(2)
subplot(3,1,1)
plot(t, x(:,1), 'g'); 
xlabel('Tiempo');
ylabel('x(t)');
xlim([0 200]);
title('Trayectorias Temporales de las soluciones ODE - No Lineal');
grid on
subplot(3,1,2)
plot(t, x(:,2), 'b'); %
xlabel('Tiempo');
ylabel('y(t)');
xlim([0 200]);
grid on
subplot(3,1,3)
plot(t, x(:,3), 'r'); 
ylabel('z(t)');
xlim([0 200]);
grid on;
%% Puntos de equilibrio
F = @(x) [-x(2) - x(3);
           x(1) + a*x(2);
           b + x(3)*(x(1) - c)];

f=@(z) a*z^2-c*z+b; %3ra ecuacion despejada

z_1=fzero(f,0);z_2=fzero(f,20);

xeq1= [z_1*a;-z_1;z_1]; xeq2= [z_2*a;-z_2;z_2];

disp('Puntos de equilibrio:')
disp(xeq1')
disp(xeq2')
%% Matriz Jacobiana

syms x y z

f1 = -y - z;
f2 = x + a*y;
f3 = b + z*(x - c);

Fv = [f1; f2; f3];

% Jacobiano simbólico
J = jacobian(Fv, [x y z]);
disp('Matriz Jacobiana simbólica:')
disp(J)
%% Linearizacion del sistema
disp('Jacobiana en el punto de equilibrio x1:')
J_eval1 = double(subs(J, [x y z], xeq1'));
disp(J_eval1)

f_lin1 = @(t, X) J_eval1 * X;

disp('Jacobiana en el punto de equilibrio x2:')
J_eval2 = double(subs(J, [x y z], xeq2'));
disp(J_eval2)

f_lin2 = @(t, X) J_eval2 * X;

% Perturbación inicial
X01 = xeq1+0.1;
X02 = xeq2+0.1;
t_span = 0:0.01:50;

% Simulación
[t, X1] = ode45(f_lin1, t_span, X01);
[~, X2] = ode45(f_lin2, t_span, X02);

% Gráfico
figure(3);
subplot(2,1,1)
plot(t, X1);
xlabel('Tiempo'); ylabel('Estados');
title('Comportamiento del sistema linealizado al rededor de x1');
legend('x', 'y', 'z');
grid on;
subplot(2,1,2)
plot(t, X2);
xlabel('Tiempo'); ylabel('Estados');
title('Comportamiento del sistema linealizado al rededor de x2');
legend('x', 'y', 'z');
grid on;

%% Autovalores
lambda1 = eig(J_eval1);
lambda2 = eig(J_eval2);
disp('Autovalores calculados en MATLAB:')
disp('Para el primer punto de equilibrio')
disp(lambda1)
disp('Para el segundo punto de equilibrio')
disp(lambda2)
