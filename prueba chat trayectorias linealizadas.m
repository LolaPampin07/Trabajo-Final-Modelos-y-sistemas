
clc; clear; close all;

%% Parámetros del sistema
a = 0.2; b = 0.2; c = 5.7;

% Condiciones iniciales
xi_0 = [0.01, -0.03, 0.04];
xi_1 = [6, -28, 28];
xi_2 = [1, 1, 1];

% Ventana temporal
t_span = linspace(0, 200, 5000);

% Sistema de Rössler
rossler = @(t,x) [-x(2) - x(3);
                  x(1) + a*x(2);
                  b + x(3)*(x(1) - c)];

% Resolver sistema no lineal
[t, cond1] = ode45(rossler, t_span, xi_0);
[~, cond2] = ode45(rossler, t_span, xi_1);
[~, cond3] = ode45(rossler, t_span, xi_2);

%% Calcular puntos de equilibrio
z_eq = roots([a, -c, b]); % ecuación cuadrática
xeq1 = [a*z_eq(1); -z_eq(1); z_eq(1)];
xeq2 = [a*z_eq(2); -z_eq(2); z_eq(2)];

%% Definir sistemas linealizados (Jacobiano en cada equilibrio)
J = @(xeq) [0, -1, -1;
            1, a, 0;
            xeq(3), 0, xeq(1)-c];

f_lin1 = @(t,X) J(xeq1)*X;
f_lin2 = @(t,X) J(xeq2)*X;

% Condiciones iniciales como perturbaciones
X0_lin1_a = xi_0' - xeq1;
X0_lin1_b = xi_1' - xeq1;
X0_lin1_c = xi_2' - xeq1;

X0_lin2_a = xi_0' - xeq2;
X0_lin2_b = xi_1' - xeq2;
X0_lin2_c = xi_2' - xeq2;

% Simulación linealizada
[t1a, X1a] = ode45(f_lin1, t_span, X0_lin1_a);
[~, X1b] = ode45(f_lin1, t_span, X0_lin1_b);
[~, X1c] = ode45(f_lin1, t_span, X0_lin1_c);

[t2a, X2a] = ode45(f_lin2, t_span, X0_lin2_a);
[~, X2b] = ode45(f_lin2, t_span, X0_lin2_b);
[~, X2c] = ode45(f_lin2, t_span, X0_lin2_c);

%% Convertir soluciones linealizadas a coordenadas originales
X1a_orig = X1a + xeq1';
X1b_orig = X1b + xeq1';
X1c_orig = X1c + xeq1';

X2a_orig = X2a + xeq2';
X2b_orig = X2b + xeq2';
X2c_orig = X2c + xeq2';

%% Gráfico combinado
figure;
hold on;
% Sistema no lineal
plot3(cond1(:,1), cond1(:,2), cond1(:,3), 'g-', 'LineWidth', 1.5);
plot3(cond2(:,1), cond2(:,2), cond2(:,3), 'b-', 'LineWidth', 1.5);
plot3(cond3(:,1), cond3(:,2), cond3(:,3), 'r-', 'LineWidth', 1.5);

% Sistema linealizado (convertido a original)
plot3(X1a_orig(:,1), X1a_orig(:,2), X1a_orig(:,3), 'g--', 'LineWidth', 1.2);
plot3(X1b_orig(:,1), X1b_orig(:,2), X1b_orig(:,3), 'b--', 'LineWidth', 1.2);
plot3(X1c_orig(:,1), X1c_orig(:,2), X1c_orig(:,3), 'r--', 'LineWidth', 1.2);

plot3(X2a_orig(:,1), X2a_orig(:,2), X2a_orig(:,3), 'g:', 'LineWidth', 1.2);
plot3(X2b_orig(:,1), X2b_orig(:,2), X2b_orig(:,3), 'b:', 'LineWidth', 1.2);
plot3(X2c_orig(:,1), X2c_orig(:,2), X2c_orig(:,3), 'r:', 'LineWidth', 1.2);

% Puntos de equilibrio
plot3(xeq1(1), xeq1(2), xeq1(3), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
plot3(xeq2(1), xeq2(2), xeq2(3), 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'k');

xlabel('x'); ylabel('y'); zlabel('z');
title('Comparación: Sistema No Lineal vs Linealizado (en coordenadas originales)');
legend('No lineal xi_0','No lineal xi_1','No lineal xi_2',...
       'Linealizado Xeq1','Linealizado Xeq1','Linealizado Xeq1',...
       'Linealizado Xeq2','Linealizado Xeq2','Linealizado Xeq2',...
       'Xeq1','Xeq2');
