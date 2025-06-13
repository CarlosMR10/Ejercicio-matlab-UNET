% Ejemplo 4.1.1: Transformada de Laplace de f(t) = 1 - cos(3t)
syms t positive;          % Declara t como variable simbólica positiva
f = 1 - cos(3*t);         % Define la función f(t)
F = laplace(f);           % Calcula la transformada de Laplace
disp('Transformada de Laplace de 1 - cos(3t):');
pretty(F);                % Muestra el resultado de forma legible



% Ejemplo 5.3: Descomposición en fracciones parciales de (2s³ + 5s² + 3s + 6)/(s³ + 6s² + 11s + 6)
num = [2 5 3 6];          % Coeficientes del numerador
den = [1 6 11 6];         % Coeficientes del denominador
[r, p, k] = residue(num, den); % Calcula residuos, polos y términos directos

% Muestra los resultados
disp('Residuos (r):'); disp(r);
disp('Polos (p):'); disp(p);
disp('Términos directos (k):'); disp(k);

% Reconstruye la fracción parcial
disp('Expresión en fracciones parciales:');
fprintf('%.4f/(s - %.4f) + %.4f/(s - %.4f) + %.4f/(s - %.4f) + %.4f\n', r(1), p(1), r(2), p(2), r(3), p(3), k);



% Ejemplo 5.6: Transformada inversa de X(s) = 2/(s(s³ + 4s² + 5s + 2))
syms s;
X = 2/(s^4 + 4*s^3 + 5*s^2 + 2*s); % Define X(s)
x = ilaplace(X);                    % Calcula la transformada inversa
disp('Solución x(t):');
pretty(simplify(x));                % Simplifica y muestra el resultado



% Ejemplo 5.9: Resolver x' + 3x = 0, x(0) = 2
syms x(t) s;
eqn = diff(x, t) + 3*x == 0;       % Define la EDO
cond = x(0) == 2;                   % Condición inicial
xSol = dsolve(eqn, cond);           % Resuelve la EDO
disp('Solución x(t):');
pretty(xSol);                       % Muestra la solución



% Ejemplo 5.11: Resolver x'' + 3x' + 2x = 0, x(0) = a, x'(0) = b
syms x(t) a b;
eqn = diff(x, t, 2) + 3*diff(x, t) + 2*x == 0; % EDO
cond = [x(0) == a, subs(diff(x), t, 0) == b];  % Condiciones iniciales
xSol = dsolve(eqn, cond);           % Resuelve la EDO
disp('Solución x(t):');
pretty(simplify(xSol));             % Muestra la solución simplificada



% Ejemplo 5.12: Resolver x'' + 2x' + 5x = 3, x(0) = 0, x'(0) = 0
syms x(t);
eqn = diff(x, t, 2) + 2*diff(x, t) + 5*x == 3; % EDO no homogénea
cond = [x(0) == 0, subs(diff(x), t, 0) == 0];  % Condiciones iniciales
xSol = dsolve(eqn, cond);           % Resuelve la EDO
disp('Solución x(t):');
pretty(expand(xSol));               % Muestra la solución expandida

% Ejemplo 5.13: Resolver y'' + 3y' + 2y = 5, y(0) = -1, y'(0) = 2
syms y(t);
eqn = diff(y, t, 2) + 3*diff(y, t) + 2*y == 5; % EDO
cond = [y(0) == -1, subs(diff(y), t, 0) == 2]; % Condiciones iniciales
ySol = dsolve(eqn, cond);           % Resuelve la EDO
disp('Solución y(t):');
pretty(ySol);                       % Muestra la solución



% Ejemplo 6.6.2.2 - Especificación de funciones de transferencia

% Método 1: Usando la variable compleja 's'
s = tf('s');
sistema1 = (s^2 + 3*s + 2) / (s^3 + 4.5*s^2 + 6.5*s + 3);
disp('Método 1 - Usando variable s:');
disp(sistema1);

% Método 2: Usando coeficientes del numerador y denominador
numerador = [1 3 2];
denominador = [1 4.5 6.5 3];
sistema2 = tf(numerador, denominador);
disp('Método 2 - Usando coeficientes:');
disp(sistema2);

% Método 3: Usando polos, ceros y ganancia
ceros = [-2; -1];
polos = [-2; -1; -1.5];
ganancia = 1;
sistema3 = zpk(ceros, polos, ganancia);
disp('Método 3 - Usando polos y ceros:');
disp(sistema3);

% Conversión entre formatos
[num, den] = zp2tf(ceros, polos, ganancia);
disp('Numerador convertido:');
disp(num);
disp('Denominador convertido:');
disp(den);

[polos2, ceros2, ganancia2] = tf2zp(numerador, denominador);
disp('Polos convertidos:');
disp(polos2);
disp('Ceros convertidos:');
disp(ceros2);


% Ejemplo 6.6.2.3 - Respuesta temporal de un sistema de primer orden

% Definir el sistema (termómetro de mercurio)
tau = 5; % Constante de tiempo de 5 segundos
s = tf('s');
sistema = 1 / (tau*s + 1);

% Configurar tiempo de simulación
tiempo = 0:0.01:20;

% 1. Respuesta al impulso
figure;
subplot(2,1,1);
impulse(sistema, tiempo);
title('Respuesta al impulso del termómetro');
ylabel('Temperatura');
grid on;

subplot(2,1,2);
[y_imp, t_imp] = impulse(sistema, tiempo);
plot(t_imp, 1 - y_imp, 'r', 'LineWidth', 2);
title('Error en respuesta al impulso');
xlabel('Tiempo (s)');
ylabel('Error');
grid on;

% 2. Respuesta al escalón
figure;
subplot(2,1,1);
step(sistema, tiempo);
hold on;
yline(1, '--k', 'Valor Esperado');
title('Respuesta al escalón del termómetro');
ylabel('Temperatura');
legend('Respuesta', 'Valor Esperado');
grid on;

subplot(2,1,2);
[y_step, t_step] = step(sistema, tiempo);
plot(t_step, 1 - y_step, 'r', 'LineWidth', 2);
yline(0, '--k');
title('Error en respuesta al escalón');
xlabel('Tiempo (s)');
ylabel('Error');
grid on;

% 3. Respuesta a onda cuadrada
figure;
[referencia, tiempo_cuad] = gensig('square', 10, 20, 0.01);
respuesta_cuad = lsim(sistema, referencia, tiempo_cuad);

subplot(2,1,1);
plot(tiempo_cuad, respuesta_cuad, 'b', tiempo_cuad, referencia, 'm--');
title('Respuesta a onda cuadrada');
ylabel('Temperatura');
legend('Respuesta', 'Referencia');
grid on;

subplot(2,1,2);
plot(tiempo_cuad, referencia - respuesta_cuad, 'b');
yline(0, '--k');
title('Error en respuesta a onda cuadrada');
xlabel('Tiempo (s)');
ylabel('Error');
grid on;

% 4. Respuesta a rampa
figure;
entrada_rampa = 1/s^2;
respuesta_rampa = impulse(entrada_rampa*sistema, tiempo);

subplot(2,1,1);
plot(tiempo, respuesta_rampa, 'b', tiempo, tiempo, 'm--');
title('Respuesta a rampa');
ylabel('Temperatura');
legend('Respuesta', 'Referencia');
grid on;

subplot(2,1,2);
error_rampa = tiempo' - respuesta_rampa;
plot(tiempo, error_rampa, 'b');
yline(error_rampa(end), '--g');
yline(0, '--k');
title('Error en respuesta a rampa');
xlabel('Tiempo (s)');
ylabel('Error');
legend('Error', 'Error final', 'Valor Esperado');
grid on;



% Ejemplo 6.7.1 - Diagrama de polos y ceros

% Definir el sistema
z = []; % No hay ceros
p = [-1 -1 -1]; % Polo triple en s = -1
k = 1/8; % Ganancia de 1/8

% Crear el sistema
sistema = zpk(z, p, k);
disp('Sistema definido:');
disp(sistema);

% Convertir a forma de polinomios
[num, den] = zp2tf(z, p, k);
disp('Numerador:');
disp(num);
disp('Denominador:');
disp(den);

% Crear función de transferencia
H = tf(num, den);

% Graficar polos y ceros
figure;
pzmap(H);
title('Diagrama de Polos y Ceros');
grid on;


% Ejemplo 7.6 - Reducción de diagrama de bloques (versión corregida)

% Definir las funciones de transferencia
s = tf('s');  % Solo si tienes el Toolbox instalado
G1 = (s + 1)/s;               % G1 = (s + 1)/s
G2 = s/(s^2 - s + 1);         % G2 = s/(s^2 - s + 1)
G3 = 2;                       % G3 = 2 (ganancia)

% --- Método 1: Reducción algebraica manual ---
Gaux = G1 * G2 / (1 + G2 * G3);
Gcl = Gaux / (1 + Gaux);
Gcl = minreal(Gcl, 1e-3);     % Simplificar con tolerancia
disp('Función de transferencia reducida (método manual):');
disp(Gcl);

% --- Método 2: Usando funciones series y feedback ---
Gaux_fb = feedback(G2, G3, -1);       % Retroalimentación negativa
Gcl_fb = feedback(series(G1, Gaux_fb), 1, -1);
Gcl_fb = minreal(Gcl_fb, 1e-3);       % Misma tolerancia
disp('Función de transferencia reducida (método feedback):');
disp(Gcl_fb);

% --- Comparación robusta ---
% Extraer coeficientes y asegurar misma longitud
[Gcl_num, Gcl_den] = tfdata(Gcl, 'v');
[Gcl_fb_num, Gcl_fb_den] = tfdata(Gcl_fb, 'v');

% Función para igualar longitudes
equalize = @(v1, v2) [zeros(1, max(length(v1), length(v2)) - length(v1)), v1];

Gcl_num_pad = equalize(Gcl_num, Gcl_fb_num);
Gcl_fb_num_pad = equalize(Gcl_fb_num, Gcl_num);

Gcl_den_pad = equalize(Gcl_den, Gcl_fb_den);
Gcl_fb_den_pad = equalize(Gcl_fb_den, Gcl_den);

% Calcular diferencia normalizada
diff_num = norm(Gcl_num_pad - Gcl_fb_num_pad);
diff_den = norm(Gcl_den_pad - Gcl_fb_den_pad);
total_diff = diff_num + diff_den;

disp('Diferencia entre métodos (debería ser cercana a cero):');
disp(total_diff);

%


% Ejemplo 6.16 - Sistema de lazo cerrado para dos tanques calentados (CORREGIDO)

% Definir las funciones de transferencia del proceso
s = tf('s');
GL = 1 / ((s + 1)*(5*s + 1));      % Función para perturbación
GM = (1/2160) / ((s + 1)*(5*s + 1)); % Función para variable manipulada

% Parámetros del sistema
Kc = 1; % Ganancia del controlador (proporcional)
GV = 500000 / 16; % Función de transferencia de la válvula (Btu/min/mA)
H = 16 / 100;     % Función de transferencia del sensor (mA/°F)

% --- Función auxiliar para igualar longitudes de vectores ---
function v_padded = pad_vector(v, desired_length)
    v_padded = [zeros(1, desired_length - length(v)), v];
end

% --- Calcular función de transferencia de lazo cerrado para cambios en carga ---
num_L = GL.num{1};
den_conv = conv(GL.den{1}, [1 0]);  % Multiplicación por 's' (derivada)
term2 = Kc * GV * H * GM.num{1};

% Igualar longitudes antes de sumar
max_len = max(length(den_conv), length(term2));
den_conv_padded = pad_vector(den_conv, max_len);
term2_padded = pad_vector(term2, max_len);

den_L = den_conv_padded + term2_padded;

GCL = tf(num_L, den_L);
GCL = minreal(GCL);
disp('Función de transferencia para cambios en carga:');
disp(GCL);

% --- Calcular función de transferencia de lazo cerrado para cambios en setpoint ---
num_R = Kc * GV * GM.num{1};
den_R = den_L; % Mismo denominador que para GCL

GCR = tf(num_R, den_R);
GCR = minreal(GCR);
disp('Función de transferencia para cambios en setpoint:');
disp(GCR);

% --- Calcular offset para sistema con control proporcional ---
offset = 1 / (1 + Kc * GV * H * (1/2160));
disp(['Offset con Kc = ', num2str(Kc), ': ', num2str(offset)]);

% --- Graficar respuesta al escalón para diferentes ganancias Kc ---
figure;
hold on;
Kc_values = [0.5, 1, 2, 5];
colors = ['b', 'r', 'g', 'm'];

for i = 1:length(Kc_values)
    Kc = Kc_values(i);
    
    % Recalcular denominador para cada Kc
    den_conv = conv(GL.den{1}, [1 0]);
    term2 = Kc * GV * H * GM.num{1};
    
    % Igualar longitudes
    max_len = max(length(den_conv), length(term2));
    den_conv_padded = pad_vector(den_conv, max_len);
    term2_padded = pad_vector(term2, max_len);
    
    den_R = den_conv_padded + term2_padded;
    
    GCR = tf(Kc * GV * GM.num{1}, den_R);
    GCR = minreal(GCR);
    
    % Graficar respuesta al escalón
    [y, t] = step(GCR, 50);
    plot(t, y, colors(i), 'DisplayName', ['Kc = ', num2str(Kc)]);
    
    % Calcular y mostrar línea de offset
    offset = 1 / (1 + Kc * GV * H * (1/2160));
    yline(offset, '--', 'Color', colors(i), 'HandleVisibility', 'off');
end

title('Respuesta al escalón para diferentes ganancias Kc');
xlabel('Tiempo (min)');
ylabel('Temperatura (°F)');
legend('Location', 'southeast');
grid on;
%


% --- Ejemplo 8.1: Respuesta a un escalón unitario ---
% Sistema de primer orden: Termómetro con tau = 0.1 min
% Función de transferencia: G(s) = 1 / (0.1s + 1)

% Definir numerador y denominador de la función de transferencia
num = [0 1];      % Numerador: 1
den = [0.1 1];    % Denominador: 0.1s + 1

% Generar la respuesta al escalón unitario
figure;
step_response = step(num, den);

% Graficar la respuesta
plot(step_response, 'Color', [0 0 0], 'LineWidth', 2);
title('Respuesta a un escalón unitario (Ejemplo 8.1)', 'FontSize', 12);
grid on;
xlabel('Tiempo (min)', 'FontSize', 12);
ylabel('Variación de la temperatura (°C)', 'FontSize', 12);
set(gca, 'FontSize', 12, 'LineWidth', 2);



% --- Ejemplo 8.2: Respuesta a un escalón unitario ---
% Sistema de segundo orden: G(s) = 25 / (s^2 + 4s + 25)

% Definir numerador y denominador de la función de transferencia
num = [0 0 25];    % Numerador: 25
den = [1 4 25];     % Denominador: s^2 + 4s + 25

% Generar la respuesta al escalón unitario
figure;
step_response = step(num, den);

% Graficar la respuesta
plot(step_response, 'Color', [0 0 0], 'LineWidth', 2);
title('Respuesta a un escalón unitario (Ejemplo 8.2)', 'FontSize', 12);
grid on;
xlabel('Tiempo (min)', 'FontSize', 12);
ylabel('Variación de la temperatura (°C)', 'FontSize', 12);
set(gca, 'FontSize', 12, 'LineWidth', 2);



% --- Ejemplo 8.3: Respuesta a un impulso unitario ---
% Sistema de primer orden: Termómetro con tau = 0.1 min
% Función de transferencia: G(s) = 1 / (0.1s + 1)

% Para obtener la respuesta al impulso, multiplicamos G(s) por s
num = [1 0];      % Numerador: s
den = [0.1 1];    % Denominador: 0.1s + 1

% Generar la respuesta al impulso unitario
figure;
impulse_response = step(num, den);  % Usamos step porque multiplicamos G(s) por s

% Graficar la respuesta
plot(impulse_response, 'Color', [0 0 0], 'LineWidth', 2);
title('Respuesta a un impulso unitario (Ejemplo 8.3)', 'FontSize', 12);
grid on;
xlabel('Tiempo (s)', 'FontSize', 12);
ylabel('Variación de la temperatura (°C)', 'FontSize', 12);
set(gca, 'FontSize', 12, 'LineWidth', 2);



% --- Ejemplo 8.4: Respuesta a un impulso unitario ---
% Sistema de segundo orden: G(s) = 1 / (s^2 + 0.2s + 1)

% Para obtener la respuesta al impulso, multiplicamos G(s) por s
num = [0 1 0];      % Numerador: s
den = [1 0.2 1];    % Denominador: s^2 + 0.2s + 1

% Generar la respuesta al impulso unitario
figure;
impulse_response = step(num, den);  % Usamos step porque multiplicamos G(s) por s

% Graficar la respuesta
plot(impulse_response, 'Color', [0 0 0], 'LineWidth', 2);
title('Respuesta a un impulso unitario (Ejemplo 8.4)', 'FontSize', 12);
grid on;
xlabel('Tiempo (s)', 'FontSize', 12);
ylabel('Variación de la temperatura (°C)', 'FontSize', 12);
set(gca, 'FontSize', 12, 'LineWidth', 2);



% --- Ejemplo 8.5: Respuesta a una rampa unitaria ---
% Sistema de primer orden: Termómetro con tau = 0.1 min
% Función de transferencia: G(s) = 1 / (0.1s + 1)

% Para obtener la respuesta a una rampa, dividimos G(s) por s
num = [0 0 1];      % Numerador: 1
den = [0.1 1 0];    % Denominador: 0.1s^2 + s

% Definir el vector de tiempo
t = 0:0.005:0.5;

% Generar la respuesta a la rampa unitaria
ramp_response = step(num, den, t);

% Graficar la respuesta y la rampa de referencia
figure;
plot(t, ramp_response, '-k', t, t, '-b', 'LineWidth', 2);
title('Respuesta a una rampa unitaria (Ejemplo 8.5)', 'FontSize', 12);
grid on;
xlabel('Tiempo (min)', 'FontSize', 12);
ylabel('Variación de la temperatura (°C)', 'FontSize', 12);
legend('Respuesta del sistema', 'Entrada rampa', 'FontSize', 10);
set(gca, 'FontSize', 12, 'LineWidth', 2);



% --- Ejemplo 8.6: Respuesta a una rampa unitaria ---
% Sistema de segundo orden: G(s) = 1 / (s^2 + s + 1)

% Para obtener la respuesta a una rampa, dividimos G(s) por s
num = [0 0 0 1];      % Numerador: 1
den = [1 1 1 0];      % Denominador: s^3 + s^2 + s

% Definir el vector de tiempo
t = 0:0.17:7;

% Generar la respuesta a la rampa unitaria
ramp_response = step(num, den, t);

% Graficar la respuesta y la rampa de referencia
figure;
plot(t, ramp_response, '-k', t, t, '-b', 'LineWidth', 2);
title('Respuesta a una rampa unitaria (Ejemplo 8.6)', 'FontSize', 12);
grid on;
xlabel('Tiempo (min)', 'FontSize', 12);
ylabel('Variación de la temperatura (°C)', 'FontSize', 12);
legend('Respuesta del sistema', 'Entrada rampa', 'FontSize', 10);
set(gca, 'FontSize', 12, 'LineWidth', 2);



% --- Ejemplo 8.7: Respuesta transitoria de dos tanques no interactuantes ---
% Tanque 1: tau1 = 0.5, R1 = 1
% Tanque 2: tau2 = 1, R2 = 1
% Función de transferencia total: G(s) = 1 / (0.5s + 1)(s + 1)

% Definir numerador y denominador para el sistema completo
num1 = [0 0 1];        % Numerador: 1
den1 = [0.5 1.5 1];    % Denominador: 0.5s^2 + 1.5s + 1

% Definir numerador y denominador para el tanque 2 (comparación)
num2 = [0 1];          % Numerador: 1
den2 = [1 1];          % Denominador: s + 1

% Definir el vector de tiempo
t = 0:0.1:5;

% Generar las respuestas al escalón unitario
[y1, ~, t] = step(num1, den1, t);
[y2, ~, t] = step(num2, den2, t);

% Graficar las respuestas
figure;
plot(t, y1, '-k', t, y2, '-b', 'LineWidth', 2);
title('Respuesta transitoria de nivel (Ejemplo 8.7)', 'FontSize', 12);
grid on;
xlabel('Tiempo (seg)', 'FontSize', 12);
ylabel('H2(t)', 'FontSize', 12);
legend('Dos tanques no interactuantes', 'Tanque 2 solo', 'FontSize', 10);
set(gca, 'FontSize', 12, 'LineWidth', 2);



% --- Ejemplo 8.7 Alternativo: Tanques interactuantes vs. no interactuantes ---
% Sistema no interactuante: G(s) = 1 / (0.5s + 1)(s + 1)
% Sistema interactuante: G(s) = 1 / (s^2 + 3s + 1)

% Definir numerador y denominador para el sistema no interactuante
num1 = [0 0 1];        % Numerador: 1
den1 = [0.5 1.5 1];    % Denominador: 0.5s^2 + 1.5s + 1

% Definir numerador y denominador para el sistema interactuante
num2 = [0 0 1];        % Numerador: 1
den2 = [1 3 1];        % Denominador: s^2 + 3s + 1

% Definir el vector de tiempo
t = 0:0.1:5;

% Generar las respuestas al escalón unitario
[y1, ~, t] = step(num1, den1, t);
[y2, ~, t] = step(num2, den2, t);

% Graficar las respuestas
figure;
plot(t, y1, '-k', t, y2, '-b', 'LineWidth', 2);
title('Respuesta transitoria de nivel (Ejemplo 8.7 Alternativo)', 'FontSize', 12);
grid on;
xlabel('Tiempo (seg)', 'FontSize', 12);
ylabel('H2(t)', 'FontSize', 12);
legend('No interactuantes', 'Interactuantes', 'FontSize', 10);
set(gca, 'FontSize', 12, 'LineWidth', 2);
