%% Constantes/datos
f_IF1 = 10.7 * 10^6; % 10.7 MHz, portadora inicial.
f_IF2 = 455 * 10^3; % 455 kHz, portadora tras conversión.
f_s = 8 * f_IF1; % = 8 * 10.7 MHz, frec de muestreo
f_m150 = 150; % 150 Hz, tono "1" de navegación (+ en der).
f_m90 = 90; % 90 Hz, tono "2" de navegación (+ en izq). 
A_c = 10; % 10 V, amplitud de la portadora
m = 0.2; % adim, índice de modulación.
n = 1:12000000; % adim, número de muestras (vector).
Ess = 1; % V podría generar también números aleatorios entre -1 y 1.

% Son colores legibles, fáciles de diferenciar para daltónicos
% https://personal.sron.nl/~pault/#sec:qualitative
color_base = '#CCBB44';
color_a = '#4477AA';
color_b = '#EE6677';


%% Funciones auxiliares 
function x_out = rectificador(x_in)
% Rectifica ondas
% Input: 
%       x_in: señal cualquiera.
% Output:
%       x_out: señal rectificada.
    x_out = (x_in > 0) .* x_in;
end

function x_out = bloqueadordc(x_in)
% Bloquea la componente continua en ondas
% Input: 
%       x_in: señal cualquiera.
% Output:
%       x_out: señal sin componente continua.
    x_out = x_in - mean(x_in);
end

function x_out = multiplicador(x_in1, x_in2)
% Multiplica dos ondas dadas
% Input: 
%       x_in1: señal cualquiera.
%       x_in2: otra señal cualquiera.
% Output:
%       x_out: señal obtenida como producto de x_in y x_OL
    x_out = x_in1 .* x_in2;
end

function x_out = mezclador(x_in, x_OL)
% Multiplica dos señales, luego filtra la salida.
% Input: 
%       x_in: señal cualquiera.
%       x_OL: otra señal cualquiera (oscilador local en ILS).
% Output:
%       x_out: señal obtenida tras el mezclado de x_in con x_ol
    x_aux = multiplicador(x_in, x_OL);
    [b, a] = butter(2, 0.04);
    x_out = filter(b, a, x_aux);
    x_out = x_out * max(x_aux)/max(x_out);
end

function x_out = detector(x_in, w_c)
% Rectifica, filtra y bloquea la componente continua de una señal
% Input: 
%       x_in: señal cualquiera.
%       w_c: pulsación de corte.
% Output:
%       x_out: señal de salida.
    x_rin = rectificador(x_in);
    [b, a] = butter(1, w_c); 
    x_out = filter(b, a, x_rin);
    x_out = x_out * max(x_rin)/max(x_out);
    x_out = bloqueadordc(x_out);
end

%% Variables y señales básicas
% Señales de navegación
w_d90 = 2*pi * f_m90/f_s; % pulsación señal 90 Hz
w_d150 = 2*pi * f_m150/f_s;  % pulsación señal 159 Hz
w_dp = 2*pi * f_IF1/f_s; %  % pulsación portadora
x_90 = sin(n * w_d90 ); % señal de navegación 90 Hz
x_150 = sin(n * w_d150); % señal de navegación 150 Hz

portadora = cos(n * w_dp); % portadora

% Señales mezclador
w_dOL = 2*pi * (f_IF1 + f_IF2)/f_s; % pulsación dek oscilador local
x_OL = cos(n * w_dOL); % oscilador local

% para la representación en frecuencia
f_muestreo = (f_s/length(n)) * (0:(length(n)-1));

% Para el detector
w_c = 2*pi * f_m150 / (f_s*pi); % pulsación de corte

%% Generación de las señales (BL, PBL)
PBL = A_c .* portadora .* (1 + m * (x_90 + x_150)); % señal PBL con modulación AM con m=0.2
BL = Ess .*  portadora .* (x_90 - x_150); % Señal BL
completa = PBL+BL;

% Figura de la señal PBL
figure
hold on
grid on
title('Señal PBL')
ylabel('Amplitud [V]')
xlabel('Muestras [adim]')
plot(n, PBL, 'color', color_base)
hold off
fig = gcf;
exportgraphics(fig, 'pbl.pdf', 'ContentType', 'vector'); 

% Figura de la señal PBL
figure
hold on
grid on
title('Señal BL')
ylabel('Amplitud [V]')
xlabel('Muestras [adim]')
plot(n, BL, 'color', color_base)
hold off
fig = gcf;
exportgraphics(fig, 'bl.pdf', 'ContentType', 'vector'); 

% Figura de la señal PBL+BL
figure
hold on
grid on
title('Señal PBL + BL')
ylabel('Amplitud [V]')
xlabel('Muestras [adim]')
plot(n, completa, 'color', color_a) 
plot(n, PBL,'color', color_b) 
legend('PBL+BL', 'PBL')
hold off
fig = gcf;
exportgraphics(fig, 'pbl+bl.pdf', 'ContentType', 'vector'); 

%% Mezclador 
% https://es.mathworks.com/help/signal/ref/butter.html
multiplicada = multiplicador(completa, x_OL); % esta no haría falta, ya que existe mezclador(), pero como me piden la gráfica, la uso.

amplificada = mezclador(completa, x_OL);

% Figura salida multiplicador
figure
hold on
grid on
title('Señal a la salida del multiplicador')
ylabel('Amplitud [V]')
xlabel('Muestras [adim]')
plot(n, multiplicada, 'color', color_base)
hold off
fig = gcf;
exportgraphics(fig, 'multiplicador.pdf', 'ContentType', 'vector'); 

% Figura salida mezclador
figure
hold on
grid on
title('Señal a la salida del mezclador')
ylabel('Amplitud [V]')
xlabel('Muestras [adim]')
plot(n, amplificada, 'color', color_base)
hold off
fig = gcf;
exportgraphics(fig, 'mezclador.pdf', 'ContentType', 'vector'); 

%% Representación en frecuencia
f_completa = abs(fft(completa, length(n)));
f_multiplicada = abs(fft(multiplicada, length(n)));
f_amplificada = abs(fft(amplificada, length(n)));

figure
hold on
grid on
title('Señal PBL + BL (en frecuencia)')
ylabel('Potencia [W]')
xlabel('Frecuencia [Hz]')
plot(f_muestreo, f_completa, 'color', color_base)
maxX = f_muestreo(find(f_completa == max(f_completa), 1, "first"));
xlim([maxX-250 maxX+250])
hold off  % https://es.mathworks.com/matlabcentral/answers/72396-x-value-on-y-max#answer_82526
fig = gcf;
exportgraphics(fig, 'pbl+blhz.pdf', 'ContentType', 'vector'); 

figure
hold on
grid on
title('Señal a la salida del multiplicador (en frecuencia)')
ylabel('Potencia [W]')
xlabel('Frecuencia [Hz]')
plot(f_muestreo, f_multiplicada, 'color', color_base)
maxX = f_muestreo(find(f_multiplicada == max(f_multiplicada), 1, "first"));
xlim([maxX-250 maxX+250])
hold off  % https://es.mathworks.com/matlabcentral/answers/72396-x-value-on-y-max#answer_82526
fig = gcf;
exportgraphics(fig, 'multiplicadorhz.pdf', 'ContentType', 'vector'); 

figure
hold on
grid on
title('Señal a la salida del mezclador (en frecuencia)')
ylabel('Potencia [W]')
xlabel('Frecuencia [Hz]')
plot(f_muestreo, f_amplificada, 'color', color_base)
maxX = f_muestreo(find(f_amplificada == max(f_amplificada), 1, "first"));
xlim([maxX-250 maxX+250])
hold off  % https://es.mathworks.com/matlabcentral/answers/72396-x-value-on-y-max#answer_82526
fig = gcf;
exportgraphics(fig, 'muestreo.pdf', 'ContentType', 'vector'); 

%% Detector
bloqueodc = detector(amplificada, w_c);

% Representación a la salida de le envolvente
figure
hold on
grid on
title('Señal a la salida del mezclador')
ylabel('Amplitud [V]')
xlabel('Muestras [admin]')
plot(n, bloqueodc, 'color', color_base)
hold off 
fig = gcf;
exportgraphics(fig, 'detector.pdf', 'ContentType', 'vector'); 

% diezmado, https://es.mathworks.com/help/signal/ref/downsample.html
diezmado = downsample(bloqueodc, 2400); % señal diezmada
n_2 = downsample(n, 2400); % muestras diezmadas
f_s2 = f_s/2400; % frecuencia diezmada
f_muestreo2 = (f_s2/length(n_2)) * (0:(length(n_2)-1));

% representación en frecuencia ANTES diezmado
f_antes = abs(fft(bloqueodc));
x = f_muestreo./(pi*length(n));
figure
hold on
grid on
title('Envolvente antes del diezmado')
ylabel('Potencia [W]')
xlabel('Frecuencia [Hz]')
plot(x, f_antes, 'color', color_base)
maxX = x(find(f_antes == max(f_antes), 1, "first"));
xlim([maxX-3*10^-6 maxX+3*10^-6])
hold off % se ve entre poco y nada, arreglar
fig = gcf;
exportgraphics(fig, 'envolvente_antes.pdf', 'ContentType', 'vector'); 

% representación en frecuencia DESPUÉS diezmado
f_despues = abs(fft(diezmado));
x = f_muestreo2./(pi*length(n_2));
figure
hold on
grid on
title('Envolvente después del diezmado')
ylabel('Potencia [W]')
xlabel('Frecuencia [Hz]')
plot(x, f_despues, 'color', color_base);
maxX = x(find(f_despues == max(f_despues), 1, "first"));
xlim([maxX-4*10^-3 maxX+4*10^-3])
hold off
fig = gcf;
exportgraphics(fig, 'envolvente_despues.pdf', 'ContentType', 'vector'); 

%% Filtros paso banda señales navegación
% https://es.mathworks.com/help/signal/ref/ellip.html
w_c90i = 2*pi*70/f_s2; % plsación corte inf 90Hz
w_c90s = 2*pi*100/f_s2; % plsación corte sup 90Hz
w_c150i = 2*pi*140/f_s2; % plsación corte inf 150Hz
w_c150s = 2*pi*170/f_s2; % plsación corte sup 150Hz

w_c90 = [w_c90i/pi, w_c90s/pi];
w_c150 = [w_c150i/pi, w_c150s/pi];

[b, a] = ellip(3, 0.2, 40, w_c90); % 90 Hz, Se sobreescriben coefs, no importan
f_90 =  filter(b, a, diezmado); % OJO con el orden de los coefs!!!

[d, c] = ellip(3, 0.2, 40, w_c150); % 150 Hz
f_150 = filter(d, c, diezmado);% OJO con el orden de los coefs!!!

f_90 = (max(diezmado)/max(f_90)) * f_90;
f_150 = (max(diezmado)/max(f_150)) * f_150;

% Representación de 90Hz y 150Hz
fft_f90 = abs(fft(f_90));
fft_f150 = abs(fft(f_150));

figure
hold on
grid on
title('Señales de navegación')
ylabel('Potencia [W]')
xlabel('Frecuencia [Hz]')
plot(f_muestreo2, fft_f90, 'color', color_a)
plot(f_muestreo2, fft_f150, 'color', color_b)
legend('90 Hz', '150 Hz')
xlim([50 200]) % aquí no hace falta hacer cosas raras
hold off 
fig = gcf;
exportgraphics(fig, 'senales.pdf', 'ContentType', 'vector'); 


% Representación espectro señal + filtros https://es.mathworks.com/help/signal/ref/freqz.html
[h, w] = freqz(b, a, 10000); % OJO con el orden de los coefs!!!
[i, y] = freqz(d, c, 10000);

figure
hold on
grid on
title('Espectros y funciones de transferencia')
ylabel('Potencia Normalizada H(j\omega) [adim]')
xlabel('Frecuencia [Hz]')
plot(f_muestreo2, f_despues/max(f_despues), 'color', color_base)
plot(w .* (f_s2/(2*pi)), abs(h), 'color', color_a)
plot(y .* (f_s2/(2*pi)), abs(i), 'color', color_b) 
legend('Señales de navegación', 'Filtro 90Hz', 'Filtro 150Hz')
xlim([50 200]) % aquí no hace falta hacer cosas raras
hold off
fig = gcf;
exportgraphics(fig, 'espectro.pdf', 'ContentType', 'vector'); 


% Representación envolvente + señales de navegación
figure
hold on
grid on
title('Envolvente y señales de navegación')
ylabel('Amplitud [V]')
xlabel('Muestras')
plot(diezmado(1,1000:end), 'color', color_base) % quito las primeras 1000 muestras porque son inestables
plot(f_90(1,1000:end), 'color', color_a, "LineStyle", "-.")
plot(f_150(1,1000:end), 'color', color_b, "LineStyle", "-.") 
legend('Señal diezmada', '90 Hz', '150Hz')
xlim([1000 4000])
hold off
fig = gcf;
exportgraphics(fig, 'envolvente+senales.pdf', 'ContentType', 'vector'); 


%% DDM
E_c = bloqueodc;
DDM = (max(f_90(1,3000:end)) - max(f_150(1,3000:end)))/(max(E_c)) % a partir de ~3000 ya se estabiliza