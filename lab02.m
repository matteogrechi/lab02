%% Dati del problema
% Le unità di misura utilizzate sono quelle del SI
% Come da figura 2 del pdf, coordinate dei punti e bracci
xQ = -.280;
yQ = .430;
OP = .150;
OG = .660;

% Range accelerazioni tangenziali e normali del velivolo
atmax = 9.8067;
atmin = -9.8067;
anmax = 9.8067;
anmin = -9.8067;

% Vel max/min estrazione carrello
vmax = 95.1722222;
vmin = 0;

% Dati carrello
S0 = 0.6;   % Superficie frontale
m = 280;    % Massa
Cd = 1.2;   % Coeff attrito aerod

% Dati impianto idraulico
P = 21e+6;  % Alta pressione
p = .2e+6;  % Bassa pressione

% Dati accumulatore
V = 1.8e-3;     % Volume
pPre = 14.2e+6; % Pressione precarica
TPre = 288.15;  % Temp precarica

% Altre costanti utilizzate
rho = 1.225;     % Densità aria condiz std
R = 8.314462618; % Costante universale dei gas
g = 9.80665;     % Accelerazione di gravità
gamma = 7/5;     % Indice caratteristico della politropica (gas azoto)

%% 2.2.1 Lunghezza a riposo e corsa
lR = sqrt((yQ-OP)^2+xQ^2);
lM = sqrt(yQ^2+(-xQ+OP)^2);
c = lM - lR;    % Corsa

%% 2.2.2 Sezioni
% Calcolo del carico massimo in estrazione e retrazione
Fa = @(theta,v) S0*1/2*rho*Cd*v.^2.*sin(theta);
bFa = @(theta) OG*sin(theta);
bFp = @(theta,alpha) OG*abs(cos(theta-alpha));
mRetta = @(theta) (OP*cos(theta)-yQ)./(OP*sin(theta)-xQ);
bF = @(m) abs(m*xQ-yQ)./sqrt(m.^2+1);
 
G = @(theta,alpha) ( Fa(theta,vmax).*bFa(theta) - m*g*bFp(theta,alpha) + m*OG*(atmax*sin(theta)-anmin*cos(theta)) )./bF(mRetta(theta));
H = @(theta,alpha) ( Fa(theta,vmin).*bFa(theta) - m*g*bFp(theta,alpha) + m*OG*(atmin*sin(theta)-anmax*cos(theta)) )./bF(mRetta(theta));

[theta,alpha] = meshgrid([0:0.001:pi/2],[-pi/6:0.01:pi/6]);
Fmax = max(max(G(theta,alpha)));
Fmin = min(min(H(theta,alpha)));

% Dimensionamento attuatore
Press = [P-p p; p-P P];
F = [Fmax; Fmin];
A = Press\F;    % Sezione camera (1) e stelo (2)

%% 2.3 Sovradimensionamento e condizioni di emergenza
Ac_ = A(1)*1.6;   % E' Ac* del pdf

n = pPre*V/(R*TPre);

Tmax = P/(n*R)*(V-c*Ac_);   % Temp max di funzionamento
pMinRichiesta = (Fmax+p*(Ac_-A(2)))/Ac_;
Tmin = (V*P)/(n*R)*(pMinRichiesta/P)^(1/gamma); % Temp min di funzionamento
if Tmax<Tmin 
    fprintf("ERRORE: l'accumulatore standard non funziona\n");
end

%% 2.3* Sovradimensionamento e condizioni di emergenza alterate
V = 5e-3;
Ac_ = A(1)*2.5;   % E' Ac* del pdf

n = pPre*V/(R*TPre);

Tmax = P/(n*R)*(V-c*Ac_);   % Temp max di funzionamento
pMinRichiesta = (Fmax+p*(Ac_-A(2)))/Ac_;
Tmin = (V*P)/(n*R)*(pMinRichiesta/P)^(1/gamma); % Temp min di funzionamento
if Tmax<Tmin 
    fprintf("ERRORE: l'accumulatore modificato non funziona\n");
end