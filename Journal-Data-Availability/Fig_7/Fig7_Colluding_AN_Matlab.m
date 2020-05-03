
% System Model
% Alice has one antennas directed towords BOB which has one antenna (Main lobe towords main lobe) with fixed distance r_RT 
% Interference for Bob: receives from other interfers(Is)(one antenna with PPP distribution)

% Colluding Eves: (one antenna with PPP distribution)the sum of the received signals from all Eves 
% They can cancel the interference from other Is

clear all;
close all;
clc;

tic
A = 10^6;                                 % Area of the Network (m^2)
r_max =sqrt(A/pi) ;                       % Radius of the area of the network (m)

lambda_E = (50:50:500)*10^-6;                      % Intensity of Eave per (m^2)
lambda = 50*10^-6;                        % Intensity of Tx per (m^2)
r_RT=15;                                  % Distance Alice-Bob

% Alice number of the antennas (Main-side lobe gains (First antenna (a1)) ) 
NS=16;                                     % #of TX antenna equal to Rx antenna
Theta_S=(2*pi)/sqrt(NS);                     % Angle of the main lobe sector for Tx & Rx  
GMS=NS;                                   % Main-lobe gain for Tx & Rx (dB)
GmS=1/(sin((3*pi)/(2*sqrt(NS))))^2;        % Side-lobe gain for Tx & Rx (dB)

% Alice number of the antennas (Main-side lobe gains (Second antenna (a2)) )
NA=16;                                     % #of TX antenna equal to Rx antenna
Theta_A=(2*pi)/sqrt(NA);                     % Angle of the main lobe sector for Tx & Rx  
GMA=NA;                                    % Main-lobe gain for Tx & Rx (dB)
GmA=1/(sin((3*pi)/(2*sqrt(NA))))^2;        % Side-lobe gain for Tx & Rx (dB)

% Bob number of the antennas (Main-side lobe gains ) 
NR=16;                                     % #of Alice antenna equal to Bob antenna
Theta_R=(2*pi)/sqrt(NR);                     % Angle of the main lobe sector for Tx & Rx  
GMR=NR;                                    % Main-lobe gain for Tx & Rx (dB)
GmR=1/(sin((3*pi)/(2*sqrt(NR))))^2;        % Side-lobe gain for Tx & Rx (dB)

% I number of the antennas (Main-side lobe gains) 
NI=16;                                     % #of TX antenna equal to Rx antenna
Theta_I=(2*pi)/sqrt(NI);                     % Angle of the main lobe sector for Tx & Rx  
GMI=NI;                                   % Main-lobe gain for Tx & Rx (dB)
GmI=1/(sin((3*pi)/(2*sqrt(NI))))^2;        % Side-lobe gain for Tx & Rx (dB)

% Eve number of the antennas (Main-side lobe gains)
N_E=16;                                        % #of Eav antenna 
Theta_E=(2*pi)/sqrt(N_E);                      % Angle of the main lobe sector for Eav
GM_E=N_E;                                     % Main-lobe gain for Eav (dB)
Gm_E=1/(sin((3*pi)/(2*sqrt(N_E))))^2;         % Side-lobe gain for Eav (dB)

% law of total expectation Is-Bob
pr_MM_IR=(Theta_R*Theta_I)/(2*pi)^2;              % Prob of GM*GM_E 
pr_Mm_IR=(Theta_R*(2*pi-Theta_I))/(2*pi)^2;       % Prob of GM*Gm_E 
pr_mM_IR=(Theta_I*(2*pi-Theta_R))/(2*pi)^2;       % Prob of Gm*GM_E 
pr_mm_IR=((2*pi-Theta_R)*(2*pi-Theta_I))/(2*pi)^2;  % Prob of Gm*Gm_E 

% law of total expectation Alice(a1)-Eves
pr_MM_ES=(Theta_S*Theta_E)/(2*pi)^2;              % Prob of GM*GM_E 
pr_Mm_ES=(Theta_S*(2*pi-Theta_E))/(2*pi)^2;       % Prob of GM*Gm_E 
pr_mM_ES=(Theta_E*(2*pi-Theta_S))/(2*pi)^2;       % Prob of Gm*GM_E 
pr_mm_ES=((2*pi-Theta_S)*(2*pi-Theta_E))/(2*pi)^2;  % Prob of Gm*Gm_E 

% law of total expectation Alice(a2)-Eves
pr_MM_EA=(Theta_A*Theta_E)/(2*pi)^2;              % Prob of GM*GM_E 
pr_Mm_EA=(Theta_A*(2*pi-Theta_E))/(2*pi)^2;       % Prob of GM*Gm_E 
pr_mM_EA=(Theta_E*(2*pi-Theta_A))/(2*pi)^2;       % Prob of Gm*GM_E 
pr_mm_EA=((2*pi-Theta_S)*(2*pi-Theta_A))/(2*pi)^2;  % Prob of Gm*Gm_E 

% Generate the average antenna gain between Tx & Rx

% law of total expectation Is-Bob
case1=pr_MM_IR;
case2=pr_MM_IR+pr_Mm_IR;
case3=pr_MM_IR+pr_Mm_IR+pr_mM_IR;
case4=pr_MM_IR+pr_Mm_IR+pr_mM_IR+pr_mm_IR;

% law of total expectation Alice(a1)-Eves
pr_M_E=Theta_E/(2*pi);
pr_m_E=1-pr_M_E;

% law of total expectation Alice(a1)-Eves
PrMmM = (Theta_S*Theta_E )/(2*pi)^2;
PrMmm = (Theta_S*(2*pi - Theta_E))/(2*pi)^2;
PrmMM = (Theta_A*Theta_E )/(2*pi)^2;
PrmMm = (Theta_A*(2*pi - Theta_E))/(2*pi)^2;
PrmmM = ((2*pi - Theta_S - Theta_A)*Theta_E )/(2*pi)^2;
Prmmm = ((2*pi - Theta_S - Theta_A)*(2*pi - Theta_E) )/(2*pi)^2;




case1E=PrMmM;
case2E=PrMmM+PrMmm;
case3E=PrMmM+PrMmm+PrmMM;
case4E=PrMmM+PrMmm+PrmMM+PrmMm;
case5E=PrMmM+PrMmm+PrmMM+PrmMm+PrmmM;
case6E=PrMmM+PrMmm+PrmMM+PrmMm+PrmmM+Prmmm;





Q=1/141.4;                                % (1/m)
d=1;                                      % reference distance (m)
C=3*10^8;                                 % Velocity of the light(m/s)
Fc=73*10^9;                               % Carrier frequency (Hz)
Beta=(C/(4*pi*Fc))^2;                     % freq indep const parameter of path loss
alpha_LoS=2;                              % path loss exponent LoS for 28 GHz
alpha_NLoS=3.4;                             % path loss exponent NLoS for 28 GHz

NF=10;                                    % Noise figure (dB)
BW=2*10^9;                                % mmWave bandwidth (Hz)
sgmasqr_dBm=-174+10*log10(BW)+NF;         % Noise power for Tx (dBm)
sgmasqr=(10^(sgmasqr_dBm/10))/1000;       % Noise power for Tx (Watt)
sgmasqr_E=sgmasqr; 

Pt1=30;                % Tx power (dBm)
Pt=(10.^(Pt1/10))/1000;                     % Tx power (Watt)

Nt=3;
Nta=3;


Pt=Pt/Nt;
Mu=0.75;

PS=Mu*Pt;
PA=Pt-PS;

NL=3;
NN=2;
% S=1;

SS1=[];
SS1_E=[];
itt=10^4;
%%
for h=1:length(lambda_E)
R1=[];
R1_E=[];
for p=1:itt
% Generate the Poisson number of nodes for Is
num_of_Tx =poissrnd(lambda*A);
All_r_Txx=r_max*sqrt(rand(num_of_Tx,1));     % The poler distance between each Tx & Origin
ang=2*pi*rand(num_of_Tx,1);                % The poler angle between each Tx & Origin
Loc_Tx=[All_r_Txx.*cos(ang),All_r_Txx.*sin(ang)];  % The postion of Tx in Cartesian coordinates

% % Get the position of the Bob which is 15m from Alice (at origin)
% Loc_TTx=[0,0];                                % Location of Alice (at origin)
% ang_RRx=2*pi*rand(1,1);
% Loc_RRx=[r_RT*cos(ang_RRx),r_RT*sin(ang_RRx)];  % Location of Bob (15 m from Alice)
% 
% % Distance between Bob and all Is
% Dis_Txs_RRx=Loc_Tx-Loc_RRx;
% All_r_Txs_RRx=sqrt(Dis_Txs_RRx(:,1).^2+Dis_Txs_RRx(:,2).^2);
% 
% % Calculat SINR at Bob
% % Alice-Bob Path loss and gamma coff. for two paths
% 
% if    rand(1)<exp(-Q*r_RT)
%        PL_RRx=Beta*(max(d,r_RT))^-alpha_LoS;     % Path loss of LoS (abs)
%       hh=gamrnd(NL,1/NL,[1,Nt]);
%  else  PL_RRx=Beta*(max(d,r_RT))^-alpha_NLoS;    % Path loss of NLoS (abs)
%     hh=gamrnd(NN,1/NN,[1,Nt]);
% end
% 
%  % Nakagami channel coff. (Alice-Bob)
%  W=(hh.^0.5);
% %  [U,S,V]=svd(W);
% %  Vs=V(:,1);
% %  Vn=V(:,2);
% 
% % Received power at Bob from Alice
% Pr=W*W'*PS*GMS*GMR*PL_RRx;
% 
% % Is-Bob Path loss and gamma coff. for one path
% All_PL=[];
% W_I=zeros(Nt,1);
% for u=1:length(All_r_Txs_RRx)
% % Calculate the path loss between each Tx & Rx 
%  if    rand(1)<exp(-Q*All_r_Txs_RRx(u))
%        PL_Los=Beta*(max(d,All_r_Txs_RRx(u)))^-alpha_LoS;     % Path loss of LoS (abs)
%        hh_i=gamrnd(NL,1/NL,[1,Nt]);
%  else  PL_Los=Beta*(max(d,All_r_Txs_RRx(u)))^-alpha_NLoS;    % Path loss of NLoS (abs)
%        hh_i=gamrnd(NN,1/NN,[1,Nt]);
%  end
%  W_i=(hh_i.^0.5)';
%  
%  All_PL=[All_PL PL_Los];
%   W_I=[W_I W_i];
% end
% 
% % Calculate the interference from Is to Bob
% I=0;
% 
% for k=1:length(All_r_Txs_RRx)
% uvar=rand(1);  
% if uvar<=case1
% Gi=GMI*GMR;
% elseif uvar>=case1&&uvar<=case2
% Gi=GMI*GmR;
% elseif uvar>=case2&&uvar<=case3
%  Gi=GmI*GMR; 
% else
%   Gi=GmI*GmR;
% end
% 
%     g=Pt*Gi*W_I(:,k+1)'*W_I(:,k+1)*All_PL(k);
%     
%     I=[I+g];
% 
% end
% 
% % Calculate the SINR at the Rx
% SINR1=Pr/(I+sgmasqr);
% 
% % Calculate rate of the channel between Tx & Rx
% R=log2(1+SINR1);
% 
% R1=[R1 R];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the Poisson number of nodes for Eves
num_of_E =poissrnd(lambda_E(h)*A);
All_r_E=r_max*sqrt(rand(num_of_E,1));       % The poler distance between each Tx & Eav
ang_E=2*pi*rand(num_of_E,1);                % The poler angle between each Tx & Eav
Loc_E=[All_r_E.*cos(ang_E),All_r_E.*sin(ang_E)]; % The postion of Eav in Cartesian coordinates

% Calculat SINR at Eves
% Alice-Eves Path loss and gamma coff. for two paths
Pr_E=0;
Pr_EA=0;
Prn_E=0;
for i=1:length(All_r_E)
if    rand(1)<exp(-Q*All_r_E(i))
       PL_Los_ES=Beta*(max(d,All_r_E(i)))^-alpha_LoS;     % Path loss of LoS (abs) 
       hh_i=gamrnd(NL,1/NL,[1,Nt]);
        hh_iA=gamrnd(NL,1/NL,[1,Nta]);
else   PL_Los_ES=Beta*(max(d,All_r_E(i)))^-alpha_NLoS;    % Path loss of NLoS (abs)
      hh_i=gamrnd(NN,1/NN,[1,Nt]);
      hh_iA=gamrnd(NN,1/NN,[1,Nta]);
end 

 % Nakagami channel coff. (Alice-Eves)
W_E=hh_i.^0.5;
W_EA=hh_iA.^0.5;

% Calculate the antenna gains between (Alice(a1&a2)-Eves)
uvar1=rand(1);  
if uvar1<=case1E
GeS=GMS*GM_E;
GeA=GmA*GM_E;
elseif uvar1>case1E&&uvar1<=case2E
GeS=GMS*Gm_E;
GeA=GmA*Gm_E;
elseif uvar1>case2E&&uvar1<=case3E
 GeS=GmS*GM_E;
 GeA=GMA*GM_E;
elseif uvar1>case3E&&uvar1<=case4E
  GeS=GmS*Gm_E;
  GeA=GMA*Gm_E;
  elseif uvar1>case4E&&uvar1<=case5E
  GeS=GmS*GM_E;
  GeA=GmA*GM_E;
  elseif uvar1>case5E&&uvar1<=case6E
  GeS=GmS*Gm_E;
  GeA=GmA*Gm_E;

end


Pr_E=Pr_E+W_E*W_E'*PS*GeS*PL_Los_ES;
Pr_EA=Pr_EA+W_EA*W_EA'*PA*GeA*PL_Los_ES;

end


% Calculate the SINR at the Eves
SINR_E=Pr_E/(Pr_EA+sgmasqr_E);

% Calculate rate of the channel between Tx & nearest Eav
R_E=log2(1+SINR_E);

 R1_E=[R1_E R_E];
end
S1=sum(R1)/itt;
SS1=[SS1 S1]

S1_E=sum(R1_E)/itt;
SS1_E=[SS1_E S1_E]
end
S=SS1-SS1_E
figure(1)
plot(Pt1,S)
grid on
Rx_rate=[ 3.8135    5.1816    6.4397    7.4586    8.1880    8.6238    8.8621    8.9759    9.0198    9.0340];
Eve_rate=[  0.1289    0.2907    0.5796    1.0096    1.5595    2.1245    2.5650    2.8797    3.0911    3.2102];
S=[    ];
toc                      