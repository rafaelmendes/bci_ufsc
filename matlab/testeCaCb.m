%% Teste com base nos autovalores medios

clc;clear;close all;

subj = '7';

filenameT = ['/files/eeg_data/compIV_dsIIa/dadosLotte_Normalizados_A0' subj 'T' '.mat'];
filenameE = ['/files/eeg_data/compIV_dsIIa/dadosLotte_Normalizados_A0' subj 'E' '.mat'];

% mi = 5;
lamA = 4;
lamB = 4;
% lam = .5;
% alfa = 100000;

load(filenameT)

clear Xo incremento fname t 
%% Canais de EEG, exclui 3 canais de EOG

Channels = 1:22;
F = F(Channels,:,:);
Fab = F(:,:,001:144);
Fcd = F(:,:,145:288);
[AB,indA,indB] = processaNaN(Fab);
[CD,indC,indD] = processaNaN(Fcd);
clear F
F = cat(3,AB,CD); % LH vs RH
% F = cat(3,CD,AB); % FE vs TO

clear AB CD Fab Fcd Xo  Channels

Na = length(indA);
Nb = length(indB);

% Na = length(indC);
% Nb = length(indD);



%% Seleciona a faixa de frequencia de interesse
[n,m,trials] = size(F);

faixa = [8 30];

SubBandL = (.5*m-1)/(banda(2))*faixa(1) + 1;
SubBandH = (.5*m-1)/(banda(2))*faixa(2) + 1;

faixaFreq = [SubBandL:SubBandH SubBandL+.5*m:SubBandH+.5*m];

q = length(faixaFreq);

clear SubBandL SubBandH


%% Determina o valor medio das COVs

C = zeros(n,n,Na+Nb);
for i=1:Na+Nb,
  aux = F(:,faixaFreq,i)*F(:,faixaFreq,i)';
%   C(:,:,i) = aux/trace(aux);
    C(:,:,i) = aux;
end

%% Projeto CSP
cA = mean(C(:,:,1:Na),3);
cB = mean(C(:,:,Na+1:Na+Nb),3);

%%

Q = cA+cB;
Sa = Q^(-.5)*cA*Q^(-.5);
Sb = Q^(-.5)*cB*Q^(-.5);

for i=1:Na, LamMaxA(i) = max(eig((Sa-Q^(-.5)*C(:,:,i)*Q^(-.5))^2));end
for i=1:Nb, LamMaxB(i) = max(eig((Sb-Q^(-.5)*C(:,:,i+Na)*Q^(-.5))^2));end

% for i=1:Na, LamMaxA(i) = sum(eig((Sa-Q^(-.5)*C(:,:,i)*Q^(-.5))^2));end
% for i=1:Nb, LamMaxB(i) = sum(eig((Sb-Q^(-.5)*C(:,:,i+Na)*Q^(-.5))^2));end

% Autovalor Maximo do Erro Normalizado ao Quadrado.
for i=1:Na, LamMaxAB(i) = max(eig((Sb-Q^(-.5)*C(:,:,i)*Q^(-.5))^2));end
for i=1:Nb, LamMaxBA(i) = max(eig((Sa-Q^(-.5)*C(:,:,i+Na)*Q^(-.5))^2));end

% for i=1:Na, LamMaxAB(i) = sum(eig((Sb-Q^(-.5)*C(:,:,i)*Q^(-.5))^2));end
% for i=1:Nb, LamMaxBA(i) = sum(eig((Sa-Q^(-.5)*C(:,:,i+Na)*Q^(-.5))^2));end

indexB = find(LamMaxB-LamMaxBA>0);
indexA = find(LamMaxA-LamMaxAB>0);

ComplementB = find(LamMaxB-LamMaxBA<0);
ComplementA = find(LamMaxA-LamMaxAB<0);


clear Q Sa Sb  LamMaxAB LamMaxBA

[U,d] = eig(cB,cB+cA);
Wc = U(:,[1 2 3 20 21 22]);

%% feature vectors
for i=1:Na+Nb,
    featCSP(:,i) = diag(Wc'*C(:,:,i)*Wc);
    featCSP(:,i) = log(featCSP(:,i)/sum(featCSP(:,i)));
end

%% Projeto LDA
xA = [featCSP(:,1:Na);ones(1,Na)];
xB = [featCSP(:,Na+1:Na+Nb);ones(1,Nb)];

mA = mean(xA,2);
mB = mean(xB,2);

Sa = xA*xA' - mA*mA';
Sb = xB*xB' - mB*mB';

w = pinv(Sa+Sb)*(mA - mB);
b = .5*w'*(mA+mB);

yT = w'*[xA xB] - b;

clear mA mB Sa Sb 
%% Evaluation

eA = sum(yT(1:Na)>0)/Na;

eB = sum(yT(Na+1:Na+Nb)<0)/Nb;


%%
H = zeros(q,q);
dp = sdpvar(q,1);
H = H + diag(dp);
clear dp

%%
Ca={};
Cb={};

for i=1:Na,
    Ca{i} = F(:,faixaFreq,i)*H*F(:,faixaFreq,i)'; % LH
end
for i =1:Nb,
    Cb{i} = F(:,faixaFreq,Na+i)*H*F(:,faixaFreq,Na+i)'; % RH
end
clear i
C_barA = Ca{1};
C_barB = Cb{1};
for i=2:Na,
    C_barA = C_barA + Ca{i};
end
clear i
C_barA = C_barA/Na;


for i=2:Nb,
    C_barB = C_barB + Cb{i};
end
clear i
C_barB = C_barB/Nb;


%%

lmi1 = [];
lmi2 = [];

for i = 1:Na,
    lmi1 = lmi1 + [[lamA*(C_barA+C_barB) (Ca{i}-C_barA)' (Ca{i} - C_barB)';(Ca{i}-C_barA) (C_barA+C_barB) zeros(n,n);(Ca{i} - C_barB) zeros(n,n) (C_barA+C_barB)]>=0];
end

for i = 1:Nb,
    lmi2 = lmi2 + [[lamB*(C_barA+C_barB) (Cb{i}-C_barA)' (Cb{i} - C_barB)';(Cb{i}-C_barA) (C_barA+C_barB) zeros(n,n);(Cb{i} - C_barB) zeros(n,n) (C_barA+C_barB)]>=0];
end


solvesdp([H>=0,lmi1,lmi2],trace(C_barB - C_barA),sdpsettings('solver','sedumi'))

clc
clear LamMaxA LamMaxB
%%
disp('feasible (+):')

%%

H = double(H);



h = diag(H);
f = linspace(faixa(1),faixa(2),.5*q);

% %% Projeto CSP com Ho otimo

% C_barA = double(C_barA);
% C_barB = double(C_barB);

% Q = C_barA+C_barB;
% Sa = Q^(-.5)*C_barA*Q^(-.5);
% Sb = Q^(-.5)*C_barB*Q^(-.5);

% for i=1:Na, LamMaxA(i) = max(eig((Sa-Q^(-.5)*double(Ca{i})*Q^(-.5))^2));end
% for i=1:Nb, LamMaxB(i) = max(eig((Sb-Q^(-.5)*double(Cb{i})*Q^(-.5))^2));end


% for i=1:Na, LamMaxAB(i) = max(eig((Sb-Q^(-.5)*double(Ca{i})*Q^(-.5))^2));end
% for i=1:Nb, LamMaxBA(i) = max(eig((Sa-Q^(-.5)*double(Cb{i})*Q^(-.5))^2));end

% clear Q Sa Sb LamMaxA LamMaxB LamMaxAB LamMaxBA

% %%

% [V,D] = eig(C_barB,C_barA+C_barB);
% [~,I] = sort(diag(D));
% W = V(:,I([1 2 3 20 21 22]));

% for i=1:Na+Nb,
%     featOtim(:,i) = diag(W'*F(:,faixaFreq,i)*H*F(:,faixaFreq,i)'*W);
%     featOtim(:,i) = log(featOtim(:,i)/sum(featOtim(:,i)));
% end

% %% Projeto LDA com Ho otim
% xxA = [featOtim(:,1:Na);ones(1,Na)];
% xxB = [featOtim(:,Na+1:Na+Nb);ones(1,Nb)];

% mA = mean(xxA,2);
% mB = mean(xxB,2);

% Sa = xxA*xxA' - mA*mA';
% Sb = xxB*xxB' - mB*mB';

% wo = pinv(Sa+Sb)*(mA - mB);
% bo = .5*w'*(mA+mB);

% yTo = wo'*[xxA xxB] - bo;

% clear mA mB Sa Sb 

% %% Evaluation Trainning Set

% eAo = sum(yTo(1:Na)>0)/Na;

% eBo = sum(yTo(Na+1:Na+Nb)<0)/Nb;

% %% Carrega dados de Validacao
% clear F Fab Fcd AB indA indB CD indC indD Na Nb
% load(filenameE)
% clear Xo incremento fname labels

% %% Canais de EEG, exclui 3 canais de EOG
% Channels = 1:22;
% F = F(Channels,:,:);
% Fab = F(:,:,001:144);
% Fcd = F(:,:,145:288);
% [AB,indA,indB] = processaNaN(Fab);
% [CD,indC,indD] = processaNaN(Fcd);
% clear F
% F = cat(3,AB,CD); % LH vs RH
% % F = cat(3,CD,AB); % FE vs TO
% clear AB CD Fab Fcd Xo  Channels

% Na = length(indA);
% Nb = length(indB);


% %% Normalizacao do F

% % F(:,:,86) = F(:,:,87); Sub01T error trial 86, channel 18
% %% Gera matrizes de COV para Ho ident
% Cy = zeros(n,n,Na+Nb);
% for i=1:Na+Nb,
%   aux = F(:,faixaFreq,i)*F(:,faixaFreq,i)';
% %   Cy(:,:,i) = aux/trace(aux);
% Cy(:,:,i) = aux;
% end
% clear aux
% %% Vetro de Caracteristica para Ho ident
% for i=1:Na+Nb,
%     featY(:,i) = diag(Wc'*Cy(:,:,i)*Wc);
%     featY(:,i) = log(featY(:,i)/sum(featY(:,i)));
% end

% yE = w'*[featY;ones(1,Na+Nb)]-b;

% %% Evaluation Validation Set

% eAe = sum(yE(1:Na)>0)/Na;

% eBe = sum(yE(Na+1:Na+Nb)<0)/Nb;


% Cyo = zeros(n,n,Na+Nb);
% for i=1:Na+Nb,
%   auxO = F(:,faixaFreq,i)*H*F(:,faixaFreq,i)';
% %   Cy(:,:,i) = aux/trace(aux);
% Cyo(:,:,i) = auxO;
% end
% clear auxO
% %%
% for i=1:Na+Nb,
%     featYo(:,i) = diag(W'*Cyo(:,:,i)*W); 
%     featYo(:,i) = log(featYo(:,i)/sum(featYo(:,i)));
% end

% yEo = wo'*[featYo;ones(1,Na+Nb)]-bo;

% %% Evaluation Validation Set

% eAeO = sum(yEo(1:Na)>0)/Na;
% eBeO = sum(yEo(Na+1:Na+Nb)<0)/Nb;


% %% 
% table = zeros(2,3);
% table(1,1) = .5*(eA+eB); % H = I, conjunto T
% table(1,2) = .5*(eAe+eBe); % H = I, conjunto E
% table(1,3) = .5*(.5*(eAe+eBe)+ .5*(eA+eB)); % H = I, conjunto T+E

% table(2,1) = .5*(eAo+eBo); % H = D, conjunto T
% table(2,2) = .5*(eAeO+eBeO); % H = D, conjunto E
% table(2,3) = .5*(.5*(eAeO+eBeO)+.5*(eAo+eBo)); % H = D, conjunto T+E

% %%
% metric_mean = {'T','E','T+E'};
% subj = {'Ho = I','Ho = D'};
% disp('----------------------------------------------------');
% disp('Accuracy (%) - Rows : H, Colums : Conjuntos');
% disp('----------------------------------------------------');
% displaytable(table,metric_mean,10,{'.4f'},subj)
% disp('----------------------------------------------------');