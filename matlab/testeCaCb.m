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

%% Normalizacao do F
% for i=1:trials,
%     F(:,faixaFreq,i) = F(:,faixaFreq,i)/norm(F(:,faixaFreq,i),'fro');
% end


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

% figure(11)
% subplot(212);plot(LamMaxB,'--o')
% hold on
% subplot(211);plot(LamMaxA,'--o')
% title('Eigenvalue Max Error H = I ')
% Autovalor Maximo do Erro Normalizado ao Quadrado.
for i=1:Na, LamMaxAB(i) = max(eig((Sb-Q^(-.5)*C(:,:,i)*Q^(-.5))^2));end
for i=1:Nb, LamMaxBA(i) = max(eig((Sa-Q^(-.5)*C(:,:,i+Na)*Q^(-.5))^2));end

% for i=1:Na, LamMaxAB(i) = sum(eig((Sb-Q^(-.5)*C(:,:,i)*Q^(-.5))^2));end
% for i=1:Nb, LamMaxBA(i) = sum(eig((Sa-Q^(-.5)*C(:,:,i+Na)*Q^(-.5))^2));end


% hold on
% subplot(212);plot(LamMaxBA,'--r*')
% hold on
% subplot(211);plot(LamMaxAB,'--r*')


indexB = find(LamMaxB-LamMaxBA>0);
indexA = find(LamMaxA-LamMaxAB>0);

ComplementB = find(LamMaxB-LamMaxBA<0);
ComplementA = find(LamMaxA-LamMaxAB<0);

% subplot(211);plot(indexA,LamMaxA(indexA),'ks')
% subplot(212);plot(indexB,LamMaxB(indexB),'ks')

clear Q Sa Sb  LamMaxAB LamMaxBA

%%
% figure(10)
% plot(diag(cA),'--o')
% hold on
% plot(diag(cB),'--ro')
% legend('A','B')
%%
[U,d] = eig(cB,cB+cA);
Wc = U(:,[1 2 3 20 21 22]);

% figure(1)
% plot(diag(d),'--o')


%% feature vectors
for i=1:Na+Nb,
    featCSP(:,i) = diag(Wc'*C(:,:,i)*Wc);
    featCSP(:,i) = log(featCSP(:,i)/sum(featCSP(:,i)));
end

% figure(1)
% plot(featCSP','--o')

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

% figure
% plot(yT,'--o')
% grid
clear mA mB Sa Sb 
%% Evaluation

eA = sum(yT(1:Na)>0)/Na;

eB = sum(yT(Na+1:Na+Nb)<0)/Nb;

% disp('CSP+LDA: T')
% .5*(eA+eB)

%% Otimizacao via LMIs
% H = diag(sdpvar(q,1));
% H = sdpvar(q);

%%
H = zeros(q,q);
dp = sdpvar(q,1);
H = H + diag(dp);
clear dp
%%
% for k=2:2,
%    ind = k - 1;
%    ds = sdpvar(q-ind,1);
%    H(k:q,1:q-ind) = H(k:q,1:q-ind) + diag(ds);
%    H(1:q-ind,k:q) = H(1:q-ind,k:q) + diag(ds);
%    clear ds
% end


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

% clear Ca Cb
%% LMI para o Medio

% lmi = [[lam*(C_barA+C_barB)-C_barB]>=0];

% lmi = [[mi*(C_barA+C_barB)-C_barB]>=0];
% lmi = [[mi*(C_barA+C_barB) C_barB';C_barB C_barA+C_barB]>=0];

% lmi = [[lam*C_barB - .5*(C_barA+C_barB)]>=0] + [[lam*.5*(C_barA+C_barB) - C_barA]>=0];

% lmi = [[lam*C_barB C_barA';C_barA C_barB]>=0];

% lmi = [[lam*C_barA (C_barA - C_barB)';(C_barA - C_barB) C_barA]>=0];

% lmi = [[lam*C_barA eye(n);eye(n) C_barB]>=0];

% lmi = [[lam*eye(n) (C_barA-C_barB)';(C_barA-C_barB) eye(n)]>=0];

%% Quasi-convex

% lmi1 = [];
% lmi2 = [];
% 
% for i = 1:length(indexA),
%     
%    lmi1 = lmi1 + [[lamA*(C_barA+C_barB) (Ca{indexA(i)}-C_barA)';Ca{indexA(i)}-C_barA (C_barA+C_barB)]>=0];
% end
% 
% for i = 1:length(indexB),
%    lmi2 = lmi2 + [[lamB*(C_barA+C_barB) (Cb{indexB(i)}-C_barB)';Cb{indexB(i)}-C_barB (C_barA+C_barB)]>=0]; 
% end


%%

lmi1 = [];
lmi2 = [];

for i = 1:Na,
    lmi1 = lmi1 + [[lamA*(C_barA+C_barB) (Ca{i}-C_barA)' (Ca{i} - C_barB)';(Ca{i}-C_barA) (C_barA+C_barB) zeros(n,n);(Ca{i} - C_barB) zeros(n,n) (C_barA+C_barB)]>=0];
end

for i = 1:Nb,
    lmi2 = lmi2 + [[lamB*(C_barA+C_barB) (Cb{i}-C_barA)' (Cb{i} - C_barB)';(Cb{i}-C_barA) (C_barA+C_barB) zeros(n,n);(Cb{i} - C_barB) zeros(n,n) (C_barA+C_barB)]>=0];
end


%%

% lmi1 = [];
% lmi2 = [];
% 
% for i = 1:length(indexA),
%     
%    lmi1 = lmi1 + [[lamA*trace(C_barA+C_barB) trace(Ca{indexA(i)}-C_barA)';trace(Ca{indexA(i)}-C_barA) trace(C_barA+C_barB)]>=0];
% end
% 
% for i = 1:length(indexB),
%    lmi2 = lmi2 + [[lamB*trace(C_barA+C_barB) trace(Cb{indexB(i)}-C_barB)';trace(Cb{indexB(i)}-C_barB) trace(C_barA+C_barB)]>=0]; 
% end


% 
%%
% lmi3 = [];
% lmi4 = [];
% 
% for i = 1:length(ComplementA),
%     
%    lmi3 = lmi3 + [[LamMaxA(ComplementA(i))*(C_barA+C_barB) (Ca{ComplementA(i)}-C_barA)';Ca{ComplementA(i)}-C_barA (C_barA+C_barB)]>=0];
% end
% 
% for i = 1:length(ComplementB),
%    lmi4 = lmi4 + [[LamMaxA(ComplementB(i))*(C_barA+C_barB) (Cb{ComplementB(i)}-C_barB)';Cb{ComplementB(i)}-C_barB (C_barA+C_barB)]>=0]; 
% end


%% Retrições de igualdade
% ineqs = [C_barA(1,1) - C_barB(1,1) == 0];

% %%
% auxA = zeros(n,n*Na);
% auxB = zeros(n,n*Nb);
% 
% 
% for i=1:Na,
%    auxA(:,1+(i-1)*n:n+(i-1)*n) =  C_barA-Ca{i};
% end
% %%
% for i=1:Nb,
% auxB(:,1+(i-1)*n:n+(i-1)*n) = C_barB-Cb{i};
% end
% %%
% somaA = C_barA+C_barB;
% 
% % auX = blkdiag(somaA,somaA);
% 
% auX = zeros((Na+1)*n,(Na+1)*n);
% 
% for i = 1:Na-1,
%     auX(1+(i-1)*n:n+(i-1)*n,1+(i-1)*n:n+(i-1)*n)  = somaA;
% end
% %%
% clear somaA
% 
% %%
% auxAA = auX;
% auxAA(1:n,1:n) = Na*lamA*auxAA(1:n,1:n);
% 
% auxAA(1:n,n+1:end) = auxA;
% auxAA(n+1:end,1:n) = auxA';
% %%
% auxBB = auX;
% auxBB(1:n,1:n) = Nb*lamB*auxBB(1:n,1:n);
% 
% auxBB(1:n,n+1:end) = auxB;
% auxBB(n+1:end,1:n) = auxB';
% 
% clear auX auxA auxB Ca Cb
% %%
% solvesdp([H>=0,lmi,auxAA>=0,auxBB>=0],1,sdpsettings('solver','sedumi'))
%%
%%
solvesdp([H>=0,lmi1,lmi2],trace(C_barB - C_barA),sdpsettings('solver','sedumi'))
% solvesdp([H>=0,lmi],1,sdpsettings('solver','sedumi'))
% solvesdp([H>=0,lmi],1,sdpsettings('solver','sedumi'))
% solvesdp([H>=0,lmi1,lmi2,trace(H)==q,Ma>=0,Mb>=0],trace(Ma)-trace(Mb),sdpsettings('solver','sedumi'))
% trace(C_barB-C_barA)
% C_barA = Ma;
% C_barB = Mb;
clc
clear LamMaxA LamMaxB
%%
disp('feasible (+):')
% % min(checkset(H>=0))
% min(checkset(lmi))
% min(checkset(lmi1))
% min(checkset(lmi2))



%%

H = double(H);

% figure
% subplot(211);plot(diag(H),'--o')
% subplot(212);plot(eig(H),'--o')
%%

h = diag(H);
f = linspace(faixa(1),faixa(2),.5*q);

% figure
% plot(f,h(.5*q+1:q),'--ro');hold on
% plot(f,h(1:.5*q),'--o')
% legend('COS','SEN')

%% Projeto CSP com Ho otimo

C_barA = double(C_barA);
C_barB = double(C_barB);

Q = C_barA+C_barB;
Sa = Q^(-.5)*C_barA*Q^(-.5);
Sb = Q^(-.5)*C_barB*Q^(-.5);

for i=1:Na, LamMaxA(i) = max(eig((Sa-Q^(-.5)*double(Ca{i})*Q^(-.5))^2));end
for i=1:Nb, LamMaxB(i) = max(eig((Sb-Q^(-.5)*double(Cb{i})*Q^(-.5))^2));end

% figure(12)
% subplot(212);plot(LamMaxB,'--o')
% hold on
% subplot(211);plot(LamMaxA,'--o')

for i=1:Na, LamMaxAB(i) = max(eig((Sb-Q^(-.5)*double(Ca{i})*Q^(-.5))^2));end
for i=1:Nb, LamMaxBA(i) = max(eig((Sa-Q^(-.5)*double(Cb{i})*Q^(-.5))^2));end
% hold on
% subplot(212);plot(LamMaxBA,'--r*')
% hold on
% subplot(211);plot(LamMaxAB,'--r*')




% subplot(211);plot(indexA,LamMaxA(indexA),'ks')
% subplot(212);plot(indexB,LamMaxB(indexB),'ks')

clear Q Sa Sb LamMaxA LamMaxB LamMaxAB LamMaxBA

%%

[V,D] = eig(C_barB,C_barA+C_barB);
[~,I] = sort(diag(D));
W = V(:,I([1 2 3 20 21 22]));

% figure(1)
% hold on
% plot(diag(D),'--ro')
% legend('Normal','Otim')
% axis([0 23 0 1])
%% Features vector Otim
for i=1:Na+Nb,
    featOtim(:,i) = diag(W'*F(:,faixaFreq,i)*H*F(:,faixaFreq,i)'*W);
    featOtim(:,i) = log(featOtim(:,i)/sum(featOtim(:,i)));
end
%%
% figure(1)
% plot(featOtim','--o')

%% Projeto LDA com Ho otim
xxA = [featOtim(:,1:Na);ones(1,Na)];
xxB = [featOtim(:,Na+1:Na+Nb);ones(1,Nb)];

mA = mean(xxA,2);
mB = mean(xxB,2);

Sa = xxA*xxA' - mA*mA';
Sb = xxB*xxB' - mB*mB';

wo = pinv(Sa+Sb)*(mA - mB);
bo = .5*w'*(mA+mB);

yTo = wo'*[xxA xxB] - bo;

% figure
% plot(yTo,'--o')
% grid
clear mA mB Sa Sb 

%% Evaluation Trainning Set
% disp('CSP+LDA: T')
% .5*(eA+eB)


eAo = sum(yTo(1:Na)>0)/Na;

eBo = sum(yTo(Na+1:Na+Nb)<0)/Nb;
% disp('CSP+LDA: T - Optimized')
% .5*(eAo+eBo)


%% Carrega dados de Validacao
clear F Fab Fcd AB indA indB CD indC indD Na Nb
load(filenameE)
clear Xo incremento fname labels

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


%% Normalizacao do F
% for i=1:Na+Nb,
%     F(:,faixaFreq,i) = F(:,faixaFreq,i)/norm(F(:,faixaFreq,i),'fro');
% end

% F(:,:,86) = F(:,:,87); Sub01T error trial 86, channel 18
%% Gera matrizes de COV para Ho ident
Cy = zeros(n,n,Na+Nb);
for i=1:Na+Nb,
  aux = F(:,faixaFreq,i)*F(:,faixaFreq,i)';
%   Cy(:,:,i) = aux/trace(aux);
Cy(:,:,i) = aux;
end
clear aux
%% Vetro de Caracteristica para Ho ident
for i=1:Na+Nb,
    featY(:,i) = diag(Wc'*Cy(:,:,i)*Wc);
    featY(:,i) = log(featY(:,i)/sum(featY(:,i)));
end

% figure
% plot(featY','--o')
yE = w'*[featY;ones(1,Na+Nb)]-b;
% figure
% plot(yE,'--o')

%% Evaluation Validation Set

eAe = sum(yE(1:Na)>0)/Na;

eBe = sum(yE(Na+1:Na+Nb)<0)/Nb;

% disp('CSP+LDA: E')
% .5*(eAe+eBe)

% disp('Total Normal :')
% .5*(.5*(eAe+eBe)+ .5*(eA+eB))
%% Otimaized

Cyo = zeros(n,n,Na+Nb);
for i=1:Na+Nb,
  auxO = F(:,faixaFreq,i)*H*F(:,faixaFreq,i)';
%   Cy(:,:,i) = aux/trace(aux);
Cyo(:,:,i) = auxO;
end
clear auxO
%%
for i=1:Na+Nb,
    featYo(:,i) = diag(W'*Cyo(:,:,i)*W); 
    featYo(:,i) = log(featYo(:,i)/sum(featYo(:,i)));
end

% figure
% plot(featYo','--o')
yEo = wo'*[featYo;ones(1,Na+Nb)]-bo;
% figure
% plot(yEo,'--o')

%% Evaluation Validation Set

eAeO = sum(yEo(1:Na)>0)/Na;
eBeO = sum(yEo(Na+1:Na+Nb)<0)/Nb;

% disp('CSP+LDA: E Optimized')
% .5*(eAeO+eBeO)

% disp('Total: Optimized')

% % .5*(.5*(eAeO+eBeO)+.5*(eAo+eBo))

%%
% plot(diag(C_barA),'--o')
% hold on
% plot(diag(C_barB),'--ro')
% legend('A','B')

% figure
% plot([sort(eig((C_barA^(-.5)*C_barB*C_barA^(-.5)))) sort(eig(cA^(-.5)*cB*cA^(-.5)))],'--o')

%% 
table = zeros(2,3);
table(1,1) = .5*(eA+eB); % H = I, conjunto T
table(1,2) = .5*(eAe+eBe); % H = I, conjunto E
table(1,3) = .5*(.5*(eAe+eBe)+ .5*(eA+eB)); % H = I, conjunto T+E

table(2,1) = .5*(eAo+eBo); % H = D, conjunto T
table(2,2) = .5*(eAeO+eBeO); % H = D, conjunto E
table(2,3) = .5*(.5*(eAeO+eBeO)+.5*(eAo+eBo)); % H = D, conjunto T+E

%%
metric_mean = {'T','E','T+E'};
subj = {'Ho = I','Ho = D'};
disp('----------------------------------------------------');
disp('Accuracy (%) - Rows : H, Colums : Conjuntos');
disp('----------------------------------------------------');
displaytable(table,metric_mean,10,{'.4f'},subj)
disp('----------------------------------------------------');