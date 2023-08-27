% Initialize Parameters
MD = 128;  % Number of cases
MX = 256;  % Number of controls
N = 3;    % Dimension
abs_shift = 0.05;
mean_shift = abs_shift * randn(1,N)/2;  % Create random shift

% Generate Data for Cases and Controls
D = makelr_0(MD,N,2,0.25) + repmat(mean_shift,MD,1);  % Data for cases
X = makelr_0(MX,N,2,0.25) - repmat(mean_shift,MX,1);  % Data for controls

% Prepare Matrices for Computations
% Computing transpose and ones vectors to use in matrix operations
Dt = transpose(D); eDn = ones(MD,1); eDt = ones(1,MD);
Xt = transpose(X); eXn = ones(MX,1); eXt = ones(1,MX);

% Compute A, B, C Matrices
A = ((1/MD)*Dt*eDn-(1/MX)*Xt*eXn)*((1/MD)*eDt*D-(1/MX)*eXt*X);
B = (1/MD)*Dt*(D - (1/MX)*eDn*(eXt*X)) + (1/MX)*Xt*(X - (1/MD)*eXn*(eDt*D));
C = (1/MD)*Dt*D - (1/MX)*Xt*X;

% Compute Eigenvectors
[UA,~] = eigs(A,1,'LM');
[UB,~] = eigs(B,1,'LM');
[UC,~] = eigs(C,1,'LM');

% Project Data
pDA = D*UA; pXA = X*UA;
pDB = D*UB; pXB = X*UB;
pDC = D*UC; pXC = X*UC;

% Compute Histograms
h_avg = mean([pDA;pDB;pXA;pXB]);
h_std = std([pDA;pDB;pXA;pXB]);
nbins = 32; hbins = linspace(h_avg-4*h_std,h_avg+4*h_std,nbins); hbins_c = linspace(0,h_avg+4*h_std,nbins/2);

% Histograms for different projections
hDA = hist(pDA,hbins)/MD; hXA = hist(pXA,hbins)/MX;
hDB = hist(pDB,hbins)/MD; hXB = hist(pXB,hbins)/MX;
hDC = hist(abs(pDC),hbins_c)/MD; hXC = hist(abs(pXC),hbins_c)/MX;

% Plotting
figure;clf;
markersize_use = 8;
subplot(1,2,1);hold on;
plot3(D(:,1),D(:,2),D(:,3),'ro','MarkerSize',markersize_use,'MarkerFaceColor',[0.8,0,0]);
plot3(X(:,1),X(:,2),X(:,3),'bo','MarkerSize',markersize_use,'MarkerFaceColor',[0,0,0.8]);
plot3(+h_std*[0;UA(1)],+h_std*[0;UA(2)],+h_std*[0;UA(3)],'-','Color',[0,1.0,0],'LineWidth',4);
plot3(-h_std*[0;UA(1)],-h_std*[0;UA(2)],-h_std*[0;UA(3)],'-','Color',[0,0.8,0],'LineWidth',4);
plot3(+h_std*[0;UB(1)],+h_std*[0;UB(2)],+h_std*[0;UB(3)],'-','Color',[1.0,0,1.0],'LineWidth',4);
plot3(-h_std*[0;UB(1)],-h_std*[0;UB(2)],-h_std*[0;UB(3)],'-','Color',[0.8,0,0.8],'LineWidth',4);
plot3(+h_std*[0;UC(1)],+h_std*[0;UC(2)],+h_std*[0;UC(3)],'-','Color',[1.0,0,0],'LineWidth',4);
plot3(-h_std*[0;UC(1)],-h_std*[0;UC(2)],-h_std*[0;UC(3)],'-','Color',[0.8,0,0],'LineWidth',4);
hold off; axis equal; axis vis3d;
title('scatterplot (rotate me!)');
subplot(3,2,2);hold on;
bar(hbins,+hDA,'r');
bar(hbins,-hXA,'b');
legend('case','ctrl');
hold off;
auc_A = auc_0(hDA,hXA); auc_A = max(auc_A,1-auc_A);
title(sprintf('case-ctrl dist after proj-A (auc = %0.3f)',auc_A));
subplot(3,2,4);hold on;
bar(hbins,+hDB,'r');
bar(hbins,-hXB,'b');
legend('case','ctrl');
hold off;
auc_B = auc_0(hDB,hXB); auc_B = max(auc_B,1-auc_B);
title(sprintf('case-ctrl dist after proj-B (auc = %0.3f)',auc_B));
subplot(3,2,6);hold on;
bar(hbins_c,+hDC,'r');
bar(hbins_c,-hXC,'b');
legend('case','ctrl');
hold off;
auc_C = auc_0(hDC,hXC); auc_C = max(auc_C,1-auc_C);
title(sprintf('case-ctrl dist after abs(proj-C) (auc = %0.3f)',auc_C));



