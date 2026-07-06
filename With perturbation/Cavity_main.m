clear 
close all

global dtau ddtau omega
global nt nq EG Tr Tl N_RR N_nr
global g0 g_p E_p delta_g g0_p gini
global output_u output_nl output_os

aa = load('Mini16_80.mat');
uuini = aa.uu;
giniT = aa.g0_p;  gini = giniT(end,:);

%%%%%%%%%%%%%%%%%%        Note           %%%%%%%%%%%%%%%%%
%%%% this programe is for a 4.6 GHz cavity at 1.5 um %%%%
%%%%    for studty the excessive noise for PengXiang    %%%%

N_rep=1; N_sc=80;       % Scalability of repetition rate and soliton number in QSS
N_RR=N_rep*N_sc;        % Effective N_QS              
EG=5e6;                 % gain saturation energy, 5e6,pJ
Tr=220;                 % Roundtrip time
Tl=(1e-5)*1e12;         % Original version (0.8e-5)*1e12,ps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Cavity parameter %%%

distan_1=0.022;            %length of the gain fiber, 4.6 GHz
step_num1=20;              %Constant step,20
step1=distan_1/step_num1;
g0=5;                      %6 for Nnr=5-10;        
%%%%%%%%%%%%%%%%%%%%%%%%%


%%% 划分时域频域网格 

nt=2048*8; Tmax=110;     %时间窗口 (2048*8,500)
nt1=512;                %白噪声带宽
ddtau=2*Tmax/nt;
%时域频域
dtau=(-nt/2:nt/2-1).*2*Tmax/nt;
omega=[(0:nt/2-1),(-nt/2:-1)].*pi/Tmax;

lambda_c=1563;             % Central wavelength
lambda_re=6*pi*1e5./(fftshift(omega+6*pi*1e5/lambda_c)); 
%%%%%%%%%%%%%%%%%%%%%%%%%


%%% 记录矩阵

% Eb_p=zeros(Nbst,2*step_num1);
% g0_p=zeros(Nbst,2*step_num1);

% gRec_p=zeros(N_nr,Nbst-1);
% ERec_p=zeros(N_nr,Nbst-1);
gRec_pp=[];
ERec_pp=[];

% output_nl=zeros(N_nr,nt);
% output_os=zeros(N_nr,nt);
% output_oslog=zeros(N_nr,nt);

% output_trainS=zeros(Nbst-1,N_nr*nt);
outputT_trainS=[];
outputW_trainS=[];

delta_g=zeros(1,2*step_num1);
Ebst_p=zeros(1,2*step_num1);

tic

Nbst=0;                           % Number of bursts
nq=1;                             % Mark the iteration
N_tot=10000;                      % 15000 for noise analysis

while Nbst<=N_tot
  
  N_nr=randi([16,17]);             % Number of the QSSs in burst, [9,10]
  % N_nr = 16;                     % 8      
  
  uu = uuini;
  uuini = QSS_sub(uu);
  
  gRec_pp = [gRec_pp;g_p(:,1)];
  ERec_pp = [ERec_pp;E_p(:,1)];
  
  Eb_p(nq,:)=sum(E_p,1);  
  % Ebst_p=Ebst_p+Eb_p(nq,:);
  Ebst_p=Eb_p(nq,:);
  
  % Burst-energy rate equation
  N_RR2=N_RR; Tbst=(1e-5)*1e12;
  Tb=N_nr*Tr*N_RR2; Tl1=(50e-6)*1e12;
  g0_p(nq,:)=g_p(N_nr,:);
  
  if nq==1
   
     RK_s1=-(g0_p(nq,:)-g0)./Tl-g0_p(nq,:).*Ebst_p(1,:)./(EG*Tb);
     RK_s2=-((g0_p(nq,:)+Tb.*RK_s1)-g0)./Tl-(g0_p(nq,:)+Tb.*RK_s1).*Ebst_p(1,:)./(EG*Tb);
     g0_p(nq+1,:)=g0_p(nq,:)+Tb/2.*(RK_s1+RK_s2);
     
  else
      
     RK_s1=-(g0_p(nq-1,:)-g0)./Tl-g0_p(nq-1,:).*Ebst_p(1,:)./(EG*Tb);
     RK_s2=-((g0_p(nq-1,:)+Tb.*RK_s1)-g0)./Tl-(g0_p(nq-1,:)+Tb.*RK_s1).*Ebst_p(1,:)./(EG*Tb);
     RK_s3=-((g0_p(nq-1,:)+Tb.*RK_s2)-g0)./Tl-(g0_p(nq-1,:)+Tb.*RK_s2).*Ebst_p(1,:)./(EG*Tb);
     RK_s4=-((g0_p(nq-1,:)+2*Tb.*RK_s3)-g0)./Tl-(g0_p(nq-1,:)+2*Tb.*RK_s3).*Ebst_p(1,:)./(EG*Tb);
     g0_p(nq+1,:)=g0_p(nq-1,:)+Tb/3.*(RK_s1+2*RK_s2+2*RK_s3+RK_s4);
     
  end
  
  delta_g=(g_p(N_nr,:)-g0_p(nq+1,:)).*(1/Tl+Ebst_p(1,:)./(EG*Tb));
  
  if nq==1
       outputT_trainS = output_nl;
       outputW_trainS = output_os;
  else
      outputT_trainS = [outputT_trainS;output_nl];
      outputW_trainS = [outputW_trainS;output_os];
  end
  
  disp(nq);
  nq = nq+1;
  Nbst=Nbst+N_nr;
  
end 

toc
temp1_u=abs(fftshift(ifft(output_u))).^2;  
temp1_norm=mat2gray(temp1_u);
temp1u_re=10.*log10(temp1_u);

output_train=mat2gray(reshape(outputT_trainS',[1,numel(outputT_trainS)]));

figure(1); clf;
subplot(121);
plot(dtau,abs(output_u).^2,'LineWidth',2); hold on;
set(gca,'Fontsize',15,'LineWidth',1.5) 
xlabel('Time (ps)'); ylabel('Intensity (W)');   

subplot(122);
plot(lambda_re,temp1u_re,'LineWidth',2); 
set(gca,'Fontsize',15,'LineWidth',1.5,'layer','top'); 
axis([1545 1575 -125 -25]); % set(gca,'ytick',(-75:25:-25));
xlabel('Wavelength (nm)'); ylabel('Intensity (dB)');

% %%
% figure(3); clf;
% plot(Eb_p(1:end-1,1),'LineWidth',2); hold on;
% set(gca,'LineWidth',1.5); set(gca,'Fontsize',15); 
% xlabel('Roundtrips'); ylabel('Energy(pJ)');
% 
% %%
% figure(4); clf;
% t_train=ddtau*N_RR*(1e-6).*(1:20:numel(output_train));
% plot(output_train(1:20:end),'LineWidth',1.5,'Color',[0.2 0.2 0.6]);
% set(gca,'Fontsize',17); set(gca,'LineWidth',1.5); set(gca,'Layer','top');
% ylim([0 1.1]); set(gca,'layer','top'); % xlim([-14 13]); set(gca,'xtick',(-12:6:12))
% set(gca,'ytick',(0:0.5:1)); % set(gca,'xtick',(-2:2:2));
% xlabel('Time (\mus)'); ylabel('Norm. int.');
% 
% %%


%%
figure(2); clf;
subplot(211);
plot(output_train(1:20:end),'LineWidth',1.5);
set(gca,'Fontsize',15,'LineWidth',1.5,'Layer','top');
ylim([0 1.1]); set(gca,'layer','top');
% xlabel('Time (\mus)'); ylabel('Norm. int.');

subplot(212);
plot((1:numel(gRec_pp)),gRec_pp,'linewidth',2); hold on;
set(gca,'linewidth',1.5,'fontsize',15); 
ylabel('Gain (m^-1)');


figure(3); clf;
subplot(211);
plot(output_train(1:20:end),'LineWidth',1.5);
set(gca,'Fontsize',15,'LineWidth',1.5,'Layer','top');
ylim([0 1.1]); set(gca,'layer','top'); % set(gca,'xtick',(-12:6:12))
% xlabel('Time (\mus)'); ylabel('Norm. int.');

subplot(212);
plot((1:numel(ERec_pp)),ERec_pp,'linewidth',2); hold on;
set(gca,'linewidth',1.5,'fontsize',17); ylabel('Energy');


%% Temporal and spectral evolution

% N_tot=size(outputT_trainS,1);
figure(4); clf;
pcolor(dtau,(1:10:Nbst)',mat2gray(outputT_trainS(1:10:Nbst,:))); shading interp;
set(gca,'LineWidth',1.5,'Fontsize',15,'Layer','top');
xlabel('Time (ps)'); ylabel('Roundtrips'); 
load 'mycolor'; colormap(mycolor);

%%
% figure(5); clf;
% pcolor((0:10:300),lambda_re',mat2gray(outputW_trainS(4500:10:4800,:)')); shading interp;
% set(gca,'LineWidth',2); set(gca,'Fontsize',16); set(gca,'Layer','top');
% ylabel('Wavelength (nm)'); xlabel('Roundtrips'); 
% load 'mycolor'; colormap(mycolor); caxis([0.1 1]);
% axis([0 300 1555 1565]); set(gca,'ytick',(1555:5:1565)); set(gca,'xtick',(0:150:300));


%%
% OS_ave=sum(outputW_trainS(N_start:N_end,:))./(N_end-N_start+1);
% figure(6); clf;
% plot(lambda_re,10.*log10(mat2gray(OS_ave)),'linewidth',2);
% axis([1520 1600 -100 0]);
% 
% lwf=[lambda_re',(10.*log10(mat2gray(OS_ave)))'];

%% Low-pass filter by applying Butterworth filter
nt_tot=numel(output_train);
omega_tot=[(0:nt_tot/2-1),(-nt_tot/2:-1)].*pi/(Tmax*nt_tot/nt);
t_train=ddtau*N_RR*(1e-6).*(1:1000:numel(output_train));
output_Fil=abs(ifft(fft(output_train).*sqrt(1./(1+(omega_tot./0.06).^16))));   % 0.13
output_wFil=abs(fftshift((fft(output_train).*sqrt(1./(1+(omega_tot./0.06).^16)))));

%%
figure(10); clf;
% subplot(211);
% plot(output_train(1:100:end),'LineWidth',1.5,'Color',[0.2 0.2 0.6]);
% set(gca,'Fontsize',17); set(gca,'LineWidth',1.5); set(gca,'Layer','top');
% ylim([0 1.1]); set(gca,'layer','top');
% subplot(212);
% plot(t_train-110,mat2gray(output_Fil(1:1000:end)),'LineWidth',1.5,'Color',[141/255 89/255 179/255]);
plot(t_train-110,mat2gray(output_Fil(1:1000:end)),'LineWidth',1.5,'Color',[188/255 143/255 221/255]);
set(gca,'fontsize',15,'lineWidth',2,'layer','top'); 
ylim([0 1.1]); set(gca,'ytick',(0:0.5:1));  
xlim([-60 60]); set(gca,'xtick',(-60:30:60));  
xlabel('Time (\mus)'); % xlim([20 90]);
% 
% % figure(11); clf;
% % loglog(fftshift(omega_tot),output_wFil,'linewidth',2);
% % set(gca,'Fontsize',17); set(gca,'LineWidth',1.5);

%%
% figure(11); clf;
% plot(2*Tmax*N_RR.*(1:numel(gRec_p))./1e6-12,gRec_pp,'linewidth',2,'Color',[0.9 0.4 0.4]); hold on;
% plot(2*Tmax*N_RR.*(1:N_nr:numel(gRec_p))./1e6-12,gRec_pp(1:N_nr:end),'.','Markersize',25,'Color',[247/255 170/255 0/255]); hold on;
% plot(2*Tmax*N_RR.*(1:numel(gRec_p))./1e6-12,5.757.*ones(1,numel(gRec_p)),'--','linewidth',2,'Color',[247/255 170/255 0/255])
% set(gca,'linewidth',1.5); set(gca,'fontsize',17); ylabel('Gain (m^-1)');
% xlim([0 2]);set(gca,'xtick',(0:2:2)); set(gca,'layer','top'); 
% ylim([0 8]); set(gca,'ytick',(0:4:8))
% xlabel('Time (\mus)');

%% 3D phase plot
% N_start=N_tot-200; N_end=N_tot; 
% % N_evo1=(N_end-N_start+1)-100;
% % N_delay=5;
% 
N_start=4000;
N_evo1=4500; N_delay=10;
EE_evo1=ERec_pp(N_start:N_evo1,1); EE_evo2=ERec_pp((N_start+N_delay:N_evo1+N_delay),1);

% Pk_real=max(outputT_trainS(N_start:N_end,:)');
% % Pk_evo=mat2gray(max(outputT_trainS(N_start:N_end,:)'));
% % Pk_evo1=Pk_evo(1:N_evo1);Pk_evo2=Pk_evo(1+N_delay:N_evo1+N_delay);
% figure(7); plot(Pk_real,'.-');
% filename1=['Data_' num2str(N_nr)];
% var_name=[(gRec0(end).*ones(size(Pk_real)))',Pk_real'];
% eval([filename1 '= var_name']);
%
figure(11); clf;
% scatter(EE_evo,Pk_evo,10,num_evo,'filled'); % I
% scatter(EE_evo1,EE_evo2,10,(N_start:N_evo1),'filled'); % II
% scatter(Pk_evo1,Pk_evo2,10,(N_start:N_evo1),'filled'); % III
% scatter3(Pk_evo1,Pk_evo2,gRec_pp(N_start:(N_start+N_evo1-1)),20,(N_start:(N_start+N_evo1-1)),'filled'); hold on;
plot3(EE_evo1,EE_evo2,gRec_pp(N_start:N_evo1),'linewidth',2,'Color',[0.6 0.6 0.6]); view([-20 40]);
xlim([30 33]); ylim([30 33]); zlim([0.8 1]); view([-36 53]);
% set(gca,'xscale','log'); set(gca,'yscale','log')
set(gca,'Fontsize',20); set(gca,'LineWidth',2); grid on; box on; 

%% Energy evolution 
% ERec_ppRe=ERec_pp(1:end); Npp=numel(ERec_ppRe);
% figure(7); clf;
% % plot(ERec_ppRe,'linewidth',2); hold on;
% plot((1:N_nr:Npp),ERec_ppRe(1:N_nr:Npp),'linewidth',3);
% set(gca,'linewidth',1.5); set(gca,'fontsize',17); ylabel('Energy');
% title(['Energy=' num2str(ERec_ppRe(end)) 'pJ'])