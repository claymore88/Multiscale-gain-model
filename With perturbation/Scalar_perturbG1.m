clear 
close all

%%%%%%%%%%%%%%%%%%        Note           %%%%%%%%%%%%%%%%%
%%%% this programe is for a 4.6 GHz cavity at 1.5 um %%%%
%%%%    for studty the excessive noise for PengXiang    %%%%

%%%%%%%饱和能量(提高Es的值等效于加大泵浦)%%%%%%%%

N_rep=1; N_sc=80;       % Scalability of repetition rate and strong correlation
N_RR=N_rep*N_sc;        % Repetition rate                
EG=5e6;                 % gain saturation energy, 5e6
Tr=220;                  % Roundtrip time
Tl=(1e-5)*1e12;         % Original version (0.8e-5)*1e12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SESAM parameter %%
l0=0.04;         %depth of the SESAM 
F_sat=15;        %saturation fluence [uJ/cm2]
tao_a=10;        %relaxation time
R_a_de=0.9;      %reflectivity at designed wavelength
beta2=-2000e-6;  
beta3=0;  beta4=0;  beta5=0;    
 
diameter_core=8.0e-6;         %diameter of the fiber,
Acore=pi*diameter_core^2/4;   %area of the fiber 
Ps=F_sat*1e-6*(Acore*1e4)*1e12/tao_a;  %saturation power [W] 
Ea=Ps*tao_a/2;
R_sat_de=R_a_de+l0;   %reflectivity after saturation

c=299792458;
distan_1=0.022;           %length of the gain fiber, 4.6 GHz
%%%%%%%调节光纤参数%%%%%%%%
beta_g=-10e-3;            %dispersion of gain fiber 
gama_g=3/1000;            %nonliearity of gain fiber 
BW=3*pi;                  %gain bandwidth, 3*pi for 24nm, 4pi for 32nm (3pi)
% g0=60/(EG/Tl);          %gain coefficient, 考虑除以2作为光场对应增益
g0=5;                     %6 for Nnr=5-10;                    
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%调节光路参数%%%%%%%%%
lambda_c=1569;             %central wavelength
step_num1=20;              %Constant step,20
step1=distan_1/step_num1;

kk=pi/(10*2*distan_1);     % birefringence ~0.2 for cavity length
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%划分时域频域网格 
nt=2048*8; Tmax=110;     %时间窗口 (2048*2,25)
nt1=512;                %白噪声带宽
ddtau=2*Tmax/nt;
%时域频域
dtau=(-nt/2:nt/2-1).*2*Tmax/nt;
omega=[(0:nt/2-1),(-nt/2:-1)].*pi/Tmax;
lambda_re=6*pi*1e5./(fftshift(omega+6*pi*1e5/lambda_c)); %wavelength matrix

%%%%%%%%%反射和透射率%%%%%%%%%
% film
R_f=0.99;       %assuming reflectivity of the film
beta_f=0e-6;    %dispersion of the film 

% GVD estimate
% lamda_1=xlsread('GVD_gf.xlsx','A:A'); ff=2*pi*c/1000./lamda_1;
% GVD_dat=xlsread('GVD_gf.xlsx','B:B'); GVD_1=GVD_dat./1e6;
% lamda_ref=1560; f_ref=2*pi*c/1000/lamda_ref;
% f0=ff-f_ref;        % angular frequency with respect to reference wavelength
% beta_coe=polyfit(f0,GVD_1,3);
% 
% % GVD_fit=beta_coe(4)+beta_coe(3).*f0+beta_coe(2).*f0.^2+beta_coe(1).*f0.^3;
% % figure;plot(f0,GVD_1); hold on; plot(f0,GVD_fit);
% 
% f_c=2*pi*c/1000/lambda_c;
% beta2=beta_coe(4)+beta_coe(3)*(f_c-f_ref)+beta_coe(2)*(f_c-f_ref)^2+beta_coe(1)*(f_c-f_ref)^3;
% beta3=beta_coe(3)+2*beta_coe(2)*(f_c-f_ref)+3*beta_coe(1)*(f_c-f_ref)^2;
% beta4=2*beta_coe(2)+6*beta_coe(1)*(f_c-f_ref); beta5=6*beta_coe(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%循环圈数/输出耦合%%%%%%%%%%%%%
Nbst=300;            % Number of bursts, 500
N_nr=16;             % Number of the QSSs in burst, 8
T_output=1-R_f;      % the output ratio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%输入信号%%%%%%%%%%%%%%%
%%% Soliton will collside 
% amp_temp=abs(wgn(1,nt1,1));  %White noise signal
% amp(1,1:nt/2-nt1/2)=0;       %
% amp(1,nt/2-nt1/2+1:nt/2+nt1/2)=(1).*amp_temp;
% amp(1,nt/2+nt1/2+1:nt)=0;   
amp=0.1.*sech((dtau+Tmax)./2);  %
uu=amp;                      %initial condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%记录矩阵
E_p=zeros(N_nr,2*step_num1);
Eb_p=zeros(Nbst,2*step_num1);
g0_p=zeros(Nbst,2*step_num1);
g_p=zeros(N_nr,2*step_num1);
gRec_p=zeros(N_nr,Nbst-1);
ERec_p=zeros(N_nr,Nbst-1);

output_nl=zeros(N_nr,nt);
output_os=zeros(N_nr,nt);
output_oslog=zeros(N_nr,nt);
%u_store=zeros(N_nr,nt);
output_trainS=zeros(Nbst-1,N_nr*nt);
outputT_trainS=zeros((Nbst-1)*N_nr,nt);
outputW_trainS=zeros((Nbst-1)*N_nr,nt);

qq=ones(1,nt).*l0;
delta_g=zeros(1,2*step_num1);
Ebst_p=zeros(1,2*step_num1);

tic
for nq=1:Nbst-1
  for nr=1:N_nr
    dispersion1_u=exp(1i*0.5*beta_g*omega.^2*step1/2);
    dispersion2_u=exp(1i*beta2/2*omega.^2+1i*beta3/6*(omega.^3) ...
                  +1i*beta4/24*(omega.^4)+1i*beta5/120*(omega.^5));     
    
    %%%% GF 1 %%%
    for n=1:step_num1
        E_p(nr,n)=trapz(dtau,abs(uu).^2); 
        % g_p(nr,n)=g0./(1+E_p(nr)/Es);
        
        if nq==1 && nr==1
           %g_p(nr,n)=g_p(end,n);
           g0_p(nq,n)=g0/(1+E_p(nr,n)/(EG*Tr/N_RR/Tl));
           g_p(nr,n)=g0_p(nq,n);
        elseif nq>1 && nr==1
           g_p(nr,n)=g0_p(nq,n); 
        elseif nr==2   
           RK_s1=-(g_p(nr-1,n)-g0)/Tl-g_p(nr-1,n)*E_p(nr-1,n)/(EG*Tr/N_RR)+delta_g(1,n);
           RK_s2=-((g_p(nr-1,n)+N_RR*Tr*RK_s1)-g0)/Tl-(g_p(nr-1,n)+N_RR*Tr*RK_s1)*E_p(nr,n)/(EG*Tr/N_RR)+delta_g(1,n);
           g_p(nr,n)=g_p(nr-1,n)+N_RR*Tr/2*(RK_s1+RK_s2); 
        else
           RK_s1=-(g_p(nr-2,n)-g0)/Tl-g_p(nr-2,n)*E_p(nr-2,n)/(EG*Tr/N_RR)+delta_g(1,n);
           RK_s2=-((g_p(nr-2,n)+N_RR*Tr*RK_s1)-g0)/Tl-(g_p(nr-2,n)+N_RR*Tr*RK_s1)*E_p(nr-1,n)/(EG*Tr/N_RR)+delta_g(1,n);
           RK_s3=-((g_p(nr-2,n)+N_RR*Tr*RK_s2)-g0)/Tl-(g_p(nr-2,n)+N_RR*Tr*RK_s2)*E_p(nr-1,n)/(EG*Tr/N_RR)+delta_g(1,n);
           RK_s4=-((g_p(nr-2,n)+2*N_RR*Tr*RK_s3)-g0)/Tl-(g_p(nr-2,n)+2*N_RR*Tr*RK_s3)*E_p(nr,n)/(EG*Tr/N_RR)+delta_g(1,n);
           g_p(nr,n)=g_p(nr-2,n)+N_RR*Tr/3*(RK_s1+2*RK_s2+2*RK_s3+RK_s4); 
        end
        
        uu_temp=ifft(uu).*dispersion1_u.*exp(-omega.^2*g_p(nr,n)/(BW^2)*step1/2+g_p(nr,n)*step1/2);
        uu1=fft(uu_temp);
        uu=uu1.*exp(1i*step1*gama_g.*(abs(uu1).^2));
        uu_temp=ifft(uu).*dispersion1_u.*exp(-omega.^2*g_p(nr,n)/(BW^2)*step1/2+g_p(nr,n)*step1/2);
        uu=fft(uu_temp); 
    end
   
    %%%% SESAM %%%%
    uu=fft(ifft(uu).*dispersion2_u);

    It=abs(uu).^2;
    for ii=1:nt-1
        RK1=-(qq(ii)-l0)/tao_a-qq(ii)*It(ii)/Ea;
        RK2=-((qq(ii)+ddtau/2*RK1)-l0)/tao_a-(qq(ii)+ddtau/2*RK1)*(It(ii)+It(ii+1))/2/Ea;
        RK3=-((qq(ii)+ddtau/2*RK2)-l0)/tao_a-(qq(ii)+ddtau/2*RK2)*(It(ii)+It(ii+1))/2/Ea;
        RK4=-((qq(ii)+ddtau*RK3)-l0)/tao_a-(qq(ii)+ddtau*RK3)*It(ii+1)/Ea;
        qq(ii+1)=qq(ii)+ddtau/6*(RK1+2*RK2+2*RK3+RK4);
    end
    
    TT_non=R_sat_de-qq;
    fai_npc=0.5*2.*(qq-l0)./2; % 0.5*2.*(qq-l0)./2
    uu=uu.*sqrt(TT_non).*exp(0.5*2*1i.*(qq-l0)./2);  % .*exp(0.5*2*1i.*(qq-l0)./2);
    %uu=uu.*sqrt(TT_non); %.*exp(1i.*qq./2);
    
    %loss_profile1=qq;
    
    %%%% GF 2 %%%
    for n=1+step_num1:step_num1+step_num1
        E_p(nr,n)=trapz(dtau,abs(uu).^2); 
        % g_p(nr,n)=g0./(1+E_p(nr)/Es);
        
        if nq==1 && nr==1
           %g_p(nr,n)=g_p(end,n);
           g0_p(nq,n)=g0/(1+E_p(nr,n)/(EG*Tr/N_RR/Tl));
           g_p(nr,n)=g0_p(nq,n);
        elseif nq>1 && nr==1
           g_p(nr,n)=g0_p(nq,n); 
        elseif nr==2   
           RK_s1=-(g_p(nr-1,n)-g0)/Tl-g_p(nr-1,n)*E_p(nr-1,n)/(EG*Tr/N_RR)+delta_g(1,n);
           RK_s2=-((g_p(nr-1,n)+N_RR*Tr*RK_s1)-g0)/Tl-(g_p(nr-1,n)+N_RR*Tr*RK_s1)*E_p(nr,n)/(EG*Tr/N_RR)+delta_g(1,n);
           g_p(nr,n)=g_p(nr-1,n)+N_RR*Tr/2*(RK_s1+RK_s2); 
        else
           RK_s1=-(g_p(nr-2,n)-g0)/Tl-g_p(nr-2,n)*E_p(nr-2,n)/(EG*Tr/N_RR)+delta_g(1,n);
           RK_s2=-((g_p(nr-2,n)+N_RR*Tr*RK_s1)-g0)/Tl-(g_p(nr-2,n)+N_RR*Tr*RK_s1)*E_p(nr-1,n)/(EG*Tr/N_RR)+delta_g(1,n);
           RK_s3=-((g_p(nr-2,n)+N_RR*Tr*RK_s2)-g0)/Tl-(g_p(nr-2,n)+N_RR*Tr*RK_s2)*E_p(nr-1,n)/(EG*Tr/N_RR)+delta_g(1,n);
           RK_s4=-((g_p(nr-2,n)+2*N_RR*Tr*RK_s3)-g0)/Tl-(g_p(nr-2,n)+2*N_RR*Tr*RK_s3)*E_p(nr,n)/(EG*Tr/N_RR)+delta_g(1,n);
           g_p(nr,n)=g_p(nr-2,n)+N_RR*Tr/3*(RK_s1+2*RK_s2+2*RK_s3+RK_s4); 
        end
        
        uu_temp=ifft(uu).*dispersion1_u.*exp(-omega.^2*g_p(nr,n)/(BW^2)*step1/2+g_p(nr,n)*step1/2);
        uu1=fft(uu_temp);
        uu=uu1.*exp(1i*step1*gama_g.*(abs(uu1).^2));
        uu_temp=ifft(uu).*dispersion1_u.*exp(-omega.^2*g_p(nr,n)/(BW^2)*step1/2+g_p(nr,n)*step1/2);
        uu=fft(uu_temp);
        
    end
    
    %集总损耗-耦合输出
    uu=fft(ifft(uu).*exp(1i*beta_f/2.*(omega.^2)));
    output_u=uu.*sqrt(T_output);
    uu=uu.*sqrt(1-T_output);
    
    output_nl(nr,:)=abs(uu).^2;
    output_os(nr,:)=abs(fftshift(ifft(uu))).^2;
    output_oslog(nr,:)=10.*log10(output_os(nr,:)); 
    
    % u_store(nr,:)=uu;
    
    %monitor%
    %disp(nr);
  end
  
  output_trainS(nq,:)=reshape(output_nl',[1,numel(output_nl)]);
  outputT_trainS((1+(nq-1)*N_nr):(nq*N_nr),:)=output_nl;
  outputW_trainS((1+(nq-1)*N_nr):(nq*N_nr),:)=output_os;
  
  gRec_p(:,nq)=g_p(:,1);
  ERec_p(:,nq)=E_p(:,1);
  
  Eb_p(nq,:)=sum(E_p,1);  
  % Ebst_p=Ebst_p+Eb_p(nq,:);
  Ebst_p=Eb_p(nq,:);
  
  % Burst-energy rate equation
  % Tbst=(5e-6)*1e12;
  % N_RR2=50; Tbst=(1e-5)*1e12;
  N_RR2=N_RR; Tbst=(1e-5)*1e12;
  Tb=N_nr*Tr*N_RR2; Tl1=(50e-6)*1e12;
  g0_p(nq,:)=g_p(N_nr,:);
  if nq==1
     
     % RK_s1=-(g0_p(nq,:)-g0)./Tl-g0_p(nq,:).*Ebst_p(1,:)./(EG*Tbst/N_RR2);
     % RK_s2=-((g0_p(nq,:)+Tb.*RK_s1)-g0)./Tl-(g0_p(nq,:)+Tb.*RK_s1).*Ebst_p(1,:)./(EG*Tbst/N_RR2);
     % g0_p(nq+1,:)=g0_p(nq,:)+Tb/2.*(RK_s1+RK_s2);
     
     RK_s1=-(g0_p(nq,:)-g0)./Tl-g0_p(nq,:).*Ebst_p(1,:)./(EG*Tb);
     RK_s2=-((g0_p(nq,:)+Tb.*RK_s1)-g0)./Tl-(g0_p(nq,:)+Tb.*RK_s1).*Ebst_p(1,:)./(EG*Tb);
     g0_p(nq+1,:)=g0_p(nq,:)+Tb/2.*(RK_s1+RK_s2);
     
  else
      
     % RK_s1=-(g0_p(nq-1,:)-g0)./Tl-g0_p(nq-1,:).*Ebst_p(1,:)./(EG*Tbst/N_RR2);
     % RK_s2=-((g0_p(nq-1,:)+Tb.*RK_s1)-g0)./Tl-(g0_p(nq-1,:)+Tb.*RK_s1).*Ebst_p(1,:)./(EG*Tbst/N_RR2);
     % RK_s3=-((g0_p(nq-1,:)+Tb.*RK_s2)-g0)./Tl-(g0_p(nq-1,:)+Tb.*RK_s2).*Ebst_p(1,:)./(EG*Tbst/N_RR2);
     % RK_s4=-((g0_p(nq-1,:)+2*Tb.*RK_s3)-g0)./Tl-(g0_p(nq-1,:)+2*Tb.*RK_s3).*Ebst_p(1,:)./(EG*Tbst/N_RR2);
     % g0_p(nq+1,:)=g0_p(nq-1,:)+Tb/3.*(RK_s1+2*RK_s2+2*RK_s3+RK_s4);
     
     RK_s1=-(g0_p(nq-1,:)-g0)./Tl-g0_p(nq-1,:).*Ebst_p(1,:)./(EG*Tb);
     RK_s2=-((g0_p(nq-1,:)+Tb.*RK_s1)-g0)./Tl-(g0_p(nq-1,:)+Tb.*RK_s1).*Ebst_p(1,:)./(EG*Tb);
     RK_s3=-((g0_p(nq-1,:)+Tb.*RK_s2)-g0)./Tl-(g0_p(nq-1,:)+Tb.*RK_s2).*Ebst_p(1,:)./(EG*Tb);
     RK_s4=-((g0_p(nq-1,:)+2*Tb.*RK_s3)-g0)./Tl-(g0_p(nq-1,:)+2*Tb.*RK_s3).*Ebst_p(1,:)./(EG*Tb);
     g0_p(nq+1,:)=g0_p(nq-1,:)+Tb/3.*(RK_s1+2*RK_s2+2*RK_s3+RK_s4);
     
  end
  
  delta_g=(g_p(N_nr,:)-g0_p(nq+1,:)).*(1/Tl+Ebst_p(1,:)./(EG*Tb));
  
  % Integral-like equation
  % EG/10 for relatively good result
  % EG/50 for spaced result
  % g0_p(nq+1,:)=g_p(N_nr,:);
  % g0_p(nq+1,:)=g_p(N_nr,:).*exp(-Eb_p(nq,:)./(500*EG));
  
  % g0_p(nq+1,:)=g_p(N_nr,:).*exp(-Ebst_p./(EG/30));  
  disp(nq);
  
end 

toc
temp1_u=abs(fftshift(ifft(output_u))).^2;  
temp1_norm=mat2gray(temp1_u);
temp1u_re=10.*log10(temp1_u);

% figure(1); clf;
% plot(dtau,abs(output_u).^2,'LineWidth',2,'Color',[0.2 0.2 0.2]); hold on;
% set(gca,'Fontsize',20); set(gca,'LineWidth',2) 
% xlabel('Time(ps)'); ylabel('Intensity (W)');   
% 
% figure(2); clf; 
% plot(lambda_re,temp1u_re,'LineWidth',2,'Color',[0.8 0.2 0.2]); 
% set(gca,'Fontsize',20); set(gca,'LineWidth',2) 
% axis([1545 1575 -125 -25]); % set(gca,'ytick',(-75:25:-25));
% set(gca,'Layer','top');
% xlabel('Wavelength (nm)'); ylabel('Intensity(dB)');
% 
% %%
% figure(3); clf;
% plot(Eb_p(1:end-1,1),'LineWidth',2); hold on;
% set(gca,'LineWidth',1.5); set(gca,'Fontsize',15); 
% xlabel('Roundtrips'); ylabel('Energy(pJ)');
% 
% %%
% figure(4); clf;
output_train=mat2gray(reshape(output_trainS',[1,numel(output_trainS)]));
% t_train=ddtau*N_RR*(1e-6).*(1:20:numel(output_train));
% plot(output_train(1:20:end),'LineWidth',1.5,'Color',[0.2 0.2 0.6]);
% set(gca,'Fontsize',17); set(gca,'LineWidth',1.5); set(gca,'Layer','top');
% ylim([0 1.1]); set(gca,'layer','top'); % xlim([-14 13]); set(gca,'xtick',(-12:6:12))
% set(gca,'ytick',(0:0.5:1)); % set(gca,'xtick',(-2:2:2));
% xlabel('Time (\mus)'); ylabel('Norm. int.');
% 
% %%
gRec_pp=zeros(1,numel(gRec_p));
ERec_pp=zeros(1,numel(ERec_p));
for i=1:Nbst-1
    gRec_pp((1+N_nr*(i-1)):(N_nr*i))=gRec_p(:,i);
    ERec_pp((1+N_nr*(i-1)):(N_nr*i))=ERec_p(:,i);
end

%%
figure(5); clf;
subplot(211);
plot(output_train(1:20:end),'LineWidth',1.5,'Color',[0.2 0.2 0.6]);
set(gca,'Fontsize',17); set(gca,'LineWidth',1.5); set(gca,'Layer','top');
ylim([0 1.1]); set(gca,'layer','top'); % set(gca,'xtick',(-12:6:12))
% set(gca,'ytick',(0:0.5:1)); % set(gca,'xtick',(-2:2:2));
% xlabel('Time (\mus)'); ylabel('Norm. int.');
subplot(212);
plot((1:numel(gRec_p)),gRec_pp,'linewidth',2); hold on;
plot((1:N_nr:numel(gRec_p)),gRec_pp(1:N_nr:end),'.','Markersize',15);
set(gca,'linewidth',1.5); set(gca,'fontsize',17); ylabel('Gain (m^-1)');
gRec0=gRec_pp(1:N_nr:end); title(['g0=' num2str(gRec0(end))]);

figure(6); clf;
subplot(211);
plot(output_train(1:20:end),'LineWidth',1.5,'Color',[0.2 0.2 0.6]);
set(gca,'Fontsize',17); set(gca,'LineWidth',1.5); set(gca,'Layer','top');
ylim([0 1.1]); set(gca,'layer','top'); % set(gca,'xtick',(-12:6:12))
% set(gca,'ytick',(0:0.5:1)); % set(gca,'xtick',(-2:2:2));
% xlabel('Time (\mus)'); ylabel('Norm. int.');
subplot(212);
plot((1:numel(ERec_p)),ERec_pp,'linewidth',2); hold on;
plot((1:N_nr:numel(ERec_p)),ERec_pp(1:N_nr:end),'.','Markersize',15);
set(gca,'linewidth',1.5); set(gca,'fontsize',17); ylabel('Energy');

%% Temporal and spectral evolution
N_tot=size(outputT_trainS,1);
figure(6); clf;
pcolor(dtau,(1:10:N_tot)',mat2gray(outputT_trainS(1:10:N_tot,:))); shading interp;
set(gca,'LineWidth',2); set(gca,'Fontsize',20); set(gca,'Layer','top');
xlabel('Time (ps)'); ylabel('Roundtrips'); 
load 'mycolor'; colormap(mycolor);

%%
% figure(7); clf;
% pcolor((0:10:300),lambda_re',mat2gray(outputW_trainS(4500:10:4800,:)')); shading interp;
% set(gca,'LineWidth',2); set(gca,'Fontsize',16); set(gca,'Layer','top');
% ylabel('Wavelength (nm)'); xlabel('Roundtrips'); 
% load 'mycolor'; colormap(mycolor); caxis([0.1 1]);
% axis([0 300 1555 1565]); set(gca,'ytick',(1555:5:1565)); set(gca,'xtick',(0:150:300));
%% Width calculation
% N_start=2500; N_end=N_tot;
% % widT=zeros(1,(Nbst-1)*N_nr);
% widW=zeros(1,(N_end-N_start+1));
% widEp=zeros(1,(N_end-N_start+1));
% for i=1:(N_end-N_start+1)
%     Wtemp=mat2gray(outputW_trainS(N_start+i-1,:));
%     locW=(Wtemp>=0.1);  % 0.5
%     widW(i)=max(lambda_re(locW))-min(lambda_re(locW));
%     widEp(i)=trapz(dtau,outputT_trainS(N_start+i-1,:));
% end
% %%
% figure(8); clf;
% yyaxis left;
% plot(widW(118:188),'.-','linewidth',2,'markersize',15); hold on;
% xlim([0 70]); set(gca,'xtick',(0:70:70));
% set(gca,'linewidth',1.5); set(gca,'fontsize',15); ylabel('FWHM');
% yyaxis right;
% plot(widEp(end,118:188),'.-','linewidth',2,'markersize',15);
% xlim([0 70]); set(gca,'xtick',(0:70:70));
% set(gca,'linewidth',1.5); set(gca,'fontsize',15); 
% 
% figure(9); clf;
% plot(lambda_re+5,mat2gray(outputW_trainS(4617+33-1,:)),'linewidth',2,'Color',[0.5 0 0]); hold on;
% % plot(lambda_re,mat2gray(outputW_trainS(4617+18-1,:)),'linewidth',2); hold on;
% plot(lambda_re+5,mat2gray(outputW_trainS(4617+60-1,:)),'linewidth',2.5,'Color',[0.9 0.5 0.5]); hold on;
% axis([1560 1570 0 1.1]); set(gca,'linewidth',1.5); set(gca,'fontsize',15); 
% set(gca,'xtick',(1560:5:1570)); set(gca,'ytick',(0:0.5:1));

%%
% OS_ave=sum(outputW_trainS(N_start:N_end,:))./(N_end-N_start+1);
% figure(9); clf;
% plot(lambda_re,10.*log10(mat2gray(OS_ave)),'linewidth',2);
% axis([1520 1600 -100 0]);
% 
% lwf=[lambda_re',(10.*log10(mat2gray(OS_ave)))'];

%% Low-pass filter by applying Butterworth filter
% nt_tot=numel(output_train);
% omega_tot=[(0:nt_tot/2-1),(-nt_tot/2:-1)].*pi/(Tmax*nt_tot/nt);
% t_train=ddtau*N_RR*(1e-6).*(1:100:numel(output_train));
% output_Fil=abs(ifft(fft(output_train).*sqrt(1./(1+(omega_tot./0.13).^16))));
% output_wFil=abs(fftshift((fft(output_train).*sqrt(1./(1+(omega_tot./0.13).^16)))));
% 
% %
% figure(10); clf;
% % subplot(211);
% % plot(output_train(1:100:end),'LineWidth',1.5,'Color',[0.2 0.2 0.6]);
% % set(gca,'Fontsize',17); set(gca,'LineWidth',1.5); set(gca,'Layer','top');
% % ylim([0 1.1]); set(gca,'layer','top');
% % subplot(212);
% plot(t_train-12,mat2gray(output_Fil(1:100:end)),'LineWidth',1.5,'Color',[180/255 19/255 35/255]);
% set(gca,'Fontsize',17); set(gca,'LineWidth',1.5); set(gca,'Layer','top');
% ylim([0 1.2]); set(gca,'layer','top'); 
% set(gca,'ytick',(0:1.2:1.2))
% xlabel('Time (\mus)');
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
% % N_start=1;
% % N_evo1=4000; N_delay=2;
% % E_pN=mat2gray(E_p);
% % E_pN=E_p;
% % EE_evo1=E_pN(N_start:N_evo1,1); EE_evo2=E_pN((N_start+N_delay:N_evo1+N_delay),1);
% % EE_evo=E_p(:,1); num_evo=(1:1000)'; 
% Pk_real=max(outputT_trainS(N_start:N_end,:)');
% % Pk_evo=mat2gray(max(outputT_trainS(N_start:N_end,:)'));
% % Pk_evo1=Pk_evo(1:N_evo1);Pk_evo2=Pk_evo(1+N_delay:N_evo1+N_delay);
% figure(7); plot(Pk_real,'.-');
% filename1=['Data_' num2str(N_nr)];
% var_name=[(gRec0(end).*ones(size(Pk_real)))',Pk_real'];
% eval([filename1 '= var_name']);
% 
% % figure(6); clf;
% % % scatter(EE_evo,Pk_evo,10,num_evo,'filled'); % I
% % % scatter(EE_evo1,EE_evo2,10,(N_start:N_evo1),'filled'); % II
% % % scatter(Pk_evo1,Pk_evo2,10,(N_start:N_evo1),'filled'); % III
% % % scatter3(Pk_evo1,Pk_evo2,gRec_pp(N_start:(N_start+N_evo1-1)),20,(N_start:(N_start+N_evo1-1)),'filled'); hold on;
% % plot3(Pk_evo1,Pk_evo2,gRec_pp(N_start:(N_start+N_evo1-1)),'linewidth',3,'Color',[220/255 100/255 100/255]); view([-20 40]);
% % xlim([0 1]); ylim([0 1]); zlim([0 8]);
% % % set(gca,'xscale','log'); set(gca,'yscale','log')
% % set(gca,'Fontsize',20); set(gca,'LineWidth',2); grid on; box on; 

%% Energy evolution 
ERec_ppRe=ERec_pp(1:end); Npp=numel(ERec_ppRe);
figure(7); clf;
% plot(ERec_ppRe,'linewidth',2); hold on;
plot((1:N_nr:Npp),ERec_ppRe(1:N_nr:Npp),'linewidth',3);
set(gca,'linewidth',1.5); set(gca,'fontsize',17); ylabel('Energy');
title(['Energy=' num2str(ERec_ppRe(end)) 'pJ'])