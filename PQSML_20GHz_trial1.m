clear 
% uu_record=uu;
% uu=uu_record;

close all

%%%%%%%%%%%%%%%%%%        Note           %%%%%%%%%%%%%%%%%
%%%% this programe is for a 11 GHz cavity at 1.5 um %%%%
%%%% for breathing simulation by HPX,YY,CXW    %%%%

%%%%%%%饱和能量(提高Es的值等效于加大泵浦)%%%%%%%%
% 

N_rep=1; N_sc=30;       % Scalability of repetition rate and strong correlation,50,30
N_RR=N_rep*N_sc;        % Repetition rate variation relative to 20GHz                
EG=3e6;                 % gain saturation energy, 3e6
Tr=50;                  % Roundtrip time
Tl=(1e-5)*1e12;         % Original version (0.8e-5)*1e12

N_nr=110;                % Number of the QSSs in burst, 70,150,180,200
Nbst=round(8e3/N_nr);                % Number of bursts, 70,25
Tb=N_nr*Tr*N_RR;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SESAM parameter %%
l0=0.03;     %depth of the SESAM 
F_sat=15;    %saturation fluence [uJ/cm2]
tao_a=10;    %relaxation time
R_a_de=0.93;      %reflectivity at designed wavelength
beta2=1025e-6;  beta3=0;  beta4=0;  beta5=0;    
 
diameter_core=8.0e-6;         %diameter of the fiber,
Acore=pi*diameter_core^2/4;   %area of the fiber 
Ps=F_sat*1e-6*(Acore*1e4)*1e12/tao_a;  %saturation power [W] 
Ea=Ps*tao_a/2;
R_sat_de=R_a_de+l0;       %reflectivity after saturation

c=299792458;
distan_1=0.005;           %length of the gain fiber, 0.92cm for ~11 GHz

%%%%%%%调节光纤参数%%%%%%%%
beta_g = -10e-3;            %dispersion of gain fiber 
gama_g = 3/1000;            %nonliearity of gain fiber 
BW = 3*pi;                  %gain bandwidth, 3*pi for 24nm, 4pi for 32nm (3pi)
% g0=60/(EG/Tl);            %gain coefficient, 考虑除以2作为光场对应增益
g0=85;                      % 85       
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%调节光路参数%%%%%%%%%
lambda_c=1563;             %central wavelength
step_num1=20;              %Constant step,20
step1=distan_1/step_num1;

kk=pi/(10*2*distan_1);     % birefringence ~0.2 for cavity length
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%划分时域频域网格 
nt=2048*2; Tmax=25;     %时间窗口 (2048*8,500)
nt1=512;                %白噪声带宽
ddtau=2*Tmax/nt;
%时域频域
dtau=(-nt/2:nt/2-1).*2*Tmax/nt;
omega=[(0:nt/2-1),(-nt/2:-1)].*pi/Tmax;
lambda_re=6*pi*1e5./(fftshift(omega+6*pi*1e5/lambda_c)); %wavelength matrix

%%%%%%%%%反射和透射率%%%%%%%%%
% film
R_f=0.99;    %assuming reflectivity of the film
beta_f=0e-6;    %dispersion of the film 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%循环圈数/输出耦合%%%%%%%%%%%%%
T_output=1-R_f;      % the output ratio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%输入信号%%%%%%%%%%%%%%%
amp_temp=abs(wgn(1,nt1,1));   %White noise signal
%%% Initial 1
% amp(1,1:nt/2-nt1/2)=0;        %
% amp(1,nt/2-nt1/2+1:nt/2+nt1/2)=(1).*amp_temp;
% amp(1,nt/2+nt1/2+1:nt)=0;   
%%% Initial 2
amp(1,1:nt)=0;
amp(1,1:nt1)=(1).*amp_temp;        %

uu=amp;                     %initial condition
% uu=sech((dtau+25)./5);        %initial condition
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
            
           % RK_s1=g_p(nr-1,n)*Ebst_p(1,n)./(EG*Tb)-g_p(nr-1,n)*E_p(nr-1,n)/(EG*Tr/N_RR);
           % RK_s2=(g_p(nr-1,n)+N_RR*Tr*RK_s1)*Ebst_p(1,n)./(EG*Tb)-(g_p(nr-1,n)+N_RR*Tr*RK_s1)*E_p(nr,n)/(EG*Tr/N_RR);
           % g_p(nr,n)=g_p(nr-1,n)+N_RR*Tr/2*(RK_s1+RK_s2); 
           
           RK_s1=-(g_p(nr-1,n)-g0)/Tl-g_p(nr-1,n)*E_p(nr-1,n)/(EG*Tr/N_RR);
           RK_s2=-((g_p(nr-1,n)+N_RR*Tr*RK_s1)-g0)/Tl-(g_p(nr-1,n)+N_RR*Tr*RK_s1)*E_p(nr,n)/(EG*Tr/N_RR);
           g_p(nr,n)=g_p(nr-1,n)+N_RR*Tr/2*(RK_s1+RK_s2); 
           
        else
            
           % RK_s1=g_p(nr-2,n)*Ebst_p(1,n)./(EG*Tb)-g_p(nr-2,n)*E_p(nr-2,n)/(EG*Tr/N_RR);
           % RK_s2=(g_p(nr-2,n)+N_RR*Tr*RK_s1)*Ebst_p(1,n)./(EG*Tb)-(g_p(nr-2,n)+N_RR*Tr*RK_s1)*E_p(nr-1,n)/(EG*Tr/N_RR);
           % RK_s3=(g_p(nr-2,n)+N_RR*Tr*RK_s2)*Ebst_p(1,n)./(EG*Tb)-(g_p(nr-2,n)+N_RR*Tr*RK_s2)*E_p(nr-1,n)/(EG*Tr/N_RR);
           % RK_s4=(g_p(nr-2,n)+2*N_RR*Tr*RK_s3)*Ebst_p(1,n)./(EG*Tb)-(g_p(nr-2,n)+2*N_RR*Tr*RK_s3)*E_p(nr,n)/(EG*Tr/N_RR);
           % g_p(nr,n)=g_p(nr-2,n)+N_RR*Tr/3*(RK_s1+2*RK_s2+2*RK_s3+RK_s4); 
           
           RK_s1=-(g_p(nr-2,n)-g0)/Tl-g_p(nr-2,n)*E_p(nr-2,n)/(EG*Tr/N_RR);
           RK_s2=-((g_p(nr-2,n)+N_RR*Tr*RK_s1)-g0)/Tl-(g_p(nr-2,n)+N_RR*Tr*RK_s1)*E_p(nr-1,n)/(EG*Tr/N_RR);
           RK_s3=-((g_p(nr-2,n)+N_RR*Tr*RK_s2)-g0)/Tl-(g_p(nr-2,n)+N_RR*Tr*RK_s2)*E_p(nr-1,n)/(EG*Tr/N_RR);
           RK_s4=-((g_p(nr-2,n)+2*N_RR*Tr*RK_s3)-g0)/Tl-(g_p(nr-2,n)+2*N_RR*Tr*RK_s3)*E_p(nr,n)/(EG*Tr/N_RR);
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
            
           % RK_s1=g_p(nr-1,n)*Ebst_p(1,n)./(EG*Tb)-g_p(nr-1,n)*E_p(nr-1,n)/(EG*Tr/N_RR);
           % RK_s2=(g_p(nr-1,n)+N_RR*Tr*RK_s1)*Ebst_p(1,n)./(EG*Tb)-(g_p(nr-1,n)+N_RR*Tr*RK_s1)*E_p(nr,n)/(EG*Tr/N_RR);
           % g_p(nr,n)=g_p(nr-1,n)+N_RR*Tr/2*(RK_s1+RK_s2); 
           
           RK_s1=-(g_p(nr-1,n)-g0)/Tl-g_p(nr-1,n)*E_p(nr-1,n)/(EG*Tr/N_RR);
           RK_s2=-((g_p(nr-1,n)+N_RR*Tr*RK_s1)-g0)/Tl-(g_p(nr-1,n)+N_RR*Tr*RK_s1)*E_p(nr,n)/(EG*Tr/N_RR);
           g_p(nr,n)=g_p(nr-1,n)+N_RR*Tr/2*(RK_s1+RK_s2);
           
        else
            
           % RK_s1=g_p(nr-2,n)*Ebst_p(1,n)./(EG*Tb)-g_p(nr-2,n)*E_p(nr-2,n)/(EG*Tr/N_RR);
           % RK_s2=(g_p(nr-2,n)+N_RR*Tr*RK_s1)*Ebst_p(1,n)./(EG*Tb)-(g_p(nr-2,n)+N_RR*Tr*RK_s1)*E_p(nr-1,n)/(EG*Tr/N_RR);
           % RK_s3=(g_p(nr-2,n)+N_RR*Tr*RK_s2)*Ebst_p(1,n)./(EG*Tb)-(g_p(nr-2,n)+N_RR*Tr*RK_s2)*E_p(nr-1,n)/(EG*Tr/N_RR);
           % RK_s4=(g_p(nr-2,n)+2*N_RR*Tr*RK_s3)*Ebst_p(1,n)./(EG*Tb)-(g_p(nr-2,n)+2*N_RR*Tr*RK_s3)*E_p(nr,n)/(EG*Tr/N_RR);
           % g_p(nr,n)=g_p(nr-2,n)+N_RR*Tr/3*(RK_s1+2*RK_s2+2*RK_s3+RK_s4); 
           
           RK_s1=-(g_p(nr-2,n)-g0)/Tl-g_p(nr-2,n)*E_p(nr-2,n)/(EG*Tr/N_RR);
           RK_s2=-((g_p(nr-2,n)+N_RR*Tr*RK_s1)-g0)/Tl-(g_p(nr-2,n)+N_RR*Tr*RK_s1)*E_p(nr-1,n)/(EG*Tr/N_RR);
           RK_s3=-((g_p(nr-2,n)+N_RR*Tr*RK_s2)-g0)/Tl-(g_p(nr-2,n)+N_RR*Tr*RK_s2)*E_p(nr-1,n)/(EG*Tr/N_RR);
           RK_s4=-((g_p(nr-2,n)+2*N_RR*Tr*RK_s3)-g0)/Tl-(g_p(nr-2,n)+2*N_RR*Tr*RK_s3)*E_p(nr,n)/(EG*Tr/N_RR);
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
  
  % delta_g=1.*(g_p(N_nr,:)-g0_p(nq+1,:)).*(1/Tl+Ebst_p(1,:)./(EG*Tb));  % 50
  
  disp(nq);
  
end 

toc
temp1_u=abs(fftshift(ifft(output_u))).^2;  
temp1_norm=mat2gray(temp1_u);
temp1u_re=10.*log10(temp1_u);

output_train=mat2gray(reshape(output_trainS',[1,numel(output_trainS)]));
t_train=ddtau*N_RR*(1e-6).*(1:20:numel(output_train));

%%

figure(1); clf;
subplot(121);
plot(dtau,abs(output_u).^2,'LineWidth',2,'Color',[0.2 0.2 0.2]); hold on;
set(gca,'Fontsize',15,'LineWidth',1.5);
xlabel('Time(ps)'); ylabel('Intensity (W)');   

subplot(122);
plot(lambda_re,temp1u_re,'LineWidth',2,'Color',[0.8 0.2 0.2]); 
set(gca,'Fontsize',15,'LineWidth',1.5);
axis([1545 1575 -125 -25]); % set(gca,'ytick',(-75:25:-25));
set(gca,'Layer','top');
xlabel('Wavelength (nm)'); ylabel('Intensity(dB)');

%%
% %%
% figure(3); clf;
% plot(Eb_p(1:end-1,1),'LineWidth',2); hold on;
% set(gca,'LineWidth',1.5); set(gca,'Fontsize',15); 
% xlabel('Roundtrips'); ylabel('Energy(pJ)');
% 

% %%
gRec_pp=zeros(1,numel(gRec_p));
ERec_pp=zeros(1,numel(ERec_p));
for i=1:Nbst-1
    gRec_pp((1+N_nr*(i-1)):(N_nr*i))=gRec_p(:,i);
    ERec_pp((1+N_nr*(i-1)):(N_nr*i))=ERec_p(:,i);
end

%%
figure(2); clf;
subplot(211);
plot(output_train(1:20:end),'LineWidth',1.5,'Color',[0.2 0.2 0.6]);
set(gca,'Fontsize',15,'LineWidth',1.5,'Layer','top');
ylim([0 1.1]); set(gca,'layer','top'); % set(gca,'xtick',(-12:6:12))
% set(gca,'ytick',(0:0.5:1)); % set(gca,'xtick',(-2:2:2));
% xlabel('Time (\mus)'); ylabel('Norm. int.');
subplot(212);
plot((1:numel(gRec_p)),gRec_pp,'linewidth',2); hold on;
plot((1:N_nr:numel(gRec_p)),gRec_pp(1:N_nr:end),'.','Markersize',15);
set(gca,'Fontsize',15,'LineWidth',1.5,'Layer','top'); ylabel('Gain (m^-1)');
gRec0=gRec_pp(1:N_nr:end); title(['g0=' num2str(gRec0(end))]);

%% Temporal and spectral evolution
N_tot=size(outputT_trainS,1);
figure(3); clf;
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
% plot(t_train-12,mat2gray(output_Fil(1:100:end))./0.75,'LineWidth',1.5,'Color',[180/255 19/255 35/255]);
% set(gca,'Fontsize',17); set(gca,'LineWidth',1.5); set(gca,'Layer','top');
% xlim([0 2]); ylim([0 1.2]); set(gca,'layer','top'); 
% set(gca,'xtick',(0:2:2)); set(gca,'ytick',(0:1.2:1.2))
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
N_start=N_tot-300; N_end=N_tot; 
% N_evo1=(N_end-N_start+1)-100;
% N_delay=5;

% N_start=1;
% N_evo1=4000; N_delay=2;
% E_pN=mat2gray(E_p);
% E_pN=E_p;
% EE_evo1=E_pN(N_start:N_evo1,1); EE_evo2=E_pN((N_start+N_delay:N_evo1+N_delay),1);
% EE_evo=E_p(:,1); num_evo=(1:1000)'; 
Pk_real=max(outputT_trainS(N_start:N_end,:)');
g0_real=[gRec0(end-5:end)';zeros(301-6,1)];
% Pk_evo=mat2gray(max(outputT_trainS(N_start:N_end,:)'));
% Pk_evo1=Pk_evo(1:N_evo1);Pk_evo2=Pk_evo(1+N_delay:N_evo1+N_delay);

figure(4); plot(Pk_real,'.-');
filename1=['Data_' num2str(N_nr)];
% var_name=[(gRec0(end).*ones(size(Pk_real)))',Pk_real'];
var_name=[g0_real,Pk_real'];
eval([filename1 '= var_name']);

% figure(6); clf;
% % scatter(EE_evo,Pk_evo,10,num_evo,'filled'); % I
% % scatter(EE_evo1,EE_evo2,10,(N_start:N_evo1),'filled'); % II
% % scatter(Pk_evo1,Pk_evo2,10,(N_start:N_evo1),'filled'); % III
% % scatter3(Pk_evo1,Pk_evo2,gRec_pp(N_start:(N_start+N_evo1-1)),20,(N_start:(N_start+N_evo1-1)),'filled'); hold on;
% plot3(Pk_evo1,Pk_evo2,gRec_pp(N_start:(N_start+N_evo1-1)),'linewidth',3,'Color',[220/255 100/255 100/255]); view([-20 40]);
% xlim([0 1]); ylim([0 1]); zlim([0 8]);
% % set(gca,'xscale','log'); set(gca,'yscale','log')
% set(gca,'Fontsize',20); set(gca,'LineWidth',2); grid on; box on; 

%%
figure(5); clf;
yyaxis left;
xRT=(1:N_nr:numel(gRec_p)); gRT=gRec_pp(1:N_nr:end); ERT=ERec_pp(1:N_nr:end);
plot(xRT(1:end),gRT(1:end),'.','Markersize',15);
set(gca,'linewidth',2,'fontsize',16); % xlim([0 8000]); ylim([0 5]);
ylabel('Gain');
yyaxis right;
plot(xRT(1:end),ERT(1:end),'linewidth',2);
set(gca,'linewidth',2,'fontsize',15); % xlim([0 8000]); ylim([0 110]);
xlabel('Roundtrips of QSS'); ylabel('Energy (pJ)'); 

%%
figure(6); clf;
plot(t_train-6,output_train(1:20:end),'LineWidth',1.5,'Color',[0.2 0.2 0.6]); 
set(gca,'Fontsize',16,'LineWidth',2,'layer','top'); set(gca,'Layer','top'); 
xlabel('Time (us)'); ylabel('Intensity (a.u.)');
xlim([0 inf]); ylim([0 1]); 
