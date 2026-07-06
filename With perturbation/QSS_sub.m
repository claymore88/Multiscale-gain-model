function uout=QSS_sub(uu)

global dtau ddtau omega
global nt nq EG Tr Tl N_RR N_nr
global g0 g_p E_p delta_g g0_p gini
global output_u output_nl output_os

%% Fiber parameter %%
distan_1=0.022;            %length of the gain fiber, 4.6 GHz
step_num1=20;              %Constant step,20
step1=distan_1/step_num1;
beta_g=-10e-3;            %dispersion of gain fiber 
gama_g=3/1000;            %nonliearity of gain fiber 
BW=3*pi;                  %gain bandwidth, 3*pi for 24nm, 4pi for 32nm (3pi)
% g0=60/(EG/Tl);          %gain coefficient, 考虑除以2作为光场对应增益
g0=5;                     % 6 for Nnr=5-10;  

%%%%%%%%%%%%%%%%%%%%%%%%%

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

%% DF parameter
R_f=0.99;              % DF reflectance
beta_f=0e-6;           % DF dispersion
T_output=1-R_f;        % the output ratio
qq=ones(1,nt).*l0;     % saturable loss

%% Initial condition

% amp_temp=abs(wgn(1,nt1,1));  %White noise signal
% amp(1,1:nt/2-nt1/2)=0;       %
% amp(1,nt/2-nt1/2+1:nt/2+nt1/2)=(1).*amp_temp;
% amp(1,nt/2+nt1/2+1:nt)=0;   
% amp=0.1.*sech((dtau+20)./2);  %
% uu=amp;                       %initial condition

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
E_p=zeros(N_nr,2*step_num1);
g_p=zeros(N_nr,2*step_num1);
output_nl=zeros(N_nr,nt);
output_os=zeros(N_nr,nt);
% output_oslog=zeros(N_nr,nt);
    
%% Cavity simulation via QSS

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
           % g0_p(nq,n)=g0/(1+E_p(nr,n)/(EG*Tr/N_RR/Tl));
           % g_p(nr,n)=g0_p(nq,n);
           g_p(nr,n)=gini(1,n);
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
    uu=uu.*sqrt(TT_non).*exp(0.5*2*1i.*(qq-l0)./2);  % .*exp(0.5*2*1i.*(qq-l0)./2);
  
    %%%% GF 2 %%%
    for n=1+step_num1:step_num1+step_num1
        E_p(nr,n)=trapz(dtau,abs(uu).^2); 
        % g_p(nr,n)=g0./(1+E_p(nr)/Es);
        
        if nq==1 && nr==1
           %g_p(nr,n)=g_p(end,n);
           %g0_p(nq,n)=g0/(1+E_p(nr,n)/(EG*Tr/N_RR/Tl));
           %g_p(nr,n)=g0_p(nq,n);
           g_p(nr,n)=gini(1,n);
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
    % output_oslog(nr,:)=10.*log10(output_os(nr,:)); 
    
    uout=uu;
    %disp(nr);
    
  end