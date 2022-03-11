function [dn,G_larg] =CarrierDensity(L,n,signal,step)
variable
format long
global h c lambda J dz R_p n_p a1 a2 a3 alpha_1 alpha_2 eta q d ...
    lambda0  g Iav E g_mp S_sp 

%2.9e24;%2.244523e+24;
lambda= 1350*1e-9;  % signal wavelength in m
E=(h*c)/(lambda); % signal energy   refrense dar https://www.youtube.com/watch?v=3T8T7u2-aVY
%E=(E)*(1/(1.6e-19)); % بر مبنای الکترون ولت
J=(17e-3)/((1000e-6))^3; % current density بزرگتر از این مقدار
%  dz = (L*1e-6)/Nz; % spatial step
 %%%%%%%%%%%%%%%%%%%%%%
 dn(1)=n;
 G_larg(1)=0;
  eps=n;
  countr=1;
  dnold=n;
for I =2:step*3
  
alpha_n=alpha_1+alpha_2*n; %equation (5)// با توجه به بررسی ها می بایست عدد به دست امده در لحظه اول زیر یازده هزار بدست بیاد


   %%%%%%%%%%%%%%%%%%%%%%%
   
tau = (((A+(B*n)+((C)*(n.^2)))).^-1);%equation (8) /// درسته
    
%%%%%%%%%%%%%%%%%%%from refrence 14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% R_p=db((A*n)+(B*(n.^2))+(C*(n.^3)))   ;   %equation (7) from Ref. 14   db 
 R_p=(A*n)+(B*(n.^2))+(C*(n.^3))  ;   %equation (7) from Ref. 14   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  lambda_p=lambda-a3*(n-n_p);%equation (3) ///درسته
    
      %material gain
% if (n-n0) ~=0
     
   g_mp=a1*(n-n0)-(a2*((lambda-lambda_p).^2));%equation (2)
 
% else
%     
%   g_mp=(a1*n0)-(a2*((lambda-lambda_p).^2));%equation (2)
%       
% end


  g=(confin_optional*(g_mp-alpha_n));%equation (4)


  G_s=exp(((confin_optional*g_mp-alpha_n)*(L)));%% single pass gain        equation (7)

% G_s=exp(confin_optional*g_mp-alpha_n*L*1e-6);%% single pass gain        equation (7)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%equation (10)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 F1=((beta_0*R_p)/(2*((confin_optional*g_mp)-alpha_n)));
  F2=((R1+R2).*((G_s-1).^2)+(2.*(1-R1*R2.*G_s).*(G_s-1)));
  F3=(L).*((confin_optional.*g_mp)-alpha_n).*(1-(R1.*R2).*(G_s.^2));
  S_sp=F1.*((F2./F3)-2);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%equation (10)%%%%%%%%%%%%%%%%%%%%%%%%%%

Pin_prob=1e-3*10.^(3/10);
Ain= sqrt(Pin_prob)./sqrt(E);

 Pump = sqrt(signal)./sqrt(E);
I_input=(Ain+Pump);
     Iav=I_input*((G_s-1)/(g*(L)));


    out=(((eta*J)/(q*d))-((n)/(tau))-((g*Iav)/E)-((g_mp)*(S_sp)));
    chek=n-(out/(-1/tau));
 
% fprintf("Step = %d , n root = %0.10f , eps =%d\n",I,n,b)

%  b=abs((chek- dn(I-1)));
% err=(abs(chek-n)/chek)*100;
% if chek < dn(I-1)
%     plot(linspace(1000,1,length(dn)),dn);
% hold on;
% 
%   break;
% 
%    end 

 n=n-(out/(-1/tau));
 if(n>n0)
 dn(countr)=n;
 G_larg(countr)=G_s;
 countr=countr+1;
 
 end

     
       
end
  
%    
       