function [azimuth1,ellipticity1,azimuth2,ellipticity2] = Azim_Elipt(Length,Carrier_Density,G_TE,G_TM,signal,step)
variable
format long

global N_TE0 confin_TE dN_TEbar_dn N_TM0 confin_TM dN_TMbar_dn S1 S2
lambda= 1350*1e-9;  % signal wavelength in m
% for Ipin=1:length(signal)
    lold=Length(1);
  for deZ = 1:1:length(Length)-1
%     deZ=end;

dz=Length(deZ)-lold;

        NTE_Z(:,deZ)= N_TE0+confin_TE*(Carrier_Density(:,deZ))*dN_TEbar_dn;
            NTM_Z(:,deZ) = N_TM0+confin_TM*(Carrier_Density(:,deZ))*dN_TMbar_dn;
          if deZ~= 1
            Delta_Phi(deZ-1)=((2*pi)/lambda)*sum(((NTE_Z(:,deZ)-NTM_Z(:,deZ)))*dz);
               MulerMatrix=[max(G_TE(:,deZ-1)).^2+max(G_TM(:,deZ-1)).^2 max(G_TE(:,deZ-1)).^2-max(G_TM(:,deZ-1)).^2 0 0;
                max(G_TE(:,deZ-1)).^2-max(G_TM(:,deZ-1)).^2 max(G_TE(:,deZ-1)).^2+max(G_TM(:,deZ-1)).^2 0 0;
                0 0 (2*max(G_TE(:,deZ-1))*max(G_TM(:,deZ-1))*cos(Delta_Phi(deZ-1))) (2*max(G_TE(:,deZ-1))*max(G_TM(:,deZ-1))*sin(Delta_Phi(deZ-1)));
                0 0 -(2*max(G_TE(:,deZ-1))*max(G_TM(:,deZ-1))*sin(Delta_Phi(deZ-1))) (2*max(G_TE(:,deZ-1))*max(G_TM(:,deZ-1))*cos(Delta_Phi(deZ-1)))];
             M= 0.5* MulerMatrix;
         S_prim2=M*S2';
          S_prim1=M*S1';
            azimuth2(deZ-1)=rad2deg(((S_prim2(3)*cos(Delta_Phi(deZ-1)))+(S_prim2(4)*sin(Delta_Phi(deZ-1))) )/S_prim2(2));
     ellipticity2(deZ-1)=rad2deg(((-S_prim2(3)*sin(Delta_Phi(deZ-1)))+(S_prim2(4)*cos(Delta_Phi(deZ-1))) )/S_prim2(1));
      azimuth1(deZ-1)=rad2deg((S_prim1(3)*cos(Delta_Phi(deZ-1)))+(S_prim1(4)*sin(Delta_Phi(deZ-1))) )/S_prim1(2);
     ellipticity1(deZ-1)=rad2deg((-S_prim1(3)*sin(Delta_Phi(deZ-1)))+(S_prim1(4)*cos(Delta_Phi(deZ-1))) )/S_prim1(1);
          end
          lold  =Length(deZ);
          

  end

signal=linspace(-10,10,step-2);
      figure(1)
     subplot(2,2,1)
%      Pin_dBm = linspace(-10,10,length(z)); 
     plot(signal,(azimuth1));
      title("azimuth 1");
       subplot(2,2,2)
      
     plot(signal,(ellipticity1));
     title("ellipticity 1");
     %%%%%%%%
        subplot(2,2,3)
     
     plot(signal,azimuth2);
      title("azimuth 2");
       subplot(2,2,4)
    
     plot(signal,ellipticity2);
       title("ellipticity 2");
          figure(3)
    
     plot(signal,Delta_Phi);
       title("delta Phi");


   end


