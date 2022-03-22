function [azimuth1,ellipticity1,azimuth2,ellipticity2] = Azim_Elipt(Length,Carrier_Density,G_TE,G_TM,signal,step)
variable
format long
delphicounter= 6;
global N_TE0 confin_TE dN_TEbar_dn N_TM0 confin_TM dN_TMbar_dn S1 S2
lambda= 1350*1e-9;  % signal wavelength in m

    lold=Length(1);
  for deZ = length(Length):-1:1




        NTE_Z(:,deZ)= N_TE0+confin_TE*(Carrier_Density(:,deZ))*dN_TEbar_dn;
            NTM_Z(:,deZ) = N_TM0+confin_TM*(Carrier_Density(:,deZ))*dN_TMbar_dn;


 Delta_Phi(deZ)=((2*pi)/lambda)*Length(deZ)*(((NTE_Z(delphicounter,deZ)-NTM_Z(delphicounter,deZ))));
               MulerMatrix=[max(G_TE(:,deZ)).^2+max(G_TM(:,deZ)).^2 max(G_TE(:,deZ)).^2-max(G_TM(:,deZ)).^2 0 0;
                max(G_TE(:,deZ)).^2-max(G_TM(:,deZ)).^2 max(G_TE(:,deZ)).^2+max(G_TM(:,deZ)).^2 0 0;
                0 0 (2*max(G_TE(:,deZ))*max(G_TM(:,deZ))*cos(max(Delta_Phi(:,deZ)))) (2*max(G_TE(:,deZ))*max(G_TM(:,deZ))*sin(max(Delta_Phi(:,deZ))));
                0 0 -(2*max(G_TE(:,deZ))*max(G_TM(:,deZ))*sin(max(Delta_Phi(:,deZ)))) (2*max(G_TE(:,deZ))*max(G_TM(:,deZ))*cos(max(Delta_Phi(:,deZ))))];
%              MulerMatrix=[ (G_TE(deZ)).^2+ (G_TM(deZ)).^2  (G_TE(deZ)).^2- (G_TM(deZ)).^2 0 0;
%                  (G_TE(deZ)).^2- (G_TM(deZ)).^2  (G_TE(deZ)).^2+ (G_TM(deZ)).^2 0 0;
%                 0 0 (2* (G_TE(deZ))* (G_TM(deZ))*cos( (Delta_Phi(deZ)))) (2* (G_TE(deZ))* (G_TM(deZ))*sin( (Delta_Phi(deZ))));
%                 0 0 -(2* (G_TE(deZ))* (G_TM(deZ))*sin( (Delta_Phi(deZ)))) (2* (G_TE(deZ))* (G_TM(deZ))*cos( (Delta_Phi(deZ))))];
             M= 0.5* MulerMatrix;
         S_prim2=M*S2';
          S_prim1=M*S1';
          %rad2deg
          
            azimuth2(deZ)=rad2deg((((S_prim2(3)*cos((Delta_Phi(deZ))))+(S_prim2(4)*sin((Delta_Phi(deZ)))) )/S_prim2(2)));
     ellipticity2(deZ)=rad2deg(((-S_prim2(3)*sin(Delta_Phi(deZ)))+(S_prim2(4)*cos(Delta_Phi(deZ))) )/S_prim2(1));
      azimuth1(deZ)=rad2deg((S_prim1(3)*cos(Delta_Phi(deZ)))+(S_prim1(4)*sin(Delta_Phi(deZ))) )/S_prim1(2);
     ellipticity1(deZ)=rad2deg((-S_prim1(3)*sin(Delta_Phi(deZ)))+(S_prim1(4)*cos(Delta_Phi(deZ))) )/S_prim1(1);

          

  end

signal=linspace(-10,10,step);
      figure(1)

     subplot(2,2,1)
     plot(signal,(ellipticity2));
      title("(a)-1");
     ylabel('Polarization Azimuth (degree)')
      xlabel('Pump Power (dBm)')
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       subplot(2,2,2)
       plot(signal,(azimuth2));
       y1 = yline(0,'--');
     y1.LabelHorizontalAlignment = 'center';
y1.Color = [.80 0 .40];
       title("(a)-2");
     ylabel('Ellepticity Angle (degree)')
      xlabel('Pump Power (dBm)')
     %%%%%%%%
     
        subplot(2,2,3)
         
     plot(signal,ellipticity1);
       title("(b)-1");
   ylabel('Polarization Azimuth (degree)')
      xlabel('Pump Power (dBm)')
      
      
       subplot(2,2,4)
       plot(signal,azimuth1);
     y2 = yline(0,'--');
     y2.LabelHorizontalAlignment = 'center';
        y2.Color = [.80 0 .40];
          title("(b)-2");
       ylabel('Ellepticity Angle (degree)')
      xlabel('Pump Power (dBm)')
      
          figure(4)
 
     plot(signal,Delta_Phi);
       title("\Delta\Phi");
 ylabel('Phase Shifted (degree)')
      xlabel('Pump Power (dBm)')

   end


