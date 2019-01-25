function H = MVAR(X,f1,f2)

            load('zef_MEG_measurements');
            load('zef_MEG_sensors');
            load('zef_MEG_V');
            load('zef_MEG_S');
            load('zef_MEG_L');
            load('zef_MEG_G');
            X = zef_MEG_measurements;
            V = zef_MEG_V;
            G = zef_MEG_G;
 
            L = zef_MEG_L;
            S = zef_MEG_S;
            A = corr(X);
            Af = fft(A);
            Xf = fft(X); 
            H = V'*(S-L)*G';

          [r1,c1]=find(H>7.5e06);
          [r2,c2]=find(H<7.5e06);
          plot(zef_MEG_sensors(:,1),zef_MEG_sensors(:,2),'o',zef_MEG_sensors(c1,1),zef_MEG_sensors(c1,2),'r-', zef_MEG_sensors(c2,1),zef_MEG_sensors(c2,2),'bo-');
          title('sensor coherence')
          

end

