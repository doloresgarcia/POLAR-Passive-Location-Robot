function [angle_Estim, position_est,decisi,CORR]=estimation_realtime(f,p,N_PKT)
load('bp_60_filled.mat'); %BP measurements 
%N_BP number of high snr bam patterns used
    N_BP=60;
    CORR=cell(1,N_PKT);
    P_AMP=zeros(N_PKT,N_BP); %for possible RSS estimation of angle
    angle_Estim=zeros(1,N_PKT);
    position_est=zeros(1,N_PKT);
    decisi=zeros(1,N_PKT);
    index=zeros(N_PKT,64);
    load(['../NEWSCRIPTS/REF_AP1RX2_FRAMES_EXT.mat']);
    FRAME_DATA1=FRAME_DATA;  %%reference frame
    load(['../NEWSCRIPTS/AP1_RX2_1_FRAMES_EXT.mat']);
    FRAME_DATA2=FRAME_DATA;  %% all frames on batch 
    bes_pos=zeros(1,N_PKT);
    
    for ii=1:N_PKT
        in=ii;
        SNR=zeros(1,64);
        SNR1=zeros(1,64);
        CIR1=zeros(64,128); %CIRS for all BP's for one packet Reference
        CIR2=zeros(64,128); %CIRS for all BP's for one packet Moving object
        x=1:128;
        time_reflections=zeros(1,N_BP);
          for jj=1:64
              CIRx = zeros(128,1);
              if FRAME_DATA1{1,jj}.VALID_PKT == 1
                    SNR(jj)=FRAME_DATA1{1,jj}.SNR;
                    CIRx=FRAME_DATA1{1,jj}.CIR;
                     if max(abs(CIRx))>0.001
                        [pea,ind]= findpeaks(abs(CIRx));
                        which=pea>max(abs(CIRx))*0.18; %only main path if overall power is more than 20%
                        ind=ind(which);
                        maxIND=min(ind);
                        offset = 64-maxIND; %fixed R=64 and center all main paths.
                        CIR1(jj,:)=circshift(CIRx,offset);
                     end
              end      
              %do the same for moving data, plot CIR2 to observe
              %reflections
              CIRx = zeros(128,1);
              s=size(FRAME_DATA2);
              if s(1)>0
                  if FRAME_DATA2 {in,jj}.VALID_PKT == 1 
                        CIRx=FRAME_DATA2{in,jj}.CIR;
                        SNR1(jj)=FRAME_DATA1{1,jj}.SNR;
                        if max(abs(CIRx))>0.001
                            [pea,ind]= findpeaks(abs(CIRx));
                            which=pea>max(abs(CIRx))*0.18;
                            ind=ind(which);
                            maxIND=min(ind);
                            offset = 64-maxIND;
                            CIR2(jj,:)=circshift(CIRx,offset);
                        end
                  end
              end
          end
          SNR1(isnan(SNR1))=0;
          SNRr=real(SNR1); 
          [h,index(ii,:)]=sort(SNRr,'descend');
          CIR1=CIR1(index(ii,1:N_BP),:);
          CIR2=CIR2(index(ii,1:N_BP),:);
%plot to check reflections
%           figure();
%           subplot(2,1,1);
%           for i=1:40
%                 plot(abs(CIR1(i,:)))
%                 hold on
%           end
%           subplot(2,1,2);
%           for i=1:40
%                 plot(abs(CIR2(i,:)))
%                 hold on
%           end
          dif=abs(CIR2)-abs(CIR1);
          corr=zeros(1,N_BP);
          for z=1:N_BP
                f=abs(CIR1(z,65:81));
                y=abs(CIR2(z,65:81));
                w=2;
                k=1;
                A=zeros(1,8);
                A1=zeros(1,8);
                for i=1:8
                    x=f(k:k+w);
                    x1=y(k:k+w);
                    k=k+2;
                    A(i)=trapz(x);
                    A1(i)=trapz(x1);
                end 
                corr(1,z)=sum((A.*A1)/(norm(A)*norm(A1)));
          end
          
          CORR{1,ii}=corr;
          for z=1:N_BP
                 if (max(abs(CIR1(z,:)))>0.001) & (max(abs(CIR2(z,:)))>0.0010)
                        [pea1,ind1]= findpeaks(dif(z,:)./norm(dif(z,:))); %normalized difference function
                        thres=0.3;
                        which=pea1>thres;
                        ind1=ind1(which);
                        ind1=ind1(ind1>65); %the reflection has to happen after the main path (condition1)
                        ind1=ind1(ind1<80); %only first order reflections
                        if length(ind1)>0
                            maxIND=min(ind1);
                            time_reflections(z)=maxIND;
                        end
                 end
          end
          b=time_reflections;
          b(b==0)=[];
          if length(b)>0
              pos1=mode(b);
              if sum(b==pos1)==1
                  for z=1:length(b)
                      if b(z)-1==pos1 | b(z)+1==pos1  %smoothing similar reflections
                          b(z)=pos1;
                      end
                  end
                  pos1=mode(b);
              end
              % only estimating the reflections seen by many BP's
              decisi(ii)=sum(b==pos1);
              position_est(ii)=pos1; 
              P_AMP(ii,:)=abs(CIR2(:,pos1)); %power of the path at main reflection
              bes_pos(ii)=max(P_AMP(ii,:)./abs(CIR2(:,64).'));
          end
    end
    theta_BP=-160:2:159;
    %calculate the angles 
    for z=1:N_PKT
            XBP=bp_60_filled(index(z,1:N_BP),:);
            corr2=zeros(160,1);
            for i=1:160
                 corr2(i)=P_AMP(z,:)*XBP(1:N_BP,i)./(norm(P_AMP(z,:))*norm(XBP(1:N_BP,i)));
            end
            [k,ind2]=max(corr2(51:150,1));  %Only angles between -60 and 60 are possible because of our frontend
            angle_Estim(z)=theta_BP(ind2+50);
    end
   
end