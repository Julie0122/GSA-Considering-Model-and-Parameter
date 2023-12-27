function[SM2,erM2,Vi2,pa_SM2,outputC2a]=M6con(z0,h1,ko1,D,P,Ta,Tm,Csn,SVC,A,Rn,b,K,f2,r,N,nn,t,Nx)
%%%R1G1M2 4个参数

w = b* (P - 399.8) ./ 365000;
outputC2a=zeros(51,N,N); outputC2K=zeros(51,N,N); outputC2f2=zeros(51,N,N); outputC2r=zeros(51,N,N);
for i=1:4
   if i==1
     M = f2 .* (Ta - Tm) + r .* Rn;
     Q = Csn .* M .* SVC .* A .* 0.001 / 86400;
     h2= 0.3 .* Q.^0.6 + z0;

        for i1=1:N
          for i2=1:N
                  [outputC2a(:,i1,i2), ~] = GRTM_Z1(K(i2),w(i1),h1,h2(i2),ko1,D,t);   
          end
        end
    %SM
    EM2(:,1)=mean(reshape(outputC2a(2:14,:,:),[13 N*N]),2 );
    VM2(:,1)=var(reshape(outputC2a(2:14,:,:),[13 N*N]),0,2 );
    rM2(:,1)=skewness(reshape(outputC2a(2:14,:,:),[13 N*N]),0,2);
    kM2(:,1)=kurtosis(reshape(outputC2a(2:14,:,:),[13 N*N]),0,2);
    %pa-ch
    Epa2a(:,1)=mean(abs(EM2-mean(outputC2a(2:14,:,:),3)),2);
    Vpa2a(:,1)=mean(abs(VM2-var(outputC2a(2:14,:,:),0,3)),2);
    rpa2a(:,1)=mean(abs(rM2-skewness(outputC2a(2:14,:,:),0,3)),2);
    kpa2a(:,1)=mean(abs(kM2-kurtosis(outputC2a(2:14,:,:),0,3)),2);
    Sub2a(:,1)= mean(  var(outputC2a(2:14,:,:),0,2),3 );
    Vi2a(:,1) =  var(mean(outputC2a(2:14,:,:),3),0,2);

   elseif i==2
      M = f2 .* (Ta - Tm) + r .* Rn;
     Q = Csn .* M .* SVC .* A .* 0.001 / 86400;
     h2= 0.3 .* Q.^0.6 + z0;
    for j1=1:N
          for j2=1:N
                  [outputC2K(:,j1,j2), ~] = GRTM_Z1(K(j1),w(j2),h1,h2(j2),ko1,D,t);   
          end
     end
    Epa2K(:,1)=mean(abs(EM2-mean(outputC2K(2:14,:,:),3)),2);
    Vpa2K(:,1)=mean(abs(VM2-var(outputC2K(2:14,:,:),0,3)),2);
    rpa2K(:,1)=mean(abs(rM2-skewness(outputC2K(2:14,:,:),0,3)),2);
    kpa2K(:,1)=mean(abs(kM2-kurtosis(outputC2K(2:14,:,:),0,3)),2);
    %sobol指数
    Sub2K(:,1)= mean(  var(outputC2K(2:14,:,:),0,2),3 );
    Vi2K(:,1)=  var(mean(outputC2K(2:14,:,:),3),0,2);
   
   elseif i==3

       for j1=1:N
          for j2=1:N
              
                 M(j1,j2) = f2(j1) .* (Ta - Tm) + r(j2) .* Rn;
                 Q(j1,j2) = Csn .* M(j1,j2) .* SVC .* A .* 0.001 / 86400;
                 h2(j1,j2) = 0.3 .* Q(j1,j2).^0.6 + z0;
                 [outputC2f2(:,j1,j2), ~] = GRTM_Z1(K(j2),w(j2),h1,h2(j1,j2),ko1,D,t);
              
          end
       end
    %pa-ch
    Epa2f2(:,1)=mean(abs(EM2-mean(outputC2f2(2:14,:,:),3)),2);
    Vpa2f2(:,1)=mean(abs(VM2-var(outputC2f2(2:14,:,:),0,3)),2);
    rpa2f2(:,1)=mean(abs(rM2-skewness(outputC2f2(2:14,:,:),0,3)),2);
    kpa2f2(:,1)=mean(abs(kM2-kurtosis(outputC2f2(2:14,:,:),0,3)),2);
    Sub2f2(:,1)= mean(  var(outputC2f2(2:14,:,:),0,2),3 );
    Vi2f2(:,1)= var(mean(outputC2f2(2:14,:,:),3),0,2);

   else
    for j1=1:N
          for j2=1:N
              
                 M(j1,j2) = f2(j2) .* (Ta - Tm) + r(j1) .* Rn;
                 Q(j1,j2) = Csn .* M(j1,j2) .* SVC .* A .* 0.001 / 86400;
                 h2(j1,j2) = 0.3 .* Q(j1,j2).^0.6 + z0;
                 [outputC2r(:,j1,j2), ~] = GRTM_Z1(K(j2),w(j2),h1,h2(j1,j2),ko1,D,t);
              
          end
     end
    Epa2r(:,1)=mean(abs(EM2-mean(outputC2r(2:14,:,:),3)),2);
    Vpa2r(:,1)=mean(abs(VM2-var(outputC2r(2:14,:,:),0,3)),2);
    rpa2r(:,1)=mean(abs(rM2-skewness(outputC2r(2:14,:,:),0,3)),2);
    kpa2r(:,1)=mean(abs(kM2-kurtosis(outputC2r(2:14,:,:),0,3)),2);
    Sub2r(:,1)= mean(  var(outputC2r(2:14,:,:),0,2),3 );
    Vi2r(:,1)=  var(mean(outputC2r(2:14,:,:),3),0,2);

   end
end

%结果写入矩阵
      Vi2=zeros(13,4,2); %行是参数 列是Si 2列是ST
      Vi2(:,1,1)=Vi2a; Vi2(:,2,1)=Vi2K; Vi2(:,3,1)=Vi2f2; Vi2(:,4,1)=Vi2r;
      Vi2(:,1,2)=Sub2a;  Vi2(:,2,2)=Sub2K;  Vi2(:,3,2)=Sub2f2; Vi2(:,4,2)=Sub2r;
      SM2=zeros(13,4,1);
      SM2(:,1,1)=EM2; SM2(:,2,1)=VM2; SM2(:,3,1)=rM2; SM2(:,4,1)=kM2;
      pa_SM2=zeros(13,4,4);% 列是EVrk 页是参数
       pa_SM2(:,1,1)=Epa2a; pa_SM2(:,2,1)=Vpa2a; pa_SM2(:,3,1)=rpa2a; pa_SM2(:,4,1)=kpa2a;
       pa_SM2(:,1,2)=Epa2K; pa_SM2(:,2,2)=Vpa2K; pa_SM2(:,3,2)=rpa2K; pa_SM2(:,4,2)=kpa2K;
       pa_SM2(:,1,3)=Epa2f2; pa_SM2(:,2,3)=Vpa2f2; pa_SM2(:,3,3)=rpa2f2; pa_SM2(:,4,3)=kpa2f2;
       pa_SM2(:,1,4)=Epa2r;  pa_SM2(:,2,4)=Vpa2r;  pa_SM2(:,3,4)=rpa2r;  pa_SM2(:,4,4)=kpa2r;

  erM2=zeros(13,length(nn),6);
 [erM2(:,:,:) ]=shoulian(outputC2a,nn,Nx); 
end
