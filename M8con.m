function[SM4,erM4,Vi4,pa_SM4,outputC4a]=M8con(z0,h1,ko1,D,P,Ta,Tm,Csn,SVC,A,Rn,b,K1,K2,f2,r,N,nn,t,Nx)
%%R1G2M2 5参

w = b* (P - 399.8) ./ 365000; 
outputC4a=zeros(51,N,N); outputC4f=zeros(51,N,N); outputC4K1=zeros(51,N,N); outputC4K2=zeros(51,N,N);outputC4r=zeros(51,N,N);  
for i=1:5
    if i==1
     M = f2 .* (Ta - Tm) + r .* Rn;
     Q = Csn .* M .* SVC .* A .* 0.001 / 86400;
     h2= 0.3 .* Q.^0.6 + z0;
       for i1=1:N
          for  i2=1:N
                 [outputC4a(:,i1,i2), ~] = GRTM_Z2(K1(i2),K2(i2),w(i1),h1,h2(i2),ko1,D,t);
          end
       end
    %SM
    EM4(:,1)=mean(reshape(outputC4a(2:14,:,:),[13 N*N]),2 );
    VM4(:,1)=var(reshape(outputC4a(2:14,:,:),[13 N*N]),0,2 );
    rM4(:,1)=skewness(reshape(outputC4a(2:14,:,:),[13 N*N]),0,2);
    kM4(:,1)=kurtosis(reshape(outputC4a(2:14,:,:),[13 N*N]),0,2);
    %pa-ch
    Epa4a(:,1)=mean(abs(EM4-mean(outputC4a(2:14,:,:),3)),2);
    Vpa4a(:,1)=mean(abs(VM4-var(outputC4a(2:14,:,:),0,3)),2);
    rpa4a(:,1)=mean(abs(rM4-skewness(outputC4a(2:14,:,:),0,3)),2);
    kpa4a(:,1)=mean(abs(kM4-kurtosis(outputC4a(2:14,:,:),0,3)),2);
    Sub4a(:,1)=mean(  var(outputC4a(2:14,:,:),0,2),3 );  
    Vi4a(:,1) =  var(mean(outputC4a(2:14,:,:),3),0,2);  

    elseif i==2
     M = f2 .* (Ta - Tm) + r .* Rn;
     Q = Csn .* M .* SVC .* A .* 0.001 / 86400;
     h2= 0.3 .* Q.^0.6 + z0;
       for i1=1:N
          for  i2=1:N
                 [outputC4K1(:,i1,i2), ~] = GRTM_Z2(K1(i1),K2(i2),w(i2),h1,h2(i2),ko1,D,t);
          end
       end

    Epa4K1(:,1)=mean(abs(EM4-mean(outputC4K1(2:14,:,:),3)),2);
    Vpa4K1(:,1)=mean(abs(VM4-var(outputC4K1(2:14,:,:),0,3)),2);
    rpa4K1(:,1)=mean(abs(rM4-skewness(outputC4K1(2:14,:,:),0,3)),2);
    kpa4K1(:,1)=mean(abs(kM4-kurtosis(outputC4K1(2:14,:,:),0,3)),2);
    %sobol指数    
    Sub4K1(:,1)=mean(  var(outputC4K1(2:14,:,:),0,2),3 );
    Vi4K1(:,1) =  var(mean(outputC4K1(2:14,:,:),3),0,2);  
  
    elseif i==3
     M = f2 .* (Ta - Tm) + r .* Rn;
     Q = Csn .* M .* SVC .* A .* 0.001 / 86400;
     h2= 0.3 .* Q.^0.6 + z0;
     for i1=1:N
          for  i2=1:N
                 [outputC4K2(:,i1,i2), ~] = GRTM_Z2(K1(i2),K2(i1),w(i2),h1,h2(i2),ko1,D,t);
          end
     end

    Epa4K2(:,1)=mean(abs(EM4-mean(outputC4K2(2:14,:,:),3)),2);
    Vpa4K2(:,1)=mean(abs(VM4-var(outputC4K2(2:14,:,:),0,3)),2);
    rpa4K2(:,1)=mean(abs(rM4-skewness(outputC4K2(2:14,:,:),0,3)),2);
    kpa4K2(:,1)=mean(abs(kM4-kurtosis(outputC4K2(2:14,:,:),0,3)),2);
    %sobol指数    
    Sub4K2(:,1)= mean(  var(outputC4K2(2:14,:,:),0,2),3 ); 
    Vi4K2(:,1) =  var(mean(outputC4K2(2:14,:,:),3),0,2);

    elseif i==4
       for j1=1:N
          for  j2=1:N
             
                 M(j1,j2) = f2(j1) .* (Ta - Tm) + r(j2) .* Rn;
                 Q(j1,j2) = Csn .* M(j1,j2) .* SVC .* A .* 0.001 / 86400;
                 h2(j1,j2) = 0.3 .* Q(j1,j2).^0.6 + z0;
                 [outputC4f(:,j1,j2), ~] = GRTM_Z2(K1(j2),K2(j2),w(j2),h1,h2(j1,j2),ko1,D,t);
             
          end
       end
    %pa-ch
    Epa4f2(:,1)=mean(abs(EM4-mean(outputC4f(2:14,:,:),3)),2);
    Vpa4f2(:,1)=mean(abs(VM4-var(outputC4f(2:14,:,:),0,3)),2);
    rpa4f2(:,1)=mean(abs(rM4-skewness(outputC4f(2:14,:,:),0,3)),2);
    kpa4f2(:,1)=mean(abs(kM4-kurtosis(outputC4f(2:14,:,:),0,3)),2);
    Sub4f2(:,1)= mean(  var(outputC4f(2:14,:,:),0,2),3 );
    Vi4f2(:,1) =  var(mean(outputC4f(2:14,:,:),3),0,2);

    else
    for j1=1:N
          for  j2=1:N
             
             M(j1,j2) = f2(j2) .* (Ta - Tm) + r(j1) .* Rn;
             Q(j1,j2) = Csn .* M(j1,j2) .* SVC .* A .* 0.001 / 86400;
             h2(j1,j2) = 0.3 .* Q(j1,j2).^0.6 + z0;
             [outputC4r(:,j1,j2), ~] = GRTM_Z2(K1(j2),K2(j2),w(j2),h1,h2(j1,j2),ko1,D,t);
         
          end
     end

    Epa4r(:,1)=mean(abs(EM4-mean(outputC4r(2:14,:,:),3)),2);
    Vpa4r(:,1)=mean(abs(VM4-var(outputC4r(2:14,:,:),0,3)),2);
    rpa4r(:,1)=mean(abs(rM4-skewness(outputC4r(2:14,:,:),0,3)),2);
    kpa4r(:,1)=mean(abs(kM4-kurtosis(outputC4r(2:14,:,:),0,3)),2);
    Sub4r(:,1)= mean(  var(outputC4r(2:14,:,:),0,2),3 );
    Vi4r(:,1) =  var(mean(outputC4r(2:14,:,:),3),0,2);

    end
end

%结果写入矩阵
      Vi4=zeros(13,5,2); %行是参数 列是Si 2列是ST
      Vi4(:,1,1)=Vi4a; Vi4(:,2,1)=Vi4K1; Vi4(:,3,1)=Vi4K2; Vi4(:,4,1)=Vi4f2; Vi4(:,5,1)=Vi4r;
      Vi4(:,1,2)=Sub4a;  Vi4(:,2,2)=Sub4K1;  Vi4(:,3,2)=Sub4K2; Vi4(:,4,2)=Sub4f2; Vi4(:,5,2)=Sub4r; 
      SM4=zeros(13,4,1);
      SM4(:,1,1)=EM4; SM4(:,2,1)=VM4; SM4(:,3,1)=rM4; SM4(:,4,1)=kM4;
      pa_SM4=zeros(13,4,5);% 行是EVrk 列是参数
       pa_SM4(:,1,1)=Epa4a;  pa_SM4(:,2,1)=Vpa4a;  pa_SM4(:,3,1)=rpa4a;  pa_SM4(:,4,1)=kpa4a;
       pa_SM4(:,1,2)=Epa4K1; pa_SM4(:,2,2)=Vpa4K1; pa_SM4(:,3,2)=rpa4K1; pa_SM4(:,4,2)=kpa4K1;
       pa_SM4(:,1,3)=Epa4K2; pa_SM4(:,2,3)=Vpa4K2; pa_SM4(:,3,3)=rpa4K2; pa_SM4(:,4,3)=kpa4K2;
       pa_SM4(:,1,4)=Epa4f2;  pa_SM4(:,2,4)=Vpa4f2;  pa_SM4(:,3,4)=rpa4f2;  pa_SM4(:,4,4)=kpa4f2;
       pa_SM4(:,1,5)=Epa4r;  pa_SM4(:,2,5)=Vpa4r;  pa_SM4(:,3,5)=rpa4r;  pa_SM4(:,4,5)=kpa4r;

       erM4=zeros(13,length(nn),6);
       [erM4(:,:,:) ]=shoulian(outputC4a,nn,Nx);

end









