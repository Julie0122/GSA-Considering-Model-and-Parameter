function[SM3,erM3,Vi3,pa_SM3,outputC3a]=M3con(z0,h1,ko1,D,P,Ta,Tm,Csn,SVC,A,a,K1,K2,f1,N,nn,t,Nx)
%%%R1G2M1  4参
 
w = a.* (P - 355.6).^0.5 ./ 365000; 
M = f1 .* (Ta - Tm);
Q = Csn .* M .* SVC .* A .* 0.001 / 86400;
h2 = 0.3 .* Q.^0.6 + z0;
outputC3a=zeros(51,N,N); outputC3K1=zeros(51,N,N); outputC3K2=zeros(51,N,N); outputC3f=zeros(51,N,N);
for i=1:4
    if i==1
        for i1=1:N
          for  i2=1:N
            
                 [outputC3a(:,i1,i2), ~] = GRTM_Z2(K1(i2),K2(i2),w(i1),h1,h2(i2),ko1,D,t);
             
          end
        end
    %SM
    EM3(:,1)=mean(reshape(outputC3a(2:14,:,:),[13 N*N]),2 );
    VM3(:,1)=var(reshape(outputC3a(2:14,:,:),[13 N*N]),0,2 );
    rM3(:,1)=skewness(reshape(outputC3a(2:14,:,:),[13 N*N]),0,2);
    kM3(:,1)=kurtosis(reshape(outputC3a(2:14,:,:),[13 N*N]),0,2);
    %pa-ch
    Epa3a(:,1)=mean(abs(EM3-mean(outputC3a(2:14,:,:),3)),2);
    Vpa3a(:,1)=mean(abs(VM3-var(outputC3a(2:14,:,:),0,3)),2);
    rpa3a(:,1)=mean(abs(rM3-skewness(outputC3a(2:14,:,:),0,3)),2);
    kpa3a(:,1)=mean(abs(kM3-kurtosis(outputC3a(2:14,:,:),0,3)),2);
    Sub3a(:,1)= mean(  var(outputC3a(2:14,:,:),0,2),3 );
    Vi3a(:,1) =  var(mean(outputC3a(2:14,:,:),3),0,2);

    elseif i==2
    for i1=1:N
          for  i2=1:N
            
                 [outputC3K1(:,i1,i2), ~] = GRTM_Z2(K1(i1),K2(i2),w(i2),h1,h2(i2),ko1,D,t);
             
          end
    end
    Epa3K1(:,1)=mean(abs(EM3-mean(outputC3K1(2:14,:,:),3)),2);
    Vpa3K1(:,1)=mean(abs(VM3-var(outputC3K1(2:14,:,:),0,3)),2);
    rpa3K1(:,1)=mean(abs(rM3-skewness(outputC3K1(2:14,:,:),0,3)),2);
    kpa3K1(:,1)=mean(abs(kM3-kurtosis(outputC3K1(2:14,:,:),0,3)),2);
    %sobol指数
    Sub3K1(:,1)= mean(  var(outputC3K1(2:14,:,:),0,2),3 );  
    Vi3K1(:,1) =  var(mean(outputC3K1(2:14,:,:),3),0,2); 

    elseif i==3

        for j1=1:N
          for  j2=1:N
            
                 [outputC3K2(:,j1,j2), ~] = GRTM_Z2(K1(j2),K2(j1),w(j2),h1,h2(j2),ko1,D,t);
           end
        end
        
    %pa-ch
    Epa3K2(:,1)=mean(abs(EM3-mean(outputC3K2(2:14,:,:),3)),2);
    Vpa3K2(:,1)=mean(abs(VM3-var(outputC3K2(2:14,:,:),0,3)),2);
    rpa3K2(:,1)=mean(abs(rM3-skewness(outputC3K2(2:14,:,:),0,3)),2);
    kpa3K2(:,1)=mean(abs(kM3-kurtosis(outputC3K2(2:14,:,:),0,3)),2);
    Sub3K2(:,1)= mean(  var(outputC3K2(2:14,:,:),0,2),3 );
    Vi3K2(:,1) =  var(mean(outputC3K2(2:14,:,:),3),0,2);

    else
   
          for  j2=1:N
             for j3=1:N
                 [outputC3f(:,j2,j3), ~] = GRTM_Z2(K1(j3),K2(j3),w(j3),h1,h2(j2),ko1,D,t);
             end
          end
    Epa3f1(:,1)=mean(abs(EM3-mean(outputC3f(2:14,:,:),3)),2);
    Vpa3f1(:,1)=mean(abs(VM3-var(outputC3f(2:14,:,:),0,3)),2);
    rpa3f1(:,1)=mean(abs(rM3-skewness(outputC3f(2:14,:,:),0,3)),2);
    kpa3f1(:,1)=mean(abs(kM3-kurtosis(outputC3f(2:14,:,:),0,3)),2);
    %sobol指数
    Sub3f1(:,1)= mean(  var(outputC3f(2:14,:,:),0,2),3 );
    Vi3f1(:,1) =  var(mean(outputC3f(2:14,:,:),3),0,2);

    end
end

%结果写入矩阵
      Vi3=zeros(13,4,2); %行是参数 列是Si 2列是ST
      Vi3(:,1,1)=Vi3a; Vi3(:,2,1)=Vi3K1; Vi3(:,3,1)=Vi3K2; Vi3(:,4,1)=Vi3f1;
      Vi3(:,1,2)=Sub3a;  Vi3(:,2,2)=Sub3K1;  Vi3(:,3,2)=Sub3K2; Vi3(:,4,2)=Sub3f1;
      SM3=zeros(13,4,1);
      SM3(:,1,1)=EM3; SM3(:,2,1)=VM3; SM3(:,3,1)=rM3; SM3(:,4,1)=kM3;
      pa_SM3=zeros(13,4,4);% 行是EVrk 列是参数
       pa_SM3(:,1,1)=Epa3a; pa_SM3(:,2,1)=Vpa3a; pa_SM3(:,3,1)=rpa3a; pa_SM3(:,4,1)=kpa3a;
       pa_SM3(:,1,2)=Epa3K1; pa_SM3(:,2,2)=Vpa3K1; pa_SM3(:,3,2)=rpa3K1; pa_SM3(:,4,2)=kpa3K1;
       pa_SM3(:,1,3)=Epa3K2; pa_SM3(:,2,3)=Vpa3K2; pa_SM3(:,3,3)=rpa3K2; pa_SM3(:,4,3)=kpa3K2;
       pa_SM3(:,1,4)=Epa3f1;  pa_SM3(:,2,4)=Vpa3f1;  pa_SM3(:,3,4)=rpa3f1;  pa_SM3(:,4,4)=kpa3f1;

 erM3=zeros(13,length(nn),6);
[erM3(:,:,:) ]=shoulian(outputC3a,nn,Nx);
end


 


