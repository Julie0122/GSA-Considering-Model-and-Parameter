function [SM1,erM1,Vi1,pa_SM,outputCa]=M1con(z0,h1,ko1,D,P,Ta,Tm,Csn,SVC,A,a,K,f1,N,nn,t,Nx)
%%%R1G1M1 3参

w = a.* (P - 355.6).^0.5 ./ 365000;  
M = f1 .* (Ta - Tm);
Q = Csn .* M .* SVC .* A .* 0.001 / 86400;
h2 = 0.3 .* Q.^0.6 + z0;
outputCa=zeros(51,N,N); outputCK=zeros(51,N,N); outputCf=zeros(51,N,N);
for i=1:3
 if i==1 
    for i1=1:N
      for  i2=1:N     
             [outputCa(:,i1,i2), ~] = GRTM_Z1(K(i2),w(i1),h1,h2(i2),ko1,D,t);    
      end
    end
    
    %SM
    EM1(:,1)=mean(reshape(outputCa(2:14,:,:),[13 N*N]),2 );
    VM1(:,1)=var(outputCa(2:14,:,:),0,[2 3] );
    rM1(:,1)=skewness(reshape(outputCa(2:14,:,:),[13 N*N]),0,2);
    kM1(:,1)=kurtosis(reshape(outputCa(2:14,:,:),[13 N*N]),0,2);
    %pa-ch
    Epa1a(:,1)=mean(abs(EM1-mean(outputCa(2:14,:,:),3)),2);
    Vpa1a(:,1)=mean(abs(VM1-var(outputCa(2:14,:,:),0,3)),2);    
    rpa1a(:,1)=mean(abs(rM1-skewness(outputCa(2:14,:,:),0,3)),2);
    kpa1a(:,1)=mean(abs(kM1-kurtosis(outputCa(2:14,:,:),0,3)),2);
    Suba(:,1) = mean(  var(outputCa(2:14,:,:),0,2),3 );  %mean(  var(outputCa(2:14,:,:),0,2),3 );
    Via(:,1)  = var(mean(outputCa(2:14,:,:),3),0,2);

 elseif i==2 
    for j1=1:N
      for  j2=1:N     
             [outputCK(:,j1,j2), ~] = GRTM_Z1(K(j1),w(j2),h1,h2(j2),ko1,D,t);    
      end
    end
    
    Epa1K(:,1)=mean(abs(EM1-mean(outputCK(2:14,:,:),3)),2);
    Vpa1K(:,1)=mean(abs(VM1-var(outputCK(2:14,:,:),0,3)),2);
    rpa1K(:,1)=mean(abs(rM1-skewness(outputCK(2:14,:,:),0,3)),2);
    kpa1K(:,1)=mean(abs(kM1-kurtosis(outputCK(2:14,:,:),0,3)),2);
    SubK(:,1) = mean(  var(outputCK(2:14,:,:),0,2),3 );
    ViK(:,1)  =  var(mean(outputCK(2:14,:,:),3),0,2);
    
 else  
    for g1=1:N
      for  g2=1:N     
             [outputCf(:,g1,g2), ~] = GRTM_Z1(K(g2),w(g2),h1,h2(g1),ko1,D,t);    
      end
    end

    Epa1f1(:,1)=mean(abs(EM1-mean(outputCf(2:14,:,:),3)),2);
    Vpa1f1(:,1)=mean(abs(VM1-var(outputCf(2:14,:,:),0,3)),2);
    rpa1f1(:,1)=mean(abs(rM1-skewness(outputCf(2:14,:,:),0,3)),2);
    kpa1f1(:,1)=mean(abs(kM1-kurtosis(outputCf(2:14,:,:),0,3)),2);
    Subf1(:,1)= mean(  var(outputCf(2:14,:,:),0,2),3 );
    Vif1(:,1) = var(mean(outputCf(2:14,:,:),3),0,2);

 end
end

%结果写入矩阵
      Vi1=zeros(13,3,2); %行是位置 列是参数参数 页是Si 2页是ST
      Vi1(:,1,1)=Via; Vi1(:,2,1)=ViK; Vi1(:,3,1)=Vif1;
      Vi1(:,1,2)=Suba;  Vi1(:,2,2)=SubK;  Vi1(:,3,2)=Subf1;
      SM1=zeros(13,4);
      SM1(:,1)=EM1; SM1(:,2)=VM1; SM1(:,3)=rM1; SM1(:,4)=kM1;
       pa_SM=zeros(13,4,3);% 行是位置  列EVrk 页是参数
       pa_SM(:,1,1)=Epa1a; pa_SM(:,2,1)=Vpa1a; pa_SM(:,3,1)=rpa1a; pa_SM(:,4,1)=kpa1a;
       pa_SM(:,1,2)=Epa1K; pa_SM(:,2,2)=Vpa1K; pa_SM(:,3,2)=rpa1K; pa_SM(:,4,2)=kpa1K;
       pa_SM(:,1,3)=Epa1f1; pa_SM(:,2,3)=Vpa1f1; pa_SM(:,3,3)=rpa1f1; pa_SM(:,4,3)=kpa1f1;
 
%%计算收敛情况 
 erM1=zeros(Nx,length(nn),6);
 [erM1(:,:,:) ]=shoulian(outputCa,nn,Nx);      
end