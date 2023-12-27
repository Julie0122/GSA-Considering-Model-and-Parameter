function [ ca, V] = GRTM_Z1(K,Fs,h1,h2,k1,D,t)
% K,Fs,D,k1,t均为数值，ca为1*51的数组
% Boundary conditions  
    L = 10000;
    nx = 51;
    x = linspace(0,L,nx);%%%将1000m分为50段，200一段
    x1 = linspace(0,L/2,nx);%0-5000m 50段 每段100m
    c0(1)=100.00;
    c0(2)=0.0;
    c0(3)=0.0;
    c0(4)=0.0;
    c0(5)=0.0;
%     k1=0.05;
    k2=0.03;
    k3=0.01;
    k4=0.005;
    k5=0.002;
%     D=10;
%   t=1000;%%%时间
% groundwater flow model
% Z1
    for jj = 1:nx
        h(jj) = sqrt(h1^2 - (h1^2-h2^2)*x(jj)/L + Fs*(L-x(jj))*x(jj)/K);
        %%% x=0到10000间51个节点的水头
    end
    m = 1;
    for jj = 1:ceil(nx/2)
        nn = jj + ceil(nx/2)-1;
        q(jj) = (K*(h1^2-h2^2)/(2*L)-Fs*(L/2-x(nn)))/h(nn);%%% x=5000到10000间26个节点（第26个点-第51个点）的比流量
        if q(jj)>=0     %%%只保留正方向的流速
           
            v1(m)=q(jj);%%% x=5000到10000间26个节点（第26个点-第51个点）的流速
            m = m + 1;
        end
    end
    v11= harmmean(v1);
    V= 10*v11;%%%%
    
% reactive transport model    
    s = diag(ones(5,1));
    s(2,1)=-k1/(k1-k2);
    s(3,1)= k1*k2/((k1-k2)*(k1-k3));
    s(3,2)=  -k2/(k2-k3);
    s(4,1)= -k1*k2*k3/((k1-k2)*(k1-k3)*(k1-k4));
    s(4,2)= k2*k3/((k2-k3)*(k2-k4));
    s(4,3)= -k3/(k3-k4);
    s(5,1)= k1*k2*k3*k4/((k1-k2)*(k1-k3)*(k1-k4)*(k1-k5));
    s(5,2)=  -k2*k3*k4/((k2-k3)*(k2-k4)*(k2-k5));
    s(5,3)=  k3*k4/((k3-k4)*(k3-k5));
    s(5,4)=  -k4/(k4-k5);
    %Si matrix
    %             si(1:5,1:5)=0;
    %             for i4=1:5
    %                 si(i4,i4)=1;
    %             end

    si = diag(ones(5,1));

    si(2,1)=k1/(k1-k2);
    si(3,1)=k1*k2/((k1-k3)*(k2-k3));
    si(3,2)=   k2/(         k2-k3);
    si(4,1)=k1*k2*k3/((k1-k4)*(k2-k4)*(k3-k4));
    si(4,2)=   k2*k3/((        k2-k4)*(k3-k4));
    si(4,3)=      k3/(                 k3-k4);
    si(5,1)=k1*k2*k3*k4/((k1-k5)*(k2-k5)*(k3-k5)*(k4-k5));
    si(5,2)=   k2*k3*k4/(        (k2-k5)*(k3-k5)*(k4-k5));
    si(5,3)=      k3*k4/(                (k3-k5)*(k4-k5));
    si(5,4)=         k4/(                         k4-k5);


    %deomposed vector
    lamda=[k1 k2 k3 k4 k5];

    %analytical boundary conditions
    b0=si*c0';

    bet1=sqrt(V^2/(4*D^2)+lamda(1)/D); % beta defined in Bear's solution
    bet2=sqrt(V^2/(4*D^2)+lamda(2)/D); % beta defined in Bear's solution
    bet3=sqrt(V^2/(4*D^2)+lamda(3)/D); % beta defined in Bear's solution
    bet4=sqrt(V^2/(4*D^2)+lamda(4)/D); % beta defined in Bear's solution
    bet5=sqrt(V^2/(4*D^2)+lamda(5)/D); % beta defined in Bear's solution

    b1=.5*b0(1)*exp(x1*V/(2*D)).*(exp(-bet1.*x1).*erfc((x1-sqrt(V^2+4*lamda(1)*D)*t)/(2*sqrt(D*t)))+...
        exp( bet1.*x1).*erfc((x1+sqrt(V^2+4*lamda(1)*D)*t)/(2*sqrt(D*t)))) ;
    b2=.5*b0(2)*exp(x1*V/(2*D)).*(exp(-bet2.*x1).*erfc((x1-sqrt(V^2+4*lamda(2)*D)*t)/(2*sqrt(D*t)))+...
        exp( bet2.*x1).*erfc((x1+sqrt(V^2+4*lamda(2)*D)*t)/(2*sqrt(D*t)))) ;
    b3=.5*b0(3)*exp(x1*V/(2*D)).*(exp(-bet3.*x1).*erfc((x1-sqrt(V^2+4*lamda(3)*D)*t)/(2*sqrt(D*t)))+...
        exp( bet3.*x1).*erfc((x1+sqrt(V^2+4*lamda(3)*D)*t)/(2*sqrt(D*t)))) ;
    b4=.5*b0(4)*exp(x1*V/(2*D)).*(exp(-bet4.*x1).*erfc((x1-sqrt(V^2+4*lamda(4)*D)*t)/(2*sqrt(D*t)))+...
        exp( bet4.*x1).*erfc((x1+sqrt(V^2+4*lamda(4)*D)*t)/(2*sqrt(D*t)))) ;
    b5=.5*b0(5)*exp(x1*V/(2*D)).*(exp(-bet5.*x1).*erfc((x1-sqrt(V^2+4*lamda(5)*D)*t)/(2*sqrt(D*t)))+...
        exp( bet5.*x1).*erfc((x1+sqrt(V^2+4*lamda(5)*D)*t)/(2*sqrt(D*t)))) ;
    bc(1,:)=b1;
    bc(2,:)=b2;
    bc(3,:)=b3;
    bc(4,:)=b4;
    bc(5,:)=b5;

    ca=s*bc;
    ca=ca(5,:)';
    
end


