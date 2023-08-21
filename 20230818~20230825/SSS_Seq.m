function d=SSS_Seq(N_id_1,N_id_2,SSSlength)

SSS_x=[1 0 0 0 0 0 0];

x0=zeros(1,SSSlength);
x1=zeros(1,SSSlength);

x0(1:7) = SSS_x;
x1(1:7) = SSS_x;

for i=1:SSSlength-length(SSS_x)

    temp0=x0(i);
    temp4=x0(i+4);
    x0(i+7)=mod(temp0+temp4,2);

end

for i=1:SSSlength-length(SSS_x)

    temp0=x1(i);
    temp1=x1(i+1);
    x1(i+7)=mod(temp0+temp1,2);

end

m0=3*floor(N_id_1/112)+N_id_2;

m1=mod(N_id_1,112)+m0+1;

n=0:126;

x0_ind=mod(n+m0,127);

x1_ind=mod(n+m1,127);

x_0=1-2*x0(x0_ind+1);

x_1=1-2*x1(x1_ind+1);

d(n+1)=x_0.*x_1;

end