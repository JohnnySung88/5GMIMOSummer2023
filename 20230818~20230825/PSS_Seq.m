function d=PSS_Seq(N_id_2,PSSlength)

Pss_x=[0 1 1 0 1 1 1];

PssSeq=zeros(1,PSSlength);

PssSeq(1:7) = Pss_x;

for i=1:PSSlength-length(Pss_x)

    temp0=PssSeq(i);
    temp4=PssSeq(i+4);
    PssSeq(i+7)=mod(temp0+temp4,2);

end

n=0:126;

m=mod(n+43*N_id_2,127);

d(n+1)=1-2*PssSeq(m+1);

end