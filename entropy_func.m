%Entropyº¯ÊýÐÎÊ½
function entr=entropy_func(p)
p=double(p);
[M,N]=size(p);
leng=M*N;
cal=zeros(1,256);
for i=1:leng
    tmp=p(i)+1;
    cal(1,tmp)=cal(1,tmp)+1;
end
entr=0;
for i=1:256
    pr(1,i)=cal(1,i)/leng;
    if pr(1,i)==0
        entr=entr;
    else
    entr=entr-pr(1,i)*log2(pr(1,i));
    end
end
digits(10);
vpa(entr)


