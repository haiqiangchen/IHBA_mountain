function [d,path,position]= verify(bestpath,threat,R,startX,startY,endX,endY,gridCount)
%�˺�����Ҫ�Ǳܿ��״�ļ���Լ����㺽�ߵľ���
%path         input:�滮��·��
%threat       input:�״����в
%R            input:�״�İ뾶

% % % ����·��
position(1,1)=startX;position(1,2)=startY;
for i=2:gridCount
    position(i,1)=bestpath(i*2-1);
    position(i,2)=bestpath(i*2);
end
position(gridCount+1,1)=endX;position(gridCount+1,2)=endY;
% % % ������ۺ�������
[n,m]=size(position);
sign=ones(n);
sign=sign-diag(diag(sign));
 
[n,m]=size(sign);
cost=ones(size(sign))*Inf;                                
for i=1:n
    for j=1:n
            cost(i,j)=sqrt( sum (  ( position(i,:)-position(j,:) ).^2  ) );
    end
end
 
%ɽ����в�ж�
% �жϽڵ��Ƿ�λ��ɽ��������䷶Χ��,���Ĵ��۾���
[a,b]=find(cost~=Inf);
Allowed=[a,b];
for i=1:length(Allowed)
    for j=1:length(threat)
        x_1=position(Allowed(i,1),1);
        y_1=position(Allowed(i,1),2);  
        x_2=position(Allowed(i,2),1);
        y_2=position(Allowed(i,2),2);
          x_min=min(x_1,x_2);
          x_max=max(x_1,x_2);
          y_min=min(y_1,y_2);
          y_max=max(y_1,y_2);
        A=(y_2-y_1);
        B=-(x_2-x_1);
        C=y_1*(x_2-x_1)-x_1*(y_2-y_1);
        d=abs(threat(j,1)*A+threat(j,2)*B+C)/sqrt(A^2+B^2);
        C_bar=A*threat(j,2)-B*threat(j,1);
        X=(-A*C-B*C_bar)/(A^2+B^2);
        Y=(C_bar*A-B*C)/(A^2+B^2);
        %������ص���λ����û�д�����ɽ��
        d_1=sqrt((x_1-threat(j,1))^2+(y_1-threat(j,2))^2);
        d_2=sqrt((x_2-threat(j,1))^2+(y_2-threat(j,2))^2);  
        if (d<R(j)&&(X<=x_max&&X>=x_min)&&(Y<=y_max&&Y>=y_min))||d_1<=R(j)||d_2<=R(j)
            cost(Allowed(i,1),Allowed(i,2))=Inf;
        end
    end
end
u=1;
%������ź��̣������ԭ����Ҫ�Ǵ���������������ѡȡ���ʵĵ㣬Ȼ����֮��������ɵ����
dist=cost(1,:);
s=zeros(size(dist));
s(1)=1;
dist(1)=0;
path=zeros(size(dist));
path(1,:)=1;
for num=2:n
    mindist=Inf;
    for i=1:length(dist)
        if s(i)==0
            if dist(i)<mindist
                mindist=dist(i);
                u=i;
            end
        end
    end
    s(u)=1;
    for w=1:length(dist)
        if s(i)==0
            if dist(u)+cost(u,w)<dist(w)
                dist(w)=dist(u)+cost(u,w);
                path(w)=u;
            end
        end
    end
end
 
i=n;
d=0;
count=0;
while i>1
j=path(i);
%plot([position(i,1),position(j,1)],[position(i,2),position(j,2)],'Linewidth',2)
%axis([0,700,0,700]);
d=d+sqrt((position(i,1)-position(j,1))^2+(position(i,2)-position(j,2))^2);
i=j;
count=count+1;
end
%������ô����Ҫ�Ǳ������˻�������ֻ����������
if count==1
    d=Inf;
end
d;path;position;