clear all clc
%改进的蝙蝠算法在无人机航向上的规划应用第X次实验
        startX=0;startY=0;                            %起开始坐标                          
        endX=700;endY=700;                            %结束坐标
        gridCount=30;                                 %段点

%% 蝙蝠算法相关参数
        Qmax=0.5;                                     %最大频率
        Qmin=0;                                       %最小频率
        Rmax=1;
        pop=20;                                       %种群个数
        N_gen=50;                                     %迭代的次数
        c1=2;                                         %粒子群算法系数
        c2=2;                                         %粒子群算法系数
        V=zeros(pop,2*gridCount);                     %速度的初始化
        Q=zeros(pop,gridCount);                       %频率的初始化
        S=zeros(pop,2*gridCount);                     %速度的初始化
        fang=0;                                       %早熟因子
        Std=0;                                        %标准系数
        F0=0.5;                                       %变异因子
        CR=0.5;                                       %杂交参数
        w=0.8;                                        %惯性权重
        pathMax=700;                                  %边界最大值
        pathMin=0;                                    %边界最小值
        a=0.9;
        path_bar_best=[];
        position_bar_best=[];
        %%%%%%%%%%%%%%%%%%%%初始化确立%%%%%%%%%%%%%%%%%%%        
        for i=1:pop           
            f(i)=rand();
            fang=fang+f(i);
            for j=1:2*gridCount
                R(i,j)=rand();
                A(i,j)=rand();
            end
        end
        %%%%%%%%%%%%%%%%%%%%参考差分进化算法中对早熟机制的处理%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fang=fang/pop;
        for i=1:pop
             f1(i)=abs(f(i)-fang);
        end
        if max(f1(i))>1
            ff=max(f1(i));
        else
            ff=1;
        end
        Std=Std+((f(i)-fang)/ff);
        Std=Std/pop;
%% 生成山峰
          threat=[304 400 0;404 320 0;440 500 0;279 310 0;560 220 0;172 527 0;....
                194 220 0;272 522 0;350 200 0;....
                 650 400 0;740 250 0;540 375 0;510 600 0];
            r=[45 50 55 10 70 65 55 25 50 30 40 40 35];

        for i=1:length(r)
              figure(1)
              [x,y,z]=sphere;
              mesh(threat(i,1)+r(i)*x,threat(i,2)+r(i)*y,abs(threat(i,3)+r(i)*z));
              hold on
        end
        view([-30,-30,70])
%% 初始化粒子数
% for o=1:3
%     switch o
%         case 1
%             N_gen=100;
%         case 2
%             N_gen=200;
%         case 3
%             N_gen=300;
%     end
    tic
%     for u=1:25
                    for i=1:pop
                        for j=1:gridCount
                                X(i,j)=startX+j*(endX-startX)/(gridCount+1);
                                Y(i,j)=startY+rand()*(endY-startY);
                                path(i,2*j-1)=X(i,j);
                                path(i,2*j)=Y(i,j);
                        end
                    end
                    for i=1:pop
                         [distance,pathpoint,positionPoint]=verify(path(i,:),threat,....
                                         r,startX,startY,endX,endY,gridCount);
                        fitness(i)=distance;
                    end
                    [bestFitness,bestindex]=min(fitness);
                    bestpath=path(bestindex,:);
                    T=std(fitness); 
                    BestFitness=Inf;
                    globalFitness=Inf;
                    pathRecord=zeros(1,gridCount+1); bestRecord=zeros(1,gridCount+1);
                    position=zeros(gridCount+1,2);
            %%  迭代开始
                    for t=1:N_gen
                        for i=1:pop    
                          Q(i)=Qmin+(Qmin-Qmax)*rand();                      %蝙蝠算法的核心公式
                          V(i,:)=V(i,:)+(bestpath-path(i,:))*Q(i);
                          S(i,:)=path(i,:)+V(i,:);
                          S(i,find(S(i,:)>pathMax))=pathMax;     
                          S(i,find(S(i,:)<pathMin))=pathMin;
                          aa=6^(-30*(-t/N_gen)^4);
                          if  rand>R(i,:)                                  %step3 进行随机扰动
                               A1=mean(A(i,:));
                               S(i,:)=bestpath(:)+A1*aa;
                          end
                          [distance,pathpoint,positionPoint]=verify(bestpath,threat,....
                                 r,startX,startY,endX,endY,gridCount);
                          fnew=distance;
                          if rand<A(i)&&fitness(i)<fnew;
                              A(i)=a*A(i);
                              R(i,:)=Rmax*(1-exp(-0.9*t));
                          end
                             a=1;b=pop;
                             dx=randperm(b-a+1)+a-1;
                             j=dx(1);k=dx(2);p=dx(3);
                          if j==i
                             j=dx(4);
                          elseif k==i
                             k=dx(4);
                         elseif p==i
                             p=dx(4);
                          end
            %%  这里的步骤是混入粒子群的步骤
                         if Std<T                                          %判断是否进入早熟
                       %% change by haiqiang at 20170402 
                             named=exp((1-N_gen)/(N_gen+1-t));             %引入变异因子
                             F=F0*2.^named;
                       %% end by haiqiang
                             V(i,:)=w*S(i,:)+c1*rand*(bestpath-S(i,:))+c2*F*(path(k,:)-path(p,:));     %粒子群算法操作
                             if rand>CR
                                 path(i,:)=V(i,:);%选择操作
                             else
                                 path(i,:)=S(i,:);
                             end
                             [distance,pathpoint,positionPoint]=verify(path(i,:),threat,....
                                 r,startX,startY,endX,endY,gridCount);
                             fmin11=distance;
                             [distance1,pathpoint,positionPoint]=verify(S(i,:),threat,....
                                 r,startX,startY,endX,endY,gridCount);
                             fminS=distance1;
            %                  if CacuFit(S(i,:))>fmin11;
                             if fminS>fmin11
                                 path(i,:)=path(i,:);
                             else
                                 path(i,:)=S(i,:);
                             end
                        end  
            %%
                        [distance,pathpoint,positionPoint]=verify(path(i,:),threat,r,startX,startY,endX,endY,gridCount);
                          if BestFitness>distance;
                              bestpath=path(i,:);
                              BestFitness=distance;
                              pathRecord=pathpoint;  
                              position=positionPoint;
                          end

                    end
                     Fmin(t)=BestFitness;
%                      sum(u,t)=BestFitness;
                end
  
%                 fiftyTimesBest(u)=BestFitness;
%     end
      toc

      %%   迭代平均画图的数据
%               sum50=[];
%                 sum50(1,N_gen)=0;
%                         for I=1:N_gen
%                             for Tm=1:u
%                                 sum50(I)=sum50(I)+sum(Tm,I);
%                             end
%                             sumTest(I)=sum50(I)/u;
%                         end
%       Name=[num2str(N_gen),'维数据图'];
%       save(Name,'sumTest')
      %%  综合保存U次数据
%     yybest=min(fiftyTimesBest);
%     yyworst=max(fiftyTimesBest);
%     yymean=mean(fiftyTimesBest);
%     yystd=std(fiftyTimesBest);
%     yytime=toc/u;
%     
%     result=[];
%     result(1,1)=yybest;
%     result(2,1)=yyworst;
%     result(3,1)=yymean;
%     result(4,1)=yystd;
%     result(5,1)=yytime;
%     Name=[num2str(N_gen),'次迭代进化的数据'];
%     save(Name,'result');
            %% 画出实体图
                record=[];                                                              
                count=1;
                i=gridCount+1;
                while i>1
                j=pathRecord(i);
                record(count,1)=position(i,1); record(count,2)=position(i,2);
                count=count+1;
                plot([position(i,1),position(j,1)],[position(i,2),position(j,2)],'Linewidth',2)
                axis([0,700,0,700]);
                i=j;
                end
                record(count,1)=position(i,1); record(count,2)=position(i,2);

                text(position(1,1)',position(1,2)','S');
                text(position(gridCount+1,1)',position(gridCount+1,2)','T');
                figure(2)
                plot(Fmin);
                % title(['最佳个体适应度变化趋势,最佳适应值=' num2str(BestFitness)])
                title(['最后适应值 =' num2str(min(BestFitness))]);
                xlabel('迭代次数')
                ylabel('适应度值')
% end