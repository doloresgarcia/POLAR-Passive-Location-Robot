function plot_results(angle_EST,position_EST,decisi,l,method)
    %%%%%%%%%%%%%%%%%%%%%%%  use decisi method
    if strcmp(method,'decisi')==1
        angle_EST3=(angle_EST(decisi>25))*-1;
        position_EST1=position_EST(decisi>25)-64;
    end
    %%%%%%%%%%%%%%%%%%%%%%% use correlation method
    if strcmp(method,'corr')==1
        angle_EST3=(angle_EST(l<0.72))*-1;
        position_EST1=position_EST(l<0.72)-64;
    end

    b=8.48; %distance between AP and RX
    x=zeros(1,length(angle_EST3));
    y=zeros(1,length(angle_EST3));
    %%%%%%%%%%%%%%%%%%%%%%%% simple geometrical approach to calculate position
    %%%%%%%%%%%%%%%%%%%%%%%% of object %%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:length(angle_EST3)
    %for i=1:length(angle_EST)
        d=position_EST1(i)*3e8/(1.76e9*0.4);
        alpha=angle_EST3(i);
        if alpha<0
             alpha=-angle_EST3(i)-45;
        else
            alpha=angle_EST3(i)+45;
        end
        dir=((d+b)^2-b^2)*0.5/(d+b-b*cos(deg2rad(alpha)));
        if alpha<0
            alpha=angle_EST3(i);
            x(i)=dir*sin(deg2rad(alpha));
            y(i)=15-dir*cos(deg2rad(alpha));
        else 
            alpha=angle_EST3(i);
            x(i)=+dir*sin(deg2rad(alpha));
            y(i)=15-dir*cos(deg2rad(alpha));
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xr=-5:0.5:3;
    yr=repmat(7,1,length(xr));

    figure();
    plot(x,y,'x')
    hold on
    for i=1:length(x)
        text(x(i),y(i),num2str(i))
    end 
    hold on
    text(-6,9,'RX')
    hold on
    text(0,15,'aP')
    hold on
    plot(xr,yr,'b')

end
