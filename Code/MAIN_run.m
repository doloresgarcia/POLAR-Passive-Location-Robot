N_PKT=48;
indexs=1;
g=1;
p=20;
CORR=cell(1,length(N_PKT));
angle_EST=zeros(length(indexs),N_PKT*g);
be_EST=zeros(length(indexs),N_PKT*g);
position_EST=zeros(length(indexs),N_PKT*g);
decisi=zeros(length(indexs),N_PKT*g);
snr=zeros(length(indexs),N_PKT*g);
Px=zeros(length(indexs),N_PKT*g);
for p=1:length(indexs)
    in=indexs(p);
    j=1;
%     corr=zeros(6,60);
    for f=1
        [angle_EST(p,j:j+N_PKT-1), position_EST(p,j:j+N_PKT-1),decisi(p,j:j+N_PKT-1),CORR]=estimation_realtime(f-1,in,N_PKT);  
        j=j+N_PKT;
    end

end

%This is a measure to check which packets to check without calculating
%relections.
l=zeros(1,N_PKT);
for i=1:N_PKT
    b=CORR{1,i}<0.6;
    c=CORR{1,i};
    c(b)=[];
    a=sort(c);
    l(i)=mean(a(1:p)); %check the 20 BP's for each packets with lowest corr 
                        %yet not so low they dont see the main path
end

%once reflections have been calculated, variable decisi is a more precise
%measure of when a frame shows an accurate value:
method='decisi'; %method='corr'
plot_results(angle_EST,position_EST,decisi,l,method) %results are plotted in units of floor tiles