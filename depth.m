close all
clear all
clc

data=load('Depth_particles');


%% %%%%%%%%%% select for each particle%%%%
for i=1:10
index=find(data(:,2)==i) ;
Z(:,i)=data(index,5);
end
%%
time_step=length(Z);
time_seconds=time_step*10;
days=time_seconds/(60*60*24);
x=linspace(0,days,time_step);
%%
figure(1)
for i=1:10
plot(x,Z(:,i))
hold on
end
%%

D_star=225273.484;
W_star=exp(-3.76715 + 1.92944 * log(D_star)-0.09815*log(D_star)^2 - 0.00575*log(D_star)^3+0.00056*log(D_star)^4)
