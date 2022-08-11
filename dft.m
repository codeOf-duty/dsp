% P controller
num=1;
den=[1 10 20];
plant=tf(num,den);
step(plant)

%p controller for closed loop

Kp=300;
contr=Kp; 
sys=feedback(contr*plant,1);
t=0:0.01:2;
step(sys,t)

%PD controller

Kp=300;
Kd=10;
contr=tf([Kd Kp],1);
sys=feedback(contr*plant,1);
t=0:0.01:2;
step(sys,t)

% Pi controller

Kp=30;
Ki=70;
contr=tf([Kp Ki],[1 0]); 
sys=feedback(contr*plant,1); 
t=0:0.01:2;
step(sys,t)


%PID contrpller

Kp=350; Ki=300; Kd=50;
contr=tf([Kd Kp Ki],[1 0]);
sys=feedback(contr*plant,1); 
t=0:0.01:2;
step(sys,t)



