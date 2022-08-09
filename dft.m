Matrix Method:


clc
clear
x = [1 1 0 3];
h = [1 1 8 0];
N=4;
for i =1 :N-1
y(:,i) = [x(N-i+1: N) x(1:N-i)];
end
z = [x' y]
y = z*h'



Tabular array method


clc
clear
x = [1 1 0 3];
h = [1 1 8 0];
N=4;
for i =1 :N-1
y(i) = [h(i:-1: 1) h(N:-1:i+1)]*x';
end
y(N) = h(N:-1:1)*x'


graphical method

clc
clear
x = [1 1 0 3];
h = [1 1 8 0];
N=4;
n = 0:1:N-1;
subplot(3 ,2,1);
stem(n,x);
for i =1:N-1
z(i,:) = [h(i:-1: 1) h(N:-1:i+1)];
subplot(3 , 2 ,i+1)
stem(n,z(i,:))
end
subplot(3,2,5);
stem(n, h(N:-1:1));
z = [z ; h(N:-1:1)];
y = z*x'
subplot(3,2,6)
stem(x,y)

linear convolution

clc
clear
x = [1 1 0 3];
h = [1 1 8 0];
N= length(x)+length(h)-1;
x = [x zeros(1,N-length(x))];
h = [h zeros(1,N-length(h))];
y= zeros(1,N);
for n=1:N
for k = 1:n
y(n) = y(n) + x(n)*h(n-k+1);
end
end
y






DFT

clc;
clear all;
N = input('Number of DFT points = ');
xn = input('Enter the sequence xn = '); 
ln = length(xn); %find the length of the sequence
xn = [xn zeros(1,N-ln)];
xk = zeros(1,N); 
ixk = zeros(1,N); 
%code block to find the DFT of the sequence
for k = 0:N-1
for n = 0:N-1
xk(k+1) = xk(k+1)+(xn(n+1)*exp((-1i)*2*pi*k*n/N));
end
end
t = 0:N-1;
subplot(2,2,1);
stem(t,xn);
ylabel('Amplitude ---->');
xlabel('n ---->');
title('Input Sequence ---->');
grid on;
magnitude = abs(xk); 
disp('DFT Sequence = ');
disp(magnitude);
%code block to plot the DFT sequence
t = 0:N-1;
subplot(2,2,2);
stem(t,magnitude);
ylabel('Amplitude ---->');
xlabel('K ---->');
title('DFT Sequence ---->');
grid on;
phase = angle(xk); 
disp('Phase = ');
disp(phase);
t = 0:N-1;
subplot(2,2,3);
stem(t,phase);
ylabel('Phase ---->');
xlabel('K ---->');
title('Phase Response');
grid on;


IDFT

for n = 0:N-1
for k = 0:N-1
ixk(n+1) = ixk(n+1)+(xk(k+1)*exp(1i*2*pi*k*n/N));
end
end
ixk = ixk./N;
%code block to plot the IDFT sequence
t = 0:N-1;
subplot(2,2,4);
stem(t,ixk);
disp('IDFT Sequence = ');
disp(ixk);
ylabel('Amplitude ---->');
xlabel('n ---->');
title('IDFT sequence ---->');
grid on;



Overall Save


clc;
clear all;
close all;
x = [2 -2 8 -2 -2 -3 -2 1 -1 9 1 3];
h= [1 2 3];
N=4
% Code to plot X(n)
subplot(2,2,1);
stem(x,'blue');
xlabel ('n---->');
ylabel ('Amplitude ---->');
title('X(n)');
%Code to plot H(n)
subplot(2,2,2);
stem(h,'black');
xlabel ('n---->');
ylabel ('Amplitude ---->');
title(' H(n)');
% Code to perform Convolution using Overlap Save Method
L=7; M=length(h);
N=L-M+1;
lx=length(x);
h=[h zeros(1,(L-M))];
x=[x zeros(1,(N-rem(lx,N)))];
x1=[zeros(1,(M-1)) x(1:N)];
for i=1:length(x)/N-1
  z(i,:)=x(i*N-1:(i+1)*N);
  y(i,:)= cconv(z(i,:),h,L);end
y=[cconv(x1,h,L);y];
P=y(:,(M:L));
P=P';
final=P(:)
% Representation of the Convoled Signal
subplot(2,2,3:4);
stem(final,'red');
xlabel ('n---->');
ylabel ('Amplitude ---->');
title('Convoled Signal');


OVERALL ADD:


clc
clear
x=[2 -2 8 -2 -2 -3 -2 1 -1 9 1 3];
h=[1 2 3];
L=7; M=length(h);
N=L-M+1;
lx=length(x);
h=[h zeros(1,(L-M))];
x=[x zeros(1,(N-rem(lx,N)))];
for i=1:length(x)/N
    z(i,:)=[x(((i-1)*N+1):i*N) zeros(1,(M-1))];
    y(i,:)=cconv(z(i,:),h,L);
end
for i=1:length(x)/N-1
  P(i,:)=[zeros(1,(i-1)*N) y(i,:) zeros(1,length(x)-L-(i-1)*N)];
end
P(i+1,:)=[zeros(1,i*N) y(i+1, (1:N))];
final=sum(P)




Sampling:

clc;
clear;
t = -10:0.1:10;
T=4;
fm=1/T;
x= cos(2*pi*fm*t);
subplot(2,2,1);
plot(t,x,'linewidth',3);
xlabel('time');
ylabel('amplitude');
grid;
title('input signal');
n1=(-10:1:10);
fs=[fm 2*fm 8*fm];
label1 = {'under sampling','uniform sampling' ,'over sampling' };
for i=1:length(fs)
z(i,:)= cos(2*pi*fm*n1/fs(i));
subplot(2,2,i+1);
stem(n1,z(i,:),'linewidth' ,3);
xlabel('number of samples');
ylabel('amplitude');
hold on;
subplot(2,2,i+1);
plot(n1,z(i,:),'linewidth',3);
xlabel('time');
ylabel('amplitude');
grid;
title(label1(i));
end

DTMF:

clc;
clear;
fs=4e3;
t=0:1/fs:0.5-1/fs;
ver=[697 770 852 941];
hor=[1209 1336 1477];
tones=[];
for k=1:4
    for m=1:3
        tone=sum(sin(2*pi*[ver(k);hor(m)].*t))';
        tones= [tones;tone;zeros(size(tone))];
    end
end
%soundsc(tones,fs)
S = timetable(seconds(0:length(tones)-1)'/fs ,tones);



FIR FILTER:

N=50;
fs=8000;
fc1=2000;
fc2=3000;
wc1= fc1/(fs/2);
wc2=fc2/(fs/2);
h = fir1(N,[wc1,wc2],'stop' ,hamming(N+1));
freqz(h,1,1000,fs)



low pass:

N=30;
fs=8000;
fc=2000;
wc= fc/(fs/2);
h = fir1(N,wc,'low' ,hamming(N+1));
freqz(h,1,1000,fs)





IIR filter:

clc; 
clear;
%wp = input("enter in rad/sec:")
%ws = input("enter rad/sec:")
%rp = input("enter ripple:")
%rs = input("enter ripple:")
fs = 1000;
fp1 = 200;
fp2 = 300;
fs1=100;
fs2=400;
wp1 = fp1/(fs/2);
wp2 = fp2/(fs/2);
ws1 = fs1/(fs/2);
ws2 = fs2/(fs/2);
kp =3;
ks = 60;

[N,wc] = buttord([wp1 wp2],[ws1 ws2],kp,ks);
[b,a] = butter(N ,wc , 'bandpass');
freqz(b,a,1000,fs);










#hardware:

DFT:
#include<stdio.h>
#include<math.h>
#define pi 3.1415926
int main(void) {
  int N = 4;
  int k, n;
  int x[4] = {
    1,
    1,
    0,
    3
  };
  int xi[4], xr[4];
  for (k = 0; k < N; k++) {
    xr[k] = 0;
    xi[k] = 0;
    for (n = 0; n < N; n++) {
      xr[k] += (x[n] * cos(2 * pi * k * n / N)) * 1000;
      xi[k] -= (x[n] * sin(2 * pi * k * n / N)) * 1000;
    }
    printf("x[%d]=%d+(%d)j\n", k, xr[k], xi[k]);
  }
}
IDFT: 
#include<stdio.h>
#include<math.h>
#define pi 3.1415926

int main(void) {
  int N = 4;
  int k, n;
  int Xr[] = {
    10,
    7,
    8,
    7
  }, Xi[] = {
    0,
    0,
    0,
    0
  };
  int xi[4], xr[4], x[4];
  for (n = 0; n < N; n++) {
    xr[n] = 0;
    xi[n] = 0;
    for (k = 0; k < N; k++) {
      xr[n] += Xr[k] * (cos(2 * pi * k * n / N)) - Xi[k] *
        (sin(2 * pi * k * n / N));
      xi[n] += Xr[k] * (sin(2 * pi * k * n / N)) + Xi[k] *
        (cos(2 * pi * k * n / N));
    }
    xr[n] = xr[n] / N;
    xi[n] = xi[n] / N;
    printf("x[%d]=%d+(%d)j\n", n, xr[n], xi[n]);
  }
}


Linear Convolution: 
#include<stdio.h>
int main(void) {
    int n, k, N1 = 4, N2 = 4;
    int N = N1 + N2 - 1;
    int temp;
    int x[] = {
      1,
      2,
      3,
      4,
      0,
      0,
      0
    }, h[] = {
      1,
      2,
      3,
      4,
      0,
      0,
      0
    }, y[7];
    for (k = 0; k < N; k++) {
      if (k < N / 2) {
        temp = h[k];
        h[k] = h[N - k - 1];
        h[N - k - 1] = temp;
      } else
        break;
    }
    for (n = 0; n < N; n++) {
      y[n] = 0;
      for (k = 0; k <= n; k++)
        y[n] = y[n] + (x[k] * h[N - 1 - n + k]);
      printf("y[%d] = %d\n", n, y[n]);
    }
    return 0;
  }
  
Circular Convolution:
  //Circular Convolution
  #include<stdio.h>

  int main(void) {
    int N = 4;
    int n, k;
    int x[] = {
      1,
      1,
      3,
      4
    }, h[] = {
      1,
      1,
      2,
      1
    }, y[4];
    for (n = 0; n < N; n++) {
      y[n] = 0;
      for (k = 0; k < N; k++) {
        if (n - k < 0)
          y[n] += x[k] * h[N + n - k];
        else
          y[n] += x[k] * h[n - k];
      }
      printf("y[%d] = %d\n", n, y[n]);
    }
    return 0;
  }
  
FIR FILTER: 
#include<stdio.h>
#include<math.h>
#define pi 3.1416
#define wc 2 * pi * 0.4

int main(void) {
  int tao = 10, n;
  int N = 21; //2 * tao + 1;
  //int wc = 2;
  int w[21], hd[21], h[21];
  for (n = 0; n < N; n++) {
    w[n] = 1000 * 0.5 * (1 - cos((2 * pi * n) / (N - 1)));
    if ((n - tao) != 0)
      hd[n] = 1000 * sin(wc * (n - tao)) / pi * (n - tao);
    else
      hd[n] = 1000 * wc / pi;
    h[n] = (hd[n] * w[n]);
    printf("w[%d]=%d, hd[%d]=%d, h[%d]=%d\n", n, w[n], n,
      hd[n], n, h[n]);
  }
  return 0;
}



IIR filter:

#include<stdio.h>
#include<math.h>
#define pi 3.1415926
int fs,T,t,x[256],y[256];
int i;
int main()
{
	fs=1;
	T=1/fs;
	for(i=0,t=0;i<25;i++,t+=T)
	{
		x[i]=(sin(2*pi*117.1875*t)+sin(2*pi*3906.25*t)+sin(2*pi*585.9375*t))*1000;
	}
	y[0]=0;
	for(i=1;i<25;i++)
	{
		y[i]=0.15835*y[i-1]+0.4208*(x[i]+x[i-1]);
		printf("y[%d]=%d\n",i,y[i]);
	}
	return 0;
}

