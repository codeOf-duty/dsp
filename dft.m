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
lx=length(x);
lh=length(h);
m=lh-1;
x=[zeros(1,m) x zeros(1,N)];
h=[h zeros(1,N-lh)];
L=N-lh+1;
k=floor((length(h))/L);
for i=0:k
y=x(1,i*L+1:i*L+N);
q=mycirconv1(y,h); %Call the mycirconv1 function.
p(i+1,:)=q;
end
p1=p(:,lh:N)';
p=p1(:)'
% Representation of the Convoled Signal
subplot(2,2,3:4);
stem(p,'red');
xlabel ('n---->');
ylabel ('Amplitude ---->');
title('Convoled Signal');


OVERALL ADD:


clc;
clear all;
Xn = [2 -2 8 -2 -2 -3 -2 1 -1 9 1 3];
Hn = [1 2 3];
L=4;
% Code to plot X(n)
subplot (2,2,1);
stem(Xn);
xlabel ('n---->');
ylabel ('Amplitude ---->');
title(' X(n)');
%Code to plot H(n)
subplot (2,2,2);
stem(Hn,'red');
xlabel ('n---->');
ylabel ('Amplitude ---->');
title(' H(n)');
% Code to perform Convolution using Overlap Add Method
NXn=length(Xn);
M=length(Hn);
M1=M-1;
R=rem(NXn,L);
N=L+M1;
Xn=[Xn zeros(1,L-R)];
Hn=[Hn zeros(1,N-M)];
K=floor(NXn/L);
y=zeros(K+1,N);
z=zeros(1,M1);
for k=0:K
Xnp=Xn(L*k+1:L*k+L);
Xnk=[Xnp z];
y(k+1,:)=mycirconv1(Xnk,Hn); %Call the mycirconv function.
end
p=L+M1;
for i=1:K
y(i+1,1:M-1)=y(i,p-M1+1:p)+y(i+1,1:M-1);
end
z1=y(:,1:L)';
y=(z1(:))'
%Code to plot the Convolved Signal
subplot (2,2,3:4);
stem(y,'black');
xlabel ('n---->');
ylabel ('Amplitude ---->');
title('Convolved Signal');


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



