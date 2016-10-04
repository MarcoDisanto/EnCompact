clear, clc, clf
div  = load('div50.dat');
grad = load('grad50.dat');
lapl = load('lapl50.dat');
n = size(div,1);
dD = size(div,2);
Div = spdiags(div,1-dD/2:dD/2,n,n+1);
dG = size(grad,2);
Grad = spdiags(grad,-dG/2:dG/2-1,n+1,n+10);
Grad(:,n+1:end) = [];
dL = size(lapl,2);
Lapl = spdiags(lapl,-(dL-1)/2:(dL-1)/2,n,n+10);
Lapl(:,n+1:end) = [];

spy(Div)
pause
spy(Grad)
pause
spy(Lapl)
pause
spy(Div*Grad-Lapl)

L = 31;
xf = linspace(0,L,n+1)';
xc = (xf(1:end-1) + xf(2:end))/2;
yc = sin(xc);
clf
plot(xc,yc,'r.-'), hold on
plot(xf,Grad*yc,'gs-')
plot(xc,Lapl*yc,'bd-')