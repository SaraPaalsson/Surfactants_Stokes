function [X,Y] = GenericPolarGrid(xyi,xyo,Nj)

Ni=length(xyi(:,1));

X=zeros(Ni,Nj);
Y=X;

% r = linspace(0.01,1,Nj)';
% r = logspace(-4,0,Nj)';
r = logspace(-2.25,0,Nj)';

Fxi=fft(xyi(:,1));
Fxo=fft(xyo(:,1));

Fyi=fft(xyi(:,2));
Fyo=fft(xyo(:,2));

    for j=1:Nj
%     Fxj=Fxi*(Nj-j)/(Nj-1)+Fxo*(j-1)/(Nj-1);
%     Fyj=Fyi*(Nj-j)/(Nj-1)+Fyo*(j-1)/(Nj-1);
    Fxj=Fxi*(1-r(j))+Fxo*r(j);
    Fyj=Fyi*(1-r(j))+Fyo*r(j);
    X(:,j)=ifft(Fxj);
    Y(:,j)=ifft(Fyj);
    end

end