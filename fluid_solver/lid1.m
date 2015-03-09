clear all

%Problem Parameters
Uw=1;%velocity of moving lid (m/s)
W=1;%Width of cavity (m)
H=1;%Height of cavity (m)
Re=100;%Reynolds No

%Convergence array
convergence=zeros(1,2000);
k=0;
%Grid Properties
np=80;
nu=np-1;
nv=np-1;

% x is alligned along the width (index i) origin at bottom left
dx=1/(np-1);
% y is alligned along the height (index i) origin at bottom left
dy=1/(np-1);
%artificial time
dt=0.01;

%Initialization
Pstr=zeros(np,np);
Ustr=zeros(np,nu);
Vstr=zeros(nv,np);
%A(j,i)
Ustr(np,:)=Uw;
qsum=1;
%------------------------------------------------------------------------------
while sqrt(qsum)>10^-3

  k=k+1;
  U=zeros(np,nu);
  U(np,:)=Uw; % top boundary moving, all others at 0
  V=zeros(nv,np);
  relax=0.7;
  %Calculation of U field from Pstr and Ustr and Vstr

  residual1 = 1;
  residual2 = 1;

  while (residual1>10^(-8))
    if(residual1>residual2&&residual1>100)
      display('U cal failed')
      display(residual1)
      break;
    end
    residual2 = residual1;
    for i= 2:nu-1
      for j=2:np-1
        Vnstr= (Vstr(j,i+1)+Vstr(j,i))/2;
        Vsstr=(Vstr(j-1,i+1)+Vstr(j-1,i))/2;
        U(j,i)=(1-relax)*U(j,i)+relax*((1/(1/dt+2*(1/Re)/dx^2+2*(1/Re)/(dy*H/W)^2))*(Ustr(j,i)/dt+U(j,i+1)*(-Ustr(j,i+1)/(2*dx)+(1/Re)/dx^2)+U(j,i-1)*(Ustr(j,i-1)/(2*dx)+(1/Re)/dx^2)+U(j+1,i)*(-Vnstr/(2*(dy*H/W))+(1/Re)/(dy*H/W)^2)+U(j-1,i)*(Vsstr/(2*(dy*H/W))+(1/Re)/(dy*H/W)^2)-1/(dx)*(Pstr(j,i+1)-Pstr(j,i))));
      end
    end

    delF = ones(1,(nu-2)*(np-2));
    c=1;
    for i=2:nu-1
      for j=2:np-1
        Vnstr=(Vstr(j,i+1)+Vstr(j,i))/2;
        Vsstr=(Vstr(j-1,i+1)+Vstr(j-1,i))/2;
        delF(c)=U(j,i)-(1/(1/dt+2*(1/Re)/dx^2+2*(1/Re)/(dy*H/W)^2))*(Ustr(j,i)/dt+U(j,i+1)*(-Ustr(j,i+1)/(2*dx)+(1/Re)/dx^2)+U(j,i-1)*(Ustr(j,i-1)/(2*dx)+(1/Re)/dx^2)+U(j+1,i)*(-Vnstr/(2*(dy*H/W))+(1/Re)/(dy*H/W)^2)+U(j-1,i)*(Vsstr/(2*(dy*H/W))+(1/Re)/(dy*H/W)^2)-1/(dx)*(Pstr(j,i+1)-Pstr(j,i)));
        c=c+1;
      end
    end

    sum = 0;
    for c=1:length(delF)
      sum = sum + delF(c)^2;
    end

    residual1= sqrt(sum/length(delF));
  end

  residual1=1;
  residual2=1;
  relax=0.7;
  while residual1>10^(-8)
    if(residual1>residual2&&residual1>100)
      display('V cal failed')
      display(residual1)
      break;
    end
    residual2=residual1;
    for i=2:np-1
      for j=2:nv-1
        Ufstr=(Ustr(j+1,i)+Ustr(j,i))/2;
        Ubstr=(Ustr(j+1,i-1)+Ustr(j,i-1))/2;
        V(j,i)=(1-relax)*V(j,i)+relax*((1/(1/dt+2*(1/Re)/dx^2+2*(1/Re)/(dy*H/W)^2))*(Vstr(j,i)/dt+V(j,i+1)*(-Ufstr/(2*dx)+(1/Re)/dx^2)+V(j,i-1)*(Ubstr/(2*dx)+(1/Re)/dx^2)+V(j+1,i)*(-Vstr(j+1,i)/(2*(dy*H/W))+(1/Re)/(dy*H/W)^2)+V(j-1,i)*(Vstr(j-1,i)/(2*(dy*H/W))+(1/Re)/(dy*H/W)^2)-1/((dy*H/W))*(Pstr(j+1,i)-Pstr(j,i))));
      end
    end
    delF=ones(1,(nv-2)*(np-2));
    c=1;
    for i=2:np-1
      for j=2:nv-1
        Ufstr=(Ustr(j+1,i)+Ustr(j,i))/2;
        Ubstr=(Ustr(j+1,i-1)+Ustr(j,i-1))/2;
        delF(c)=V(j,i)-(1/(1/dt+2*(1/Re)/dx^2+2*(1/Re)/(dy*H/W)^2))*(Vstr(j,i)/dt+V(j,i+1)*(-Ufstr/(2*dx)+(1/Re)/dx^2)+V(j,i-1)*(Ubstr/(2*dx)+(1/Re)/dx^2)+V(j+1,i)*(-Vstr(j+1,i)/(2*(dy*H/W))+(1/Re)/(dy*H/W)^2)+V(j-1,i)*(Vstr(j-1,i)/(2*(dy*H/W))+(1/Re)/(dy*H/W)^2)-1/((dy*H/W))*(Pstr(j+1,i)-Pstr(j,i)));
        c=c+1;
      end
    end

    sum=0;

    for c=1:length(delF)
      sum=sum+delF(c)^2;
    end

    residual1=sqrt(sum/length(delF));

  end

  %updating velocity fields
  Ustr = U;
  Vstr = V;

  %calculation of error in continuity from Ustr and Vstr
  Q = zeros(np,np);
  qsum = 0 ;
  for i=2:np-1
    for j=2:np-1
      Q(j,i)=(Ustr(j,i)-Ustr(j,i-1))/dx+(Vstr(j,i)-Vstr(j-1,i))/(dy*H/W);
      qsum=qsum+Q(j,i)^2;
    end
  end
  display(k)
  sqrt(qsum)
  cnvergence(k) =  sqrt(qsum);

  if (k>2&& convergence(k)>1.5*convergence(k-1))
    display('Iterations are diverging')
    break
  end
  save('Re_V4_100_np_80_conv_o3','Ustr','Vstr','Pstr','convergence','k','Re','np','dx','dy')
  %Calculation of P correction from Pstr, Ustr and Vstr
  Pcor = zeros(np,np);
  relax = 0.7;
  residual1 = 1;
  residual2  = 1;

  while residual1>10^(-8)

    if(residual1>residual2 && residual1>100)
      display('Pcor calculation failed')
      display(residual1)
      break;
    end

    residual2 = residual1;
    A=-(1/(dx))*(1/(1/dt+2*(1/Re)/dx^2+2*(1/Re)/(dy*H/W)^2));
    B=-(1/((dy*H/W)))*(1/(1/dt+2*(1/Re)/dx^2+2*(1/Re)/(dy*H/W)^2));
    for i = 2:np-1
      for j = 2:np-1
        if(i==2&&j==2)
          Pcor(i,j)=(1-relax)*Pcor(i,j)+(relax)*(1/(A/dx+B/(dy*H/W))*(A/dx*Pcor(i+1,j)+B/(dy*H/W)*Pcor(i,j+1)+Q(i,j)));
        elseif(i==2&&j==np-1)
          Pcor(i,j)=(1-relax)*Pcor(i,j)+(relax)*(1/(A/dx+B/(dy*H/W))*(A/dx*Pcor(i+1,j)+B/(dy*H/W)*Pcor(i,j-1)+Q(i,j)));
        elseif(i==np-1&&j==2)
          Pcor(i,j)=(1-relax)*Pcor(i,j)+(relax)*(1/(A/dx+B/(dy*H/W))*(A/dx*Pcor(i-1,j)+B/(dy*H/W)*Pcor(i,j+1)+Q(i,j)));
        elseif(i==np-1&&j==np-1)
          Pcor(i,j)=(1-relax)*Pcor(i,j)+(relax)*(1/(A/dx+B/(dy*H/W))*(+A/dx*Pcor(i-1,j)+B/(dy*H/W)*Pcor(i,j-1)+Q(i,j)));
        elseif(i==2)
          Pcor(i,j)=(1-relax)*Pcor(i,j)+(relax)*(1/(A/dx+2*B/(dy*H/W))*(A/dx*Pcor(i+1,j)+B/(dy*H/W)*Pcor(i,j+1)+B/(dy*H/W)*Pcor(i,j-1)+Q(i,j)));
        elseif(i==np-1)
          Pcor(i,j)=(1-relax)*Pcor(i,j)+(relax)*(1/(A/dx+2*B/(dy*H/W))*(A/dx*Pcor(i-1,j)+B/(dy*H/W)*Pcor(i,j+1)+B/(dy*H/W)*Pcor(i,j-1)+Q(i,j)));
        elseif(j==2)
          Pcor(i,j)=(1-relax)*Pcor(i,j)+(relax)*(1/(2*A/dx+B/(dy*H/W))*(A/dx*Pcor(i+1,j)+A/dx*Pcor(i-1,j)+B/(dy*H/W)*Pcor(i,j+1)+Q(i,j)));
        elseif(j==np-1)
          Pcor(i,j)=(1-relax)*Pcor(i,j)+(relax)*(1/(2*A/dx+B/(dy*H/W))*(A/dx*Pcor(i+1,j)+A/dx*Pcor(i-1,j)+B/(dy*H/W)*Pcor(i,j-1)+Q(i,j)));
        else
          Pcor(i,j)=(1-relax)*Pcor(i,j)+(relax)*(1/(2*A/dx+2*B/(dy*H/W))*(A/dx*Pcor(i+1,j)+A/dx*Pcor(i-1,j)+B/(dy*H/W)*Pcor(i,j+1)+B/(dy*H/W)*Pcor(i,j-1)+Q(i,j)));

        end

      end

    end

    delF=ones(1,(np-2)*(np-2));
    c=1;
    for i=2:np-1
      for j=2:np-1
        if(i==2&&j==2)
          delF(c)=Pcor(i,j)-1/(A/dx+B/(dy*H/W))*(A/dx*Pcor(i+1,j)+B/(dy*H/W)*Pcor(i,j+1)+Q(i,j));
          elseif(i==2&&j==np-1)
            delF(c)=Pcor(i,j)-1/(A/dx+B/(dy*H/W))*(A/dx*Pcor(i+1,j)+B/(dy*H/W)*Pcor(i,j-1)+Q(i,j));
          elseif(i==np-1&&j==2)
            delF(c)=Pcor(i,j)-1/(A/dx+B/(dy*H/W))*(A/dx*Pcor(i-1,j)+B/(dy*H/W)*Pcor(i,j+1)+Q(i,j));
          elseif(i==np-1&&j==np-1)
            delF(c)=Pcor(i,j)-1/(A/dx+B/(dy*H/W))*(+A/dx*Pcor(i-1,j)+B/(dy*H/W)*Pcor(i,j-1)+Q(i,j));
          elseif(i==2)
            delF(c)=Pcor(i,j)-1/(A/dx+2*B/(dy*H/W))*(A/dx*Pcor(i+1,j)+B/(dy*H/W)*Pcor(i,j+1)+B/(dy*H/W)*Pcor(i,j-1)+Q(i,j));
          elseif(i==np-1)
            delF(c)=Pcor(i,j)-1/(A/dx+2*B/(dy*H/W))*(A/dx*Pcor(i-1,j)+B/(dy*H/W)*Pcor(i,j+1)+B/(dy*H/W)*Pcor(i,j-1)+Q(i,j));
          elseif(j==2)
            delF(c)=Pcor(i,j)-1/(2*A/dx+B/(dy*H/W))*(A/dx*Pcor(i+1,j)+A/dx*Pcor(i-1,j)+B/(dy*H/W)*Pcor(i,j+1)+Q(i,j));
          elseif(j==np-1)
            delF(c)=Pcor(i,j)-1/(2*A/dx+B/(dy*H/W))*(A/dx*Pcor(i+1,j)+A/dx*Pcor(i-1,j)+B/(dy*H/W)*Pcor(i,j-1)+Q(i,j));
          else
            delF(c)=Pcor(i,j)-1/(2*A/dx+2*B/(dy*H/W))*(A/dx*Pcor(i+1,j)+A/dx*Pcor(i-1,j)+B/(dy*H/W)*Pcor(i,j+1)+B/(dy*H/W)*Pcor(i,j-1)+Q(i,j));
        end
        c=c+1;
      end
    end
    sum=0;

    for c=1:length(delF)
      sum=sum+delF(c)^2;
    end
    residual1=sqrt(sum/length(delF));
  end
  Pcor(1,:)=Pcor(2,:);
  Pcor(np,:)=Pcor(np-1,:);
  Pcor(:,np)=Pcor(:,np-1);
  Pcor(:,1)=Pcor(:,2);
  %display('Pcor calculation')
  %display(residual);
  Pcor;
  Ucor=zeros(np,nu);
  Vcor=zeros(nv,np);
  for i=2:nu-1
    for j=2:np-1
      Ucor(j,i)=(1/(1/dt+2*(1/Re)/dx^2+2*(1/Re)/(dy*H/W)^2))*(-1/(dx)*(Pcor(j,i+1)-Pcor(j,i)));
    end
  end
  for i=2:np-1
    for j=2:nv-1
      Vcor(j,i)=(1/(1/dt+2*(1/Re)/dx^2+2*(1/Re)/(dy*H/W)^2))*(-1/((dy*H/W))*(Pcor(j+1,i)-Pcor(j,i)));
    end
  end
  alp=0.1;
  al=1;
  Ustr=Ustr+al*Ucor;
  Vstr=Vstr+al*Vcor;
  Pstr=Pstr+alp*Pcor;

%calculation of error in continuity from Ustr and Vstr-test

% Q=zeros(np,np);

% qsum=0;

% for i=2:np-1

% for j=2:np-1

% Q(j,i)=(Ustr(j,i)-Ustr(j,i-1))/dx+(Vstr(j,i)-Vstr(j-1,i))/dy;

% qsum=qsum+Q(j,i)^2;

% end

% end

% qsum;
end

save('Re_V4_100_np_80_conv_o3','Ustr','Vstr','Pstr','convergence','k','Re', 'np','dx','dy')
