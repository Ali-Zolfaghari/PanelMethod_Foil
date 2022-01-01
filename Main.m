%***************************************************************************************************
%*   Simulate pressure distribution on foil, using panel method (source and vortex), by presented code.
%*   I take no responsibilities for any errors in the code or damage thereby.
%*   Please notify me at zolfaghari1992iut@gmail.com if the code is used in any type of application.
%***************************************************************************************************
%*   Developer   : Ali Zolfaghari Sichani (06-05-2015)
%***************************************************************************************************
%*   References  : 
%*   'Foundations of Aerodynamics: Bases of Aerodynamic Design',Fifth edition, John Wiley & Sons, Inc.
%*   Kuethe, Arnold M. and Chow, Chuen-Yen.
%*   'Theory of Wing Sections', Abbott, I.H. and Von Doenhoff, A.E.
%*   
%***************************************************************************************************
%*   Inputs      :
%*   get input when run the code
%*   Outputs      :
%*   all data export in Tecplot format
%***************************************************************************************************



clear,clc
close all
format compact
format long



%   Get input of airfoil
%   The profile specified is a four-digit unmodified NACA airfoil.
%   Automatically generate points on the contour of the corresponding NACA airfoil.
%   The non-dimensional velocity parallel to each panel (v/vinf) 
check=0;
tot=0;
while(check==0)
disp('NACA AIRFOIL ');
foil=input(' --> : \n','s');
NC=str2double(foil);
dn=size(foil,2);
if ( dn==5 || dn==4 );
for i=1:dn
    naca(dn-i+1)= mod(NC,10);
    tot=tot+naca(dn-i+1);
    NC=0.1*(NC-mod(NC,10));
end
check=1;
if(tot==0)
    check=0;
    disp('Please Enter Correct Input');
end
else
    check=0;
    disp('Please Enter Correct Input');
end
end

%	Get input of Method , Angle of attack , Panel number
inchck=0;
while(inchck==0)
    disp('" Source Method = 0 " " Vortex Method = 1 "');
    METHOD=input(' --> : \n');
    if(METHOD==0 || METHOD==1)
        inchck=1;
    else
        inchck=0;
        disp('Please Enter Correct Input');
    end
end

inchck=0;
while(inchck==0)
    disp('Angle of attack degrees ');
    angle=input(' --> : \n');
    if(angle>=0 || angle <30)
        alpha=angle*pi/180;
        inchck=1;
    else
        inchck=0;
        disp('Please Enter Correct Input');
    end
end

inchck=0;
while(inchck==0)
    disp('Total number of panels ');
    Ne=input(' --> : \n');
    if(Ne>3)
        Ne=double((uint32(Ne/2)))*2;
        Np=Ne+1;
        inchck=1;
    else
        inchck=0;
        disp('Please Enter Correct Input');
    end
end

%   Genarate airfoil coordinates
if dn==4;
    m=naca(1)/100;
    p=naca(2)/10;
    thk=(naca(3)*10+naca(4))/100;
elseif dn==5;
    thk=(naca(4)*10+naca(5))/100;
    des=naca(1)*100+naca(2)*10+naca(3);
    pst=uint32((naca(2)*10+naca(3))/2);
    p=double(pst)/100;
    switch des;
        case(210), m=0.0580;
            k1=361.4;
        case(220), m=0.1260;
            k1=51.64;
        case(230), m=0.2025;
            k1=15.957;
        case(240), m=0.2900;
            k1=6.643;
        case(250), m=0.3910;
            k1=3.230;
        case(221), m=0.1300;
            k1=51.990;
            kr=0.000764;
        case(231), m=0.2170;
            k1=15.793;
            kr=0.00677;
        case(241), m=0.3180;
            k1=6.520;
            kr=0.0303;
        case(251), m=0.4410;
            k1=3.191;
            kr=0.1355;
    end
end
cang=0;
astp=2*pi/Ne;
for i=1:Np;
    xs(i)=0.5+0.5*cos((i-1)*astp);
end
if dn==4;
    if p~=0;
        med=(Ne/2)+1;
        x(med)=0;
        y(med)=0;
        yl=y(med);
        ycr(med)=yl;
        for i=(med+1):Np;
            if xs(i)>=p;
                ycc=(m/(1-p)^2)*((1-2*p)+2*p*xs(i)-xs(i)^2);
            else
                ycc=(m/(p^2))*(2*p*xs(i)-xs(i)^2);
            end
            ycr(i)=ycc;
            yt=(thk/0.2)*(0.2969*xs(i)^0.5-0.126*xs(i)-0.3516*xs(i)^2+0.2843*xs(i)^3-0.1015*xs(i)^4-0.0021*xs(i)^8);
            dyc=ycc-yl;
            dxs=xs(i)-xs(i-1);
            tc=atan2(dyc,dxs);
            x(i)=xs(i)-yt*sin(tc);
            y(i)=ycc+yt*cos(tc);
            x(2*med-i)=xs(2*med-i)+yt*sin(tc);
            y(2*med-i)=ycc-yt*cos(tc);
            yl=ycc;
        end;
        y(1)=0;
        x(1)=1;
        y(Np)=0;
        x(Np)=1;
    else
        med=(Ne/2)+1;
        y(med)=0;
        x=xs;
        for i=(med+1):Np;
            yt=(thk/0.2)*(0.2969*xs(i)^0.5-0.126*xs(i)-0.3516*xs(i)^2+0.2843*xs(i)^3-0.1015*xs(i)^4-0.0021*xs(i)^8);
            y(i)=yt;
            y(2*med-i)=-yt;
        end;
        y(1)=0;
        x(1)=1;
        y(Np)=0;
        x(Np)=1;
    end
else
    if naca(3)==0;
        med=(Ne/2)+1;
        x(med)=0;
        y(med)=0;
        yl=y(med);
        ycr(med)=yl;
        for i=(med+1):Np;
            if xs(i)>=m;
                ycc=((k1*m^3)/6)*(1-xs(i));
            else
                ycc=(k1/6)*(xs(i)^3-3*m*xs(i)^2+m^2*(3-m)*xs(i));
            end
            ycr(i)=ycc;
            yt=(thk/0.2)*(0.2969*xs(i)^0.5-0.126*xs(i)-0.3516*xs(i)^2+0.2843*xs(i)^3-0.1015*xs(i)^4-0.0021*xs(i)^8);
            dyc=ycc-yl;
            dxs=xs(i)-xs(i-1);
            tc=atan2(dyc,dxs);
            x(i)=xs(i)-yt*sin(tc);
            y(i)=ycc+yt*cos(tc);
            x(2*med-i)=xs(2*med-i)+yt*sin(tc);
            y(2*med-i)=ycc-yt*cos(tc);
            yl=ycc;
        end
        y(1)=0;
        x(1)=1;
        y(Np)=0;
        x(Np)=1;
    else
        med=(Ne/2)+1;
        x(med)=0;
        y(med)=0;
        yl=y(med);
        ycr(med)=yl;
        for i=(med+1):Np;
            if xs(i)>=m;
                ycc=(k1/6)*(kr*(xs(i)-m)^3-kr*xs(i)*(1-m)^3-xs(i)*m^3+m^3);
            else
                ycc=(k1/6)*((xs(i)-m)^3-kr*xs(i)*(1-m)^3-xs(i)*m^3+m^3);
            end;
            ycr(i)=ycc;
            yt=(thk/0.2)*(0.2969*xs(i)^0.5-0.126*xs(i)-0.3516*xs(i)^2+0.2843*xs(i)^3-0.1015*xs(i)^4-0.0021*xs(i)^8);
            dyc=ycc-yl;
            dxs=xs(i)-xs(i-1);
            tc=atan2(dyc,dxs);
            x(i)=xs(i)-yt*sin(tc);
            y(i)=ycc+yt*cos(tc);
            x(2*med-i)=xs(2*med-i)+yt*sin(tc);
            y(2*med-i)=ycc-yt*cos(tc);
            yl=ycc;
        end
        y(1)=0;
        x(1)=1;
        y(Np)=0;
        x(Np)=1;
    end
end

XR=x;YR=y;

%   x  = airfoil x-location (control panel )
%   y  = airfoil y-location (control panel )
%   xc = control point x-location (control point)
%   yc = control point y-location (control point)

ch1ck=0;
while(ch1ck==0)
    if(METHOD==0)
    disp('" Point on airfoil & Panel on airfoil = 0 " ');
    disp('" Point on camber  & Panel on airfoil = 1 "');
    ctrl=input(' --> : \n');
    else
    ctrl=0;
    end
    if (ctrl==0)
        STR=['Point on airfoil & Panel on airfoil NACA ',foil];
        for i=1:Ne
            xc(i)=(x(i)+x(i+1))/2;
            yc(i)=(y(i)+y(i+1))/2;
        end
        XC=x;
        YC=y;
        ch1ck=1;
    elseif (ctrl==1)
        STR=['Point on camber & Panel on airfoil NACA ',foil];
        ch2ck=0;
        while(ch2ck==0)
            disp('" Distance LE from starting control point on camber : " ');
            disp('" Half of thickness = 0 " ');
            disp('" Choose a Number = 1 " ');
            ST=input(' --> : \n');
            if (ST==0)
                stp=thk/2;
                ch2ck=1;
            elseif (ST==1)
                disp('" Enter start Distance between 0 to 1 " : ');
                stp=input(' --> : \n');
                if (stp==0)
                    stp=0.00001;
                elseif (stp==1)
                    stp=0.99999;
                end
                ch2ck=1;
            else
                ch2ck=0;
                disp('Please Enter Correct Input');
            end
        end
        ch3ck=0;
        while(ch3ck==0)
            disp('" Distance TE from ending control point on camber : " ');
            disp('" Half of thickness = 0 " ');
            disp('" Choose a Number = 1 " ');
            ED=input(' --> : \n');
            if (ED==0)
                edp=thk/2;
                ch3ck=1;
            elseif (ED==1)
                disp(['" Enter end Distance between 0 to ',num2str(1-stp),' " : ']);
                edp=input(' --> : \n');
                if (edp==0)
                    edp=0.00001;
                elseif (edp==1)
                    edp=0.99999;
                end
                ch3ck=1;
            else
                ch3ck=0;
                disp('Please Enter Correct Input');
            end
        end
        for i=1:(Np+1)/2
            XC(i)=(x(i)+x(Np+1-i))/2;
            YC(i)=(y(i)+y(Np+1-i))/2;
        end
        clear x y
        for i=1:Ne/2
            x(i)=((1+stp-edp)/2)+((1-stp-edp)/2)*cos((i-1)*(astp));
        end
        y=pchip(XC,YC,x);
        for i=1:size(x,2)-1
            xc(i)=(x(i)+x(i+1))/2;
            yc(i)=(y(i)+y(i+1))/2;
        end
        XX=[x,xc(end:-1:1)];
        YY=[y,yc(end:-1:1)];
        XX(Ne)=XX(1);
        YY(Ne)=YY(1);
        clear x y xc yc
        x=XR;
        y=YR;
        xc=XX;
        yc=YY;
        ch1ck=1;
    else
        ch1ck=0;
        disp('Please Enter Correct Input');
    end
end

%	Calculate s , t , rhs
%   t  = angle of panel relative to horizontal (in radians)
%   vd = non-dimensional (v/vinf) velocity parallel to panel
%   cp =pressure coefficient at panel
for i=1:Ne 
    s(i)=((x(i+1)-x(i))^2+(y(i+1)-y(i))^2)^(0.5);
    t(i)=atan2((y(i+1)-y(i)),(x(i+1)-x(i)));
    rhs(i)=sin(t(i)-alpha);
end

if (METHOD==0)
    
% Source Method
    for i=1:Ne
        for j=1:Ne
            A(i,j)=-((xc(i)-x(j))*cos(t(j)))-((yc(i)-y(j))*sin(t(j)));
            B(i,j)=((xc(i)-x(j))^2)+((yc(i)-y(j))^2);
            C(i,j)=sin(t(i)-t(j));
            D(i,j)=((yc(i)-y(j))*cos(t(i)))-((xc(i)-x(j))*sin(t(i)));
            E(i,j)=((xc(i)-x(j))*sin(t(j)))-((yc(i)-y(j))*cos(t(j)));
        end
    end
    
    for i=1:Ne
        for j=1:Ne
            if(i~=j)
                I(i,j)=((C(i,j)/2)*log(((s(j)^2)+2*A(i,j)*s(j)+B(i,j))/(B(i,j))))+...
                    ((D(i,j)-A(i,j)*C(i,j))/E(i,j))*((atan((s(j)+A(i,j))/E(i,j)))-atan(A(i,j)/E(i,j)));
            end
        end
    end
    
    for i=1:Ne
        I(i,i)=pi;
    end
    
    LANDA=(I)\(rhs');
    
    for i=1:Ne
        for j=1:Ne
            if(i~=j)
                Intgral(i,j)=((D(i,j)-A(i,j)*C(i,j))/(2*E(i,j)))*log(((s(j)^2)+2*A(i,j)...
                    *s(j)+B(i,j))/(B(i,j)))-(C(i,j)*((atan((s(j)+A(i,j))/E(i,j)))-atan(A(i,j)/E(i,j))));
            end
        end
    end
    
    for i=1:Ne
        Intgral(i,i)=0;
    end
    
    for i=1:Ne
        for j=1:Ne
            sigma(i,j)=LANDA(j)*Intgral(i,j);
        end
    end
    zigma=sum(sigma');
    cx=0;
    for i=1:Ne
        vd(i)=cos(t(i)-alpha)+zigma(i);
        cp(i)=1-(vd(i)^2);
        cx=cx-cp(i)*s(i)*cos(t(i)-alpha);
    end
    
    Error=0;
    for k=1:Ne
        Error=Error+LANDA(k)*s(k);
    end
    
    
elseif(METHOD==1)
    
% Vortex Method
    for i=1:Ne;
        for j=1:Ne;
            if i==j;
                cn1(i,j)=-1;
                cn2(i,j)=1;
                ct1(i,j)=pi/2;
                ct2(i,j)=pi/2;
            else
                A=-(xc(i)-x(j))*cos(t(j))-(yc(i)-y(j))*sin(t(j));
                B=(xc(i)-x(j))^2+(yc(i)-y(j))^2;
                C=sin(t(i)-t(j));
                D=cos(t(i)-t(j));
                E=(xc(i)-x(j))*sin(t(j))-(yc(i)-y(j))*cos(t(j));
                F=log(1+(s(j)^2+2*A*s(j))/B);
                G=atan2((E*s(j)),(B+A*s(j)));
                P=(xc(i)-x(j))*sin(t(i)-2*t(j))+(yc(i)-y(j))*cos(t(i)-2*t(j));
                Q=(xc(i)-x(j))*cos(t(i)-2*t(j))-(yc(i)-y(j))*sin(t(i)-2*t(j));
                cn2(i,j)=D+0.5*Q*F/s(j)-(A*C+D*E)*G/s(j);
                cn1(i,j)=0.5*D*F+C*G-cn2(i,j);
                ct2(i,j)=C+0.5*P*F/s(j)+(A*D-C*E)*G/s(j);
                ct1(i,j)=0.5*C*F-D*G-ct2(i,j);
            end
        end
    end
    
    for i=1:Ne;
        an(i,1)=cn1(i,1);
        an(i,Np)=cn2(i,Ne);
        at(i,1)=ct1(i,1);
        at(i,Np)=ct2(i,Ne);
        for j=2:Ne;
            an(i,j)=cn1(i,j)+cn2(i,(j-1));
            at(i,j)=ct1(i,j)+ct2(i,(j-1));
        end
    end
    an(Np,1)=1;
    an(Np,Np)=1;
    rhs(Np)=0;
    for j=2:Ne;
        an(Np,j)=0;
    end
    
    gamma=an\rhs';
    %   The actual circulation density, GAMMA=gamma*2*pi*Vinf.)
    cx=0;
    cy=0;
    for i=1:Ne;
        smt=0;
        for j=1:Np;
            sm=at(i,j)*gamma(j);
            smt=smt+sm;
        end;
        vd(i)=cos(t(i)-alpha)+smt;
        cp(i)=1-(vd(i))^2;
        cx=cx+cp(i)*s(i)*sin(t(i));
        cy=cy-cp(i)*s(i)*cos(t(i));
    end
    clift=cy*cos(alpha)-cx*sin(alpha);
    cdrag=cy*sin(alpha)+cx*cos(alpha);
    
end

%   Plot
figure(1);
subplot(2,2,[1 2]),plot(x,y,'*-',xc,yc,'.r',XR,YR,'g','LineWidth',2);
grid on;xlabel('X');ylabel('Y');legend('Control Panel','Control Point','Airfoil');title(STR);
subplot(2,2,[3 4]),plot(x(1:end-1),-cp,'.-','LineWidth',2);
grid on;xlabel('X/C');ylabel('C_P');title('Pressure Coeffiecient');

exp=x(1:end-1);
fileID=fopen('X.dat','w');
fprintf(fileID,'%8.6f \r\n',exp);
fclose(fileID);

fileID=fopen('CP.dat','w');
fprintf(fileID,'%8.6f \r\n',-cp);
fclose(fileID);

