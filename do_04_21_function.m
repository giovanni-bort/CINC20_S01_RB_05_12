% Processing
    %-----------------------------------------------------------------------
    ECG=ECG/1000;
    I=ECG(1,:); II=ECG(2,:); III=ECG(3,:);
    aVR=ECG(4,:); aVL=ECG(5,:); aVF=ECG(6,:);
    V1=ECG(7,:); V2=ECG(8,:); V3=ECG(9,:);
    V4=ECG(10,:); V5=ECG(11,:); V6=ECG(12,:);

   [A,a]=find(I>2); I(a)=2; [A,a]=find(II>2); II(a)=2;
   [A,a]=find(III>2); III(a)=2;[A,a]=find(V1>2); V1(a)=2;
   [A,a]=find(V2>2); V2(a)=2;[A,a]=find(V3>2); V3(a)=2;
   [A,a]=find(V4>2); V4(a)=2; [A,a]=find(V5>2); V5(a)=2;
   [A,a]=find(V6>2); V6(a)=2;

    [A,a]=find(I<-2); I(a)=-2; [A,a]=find(II<-2); II(a)=-2;
    [A,a]=find(III<-2); III(a)=-2;[A,a]=find(V1<-2); V1(a)=-2;
    [A,a]=find(V2<-2); V2(a)=-2;[A,a]=find(V3<-2); V3(a)=-2;
    [A,a]=find(V4<-2); V4(a)=-2; [A,a]=find(V5<-2); V5(a)=-2;
    [A,a]=find(V6<-2); V6(a)=-2;


PVC(j)=0;

% QRS detectin
    d1=I(1:2:length(II));      % 250 Hz
        if(sum(abs(d1))==0), d1=II (1:2:length(II)); end       %MODFIED 22.04.20
        if(sum(abs(d1))==0), d1=III(1:2:length(II)); end       %MODFIED 22.04.20
        if(sum(abs(d1))==0), d1=aVR(1:2:length(II)); end        %MODFIED 22.04.20
    d2=V2(1:2:length(II));
        if(sum(abs(d2))==0), d2=V3(1:2:length(II)); end       %MODFIED 22.04.20
        if(sum(abs(d2))==0), d2=V4(1:2:length(II)); end       %MODFIED 22.04.20
        if(sum(abs(d2))==0), d2=V5(1:2:length(II)); end       %MODFIED 22.04.20
  
    
    
    [QRS]=QRS_det(d1,d2);      % go to function
    QRS=QRS*2;  % 500Hz
     figure(1); clf
      x=(1:1:length(I))/Hz;
      plot(x,I, QRS/Hz, I(QRS), 'ro')
%      figure(3)
 %    plot(x,V2)

    if length(QRS)<=4;
        I=I/2; V2=V2/2;
        [QRS]=QRS_det(d1,d2);      % go to function
        QRS=QRS*2;                  % 500 Hz
    end    
    if QRS(1)==QRS(2); QRS=QRS(2:length(QRS)); end
   
 
 if length(QRS)<5 | length(I) - QRS(length(QRS)) >2000
ecg=II;
if(sum(abs(ecg(:)))==0),ecg=III;end    %MODFIED 22.04.20
if(sum(abs(ecg(:)))==0),ecg=V1;end    %MODFIED 22.04.20
if(sum(abs(ecg(:)))==0),ecg=V2;end    %MODFIED 22.04.20
if(sum(abs(ecg(:)))==0),ecg=V3;end    %MODFIED 22.04.20



fs=1000;
QRS=0;
[QRS,sign,en_thres] = qrs_detect2(ecg',0.25,0.6,fs);
figure(1); clf
x=(1:1:length(I))/Hz;
plot(x,II, QRS/Hz, II(QRS), 'ro')
end

I_AVB(j)=0; AF(j)=0; PVC(j)=PVC(j); PAC(j)=0; LBBB(j)=0; RBBB(j)=0; Normal(j)=0; STD(j)=0; STE(j)=0;


%======================================================================
%		Preprocessing
%======================================================================n=7;					% 7*2 = 14 ms smoothing interval on each side
T=1/Hz;					% [s] sampling period
Fc=0.64;				% [Hz]
a=1; b=[1 1 1 1 1 1 1 1 1 1]/10;
n=7;

X=I;        Y=filter(b,a,X);        % 50 Hz
X=Y;        [Y]=smooth(X,n);			% Smoothing
X=Y;		[Y]=drift_Ivo(X,Fc,T);	I=Y;	% Drift filtration

X=II;       Y=filter(b,a,X);        % 50 Hz
X=Y;        [Y]=smooth(X,n);			% Smoothing
X=Y;		[Y]=drift_Ivo(X,Fc,T);	II=Y;	% Drift filtration

X=V1;       Y=filter(b,a,X);        % 50 Hz
X=Y;        [Y]=smooth(X,n);			% Smoothing
X=Y;		[Y]=drift_Ivo(X,Fc,T);	V1=Y;	% Drift filtration

X=V2;       Y=filter(b,a,X);        % 50 Hz
X=Y;        [Y]=smooth(X,n);			% Smoothing
X=Y;		[Y]=drift_Ivo(X,Fc,T);	V2=Y;	% Drift filtration

X=V3;       Y=filter(b,a,X);        % 50 Hz
X=Y;        [Y]=smooth(X,n);			% Smoothing
X=Y;		[Y]=drift_Ivo(X,Fc,T);	V3=Y;	% Drift filtration

X=V4;       Y=filter(b,a,X);        % 50 Hz
X=Y;        [Y]=smooth(X,n);			% Smoothing
X=Y;		[Y]=drift_Ivo(X,Fc,T);	V4=Y;	% Drift filtration

X=V5;       Y=filter(b,a,X);        % 50 Hz
X=Y;        [Y]=smooth(X,n);			% Smoothing
X=Y;		[Y]=drift_Ivo(X,Fc,T);	V5=Y;	% Drift filtration

X=V6;       Y=filter(b,a,X);        % 50 Hz
X=Y;        [Y]=smooth(X,n);			% Smoothing
X=Y;		[Y]=drift_Ivo(X,Fc,T);	V6=Y;	% Drift filtration


% Fiducial point = at the max peak
FP=0;
left=80*Hz/1000;    %180 ms
right=100*Hz/1000;   %180 ms

if QRS(1)-left<=0; QRS=QRS(2:length(QRS)); end
if QRS(1)-left<=0; QRS=QRS(2:length(QRS)); end

Max=max(I(QRS(1)-left:QRS(1)+right)); Min=abs(min(I(QRS(1)-left:QRS(1)+right)));
if QRS(1)- left <= 0; QRS=QRS(2:length(QRS)); end
if Max>=Min;
    for i=1:length(QRS)-1;
        [M,m]=max(I(QRS(i)-left:QRS(i)+right));
        FP(i)=QRS(i)-left+m;
    end    
else
    for i=1:length(QRS)-1;
        [M,m]=min(I(QRS(i)-left:QRS(i)+right));
        FP(i)=QRS(i)-left+m ;
    end    
end


x1=(1:1:length(I))/Hz;
figure(1); clf
plot(x1,I, FP/Hz, I(FP),'ro') 
title(filename)

% PVC(j) & PAC(j)
PAC(j)=0; 
if FP(1)==FP(2); FP=FP(2:length(FP)); end
RR=FP(2:length(FP))-FP(1:length(FP)-1);
RRm=mean(RR);
RRm_20=RRm*20/100;  % 20%

RR1=FP(2:2:length(FP))-FP(1:2:length(FP)-1);
RRm1=mean(RR1);
if abs(RRm-RRm1) >70; PAC(j)=round(length(FP)/2); end
for i=1:length(FP)-2; 
    if RRm-RR(i)>0.05*RRm && 1.4*RR(i)<RR(i+1) &&...
            RRm-RR(i)>0.05*RRm && 1.3*RR(i)<RR(i+1) &&...
            RR(i)+RR(i+1)<2.2*RRm && RR(i+1)>RRm
         PAC(j)=PAC(j)+1;
    end
end

% PVC
Lead=I;

if(sum(abs(Lead(:)))==0),Lead=II;end   %*** MODIFIED 22.4.20
if(sum(abs(Lead(:)))==0),Lead=III;end  %*** MODIFIED 22.4.20
if(sum(abs(Lead(:)))==0),Lead=V1;end   %*** MODIFIED 22.4.20
if(sum(abs(Lead(:)))==0),Lead=V2;end   %*** MODIFIED 22.4.20



if mean(Lead(FP))>0
    if max(Lead(FP))*0.8>mean(Lead(FP))
    PVC(j)=1; PAC(j)=0;
    end
    if min(Lead(FP))<0.7*mean(Lead(FP))
    PVC(j)=1; PAC(j)=0;
    end
else
    if min(Lead(FP))*0.8<mean(Lead(FP))
    PVC(j)=1; PAC(j)=0;
    end
    if max(Lead(FP))>0.7*mean(Lead(FP))
    PVC(j)=1; PAC(j)=0;
    end
end

Left=650*Hz/1000;    %650 ms left
Right=round(RRm)+150*Hz/1000;
Right=round(Right/2)*2; % Right must be even

if RRm<250
    Left=600*Hz/1000;
    Right=round(RRm)+50*Hz/1000;
    Right=round(Right/2)*2;
end


% Average of all leads
A=find(RR<150);
FP(A)=0;

for i=length(A):-1:1;
    n=A(i);
    FP=[FP(1:n-1), FP(n+1:length(FP))];
end
if FP(length(FP))==0; FP=FP((1):length(FP)-1); end


[FP_,matrix]=cross_corr(I,FP,Hz,Left,Right);     I_matrix=matrix;  I_mean=mean(I_matrix);
[FP_,matrix]=cross_corr(II,FP,Hz,Left,Right);    II_matrix=matrix; II_mean=mean(II_matrix);
[FP_,matrix]=cross_corr(III,FP,Hz,Left,Right);   III_matrix=matrix; III_mean=mean(III_matrix);
if mean(V1)~=0      %Peripheral leads only
[FP_,matrix]=cross_corr(V1,FP,Hz,Left,Right);    V1_matrix=matrix; V1_mean=mean(V1_matrix);
[FP_,matrix]=cross_corr(V2,FP,Hz,Left,Right);    V2_matrix=matrix; V2_mean=mean(V2_matrix);
[FP_,matrix]=cross_corr(V3,FP,Hz,Left,Right);    V3_matrix=matrix; V3_mean=mean(V3_matrix);
[FP_,matrix]=cross_corr(V4,FP,Hz,Left,Right);    V4_matrix=matrix; V4_mean=mean(V4_matrix);
[FP_,matrix]=cross_corr(V5,FP,Hz,Left,Right);    V5_matrix=matrix; V5_mean=mean(V5_matrix);
if mean(V6)~=0
[FP_,matrix]=cross_corr(V6,FP,Hz,Left,Right);    V6_matrix=matrix; V6_mean=mean(V6_matrix);
end
end
x1=1/Hz:1/Hz:length(I_mean)/Hz;
figure(2); clf
subplot(3,3,1); plot(x1,I_mean); title('I')
subplot(3,3,2); plot(x1,II_mean); title('II')
subplot(3,3,3); plot(x1,III_mean); title('III')
if mean(V1)==0; 
    V1_mean(1:length(I_mean))=0;
    V2_mean(1:length(I_mean))=0;
    V3_mean(1:length(I_mean))=0;
    V4_mean(1:length(I_mean))=0;
    V5_mean(1:length(I_mean))=0;
    V6_mean(1:length(I_mean))=0;
end
    if mean(V6)==0;     V6_mean(1:length(I_mean))=0; end

x2=(1:1:length(V1_mean))/Hz;
subplot(3,3,4); plot(x2,V1_mean); title('V1')
subplot(3,3,5); plot(x2,V2_mean); title('V2')
subplot(3,3,6); plot(x2,V3_mean); title('V3')
subplot(3,3,7); plot(x2,V4_mean); title('V4')
subplot(3,3,8); plot(x2,V5_mean); title('V5')
x3=(1:1:length(V6_mean))/Hz;
subplot(3,3,9); plot(x3,V6_mean); title('V6')

% Find the izoelectric point in II or V1
a=90;  b=200;  %  Offset
Flat=20*Hz/1000; 
Lead=I_mean;
[Max,t1]=max(Lead(a:b));
[Min,t2]=min(Lead(a:b));
Min=abs(Min);
if Max>=Min; 
    Crit=0.02*Max;	% make Crit proportional to the QRSampl;
    From=t1+a;
else Crit=0.02*Min;	% make Crit proportional to the QRSampl;j
    From=t2+a;
end

To=max(From-60, 12);
[Pnt,Crit]=Flat_l(Lead,Crit,From,To,Flat);  % go to function
Iz_1=Pnt-Flat;

Flat=20*Hz/1000;        % 20 ms
From=From+30;            			% from the point of the peak
To=From+80*Hz/1000;                % to 100 ms right ...
[Pnt,Crit]=flat_right(Lead,Crit,From,To,Flat);  % go to function
J_1=Pnt;

x1=(1:1:length(I_mean))/Hz;
subplot (3,3,1);
plot(x1,Lead), hold on
plot(Iz_1/Hz, Lead(Iz_1), 'ro', J_1/Hz, Lead(J_1), 'ro')
title ('I')

Lead=V1_mean;
Iz2=Iz_1;
if V1_mean~=0;
[Max,t1]=max(Lead(a:b));
[Min,t2]=min(Lead(a:b));
Min=abs(Min);
if Max>=Min; 
    Crit=0.02*Max;	% make Crit proportional to the QRSampl;
    From=t1+a;
else Crit=0.02*Min;	% make Crit proportional to the QRSampl;j
    From=t2+a;
end

To=max(From-60, 12);
[Pnt,Crit]=Flat_l(Lead,Crit,From,To,Flat);  % go to function

Iz_2=Pnt-Flat;       % 20 ms
From=From+30;            			% from the point of the peak
To=From+80*Hz/1000;                % to 100 ms right ...
[Pnt,Crit]=flat_right(Lead,Crit,From,To,Flat);  % go to function
J_2=Pnt;

x1=(1:1:length(Lead))/Hz;
subplot (3,3,4);
plot(x1,Lead, Iz_2/Hz, Lead(Iz_2), 'ro', J_2/Hz, Lead(J_2), 'ro')

Iz=max(Iz_1, Iz_2);

Iz2=Iz;
if abs(J_1-J_2) <20; J=max(J_1, J_2);
    else J=round((J_1+J_2)/2);
end


% RBBB(j) criteria 
% 'M' shaped  QRS in V1; wide S in V6 & I
%	Searching for QRS segmentation 
RBBB1=0; RBBB2=0; RBBB3=0; RBBB4=0; RBBB5=0; RBBB(j)=0; LBBB(j)=0;
a=15;
Lead=V2_mean;

    Iz=Iz+a;    J=J-a;
    [Zero1]=Deriv(Lead,Iz,J); % go to function
    Iz=Iz-a;   J=J+a;   % restore J & Iz
      
    
    x1=(1:1:length(V1_mean))/Hz;
    subplot(3,3,5);
    plot(x1,Lead, Iz/Hz,Lead (Iz), 'ko',J/Hz,Lead (J), 'ko'); hold on
    plot(Zero1/Hz, Lead    (Zero1), 'go')   
    if length(Zero1)>=3; title ('RBBB'); RBBB1=1; 
    else title ('V2'); RBBB1=0; end


    Lead=V1_mean;
    Iz=Iz+a;    J=J-a;
    [Zero1]=Deriv(Lead,Iz,J); % go to function
    Iz=Iz-a;   J=J+a;   % restore J & Iz

    m=1; A=Lead(Iz:J);
    D1(1+m:length(A)-m) = Lead(1:length(A)-m-m)- Lead(1+m:length(A)-m);

    x1=(1:1:length(V1_mean))/Hz;
    subplot(3,3,4);
    plot(x1,Lead, Iz/Hz,Lead (Iz), 'ko',J/Hz,Lead (J), 'ko'); hold on
    hold on
    plot(Zero1/Hz, Lead(Zero1), 'go')   
    if length(Zero1)>=3; title('RBBB'); RBBB2=1; 
    else title ('V1'); RBBB2=0; end
    
     Lead=V3_mean;
    
    Iz=Iz+a;    J=J-a;
    [Zero1]=Deriv(Lead,Iz,J); % go to function
    m=1; A=Lead(Iz:J);
    D1(1+m:length(A)-m) = Lead(1:length(A)-m-m)- Lead(1+m:length(A)-m);
      A=mean(D1(8:14));


    x1=(1:1:length(V1_mean))/Hz;
    Iz=Iz-a;   J=J+a;       % restore J & Iz
    subplot(3,3,6);
    plot(x1,Lead, Iz/Hz,Lead (Iz), 'ko',J/Hz,Lead (J), 'ko'); hold on
    plot(Zero1/Hz, Lead(Zero1), 'go')   
    if length(Zero1)>=3; title ('RBBB'); RBBB3=1; 
    else title ('V3'); RBBB3=0; end
    
    
    Lead=V4_mean;
    Iz=Iz+a;    J=J-a;
    [Zero1]=Deriv(Lead,Iz,J); % go to function
    Iz=Iz-a;   J=J+a;   % restore J & Iz
    
    m=1; A=Lead(Iz:J);
    D1(1+m:length(A)-m) = Lead(1:length(A)-m-m)- Lead(1+m:length(A)-m);
    
    x1=(1:1:length(V1_mean))/Hz;
    subplot(3,3,7);
    plot(x1,Lead, Iz/Hz,Lead (Iz), 'ko',J/Hz,Lead (J), 'ko'); hold on
    plot(Zero1/Hz, Lead(Zero1), 'go')   
    if length(Zero1)>=3; title ('RBBB'); RBBB4=1; 
    else title ('V4'); RBBB4=0; end

%     Lead=III_mean
% 
%     Iz=Iz+a;    J=J-a;
%     [Zero1]=Deriv(Lead,Iz,J); % go to function
%     Iz=Iz-a;   J=J+a;   % restore J & Iz
%       
%     
%     x1=(1:1:length(V1_mean))/Hz;
%     subplot(3,3,3);
%     plot(x1,Lead, Iz/Hz,Lead (Iz), 'ko',J/Hz,Lead (J), 'ko'); hold on
%     plot(Zero1/Hz, Lead    (Zero1), 'go')   
%     if length(Zero1)>=3; title ('RBBB'); RBBB5=1; 
%     else title ('III'); RBBB5=0; end

if RBBB1|RBBB2|RBBB3|RBBB4; RBBB(j)=1; end

% LBBB(j) criteria 
% monophasic R, no Q in I & V6; QS in V1
LBBB(j)=0;

if J-Iz>120*Hz/1000    
    Lead=I_mean; Minus=0; Plus=0;
    Lead=Lead-Lead(Iz);
    for i=Iz:J
        if Lead(i)>0; Plus=Plus+Lead(i);
        else Minus=Minus+Lead(i);
        end
    end

    if Plus > 15*abs(Minus); LBBB(j)=1;
    %elseif Minus> 5*Plus; LBBB(j)=0;
    end
    
    Lead=V6_mean; 
    a=90;  b=200;  %  Offset
    [Max,t1]=max(Lead(a:b));
    [Min,t2]=min(Lead(a:b));
    Min=abs(Min);
if Max>=Min; 
    Crit=0.02*Max;	% make Crit proportional to the QRSampl;
    From=t1+a;
else Crit=0.02*Min;	% make Crit proportional to the QRSampl;j
    From=t2+a;
end

    prag=0.8*Lead(From);
    for i=From:-1:From-40
        if Lead(i)<=prag; break
        end
    end
    pl=i;
    
    for i=From:From+40
        if Lead(i)<=prag
            0.8*Lead(From); break
        end
    end
    pr=i;
    
    if pr-pl>70*Hz/1000; LBBB(j)=1;  end
  
    if pr-pl<30*Hz/1000;LBBB(j)=0;   end
  
   figure(2);
    subplot(3,3,9);
    
    x1=(1:1:length(Lead))/Hz;
    plot(x1,Lead, From/Hz, Lead(From),'ko',...
        pl/Hz, Lead(pl),'ko', pr/Hz, Lead(pr),'ko')
    title ('V6')
end

Iz2=Iz; J2=J;

% STD(j) STE(j)
ST1=(V1_mean(J2)-V1_mean(Iz2));
ST2=(V2_mean(J2)-V2_mean(Iz2));
ST3=(V3_mean(J2)-V3_mean(Iz2));
ST4=(V4_mean(J2)-V4_mean(Iz2));
ST5=(V5_mean(J2)-V5_mean(Iz2));
ST6=(V6_mean(J2)-V6_mean(Iz2));

ma=[ST1,ST2,ST3,ST4,ST5,ST6];
[Max,m]=max(abs(ma));
ST=ma(m);
if ST>0.10; STE(j)=1; else STE(j)=0; end
if ST<-0.10; STD(j)=1; else STD(j)=0; end

end

%---------
% I_AVB(j)
% Find P wave in Average ecg

% Positive P-wave


prag2=0.0050;     % 5 uV threshold

P_wave=0; 
Lead=II_mean;

for i=Iz2-20:-1:1;
    if Lead(i+20)-Lead(i)>4*prag2 && Lead(i+20)-Lead(i+40)>4*prag2 &&...
            Lead(i+20)-Lead(i+10)>prag2 && Lead(i+20)-Lead(i+30)>prag2;
        P_wave=1;
        break
    end
end

P_x=i+20;
    figure(2);
    x=(1:1:length(I_mean))/Hz;
    subplot(3,3,2);
    
    
    Iz=Iz2;
    if P_wave==1
    plot(x, Lead, P_x/Hz, Lead(P_x), 'go', Iz/Hz, Lead(Iz), 'go');
    else plot(x, Lead, Iz/Hz, Lead(Iz), 'go');
    end
    title('II')    
    
if P_wave==1 && Iz-P_x > 125*Hz/1000; I_AVB(j)=1;
else I_AVB(j)=0;
end

% Double check for AF(j)
Lead=II;
AF(j)=0;
B=0;
for i=1:length(FP)-1
    RR(i)=FP(i+1) - FP(i);
end
[R,m]=max(RR(1:length(RR)-2));
    From=FP(m+1)-20;
    To=From-4;
    [Pnt,Crit]=Flat_l(Lead,Crit,From,To,Flat);  % go to function
    Iz_=Pnt-Flat;
    From=FP(m);
    From=From+20;            			% from the point of the peak
    To=From+100*Hz/1000;                % to 100 ms right ...
    [Pnt,Crit]=flat_right(Lead,Crit,From,To,Flat);  % go to function
    J_=Pnt;
   
    Iz=J_; J=Iz_;
    
     figure(1); hold on
  %      plot(J/Hz, Lead(J), 'ro', Iz/Hz, Lead(Iz), 'ro')
    [Zero1]=Deriv(Lead,Iz,J); % go to function
    B=[Iz,Zero1,J]; 
           
    if length(B)>=12 && P_wave==0; AF(j)=1;
        figure(1);
        x1=(1:1:length(Lead))/Hz;
        plot(x1, Lead, B/Hz, Lead(B), 'go')
        title(II);
        else AF(j)=0;
    end
figure(2)

if AF(j)==1; I_AVB(j)=0; end

if PAC(j)~0; PAC(j)=1; end
if I_AVB(j) + AF(j) + PVC(j) + PAC(j) + LBBB(j) + RBBB(j) +...
        STD(j) + STE(j) ==0; Normal(j) = 1;
else Normal(j)=0;
end

% Conflictig are AF(j)-AVB, AF(j)-PAC(j), AF(j)-N, LBBB(j)-RBBB(j), AVB-RBBB(j) 
if AF(j)==1; PAC(j)=0; I_AVB(j)=0; end
if LBBB(j)==1; RBBB(j)=0; end 
if RBBB(j)==1; LBBB(j)=0; end

Text=['AF(j)=', num2str(AF(j)), '  I-AVB(j)=', num2str(I_AVB(j))', '  LBBB(j)=', num2str(LBBB(j)),...
     ' Nor(j)=', num2str(Normal(j)), '  PAC(j)=', num2str(PAC(j)), '  PVC(j)=', num2str(PVC(j)),...
     '  RBBB(j)=', num2str(RBBB(j)), '  STD(j)=', num2str(STD(j)), '  STE(j)=', num2str(STE(j))];
 fprintf('DIAGN: %s\n',Text);

% % %  save_cinc_ivaylo;
 



function plot(a,b,c,d,e,f,g,h,i,j,k,l,m)
end
function subplot(a,b,c,d)
end
function figure(a,b,c)
end
function title(a,b,c,d)
end
function clf
end
function hold (a,b,c)
end
