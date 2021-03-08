%================================test3===================
%������˳���Ϊ2km,�����ʱ�ӷ�ΧΪ0-10��s
%α��㶨λ��λ����ȡ�����źŷ�����ʱ�估ʱ�䷶Χ�ڣ����źŵļ��
%Ϊ��֤�źż���ڼ�϶��Ҫôԭʼ�ź�����ǿ��С�Ĳ���ֱ���˳���Ҫô���˱��������ʹ������ź�ֻ��ת��������ź�
period=200e-6;
occupy=0.99;
pulseLength=period*occupy;
delta_t=0.1e-6;%������Ϊ1MHz
num=80000;
t=delta_t:delta_t:num*delta_t;
alpha1=1e-3;
GaussinCenter1=2e-3;
interregnum=4e-3;%���źŷ���ǰ����Ϊ3ms
GaussinCenter2=GaussinCenter1+interregnum;
GaussinWindow1=exp(-(t-GaussinCenter1).^2./alpha1^2);
GaussinWindow2=exp(-(t-GaussinCenter2).^2./alpha1^2);
delay1=5e-6;
delay2=12e-6;
delayRemain=10*max(abs(delay1),abs(delay2));
delayRemain=(delayRemain/delta_t);
tt=delta_t:delta_t:(num+delayRemain)*delta_t;
lamda1=1550e-9;
lamda2=1551e-9;
v1=3e8/lamda1;
v2=3e8/lamda2;
f1=10000; %����Ƶ��20kHz,����Ƶ��1MHz�����£��ź�Ƶ������Ϊ5kHz,���ܼ���λ����Ϣ
f2=10000; %Ƶ�ʹ��ͣ������ڼ�ⲻ���㹻����Ϣ��Ƶ�ʹ�����λ����ܳ�Խ�����࣬�����׼ȷ
modulation=pulse(pulseLength,tt,period);
intrusion1=1e-7*sin(2*pi*f1*t);
intrusion2=1e-7*sin(2*pi*f2*t);
intrusion1=exp(1i*2*pi*v1*t).*exp(1i*(2*pi/lamda1.*(intrusion1+2000)))+exp(1i*2*pi*v1*t).*exp(1i*(2*pi/lamda1*2000-pi));
intrusion2=exp(1i*2*pi*v2*t).*exp(1i*(2*pi/lamda2.*(intrusion2+2000)))+exp(1i*2*pi*v2*t).*exp(1i*(2*pi/lamda2*2000-pi));
intrusion1=intrusion1.^2;
intrusion2=intrusion2.^2;
y11=[intrusion1.*GaussinWindow1,zeros(1,int16(delayRemain))];
y12=[zeros(1,int16(delay1/delta_t)),intrusion1.*GaussinWindow1,zeros(1,int16(delayRemain)-int16(delay1/delta_t))];%�����ʱ���Ϊ0-10us����0-10delta_t
y21=[intrusion2.*GaussinWindow2,zeros(1,int16(delayRemain))];
y22=[zeros(1,int16(abs(delay2)/delta_t)),intrusion2.*GaussinWindow2,zeros(1,int16(delayRemain)-int16(abs(delay2)/delta_t))];%�����ʱ���Ϊ0-10us����0-10delta_t
if delay1<0&&delay2<0
signal1=abs(y11+y21);
signal2=abs(y12+y22);
end
if delay1<0&&delay2>0
signal1=abs(y11+y22);
signal2=abs(y12+y21);
end
if delay1>0&&delay2<0
signal1=abs(y12+y21);
signal2=abs(y11+y22);
end
if delay1>0&&delay2>0
signal1=abs(y12+y22);
signal2=abs(y11+y21);
end
modulatedSignal1=signal1.*modulation;
modulatedSignal2=signal2.*modulation;
figure(1)
subplot(211)
plot(tt,signal1,tt,signal2)
xlabel("t/s")
title("ԭʼ�ź�")
subplot(212)
plot(tt,modulatedSignal1,tt,modulatedSignal2)
xlabel("t/s")
title("�����ź�")
figure(2)
tcorr=-(num+delayRemain-1)*delta_t:delta_t:(num+delayRemain-1)*delta_t;
subplot(211)
plot(tcorr,xcorr(signal1,signal2))
title("ԭʼ�źŻ���ؽ��")
xlabel("t/s")
subplot(212)
tic
plot(tcorr,xcorr(modulatedSignal1,modulatedSignal2))
xlabel("t/s")
title("ģ���źŻ���ؽ��")
toc
%��ѡ����
n=period/delta_t;
m=floor(num*delta_t/period);
partT=delta_t:delta_t:n*delta_t;
partTcorr=-(n-1)*delta_t:delta_t:(n-1)*delta_t;
n=int32(n);
A1=reshape(modulatedSignal1(1:m*n),n,m);
A2=reshape(modulatedSignal2(1:m*n),n,m);
recordTime=period:period:m*period;
record=zeros(1,m);
for i=1:m
%     figure(4)
    part1=xcorr(A1(:,i),A2(:,i));
%     plot(partTcorr,part1)
    value=max(part1);
    amplitude=max(max(A1(:,i)),max(A2(:,i)));
    result=partTcorr((part1==value));
    result=result(1);
    if amplitude>0.05
    record(i)=result;
    end
end

figure(3)
subplot(211)
plot(tt,modulatedSignal1,tt,modulatedSignal2)
axis([0,(num+delayRemain)*delta_t,min(modulatedSignal1),max(modulatedSignal2)])
title("ԭʼ�ź�")
xlabel("t/s")
subplot(212)
stem(recordTime,record)
title("�������廥��ؽ��")
xlabel("t/s")
%����ն�Ϻ������ƽ��
%����������źų����쳣�����̶Ժ��������źŽ��в�����أ������ֽϴ�����Ϊ�ö�ʱ��ͬʱ�����ĵڶ�����

%�ƶ����崰����ʱ������������ʱ����=================================��
[m,n]=size(signal1);
period_n=int32(period/delta_t);
newrecord=100e-6*ones(1,n-period_n-1);
for i=1:n-period_n-1
    part1=signal1(i:i+period_n-1);
    part2=signal2(i:i+period_n-1);
    result=xcorr(part1,part2);
    value=max(result);
    delay=partTcorr((result==value));
    delay=delay(1);
    newrecord(i)=delay;
end
figure(5)
chart=-40*delta_t:0.2*delta_t:40*delta_t-0.2*delta_t;
table=zeros(1,400);
[~,n]=size(newrecord);
for i=1:1:400
    count=0;
    for j=1:n
        if int16(newrecord(j)*1e6)==int16(chart(i)*1e6)
            count=count+1;
        end
    end
    table(i)=count;
end
subplot(121)
plot(chart,table)
xlabel("t/us")
title("��λ1")
[~,Table]=meshgrid(1:5,table);
subplot(122)
imagesc(1:5,chart,Table)
xlabel("al")
ylabel("t/us")
title("��λ2")



function result=pulse(pulseLength,time,period)
result=zeros(size(time));
for t=time
    shadow=mod(t,period);
    if shadow<pulseLength
        result((time==t))=1;
    else
        result((time==t))=0;
    end
end
end