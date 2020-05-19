function[]=bilateral_vid1(vid,sigd,sigr)
nframes=vid.Duration*vid.FrameRate;
I=read(vid,nframes);
I=rgb2gray(I);
I=im2double(I);
F=zeros([size(I,1) size(I,2) 1 nframes],class(I));
sum=0;
for i=1:nframes
    img=read(vid,i);
    img=rgb2gray(img);
    img=im2double(img);
    sum=sum+img;
end
avg=sum/nframes;
s=0;    
for n=1:nframes
    f=read(vid,n);
f=rgb2gray(f);
f=im2double(f);
[r,c]=size(f);
f=f-avg;
K=ceil(6*sigd+1);
s=ceil((K/2));
gx=zeros(K,K);
gy=zeros(K,K);
for i=1:K
    for j=1:K
        gx(i,j)=-(j-s)*exp((-(i-s)^2-(j-s)^2)/(2*sigd^2));
        gy(i,j)=-(i-s)*exp((-(i-s)^2-(j-s)^2)/(2*sigd^2));
    end
end
pad_fig=zeros(r+K-1,c+K-1);
gradx=pad_fig;
grady=pad_fig;
for i=1:(r+K-1)
    for j=1:(c+K-1)
        if ((i<=(s-1))&&(j<=(s-1)))
            pad_fig(i,j)=f(1,1);
        elseif ((i<=(s-1))&&(j>(c+(s-1))))
            pad_fig(i,j)=f(1,c);
        elseif ((i>(r+s-1))&&(j<=s-1))
            pad_fig(i,j)=f(r,1);
        elseif ((i>(r+s-1))&&(j>(c+s-1)))
            pad_fig(i,j)=f(r,c);
        elseif ((i<=s-1)&&(j>s-1)&&(j<=(c+s-1)))
            pad_fig(i,j)=f(1,j-s+1);
        elseif ((i>(r+s-1))&&(j>s-1)&&(j<=(c+s-1)))
            pad_fig(i,j)=f(r,j-s+1);
        elseif ((i>s-1)&&(i<=(r+s-1))&&(j<=s-1))
            pad_fig(i,j)=f(i-s+1,1);
        elseif ((i>s-1)&&(i<=(r+s-1))&&(j>(c+s-1)))
            pad_fig(i,j)=f(i-s+1,c);
        else pad_fig(i,j)=f(i-s+1,j-s+1);
        end
    end
end
for i=1:r
    for j=1:c
        for k=1:K
            for l=1:K
                gradx(i+s-1,j+s-1)=gradx(i+s-1,j+s-1)+pad_fig(i+k-1,j+l-1)*gx(k,l);
                grady(i+s-1,j+s-1)=grady(i+s-1,j+s-1)+pad_fig(i+k-1,j+l-1)*gy(k,l);
            end
        end
    end
end
Gx=zeros(r,c);
Gy=Gx;
for i=1:r
    for j=1:c
        Gx(i,j)=gradx(i+s-1,j+s-1);
        Gy(i,j)=grady(i+s-1,j+s-1);
    end
end
        
        
g1=Gx.^2;
g2=Gy.^2;
g3=g1+g2;
grad_img=sqrt(g3);



new_img=zeros(r,c);
for i=2:r-1
    for j=2:c-1
        for k=1:3
            for l=1:3
                new_img(i,j)=new_img(i,j)+grad_img(i,j)*(1.1-exp(-((f(i,j)-f(i+k-2,j+l-2))^2)/(2*sigr^2)));
            end
        end
    end
end

     
non_max=new_img;
for i=2:r-1
    for j=2:c-1
        if (Gx==0)
            tangent=5;
        else
            tangent=(Gy(i,j)/Gx(i,j));
        
        end
        if(-0.4142<tangent && tangent<=0.4142)%horizontal
            if(new_img(i,j)<new_img(i,j+1) || new_img(i,j)<new_img(i,j-1))
                non_max(i,j)=0;
            end
        end
        if(0.4142<tangent && tangent<=2.4142)%45deg
            if(new_img(i,j)<new_img(i+1,j+1) || new_img(i,j)<new_img(i-1,j-1))
                non_max(i,j)=0;
            end
        end
        if( abs(tangent) >2.4142)%vertical
            if(new_img(i,j)<new_img(i-1,j) || new_img(i,j)<new_img(i+1,j))
                non_max(i,j)=0;
            end
        end
        if(-2.4142<tangent && tangent<= -0.4142)%-45deg
            if(new_img(i,j)<new_img(i-1,j+1) || new_img(i,j)<new_img(i+1,j-1))
                non_max(i,j)=0;
            end
        end
        
    end
end

final=non_max;
[r,c]=size(non_max);
l=r*c;
v=zeros(l,1);
k=0;
for i=1:r
    for j=1:c
        k=k+1;
        v(k,1)=non_max(i,j);
    end
end
g=hist(v,0:255);
p=g/numel(non_max);
P1=cumsum(p(1:254));
P1=P1(2:254);
P2=zeros(253);
P3=zeros(1,253);
a=256;
for i=3:255
    a=a-1;
    for j=3:255
        if(j<=a)
            for k=j:a
            P2(i-2,j-2)=P2(i-2,j-2)+p(1,k);
            end
        end
        if(j>a)
            break;
        end
    end
end
mg=0;
for i=1:256
    mg=mg+(i-1)*p(1,i);
end
for i=4:256
    for j=i:256
    P3(1,i-3)=P3(1,i-3)+p(1,j);
    end
end
    
m1=zeros(1,253);
for i=2:254
    for j=1:i
    m1(1,i-1)=m1(1,i-1)+(j-1)*p(1,j);
    end
    m1(1,i-1)=m1(1,i-1)/P1(1,i-1);
end
a=256;
m2=zeros(253);
m3=zeros(1,253);
for i=3:255
    a=a-1;
    for j=3:255
        if(j<=a)
            for k=j:a
                m2(i-2,j-2)=m2(i-2,j-2)+(k-1)*p(1,k);
            end
             m2(i-2,j-2)=m2(i-2,j-2)/P2(i-2,j-2);
        end
        if(j>a)
            break;
        end
    end
end
for i=4:256
    for j=i:256
        m3(1,i-3)=m3(1,i-3)+(j-1)*p(1,j);
    end
    m3(1,i-3)=m3(1,i-3)/P3(1,i-3);
end
sigb=zeros(253);
max=sigb(1,1);
k1=0;
for i=2:254
    k1=k1+1;
    k2=0;
    for j=i+1:255
        k2=k2+1;
        sigb(k1,k2)=P1(1,k1)*(m1(1,k1)-mg)^2+P2(k2,k1)*(m2(k2,k1)-mg)^2+P3(1,254-k2)*(m3(1,254-k2)-mg)^2;
        if(sigb(k1,k2)>max)
            max=sigb(k1,k2);
            K1=k1-1;
            K2=k2-1;
        end
    end
end
for i=1:r
    for j=1:c
    if(non_max(i,j)>K2/255)
        final(i,j)=1;
    elseif(non_max(i,j)>K1/255)
        final(i,j)=0.5;
    else final(i,j)=0;
    end
    end
end

flag=1;
while (flag==1)
    flag=0;
for i=2:r-1
    for j=2:c-1
        if (final(i,j)==0.5)
            if (final(i-1,j-1)==1||final(i-1,j)==1||final(i-1,j+1)==1||final(i,j-1)==1||final(i,j+1)==1||final(i+1,j-1)==1||final(i+1,j)==1||final(i+1,j+1)==1)
                final(i,j)=1;
                flag=1;
            else final(i,j)=0;
            end
        end
    end
end
end
F(:,:,:,n)=final;
end

framerate=vid.Framerate;
figure,implay(F,framerate);
imshow(bg,[]);
disp(K1);
disp(K2);
end

