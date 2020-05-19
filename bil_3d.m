function []=bil_3d(vid,sigd,sigr)
nframes=vid.Duration*vid.FrameRate; %no. of frames
I=read(vid,nframes); %convert the last frame into image
I=rgb2gray(I);
I=im2double(I);
F=zeros([size(I,1) size(I,2) 1 nframes],class(I));%object to hold processed frames
F1=zeros([size(I,1) size(I,2) 1 nframes],class(I));
sum=0;
for i=1:nframes
    img=read(vid,i); %convert each frame to image
    img=rgb2gray(img);
    img=im2double(img);
    sum=sum+img;
end
avg=sum/nframes; %find average of frames
imshow(avg,[]);    
for n=1:nframes
    f=read(vid,n); 
f=rgb2gray(f);
f=im2double(f);
[r,c]=size(f);
f=f-avg; %subtract bakground from each frame
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
F1(:,:,:,n)=grad_img;
end

gt=zeros(1,K);
for i=1:K
    gt(i)=-(i-s)*exp(-(i-s)^2/(2*sigd^2));
end
for n=1:nframes-K+1 %(1st+s-1 th frame to nth frame -(s-1))
    img=zeros(r,c);
    for i=1:r
        for j=1:c
            for k=1:K
                a=F1(:,:,:,n+k-1);
                img(i,j)=img(i,j)+a(i,j)*gt(1,k);
            end
        end
    end
    F1(:,:,:,n+s-1)=img;
end

for n=1:nframes
grad_img1=F1(:,:,:,n);
new_img=zeros(r,c);
for i=2:r-1
    for j=2:c-1
        for k=1:3
            for l=1:3
                new_img(i,j)=new_img(i,j)+grad_img1(i+k-2,j+l-2)*(1.1-exp(-((f(i,j)-f(i+k-2,j+l-2))^2)/(2*sigr^2)));
            end
        end
    end
end
F1(:,:,:,n)=new_img;
end

for n=1:nframes-2
    img1=zeros(r,c);
    for i=1:r
        for j=1:c
            for k=1:3
                a1=F1(:,:,:,n+k-1);
                img1(i,j)=img1(i,j)+a1(i,j)*(1.1-exp(-(f(i,j)-a1(i,j))^2/(2*sigr^2)));
            end
        end
    end
    F1(:,:,n+1)=img1;
end

for n=1:nframes
    
non_max=F1(:,:,:,n);
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
for k=1:l
        v(k,1)=non_max(k);
end

g=hist(v,0:255);
p=g/numel(non_max);
p=im2double(p);
P=zeros(1,256);
for i=1:256
    for j=1:i
        P(1,i)=P(1,i)+p(1,j);
    end
end
m=zeros(1,256);
for i=1:256
    for j=1:i
        m(1,i)=m(1,i)+(j-1)*p(1,j);
    end
end
mg=0;
for i=1:256
    mg=mg+(i-1)*p(1,i);
end
sigb=zeros(1,256);
max=sigb(1,1);
for i=1:256
    sigb(1,i)=((mg*P(1,i)-m(1,i))^2)/(P(1,i)*(1-P(1,i)));
    
        if (sigb(1,i)>max)
            max=sigb(1,i);
            k=i-1;
        end
end
for i=1:r
    for j=1:c
        if(non_max(i,j)>k)
            final(i,j)=1;
        else
            final(i,j)=0;
        end
    end
end
F(:,:,:,n)=final;
end
framerate=vid.FrameRate;
implay(F,framerate);

end


