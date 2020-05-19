function [final] = bilateral_3(g,sigd,sigr)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if size(g,3)==1
    f=g; % if image is grayscale
elseif size(g,3)==3 %if image is rgb
f=rgb2gray(g); 
end
f=im2double(f);
[r,c]=size(f);
K=ceil(3*sigd)*2+1; %size of DoG filter
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
                new_img(i,j)=new_img(i+k-2,j+l-2)+grad_img(i+k-2,j+l-2)*(1.1-exp(-((f(i,j)-f(i+k-2,j+l-2))^2)/(2*sigr^2)));
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


final=otsu(non_max);


%figure,imshow(grad_img,[]);
%title('gradient');
%figure,imshow(new_img,[]);
%title('after bilateral');
%figure,imshow(non_max,[]);
%title('after nonmax suppression');
%figure,imshow(final,[]);
%title('edge map');
%final=im2double(final);
%g=im2double(g);
%g(:,:,1)=f;
%g(:,:,2)=g(:,:,1);
%g(:,:,3)=g(:,:,1);
%g(:,:,1)=final+g(:,:,1);
%g1=cat(3,g(:,:,1),g(:,:,2),g(:,:,3));


%figure,imshow(g1,[]);
%title('Superimposed');

%disp(k);





                
            
        
        

end