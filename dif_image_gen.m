%对两个图像计算UACI及NPCR
clear all;
clc;
pix1=imread('en_CT_Head.bmp');
pix2=imread('en_dk4.bmp');
[M,N]=size(pix1);
pix1=double(pix1);
pix2=double(pix2);
pixd=pix1-pix2;
NPCR=nnz(pixd)/(M*N);
counter = 0;
  for i=1:1:M
     for j=1:1:N
       diff = abs(pix1(i,j) - pix2(i,j))/255;
       counter = counter + diff;
    end    
 end
uaci_ratio = counter / (M*N);
imwrite(pixd,'diffimage-dk4.tif');

