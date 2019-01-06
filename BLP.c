// 速度测试，采用Logistic映射产生diffusion密钥，每轮confuse使用3轮cat

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"

// 加密标准图像 lenna.raw，定义图像大小和灰度级别
#define WIDTH  512
#define LENGTH 512
#define GREY_LEVEL 256

// 混沌密钥流长度等于图像大小
#define KEY_STREAM_LENGTH WIDTH*LENGTH
//宏定义要加密的图像，以便于多轮加密
#define cipher_count 2

#define pic_o  "lenna.raw"
#define pic_en  "lenna_en.raw"
#define pic_de  "lenna_de.raw"


// arnold_trans 函数实现对图像的 arnold cat变换，img 为原始图像，a，b 为 arnold 变换参数，count 为变换次数
void arnold_trans(unsigned char (*img)[LENGTH], int a, int b, int count) {
    int loop, x, y, x_new, y_new, i, j;
    unsigned char img_temp[WIDTH][LENGTH];
    unsigned char *img_pointer, *img_temp_pointer;

    for(loop=0; loop<count; loop++) {
        for(x=0; x<WIDTH; x++)
            for(y=0; y<LENGTH; y++) {
                x_new = (x+a*y) % WIDTH;
                y_new = (b*x+(a*b+1)*y) % LENGTH;
                *(*(img_temp + x_new) + y_new) = *(*(img + x) + y);
            }
        img_pointer = img[0];
        img_temp_pointer = img_temp[0];
        for(i=0; i<WIDTH; i++)
            for(j=0; j<LENGTH; j++) {
                *img_pointer = *img_temp_pointer;
                img_pointer++;
                img_temp_pointer++;
        }
    }
}

// mod 函数用于提供逆 arnold 变换中使用的模函数，逆变换时由于存在负数，不能直接使用 C 提供的 % 模运算符
int mod(int a, int b) {
    if(a>=0)
        return a%b;
    else {
        while(a<0) {
            a=a+b;
        }
    }
    return a%b;
}

// 实现对图像的逆 arnold 变换 (inverse transform)，img 为待变换图像，a，b为 inv_arnold 变换参数，count 为变换次数
void inv_arnold_trans(unsigned char (*img)[LENGTH], int a, int b, int count) {
    int loop, x, y, x_new, y_new, i, j;
    unsigned char img_temp[WIDTH][LENGTH];
    unsigned char *img_pointer, *img_temp_pointer;

    for(loop=0; loop<count; loop++) {
        for(x=0; x<WIDTH; x++)
            for(y=0; y<LENGTH; y++) {
                x_new = mod(((a*b+1)*x-a*y), WIDTH);
                y_new = mod((-b*x+y), LENGTH);
                *(*(img_temp + x_new) + y_new) = *(*(img + x) + y);
            }
        img_pointer = img[0];
        img_temp_pointer = img_temp[0];
        for(i=0; i<WIDTH; i++)
            for(j=0; j<LENGTH; j++) {
                *img_pointer = *img_temp_pointer;
                img_pointer++;
                img_temp_pointer++;
        }
    }
}

void permutrans(unsigned char (*img)[LENGTH],double x0,double y0)
{
    unsigned char bitplane8[WIDTH][LENGTH];
    unsigned char bitplane7[WIDTH][LENGTH];
    unsigned char bitplane6[WIDTH][LENGTH];
    unsigned char bitplane5[WIDTH][LENGTH];
    unsigned char bitplane14[WIDTH][LENGTH];
    double miu=3.99999;
    int i,j,p14,q14,p5,q5,p6,q6,p7,q7,p8,q8;
    for(i=0;i<WIDTH;i++)
    for(j=0;j<LENGTH;j++)
    {
            bitplane14[i][j]=img[i][j]&15;
            bitplane5[i][j]=img[i][j]&16;
            bitplane6[i][j]=img[i][j]&32;
            bitplane7[i][j]=img[i][j]&64;
            bitplane8[i][j]=img[i][j]&128;
    }
    //首先计算比特级置乱的参数
    for (i=1;i<109;i++)
    {
        x0=miu*x0*(1-x0);
        y0=miu*y0*(1-y0);
        if(i==101)
        {
            p14=((long long)(x0*1.0e+14)) % 512;
            q14=((long long)(y0*1.0e+14)) % 512;
        }
        if(i==105)
        {
            p5=((long long)(x0*1.0e+14)) % 512;
            q5=((long long)(y0*1.0e+14)) % 512;
        }
        if(i==106)
        {
            p6=((long long)(x0*1.0e+14)) % 512;
            q6=((long long)(y0*1.0e+14)) % 512;
        }
        if(i==107)
        {
            p7=((long long)(x0*1.0e+14)) % 512;
            q7=((long long)(y0*1.0e+14)) % 512;
        }
        if(i==108)
        {
            p8=((long long)(x0*1.0e+14)) % 512;
            q8=((long long)(y0*1.0e+14)) % 512;
        }

    }

    arnold_trans(bitplane14, p14, q14, 1);
    arnold_trans(bitplane5, p5, q5, 1);
    arnold_trans(bitplane6, p6, q6, 1);
    arnold_trans(bitplane7, p7, q7, 1);
    arnold_trans(bitplane8, p8, q8, 1);
    for(i=0;i<WIDTH;i++)
    for(j=0;j<LENGTH;j++)
    {
        img[i][j]=bitplane14[i][j]+bitplane5[i][j]+bitplane6[i][j]+bitplane7[i][j]+bitplane8[i][j];
    }
}

void inv_permutrans(unsigned char (*img)[LENGTH],double x0,double y0)
{
    unsigned char bitplane8[WIDTH][LENGTH];
    unsigned char bitplane7[WIDTH][LENGTH];
    unsigned char bitplane6[WIDTH][LENGTH];
    unsigned char bitplane5[WIDTH][LENGTH];
    unsigned char bitplane14[WIDTH][LENGTH];
    double miu=3.99999;
    int i,j,p14,q14,p5,q5,p6,q6,p7,q7,p8,q8;
    for(i=0;i<WIDTH;i++)
    for(j=0;j<LENGTH;j++)
    {
            bitplane14[i][j]=img[i][j]&15;
            bitplane5[i][j]=img[i][j]&16;
            bitplane6[i][j]=img[i][j]&32;
            bitplane7[i][j]=img[i][j]&64;
            bitplane8[i][j]=img[i][j]&128;
    }
    //首先计算比特级置乱的参数
    for (i=1;i<109;i++)
    {
        x0=miu*x0*(1-x0);
        y0=miu*y0*(1-y0);
        if(i==101)
        {
            p14=((long long)(x0*1.0e+14)) % 512;
            q14=((long long)(y0*1.0e+14)) % 512;
        }
        if(i==105)
        {
            p5=((long long)(x0*1.0e+14)) % 512;
            q5=((long long)(y0*1.0e+14)) % 512;
        }
        if(i==106)
        {
            p6=((long long)(x0*1.0e+14)) % 512;
            q6=((long long)(y0*1.0e+14)) % 512;
        }
        if(i==107)
        {
            p7=((long long)(x0*1.0e+14)) % 512;
            q7=((long long)(y0*1.0e+14)) % 512;
        }
        if(i==108)
        {
            p8=((long long)(x0*1.0e+14)) % 512;
            q8=((long long)(y0*1.0e+14)) % 512;
        }

    }

    inv_arnold_trans(bitplane14, p14, q14, 1);
    inv_arnold_trans(bitplane5, p5, q5, 1);
    inv_arnold_trans(bitplane6, p6, q6, 1);
    inv_arnold_trans(bitplane7, p7, q7, 1);
    inv_arnold_trans(bitplane8, p8, q8, 1);
    for(i=0;i<WIDTH;i++)
    for(j=0;j<LENGTH;j++)
    {
        img[i][j]=bitplane14[i][j]+bitplane5[i][j]+bitplane6[i][j]+bitplane7[i][j]+bitplane8[i][j];
    }
}

void diffusion(unsigned char (*img)[LENGTH],double keyd)
 {
    int i;
    double tmp1,tmp2,xn;
    // 由于加密时某个点除了与密钥流元素相关，还与前一个点相关，对于第一个点，seed 用做在加密时的初始种子 c(-1), 由其生成 c(0)
    // key_stream 数组存储密钥流
    unsigned char *img_pointer,mask,previous_cipher;
    img_pointer=img[0];
    xn=keyd;
    // 以下开始加密操作
    tmp1=4*keyd*(1-keyd);
    previous_cipher=((long long)(tmp1*1.0e+14)) % 256;
    for(i=0;i<KEY_STREAM_LENGTH;i++)
    {
        xn=4*xn*(1-xn);
        tmp2=previous_cipher/1000;
        tmp1=4*tmp2*(1-tmp2);
        mask=((unsigned int) (tmp1*1000)+(long long)(xn*1.0e+10) % 256)%256;
        *img_pointer=(*img_pointer)^mask;
        previous_cipher=*img_pointer;
        img_pointer++;
    }

}
void inv_diffusion(unsigned char (*img)[LENGTH],double keyd)
 {
    int i;
    double tmp1,tmp2,xn=keyd;
    // 由于加密时某个点除了与密钥流元素相关，还与前一个点相关，对于第一个点，seed 用做在加密时的初始种子 c(-1), 由其生成 c(0)
    // key_stream 数组存储密钥流
    unsigned char *img_pointer,mask,previous_cipher,current_cipher;
    img_pointer=img[0];
    // 以下开始加密操作
    tmp1=4*keyd*(1-keyd);
    previous_cipher=((long long)(tmp1*1.0e+14)) % 256;
    for(i=0;i<KEY_STREAM_LENGTH;i++)
    {
        xn=4*xn*(1-xn);
        tmp2=previous_cipher/1000;
        tmp1=4*tmp2*(1-tmp2);
        mask=((unsigned int)(tmp1*1000)+(long long)(xn*1.0e+10) % 256)%256;
        current_cipher=*img_pointer;
        *img_pointer=(*img_pointer)^mask;
        previous_cipher=current_cipher;
        img_pointer++;
    }
}

int encryption(char *file) {

    // 读取明文图像
    FILE *img_fp;
    unsigned char img[WIDTH][LENGTH];
    int i;
    double x0 = 0.12345678912342,y0 = 0.87865765433212,keyd=0.34565487923280;
    if ((img_fp = fopen(file, "rb")) == NULL) {
        printf("cannot open file\n");
        return -1;
    }
    fread(img, sizeof(unsigned char), WIDTH*LENGTH, img_fp);

    // 关闭文件
    fclose(img_fp);

    // 对 img 进行 arnold 变换, 3轮 40, 8 为密钥
    for (i=0;i<cipher_count;i++)
    {
        permutrans(img,x0,y0);
        diffusion(img,keyd);
    }

    // 将加密图像写入磁盘
    if ((img_fp = fopen(pic_en, "wb")) == NULL) {
        printf("cannot open file to write\n");
        return -1;
    }
    fwrite(img, sizeof(unsigned char), WIDTH*LENGTH, img_fp);
    fclose(img_fp);
}

int decryption(char *file) {

    // 读取密文图像
    FILE *img_fp;
    unsigned char img[WIDTH][LENGTH];
    int i;
    double x0 = 0.12345678912342,y0 = 0.87865765433212,keyd=0.34565487923280;
    if ((img_fp = fopen(file, "rb")) == NULL) {
        printf("cannot open file\n");
        return -1;
    }
    fread(img, sizeof(unsigned char), WIDTH*LENGTH, img_fp);
    fclose(img_fp);

    // 对 confusion 后的图像进行逆 arnold 变换解密
    for (i=0;i<cipher_count;i++)
    {

        inv_diffusion(img,keyd);
        inv_permutrans(img,x0,y0);

    }
    // 将解密图像写入磁盘
    if ((img_fp = fopen(pic_de, "wb")) == NULL) {
        printf("cannot open file to write\n");
        return -1;
    }
    fwrite(img, sizeof(unsigned char), WIDTH*LENGTH, img_fp);
    fclose(img_fp);
}

int main(void) {

    long starter, ender;
    // 测试整个加密过程花费的时间
    int i;
    starter = clock();
    encryption(pic_o);
    ender = clock();
    printf("加密消耗%f秒\n", (double)(ender - starter)/CLOCKS_PER_SEC);

    starter = clock();
    decryption(pic_en);
    ender = clock();
    printf("解密消耗%f秒\n", (double)(ender - starter)/CLOCKS_PER_SEC);
    system("pause");

}
