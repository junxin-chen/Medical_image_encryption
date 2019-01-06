// �ٶȲ��ԣ�����Logisticӳ�����diffusion��Կ��ÿ��confuseʹ��3��cat

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"

// ���ܱ�׼ͼ�� lenna.raw������ͼ���С�ͻҶȼ���
#define WIDTH  512
#define LENGTH 512
#define GREY_LEVEL 256

// ������Կ�����ȵ���ͼ���С
#define KEY_STREAM_LENGTH WIDTH*LENGTH
//�궨��Ҫ���ܵ�ͼ���Ա��ڶ��ּ���
#define cipher_count 2

#define pic_o  "lenna.raw"
#define pic_en  "lenna_en.raw"
#define pic_de  "lenna_de.raw"


// arnold_trans ����ʵ�ֶ�ͼ��� arnold cat�任��img Ϊԭʼͼ��a��b Ϊ arnold �任������count Ϊ�任����
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

// mod ���������ṩ�� arnold �任��ʹ�õ�ģ��������任ʱ���ڴ��ڸ���������ֱ��ʹ�� C �ṩ�� % ģ�����
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

// ʵ�ֶ�ͼ����� arnold �任 (inverse transform)��img Ϊ���任ͼ��a��bΪ inv_arnold �任������count Ϊ�任����
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
    //���ȼ�����ؼ����ҵĲ���
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
    //���ȼ�����ؼ����ҵĲ���
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
    // ���ڼ���ʱĳ�����������Կ��Ԫ����أ�����ǰһ������أ����ڵ�һ���㣬seed �����ڼ���ʱ�ĳ�ʼ���� c(-1), �������� c(0)
    // key_stream ����洢��Կ��
    unsigned char *img_pointer,mask,previous_cipher;
    img_pointer=img[0];
    xn=keyd;
    // ���¿�ʼ���ܲ���
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
    // ���ڼ���ʱĳ�����������Կ��Ԫ����أ�����ǰһ������أ����ڵ�һ���㣬seed �����ڼ���ʱ�ĳ�ʼ���� c(-1), �������� c(0)
    // key_stream ����洢��Կ��
    unsigned char *img_pointer,mask,previous_cipher,current_cipher;
    img_pointer=img[0];
    // ���¿�ʼ���ܲ���
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

    // ��ȡ����ͼ��
    FILE *img_fp;
    unsigned char img[WIDTH][LENGTH];
    int i;
    double x0 = 0.12345678912342,y0 = 0.87865765433212,keyd=0.34565487923280;
    if ((img_fp = fopen(file, "rb")) == NULL) {
        printf("cannot open file\n");
        return -1;
    }
    fread(img, sizeof(unsigned char), WIDTH*LENGTH, img_fp);

    // �ر��ļ�
    fclose(img_fp);

    // �� img ���� arnold �任, 3�� 40, 8 Ϊ��Կ
    for (i=0;i<cipher_count;i++)
    {
        permutrans(img,x0,y0);
        diffusion(img,keyd);
    }

    // ������ͼ��д�����
    if ((img_fp = fopen(pic_en, "wb")) == NULL) {
        printf("cannot open file to write\n");
        return -1;
    }
    fwrite(img, sizeof(unsigned char), WIDTH*LENGTH, img_fp);
    fclose(img_fp);
}

int decryption(char *file) {

    // ��ȡ����ͼ��
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

    // �� confusion ���ͼ������� arnold �任����
    for (i=0;i<cipher_count;i++)
    {

        inv_diffusion(img,keyd);
        inv_permutrans(img,x0,y0);

    }
    // ������ͼ��д�����
    if ((img_fp = fopen(pic_de, "wb")) == NULL) {
        printf("cannot open file to write\n");
        return -1;
    }
    fwrite(img, sizeof(unsigned char), WIDTH*LENGTH, img_fp);
    fclose(img_fp);
}

int main(void) {

    long starter, ender;
    // �����������ܹ��̻��ѵ�ʱ��
    int i;
    starter = clock();
    encryption(pic_o);
    ender = clock();
    printf("��������%f��\n", (double)(ender - starter)/CLOCKS_PER_SEC);

    starter = clock();
    decryption(pic_en);
    ender = clock();
    printf("��������%f��\n", (double)(ender - starter)/CLOCKS_PER_SEC);
    system("pause");

}
