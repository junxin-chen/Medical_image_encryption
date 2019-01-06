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
#define cipher_count 1
#define permut_count 2

#define pic_o  "lenna.raw"
#define pic_en  "lenna_en.raw"
#define pic_de  "lenna_de.raw"

unsigned char bitplane8[WIDTH][LENGTH];
unsigned char bitplane7[WIDTH][LENGTH];
unsigned char bitplane6[WIDTH][LENGTH];
unsigned char bitplane5[WIDTH][LENGTH];
unsigned char bitplane4[WIDTH][LENGTH];
unsigned char bitplane3[WIDTH][LENGTH];
unsigned char bitplane2[WIDTH][LENGTH];
unsigned char bitplane1[WIDTH][LENGTH];


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

void permutrans(unsigned char (*img)[LENGTH],double miux, double x0,double miuy,double y0)
{


    int i,j,p1,q1,p2,q2,p3,q3,p4,q4,p5,q5,p6,q6,p7,q7,p8,q8;

    for(i=0;i<WIDTH;i++)
    for(j=0;j<LENGTH;j++)
    {
            bitplane1[i][j]=img[i][j]&1;
            bitplane2[i][j]=img[i][j]&2;
            bitplane3[i][j]=img[i][j]&4;
            bitplane4[i][j]=img[i][j]&8;
            bitplane5[i][j]=img[i][j]&16;
            bitplane6[i][j]=img[i][j]&32;
            bitplane7[i][j]=img[i][j]&64;
            bitplane8[i][j]=img[i][j]&128;
    }
    //���ȼ�����ؼ����ҵĲ���
    for (i=1;i<109;i++)
    {
        x0=miux*x0*(1-x0);
        y0=miuy*y0*(1-y0);
        if(i==101)
        {
            p1=((long long)(x0*1.0e+14)) % 512;
            q1=((long long)(y0*1.0e+14)) % 512;
        }
        {
            p2=((long long)(x0*1.0e+14)) % 512;
            q2=((long long)(y0*1.0e+14)) % 512;
        }
        {
            p3=((long long)(x0*1.0e+14)) % 512;
            q3=((long long)(y0*1.0e+14)) % 512;
        }
        {
            p4=((long long)(x0*1.0e+14)) % 512;
            q4=((long long)(y0*1.0e+14)) % 512;
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

    arnold_trans(bitplane1, p1, q1, permut_count);
    arnold_trans(bitplane2, p2, q2, permut_count);
    arnold_trans(bitplane3, p3, q3, permut_count);
    arnold_trans(bitplane4, p4, q4, permut_count);
    arnold_trans(bitplane5, p5, q5, permut_count);
    arnold_trans(bitplane6, p6, q6, permut_count);
    arnold_trans(bitplane7, p7, q7, permut_count);
    arnold_trans(bitplane8, p8, q8, permut_count);
    for(i=0;i<WIDTH;i++)
    for(j=0;j<LENGTH;j++)
    {
        img[i][j]=img[i][j]=bitplane1[i][j]+bitplane2[i][j]+bitplane3[i][j]+bitplane4[i][j]+bitplane5[i][j]+bitplane6[i][j]+bitplane7[i][j]+bitplane8[i][j];
    }
}

void inv_permutrans(unsigned char (*img)[LENGTH],double miux, double x0,double miuy,double y0)
{


    int i,j,p1,q1,p2,q2,p3,q3,p4,q4,p5,q5,p6,q6,p7,q7,p8,q8;

    for(i=0;i<WIDTH;i++)
    for(j=0;j<LENGTH;j++)
    {
            bitplane1[i][j]=img[i][j]&1;
            bitplane2[i][j]=img[i][j]&2;
            bitplane3[i][j]=img[i][j]&4;
            bitplane4[i][j]=img[i][j]&8;
            bitplane5[i][j]=img[i][j]&16;
            bitplane6[i][j]=img[i][j]&32;
            bitplane7[i][j]=img[i][j]&64;
            bitplane8[i][j]=img[i][j]&128;
    }
    //���ȼ�����ؼ����ҵĲ���
    for (i=1;i<109;i++)
    {
        x0=miux*x0*(1-x0);
        y0=miuy*y0*(1-y0);
        if(i==101)
        {
            p1=((long long)(x0*1.0e+14)) % 512;
            q1=((long long)(y0*1.0e+14)) % 512;
        }
        {
            p2=((long long)(x0*1.0e+14)) % 512;
            q2=((long long)(y0*1.0e+14)) % 512;
        }
        {
            p3=((long long)(x0*1.0e+14)) % 512;
            q3=((long long)(y0*1.0e+14)) % 512;
        }
        {
            p4=((long long)(x0*1.0e+14)) % 512;
            q4=((long long)(y0*1.0e+14)) % 512;
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

    inv_arnold_trans(bitplane1, p1, q1, permut_count);
    inv_arnold_trans(bitplane2, p2, q2, permut_count);
    inv_arnold_trans(bitplane3, p3, q3, permut_count);
    inv_arnold_trans(bitplane4, p4, q4, permut_count);
    inv_arnold_trans(bitplane5, p5, q5, permut_count);
    inv_arnold_trans(bitplane6, p6, q6, permut_count);
    inv_arnold_trans(bitplane7, p7, q7, permut_count);
    inv_arnold_trans(bitplane8, p8, q8, permut_count);
    for(i=0;i<WIDTH;i++)
    for(j=0;j<LENGTH;j++)
    {
        img[i][j]=bitplane1[i][j]+bitplane2[i][j]+bitplane3[i][j]+bitplane4[i][j]+bitplane5[i][j]+bitplane6[i][j]+bitplane7[i][j]+bitplane8[i][j];
    }
}

int encryption(char *file) {

    // ��ȡ����ͼ��
    FILE *img_fp;
    unsigned char img[WIDTH][LENGTH];
    int i;
    double miux=3.9991234, x0 = 0.12345678912342,miuy=3.9995678,y0 = 0.87865765433212;
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
        permutrans(img,miux,x0,miuy,y0);
        //diffusion(img,keyd);
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
    double miux=3.9991234, x0 = 0.12345678912342,miuy=3.9995678,y0 = 0.87865765433212;
    if ((img_fp = fopen(file, "rb")) == NULL) {
        printf("cannot open file\n");
        return -1;
    }
    fread(img, sizeof(unsigned char), WIDTH*LENGTH, img_fp);
    fclose(img_fp);

    // �� confusion ���ͼ������� arnold �任����
    for (i=0;i<cipher_count;i++)
    {

        //inv_diffusion(img,keyd);
        inv_permutrans(img,miux,x0,miuy,y0);

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
