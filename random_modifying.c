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

#define pic_o  "CT_Head.raw"
#define pic_en  "CT_Head_d8.raw"


int encryption(char *file) {

    // ��ȡ����ͼ��
    FILE *img_fp;
    unsigned char img[WIDTH][LENGTH];
    unsigned char x,y,tmp,i;
    x=380;
    y=200;

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
        tmp=*(*(img + x) + y);
        tmp=tmp^1;
        *(*(img + x) + y)=tmp;
    }

    // ������ͼ��д�����
    if ((img_fp = fopen(pic_en, "wb")) == NULL) {
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
    //printf("��������%f��\n", (double)(ender - starter)/CLOCKS_PER_SEC);

}
