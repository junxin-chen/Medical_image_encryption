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
#define cipher_count 1

#define pic_o  "CT_Head.raw"
#define pic_en  "CT_Head_d8.raw"


int encryption(char *file) {

    // 读取明文图像
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

    // 关闭文件
    fclose(img_fp);

    // 对 img 进行 arnold 变换, 3轮 40, 8 为密钥
    for (i=0;i<cipher_count;i++)
    {
        tmp=*(*(img + x) + y);
        tmp=tmp^1;
        *(*(img + x) + y)=tmp;
    }

    // 将加密图像写入磁盘
    if ((img_fp = fopen(pic_en, "wb")) == NULL) {
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
    //printf("加密消耗%f秒\n", (double)(ender - starter)/CLOCKS_PER_SEC);

}
