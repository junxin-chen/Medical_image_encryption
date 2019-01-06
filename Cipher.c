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

#define pic_en1  "en_CT_Abdomen.raw"
#define pic_en2  "en_CT_Head.raw"
#define pic_en3  "en_CT_Paranasal_sinus.raw"
#define pic_en4  "en_MR_Cervical_vertebra.raw"
#define pic_en5  "en_MR_Knee.raw"
#define pic_en6  "en_MR_Prostate.raw"
#define pic_en7  "en_MR_Waist.raw"
#define pic_en8  "en_X_Lungs.raw"



double var1[KEY_STREAM_LENGTH] = {0}; //随机变量流1
double var2[KEY_STREAM_LENGTH] = {0}; //随机变量流2
double var3[KEY_STREAM_LENGTH] = {0}; //随机变量流3
double var4[KEY_STREAM_LENGTH] = {0}; //随机变量流4


unsigned int pv1[KEY_STREAM_LENGTH] = {0};//置乱向量1
unsigned int pv2[KEY_STREAM_LENGTH] = {0};//置乱向量2
unsigned int pv1_tmp[KEY_STREAM_LENGTH] = {0};//置乱向量1
unsigned int pv2_tmp[KEY_STREAM_LENGTH] = {0};//置乱向量2

unsigned char mask1[KEY_STREAM_LENGTH] = {0};//掩码1
unsigned char mask2[KEY_STREAM_LENGTH] = {0};//掩码2


void generate_variables(double x0,double y0,double z0,double w0)
{
	double h=0.0005,a = 12,b = 23,c = 1,d = 2.1,k = 6,r = 0.2;   //步进值及控制参数
	double K1,K2,K3,K4,L1,L2,L3,L4,M1,M2,M3,M4,N1,N2,N3,N4;
	int i=0;

	var1[0]=x0;var2[0]=y0;var3[0]=z0;var4[0]=w0;

	for (i=0;i<300;i++)//热身200次
	{
        K1 = a * (var2[i] - var1[i]);
        L1 = b * var1[i] - var1[i] * var3[i] - c * var2[i] + var4[i];
        M1 = var1[i] * var2[i] - d * var3[i];
        N1 = -k * var2[i] - r * var4[i];

        K2 = a * ((var2[i] + h/2 * L1) - (var1[i] + h/2 * K1));
        L2 = b * (var1[i] + h/2 * K1) - (var1[i] + h/2 * K1) * (var3[i] + h/2 * M1) - c * (var2[i] + h/2 * L1) + (var4[i] + h/2 * N1);
        M2 = (var1[i] + h/2 * K1) * (var2[i] + h/2 * L1) - d * (var3[i] + h/2 * M1);
        N2 = -k * (var2[i] + h/2 * L1) - r * (var4[i] + h/2 * N1);

        K3 = a * ((var2[i] + h/2 * L2) - (var1[i] + h/2 * K2));
        L3 = b * (var1[i] + h/2 * K2) - (var1[i] + h/2 * K2) * (var3[i] + h/2 * M2) - c * (var2[i] + h/2 * L2) + (var4[i] + h/2 * N2);
        M3 = (var1[i] + h/2 * K2) * (var2[i] + h/2 * L2) - d * (var3[i] + h/2 * M2);
        N3 = -k * (var2[i] + h/2 * L2) - r * (var4[i] + h/2 * N2);

        K4 = a * ((var2[i] + h * L3) - (var1[i] + h * K3));
        L4 = b * (var1[i] + h * K3) - (var1[i] + h * K3) * (var3[i] + h * M3) - c * (var2[i] + h * L3) + (var4[i] + h * N3);
        M4 = (var1[i] + h * K3) * (var2[i] + h * L3) - d * (var3[i] + h * M3);
        N4 = -k * (var2[i] + h * L3) - r * (var4[i] + h * N3);

        var1[i+1] = var1[i] + h/6 * (K1 + 2*K2 + 2*K3 + K4);
        var2[i+1] = var2[i] + h/6 * (L1 + 2*L2 + 2*L3 + L4);
        var3[i+1] = var3[i] + h/6 * (M1 + 2*M2 + 2*M3 + M4);
        var4[i+1] = var4[i] + h/6 * (N1 + 2*N2 + 2*N3 + N4);

	}
	var1[0]=var1[300];
	var2[0]=var2[300];
	var3[0]=var3[300];
	var4[0]=var4[300];

	for (i=0;i<KEY_STREAM_LENGTH;i++)// 生成混沌序列
	{
        K1 = a * (var2[i] - var1[i]);
        L1 = b * var1[i] - var1[i] * var3[i] - c * var2[i] + var4[i];
        M1 = var1[i] * var2[i] - d * var3[i];
        N1 = -k * var2[i] - r * var4[i];

        K2 = a * ((var2[i] + h/2 * L1) - (var1[i] + h/2 * K1));
        L2 = b * (var1[i] + h/2 * K1) - (var1[i] + h/2 * K1) * (var3[i] + h/2 * M1) - c * (var2[i] + h/2 * L1) + (var4[i] + h/2 * N1);
        M2 = (var1[i] + h/2 * K1) * (var2[i] + h/2 * L1) - d * (var3[i] + h/2 * M1);
        N2 = -k * (var2[i] + h/2 * L1) - r * (var4[i] + h/2 * N1);

        K3 = a * ((var2[i] + h/2 * L2) - (var1[i] + h/2 * K2));
        L3 = b * (var1[i] + h/2 * K2) - (var1[i] + h/2 * K2) * (var3[i] + h/2 * M2) - c * (var2[i] + h/2 * L2) + (var4[i] + h/2 * N2);
        M3 = (var1[i] + h/2 * K2) * (var2[i] + h/2 * L2) - d * (var3[i] + h/2 * M2);
        N3 = -k * (var2[i] + h/2 * L2) - r * (var4[i] + h/2 * N2);

        K4 = a * ((var2[i] + h * L3) - (var1[i] + h * K3));
        L4 = b * (var1[i] + h * K3) - (var1[i] + h * K3) * (var3[i] + h * M3) - c * (var2[i] + h * L3) + (var4[i] + h * N3);
        M4 = (var1[i] + h * K3) * (var2[i] + h * L3) - d * (var3[i] + h * M3);
        N4 = -k * (var2[i] + h * L3) - r * (var4[i] + h * N3);

        var1[i+1] = var1[i] + h/6 * (K1 + 2*K2 + 2*K3 + K4);
        var2[i+1] = var2[i] + h/6 * (L1 + 2*L2 + 2*L3 + L4);
        var3[i+1] = var3[i] + h/6 * (M1 + 2*M2 + 2*M3 + M4);
        var4[i+1] = var4[i] + h/6 * (N1 + 2*N2 + 2*N3 + N4);
	}
}

void generate_pv_and_mask()
{
    unsigned int i,pv_tmp;
    double temp1,temp2,temp3,temp4;

    //按序生成pv1和pv2，以及pixel swapping的交换矩阵
    for (i=0;i<KEY_STREAM_LENGTH;i++)
    {
        temp1=fabs(var1[i])-floor(fabs(var1[i]));
        pv1_tmp[i]=i+((long long)(temp1*1.0e+14)) % (KEY_STREAM_LENGTH-i);

        temp2=fabs(var2[i])-floor(fabs(var2[i]));
        pv2_tmp[i]=i+((long long)(temp2*1.0e+14)) % (KEY_STREAM_LENGTH-i);

        pv1[i]=i;
        pv2[i]=i;

        //if (i>262100)
        //printf("var1:%f, temp1:%f, pv1_tmp:%u,  var2:%10f, temp2:%f, pv2:%u\n",var1[i],temp1,pv1_tmp[i],var2[i],temp2,pv2_tmp[i]);
        //printf("pv1:%u, pv2:%u\n",pv1[i],pv2[i]);
    }
    //使用交换矩阵生成最红的pv
    for (i=0;i<KEY_STREAM_LENGTH;i++)
    {
        pv_tmp=pv1[i];
        pv1[i]=pv1[pv1_tmp[i]];
        pv1[pv1_tmp[i]]=pv_tmp;

        pv_tmp=pv2[i];
        pv2[i]=pv2[pv2_tmp[i]];
        pv2[pv2_tmp[i]]=pv_tmp;

        //if (i>262100)
        //printf("pv1:%u, pv2:%u\n",pv1[i],pv2[i]);
    }

    //生成用于扩散的mask
    for(i=0;i<KEY_STREAM_LENGTH;i++)
    {
    	temp3=fabs(var3[i])-floor(fabs(var3[i]));
        mask1[i]= ((long long)(temp3*1.0e+14)) % GREY_LEVEL;

        temp4=fabs(var4[i])-floor(fabs(var4[i]));
        mask2[i]= ((long long)(temp4*1.0e+14)) % GREY_LEVEL;
        //if (i>262100)
        //printf("var3:%f, temp3:%f, mask1:%d,var4:%10f, temp4:%f, mask2:%d\n",var3[i],temp3,mask1[i],var4[i],temp4,mask2[i]);
        //{tt=0.934;  printf("%lf\n",fabs(tt));}
     }
}

unsigned cirshift(unsigned char p,char direction,char loop)
{
	int i,j;
	unsigned char temp=0;//用于取出左侧的值
	//如果移位位数为0，直接返回
	if (loop==0)
	return p;
	//否则执行以下命令
	if(direction == 0)  //如果循环左移
	{
		for(i=0;i<loop;i++)//每左移移位为1个循环步骤
		{
			temp=p&128;//取出循环左移时要消失的第一位，取到temp的第一位
			p=p<<1;//P左移1位，则最高位消失，其他位左移，最低位补零
			temp=temp>>7;//temp右移7位，将取出的最高位移到最低位
			p=p|temp;//将temp的最低位补到p的最低位，由于temp的其他位均为0，故按位或的时候不影响其他位
		}
		return p;
	}//循环左移完成
	if(direction == 1)//如果循环右移
	{
		for(i=0;i<loop;i++)
		{
			temp=p&1;//取出循环右移时要消失的最低位，取到temp的最低位
			p=p>>1;//P右移1位，则最低位消失，其他位右移，最高位补零
			temp=temp<<7;//temp左移7位，将取出的p的最低位移至最高位
			p=p|temp;//将temp的最高位补到p的最高位，由于temp的其他位为0，按位或的时候不影响其他位
		}
		return p;
	}
	printf("cirshift error\n");
}

void permutation(unsigned char (*img)[LENGTH])
{
    unsigned char img_temp[WIDTH][LENGTH];
    unsigned char *img_pointer, *img_temp_pointer, current_plain,previous_plain,previous_cipher;

    char direction,loop;
    unsigned int i, pv1i,pv2i;

    img_pointer = img[0];
    img_temp_pointer = img_temp[0];

    previous_plain=123;
    previous_cipher=123;

    for (i=0;i<KEY_STREAM_LENGTH;i++)
    {
        pv1i=pv1[i];
        pv2i=pv2[i];
        //先进行置乱部分
        //先获得移位的参数
        direction=pv2i&1;//移位的方向由pv2的最后一位决定
        loop=pv1i&7;//移位的位数由pv1的后三位决定
        current_plain=*(img_pointer+pv1i);
        current_plain=cirshift(current_plain,direction,loop);
        current_plain=current_plain^previous_plain;
        *(img_temp_pointer+pv2i)=current_plain;
        previous_plain=*(img_pointer+pv1i);
        //置乱加密完成
        }

    for (i=0;i<KEY_STREAM_LENGTH;i++)
    *(img_pointer+i) = *(img_temp_pointer+i);
}

void inv_permutation(unsigned char (*img)[LENGTH])
{
    unsigned char img_temp[WIDTH][LENGTH];
    unsigned char *img_pointer, *img_temp_pointer,current_cipher,previous_plain;

    char direction,loop;
    unsigned int i, pv1i,pv2i;

    img_pointer = img[0];
    img_temp_pointer = img_temp[0];

    previous_plain=123;
    for (i=0;i<KEY_STREAM_LENGTH;i++)
    {
        pv1i=pv1[i];
        pv2i=pv2[i];

        //开始置乱解密
        direction=(pv2i^1)&1;//移位的方向由pv2的最后一位决定
        loop=pv1i&7;//移位的位数由pv1的后三位决定
        current_cipher=*(img_pointer+pv2i);
        current_cipher=current_cipher^previous_plain;
        current_cipher=cirshift(current_cipher,direction,loop);
        *(img_temp_pointer+pv1i)=current_cipher;
        previous_plain=current_cipher;
    }

    for (i=0;i<KEY_STREAM_LENGTH;i++)
    *(img_pointer+i) = *(img_temp_pointer+i);

}

void diffusion(unsigned char (*img)[LENGTH])
{
    unsigned char img_temp[WIDTH][LENGTH];
    unsigned char *img_pointer, *img_temp_pointer, current_plain,previous_cipher,current_cipher;

    unsigned int i, pv1i,pv2i;
    unsigned char mask1i,mask2i;
    unsigned char t1,t2,t3;

    img_pointer = img[0];
    img_temp_pointer = img_temp[0];

    previous_cipher=123;

    for (i=0;i<KEY_STREAM_LENGTH;i++)
    {
        pv1i=pv1[i];
        pv2i=pv2[i];
        mask1i=mask1[i];
        mask2i=mask2[i];
        current_plain=*(img_pointer+pv1i);
        t1=(mask1i+mask2i)%GREY_LEVEL;
        t2=mask2i^current_plain^previous_cipher;
        t3=(mask1i+t2)%GREY_LEVEL;
        current_cipher=t1^t3;
        previous_cipher=current_cipher;
        *(img_temp_pointer+pv2i)=current_cipher;
        }

    for (i=0;i<KEY_STREAM_LENGTH;i++)
    *(img_pointer+i) = *(img_temp_pointer+i);
}

void inv_diffusion(unsigned char (*img)[LENGTH])
{
    unsigned char img_temp[WIDTH][LENGTH];
    unsigned char *img_pointer, *img_temp_pointer, current_plain,previous_cipher,current_cipher;

    unsigned int i, pv1i,pv2i;
    unsigned char mask1i,mask2i;
    unsigned char t1,t2,t3;

    img_pointer = img[0];
    img_temp_pointer = img_temp[0];

    previous_cipher=123;

    for (i=0;i<KEY_STREAM_LENGTH;i++)
    {
        pv1i=pv1[i];
        pv2i=pv2[i];
        mask1i=mask1[i];
        mask2i=mask2[i];
        current_cipher=*(img_pointer+pv2i);
        t1=(mask1i+mask2i)%GREY_LEVEL;
        t2=current_cipher^t1;
        t3=(t2-mask1i+GREY_LEVEL)%GREY_LEVEL;
        current_plain=t3^previous_cipher^mask2i;
        previous_cipher=current_cipher;
        *(img_temp_pointer+pv1i)=current_plain;
        }

    for (i=0;i<KEY_STREAM_LENGTH;i++)
    *(img_pointer+i) = *(img_temp_pointer+i);
}

void read_file(char *file, unsigned char (*img)[LENGTH])
{
    //仅用于读取文件
    FILE *img_fp;
    if ((img_fp = fopen(file, "rb")) == NULL) {
        printf("cannot open file for open");
        printf(" %s\n",file);
        return;
    }
    fread(img, sizeof(unsigned char), WIDTH*LENGTH, img_fp);
    fclose(img_fp);
}

void write_file(char *file, unsigned char (*img)[LENGTH])
{
    //仅用于读取文件
    FILE *img_fp;
    if ((img_fp = fopen(file, "wb")) == NULL) {
        printf("cannot open file for write");
        printf(" %s\n",file);
        return;
    }
    fwrite(img, sizeof(unsigned char), WIDTH*LENGTH, img_fp);
    fclose(img_fp);
}

int encryption(unsigned char (*img)[LENGTH],double x0,double y0,double z0,double w0)
{
    generate_variables(x0,y0,z0,w0);
    generate_pv_and_mask();
    permutation(img);
    //permutation(img);
    diffusion(img);
}

int decryption(unsigned char (*img)[LENGTH],double x0,double y0,double z0,double w0)
{
    generate_variables(x0,y0,z0,w0);
    generate_pv_and_mask();
    inv_diffusion(img);
    inv_permutation(img);
    //inv_permutation(img);
}


void npcruaci(unsigned char (*img)[LENGTH],unsigned char (*img2)[LENGTH],double *npcr, double *uaci)
{
    int i,j;
    double counter=0,diff=0,p1,p2,pd;
    for (i=0;i<WIDTH;i++)
    for (j=0;j<LENGTH;j++)
    {
        p1=img[i][j];
        p2=img2[i][j];
        if(p1!=p2)
        counter++;
        pd=fabs(p1-p2)/255;
        diff=diff+pd;
    }
    *npcr=(100*counter)/(WIDTH*LENGTH);
    *uaci=(100*diff)/(WIDTH*LENGTH);
}

void entropy(unsigned char *data,double *shang,double *kafang)
{
	int i, totalpixel;
	double Cal[256] = {0},Cal2[256] = {0},tmp;
	double entr = 0,chip=0;
	unsigned char *mdata = data;

	totalpixel=WIDTH*LENGTH;

	for(i=0; i<totalpixel; i++){
		Cal[mdata[i]]+=1;
	}

    //计算信息熵
	for(i=0;i<256;i++)
	{
		Cal2[i]=Cal[i];
		Cal[i] /= totalpixel;
		if(Cal[i]==0)
			entr = entr;
		else
			entr = entr - Cal[i]*(log(Cal[i])/log(2));
	}
	//计算卡方分布
	for(i=0;i<256;i++)
	{
		//printf("%d,%f\n",i,Cal[i]);
		tmp=fabs(Cal2[i]-1024);
		chip=chip+(tmp*tmp)/1024;
	}

	*shang=entr;
	*kafang=chip;
}

void load_npcr(double k1, double k2, double k3, double k4)
{
    //计算并打印npcr和UACi情况
    int i;
    unsigned char img[WIDTH][LENGTH];
    unsigned char img_d[WIDTH][LENGTH];
    double npcr, uaci;
    //待测试的文件对
    char *file1="CT_Abdomen.raw",   *file1_d="CT_Abdomen_d.raw";
    char *file2="CT_Head.raw",       *file2_d="CT_Head_d.raw";
    char *file3="CT_Paranasal_sinus.raw",    *file3_d="CT_Paranasal_sinus_d.raw";
    char *file4="MR_Cervical_vertebra.raw",      *file4_d="MR_Cervical_vertebra_d.raw";
    char *file5="MR_Knee.raw",      *file5_d="MR_Knee_d.raw";
    char *file6="MR_Prostate.raw",      *file6_d="MR_Prostate_d.raw";
    char *file7="MR_Waist.raw",      *file7_d="MR_Waist_d.raw";
    char *file8="X_Lungs.raw",      *file8_d="X_Lungs_d.raw";
    //计算4轮内的NPCR和UACI
    printf("===========================NPCR and UACI TEST=============================\n");
    printf("%-15s,  %-17s,%-17s,%-17s,%-17s\n","Test img","Round 1","Round 2","Round 3","Round 4");
    //第一对图片
    printf("%-15s,",file1);
    //读取明文图像
    read_file(file1,img);
    read_file(file1_d,img_d);
    for (i=1;i<5;i++)
    {
        encryption(img,k1,k2,k3,k4);
        encryption(img_d,k1,k2,k3,k4);
        npcruaci(img,img_d,&npcr,&uaci);
        printf("%8.4f,%-8.4f,", npcr,uaci);
    }
    printf("\n");
    //第二对图片
    printf("%-15s,",file2);
    read_file(file2,img);
    read_file(file2_d,img_d);
    for (i=1;i<5;i++)
    {
        encryption(img,k1,k2,k3,k4);
        encryption(img_d,k1,k2,k3,k4);
        npcruaci(img,img_d,&npcr,&uaci);
        printf("%8.4f,%-8.4f,", npcr,uaci);
    }
    printf("\n");
    //第三对图片
    printf("%-15s,",file3);
    read_file(file3,img);
    read_file(file3_d,img_d);
    for (i=1;i<5;i++)
    {
        encryption(img,k1,k2,k3,k4);
        encryption(img_d,k1,k2,k3,k4);
        npcruaci(img,img_d,&npcr,&uaci);
        printf("%8.4f,%-8.4f,", npcr,uaci);
    }
    printf("\n");
    //第四对图片
    printf("%-15s,",file4);
    read_file(file4,img);
    read_file(file4_d,img_d);
    for (i=1;i<5;i++)
    {
        encryption(img,k1,k2,k3,k4);
        encryption(img_d,k1,k2,k3,k4);
        npcruaci(img,img_d,&npcr,&uaci);
        printf("%8.4f,%-8.4f,", npcr,uaci);
    }
    printf("\n");
    //第五对图片
    printf("%-15s,",file5);
    read_file(file5,img);
    read_file(file5_d,img_d);
    for (i=1;i<5;i++)
    {
        encryption(img,k1,k2,k3,k4);
        encryption(img_d,k1,k2,k3,k4);
        npcruaci(img,img_d,&npcr,&uaci);
        printf("%8.4f,%-8.4f,", npcr,uaci);
    }
    printf("\n");
    //第6对图片
    printf("%-15s,",file6);
    read_file(file6,img);
    read_file(file6_d,img_d);
    for (i=1;i<5;i++)
    {
        encryption(img,k1,k2,k3,k4);
        encryption(img_d,k1,k2,k3,k4);
        npcruaci(img,img_d,&npcr,&uaci);
        printf("%8.4f,%-8.4f,", npcr,uaci);
    }
    printf("\n");
    //第7对图片
    printf("%-15s,",file7);
    read_file(file7,img);
    read_file(file7_d,img_d);
    for (i=1;i<5;i++)
    {
        encryption(img,k1,k2,k3,k4);
        encryption(img_d,k1,k2,k3,k4);
        npcruaci(img,img_d,&npcr,&uaci);
        printf("%8.4f,%-8.4f,", npcr,uaci);
    }
    printf("\n");
    //第8对图片
    printf("%-15s,",file8);
    read_file(file8,img);
    read_file(file8_d,img_d);
    for (i=1;i<5;i++)
    {
        encryption(img,k1,k2,k3,k4);
        encryption(img_d,k1,k2,k3,k4);
        npcruaci(img,img_d,&npcr,&uaci);
        printf("%8.4f,%-8.4f,", npcr,uaci);
    }
    printf("\n===========================NPCR and UACI TEST=============================\n");
}


void load_npcr2(double k1, double k2, double k3, double k4)
{
    //计算并打印npcr和UACi情况
    int i;
    unsigned char img[WIDTH][LENGTH];
    unsigned char img_d[WIDTH][LENGTH];
    double npcr, uaci;
    //待测试的文件对
    char *file1="CT_Head.raw",   *file1_d="CT_Head_d1.raw";
    char *file2="CT_Head.raw",       *file2_d="CT_Head_d2.raw";
    char *file3="CT_Head.raw",    *file3_d="CT_Head_d3.raw";
    char *file4="CT_Head.raw",      *file4_d="CT_Head_d4.raw";
    char *file5="CT_Head.raw",      *file5_d="CT_Head_d5.raw";
    char *file6="CT_Head.raw",      *file6_d="CT_Head_d6.raw";
    char *file7="CT_Head.raw",      *file7_d="CT_Head_d7.raw";
    char *file8="CT_Head.raw",      *file8_d="CT_Head_d8.raw";
    //计算4轮内的NPCR和UACI
    printf("===========================NPCR and UACI TEST22222=============================\n");
    printf("%-15s,  %-17s,%-17s,%-17s,%-17s\n","Test img","Round 1","Round 2","Round 3","Round 4");
    //第一对图片
    printf("%-15s,",file1);
    //读取明文图像
    read_file(file1,img);
    read_file(file1_d,img_d);
    for (i=1;i<5;i++)
    {
        encryption(img,k1,k2,k3,k4);
        encryption(img_d,k1,k2,k3,k4);
        npcruaci(img,img_d,&npcr,&uaci);
        printf("%8.4f,%-8.4f,", npcr,uaci);
    }
    printf("\n");
    //第二对图片
    printf("%-15s,",file2);
    read_file(file2,img);
    read_file(file2_d,img_d);
    for (i=1;i<5;i++)
    {
        encryption(img,k1,k2,k3,k4);
        encryption(img_d,k1,k2,k3,k4);
        npcruaci(img,img_d,&npcr,&uaci);
        printf("%8.4f,%-8.4f,", npcr,uaci);
    }
    printf("\n");
    //第三对图片
    printf("%-15s,",file3);
    read_file(file3,img);
    read_file(file3_d,img_d);
    for (i=1;i<5;i++)
    {
        encryption(img,k1,k2,k3,k4);
        encryption(img_d,k1,k2,k3,k4);
        npcruaci(img,img_d,&npcr,&uaci);
        printf("%8.4f,%-8.4f,", npcr,uaci);
    }
    printf("\n");
    //第四对图片
    printf("%-15s,",file4);
    read_file(file4,img);
    read_file(file4_d,img_d);
    for (i=1;i<5;i++)
    {
        encryption(img,k1,k2,k3,k4);
        encryption(img_d,k1,k2,k3,k4);
        npcruaci(img,img_d,&npcr,&uaci);
        printf("%8.4f,%-8.4f,", npcr,uaci);
    }
    printf("\n");
    //第五对图片
    printf("%-15s,",file5);
    read_file(file5,img);
    read_file(file5_d,img_d);
    for (i=1;i<5;i++)
    {
        encryption(img,k1,k2,k3,k4);
        encryption(img_d,k1,k2,k3,k4);
        npcruaci(img,img_d,&npcr,&uaci);
        printf("%8.4f,%-8.4f,", npcr,uaci);
    }
    printf("\n");
    //第6对图片
    printf("%-15s,",file6);
    read_file(file6,img);
    read_file(file6_d,img_d);
    for (i=1;i<5;i++)
    {
        encryption(img,k1,k2,k3,k4);
        encryption(img_d,k1,k2,k3,k4);
        npcruaci(img,img_d,&npcr,&uaci);
        printf("%8.4f,%-8.4f,", npcr,uaci);
    }
    printf("\n");
    //第7对图片
    printf("%-15s,",file7);
    read_file(file7,img);
    read_file(file7_d,img_d);
    for (i=1;i<5;i++)
    {
        encryption(img,k1,k2,k3,k4);
        encryption(img_d,k1,k2,k3,k4);
        npcruaci(img,img_d,&npcr,&uaci);
        printf("%8.4f,%-8.4f,", npcr,uaci);
    }
    printf("\n");
    //第8对图片
    printf("%-15s,",file8);
    read_file(file8,img);
    read_file(file8_d,img_d);
    for (i=1;i<5;i++)
    {
        encryption(img,k1,k2,k3,k4);
        encryption(img_d,k1,k2,k3,k4);
        npcruaci(img,img_d,&npcr,&uaci);
        printf("%8.4f,%-8.4f,", npcr,uaci);
    }
    printf("\n===========================NPCR and UACI TEST  2222=============================\n");
}





void load_entropy_chiq(double k1, double k2, double k3, double k4)
{
    //计算并打印信息熵和卡方测试结果
    int i,j,totalpixel;
    FILE *img_fp;
    unsigned char img[WIDTH][LENGTH],ciphertext1d[WIDTH*LENGTH];
    double shang,kafang;
    //待测试的文件对
    char *file1="CT_Abdomen.raw";
    char *file2="CT_Head.raw";
    char *file3="CT_Paranasal_sinus.raw";
    char *file4="MR_Cervical_vertebra.raw";
    char *file5="MR_Knee.raw";
    char *file6="MR_Prostate.raw";
    char *file7="MR_Waist.raw";
    char *file8="X_Lungs.raw";

    //打印头
    totalpixel=WIDTH*LENGTH;
    printf("\n===========Entropy and Chi-Squre TEST=============\n");
    printf("%-25s,%-17s,%-17s\n","Test img","Entropy","Chi-Square");
    //printf("%-15s\n",cipher_count);

    //第1副图像
    //读取明文图像
    read_file(file1,img);
    for (i=0;i<cipher_count;i++)
    encryption(img,k1,k2,k3,k4);
    for (i=0;i<WIDTH;i++)
    for (j=0;j<LENGTH;j++)
    ciphertext1d[i*LENGTH+j]=img[i][j];
    entropy(ciphertext1d,&shang,&kafang);
    printf("%-25s,%-17.9f,%-17.0f\n",file1,shang,kafang);
    //保存以便计算局部熵

    write_file(pic_en1,img);

    //第2幅图片
    read_file(file2,img);
    for (i=0;i<cipher_count;i++)
    encryption(img,k1,k2,k3,k4);
    for (i=0;i<WIDTH;i++)
    for (j=0;j<LENGTH;j++)
    ciphertext1d[i*LENGTH+j]=img[i][j];
    entropy(ciphertext1d,&shang,&kafang);
    printf("%-25s,%-17.9f,%-17.0f\n",file2,shang,kafang);
    //保存以便计算局部熵
    write_file(pic_en2,img);

    //第3幅图片
    read_file(file3,img);
    for (i=0;i<cipher_count;i++)
    encryption(img,k1,k2,k3,k4);
    for (i=0;i<WIDTH;i++)
    for (j=0;j<LENGTH;j++)
    ciphertext1d[i*LENGTH+j]=img[i][j];
    entropy(ciphertext1d,&shang,&kafang);
    printf("%-25s,%-17.9f,%-17.0f\n",file3,shang,kafang);
    //保存以便计算局部熵
    write_file(pic_en3,img);
    //第4幅图片
    read_file(file4,img);
    for (i=0;i<cipher_count;i++)
    encryption(img,k1,k2,k3,k4);
    for (i=0;i<WIDTH;i++)
    for (j=0;j<LENGTH;j++)
    ciphertext1d[i*LENGTH+j]=img[i][j];
    entropy(ciphertext1d,&shang,&kafang);
    printf("%-25s,%-17.9f,%-17.0f\n",file4,shang,kafang);
    //保存以便计算局部熵
    write_file(pic_en4,img);
    //第5幅图片
    read_file(file5,img);
    for (i=0;i<cipher_count;i++)
    encryption(img,k1,k2,k3,k4);
    for (i=0;i<WIDTH;i++)
    for (j=0;j<LENGTH;j++)
    ciphertext1d[i*LENGTH+j]=img[i][j];
    entropy(ciphertext1d,&shang,&kafang);
    printf("%-25s,%-17.9f,%-17.0f\n",file5,shang,kafang);
    //保存以便计算局部熵
    write_file(pic_en5,img);
    //第6幅图片
    read_file(file6,img);
    for (i=0;i<cipher_count;i++)
    encryption(img,k1,k2,k3,k4);
    for (i=0;i<WIDTH;i++)
    for (j=0;j<LENGTH;j++)
    ciphertext1d[i*LENGTH+j]=img[i][j];
    entropy(ciphertext1d,&shang,&kafang);
    printf("%-25s,%-17.9f,%-17.0f\n",file6,shang,kafang);
    //保存以便计算局部熵
    write_file(pic_en6,img);
    //第7幅图片
    read_file(file7,img);
    for (i=0;i<cipher_count;i++)
    encryption(img,k1,k2,k3,k4);
    for (i=0;i<WIDTH;i++)
    for (j=0;j<LENGTH;j++)
    ciphertext1d[i*LENGTH+j]=img[i][j];
    entropy(ciphertext1d,&shang,&kafang);
    printf("%-25s,%-17.9f,%-17.0f\n",file7,shang,kafang);
    //保存以便计算局部熵
    write_file(pic_en7,img);
    //第8幅图片
    read_file(file8,img);
    for (i=0;i<cipher_count;i++)
    encryption(img,k1,k2,k3,k4);
    for (i=0;i<WIDTH;i++)
    for (j=0;j<LENGTH;j++)
    ciphertext1d[i*LENGTH+j]=img[i][j];
    entropy(ciphertext1d,&shang,&kafang);
    printf("%-25s,%-17.9f,%-17.0f\n",file8,shang,kafang);
    //保存以便计算局部熵
    write_file(pic_en8,img);
    printf("\n===========Entropy and Chi-Squre TEST=============\n\n\n\n");
}


void load_keysensi(double k1, double dk1, double k2, double dk2, double k3, double dk3, double k4, double dk4)
{
    //计算并打印像素相关性
    int i;
    FILE *img_fp;
    unsigned char img[WIDTH][LENGTH],img_e1[WIDTH][LENGTH],img_e2[WIDTH][LENGTH],img_de[WIDTH][LENGTH];
    double npcr, uaci;
    printf("===========Key Sensitivity TEST=============\n");
    //读取待测试的文件
    //基准明文
    char *file="CT_Head.raw";
    read_file(file,img);
    //产生基准密文
    read_file(file,img_e1);
    for (i=0;i<cipher_count;i++)
    encryption(img_e1,k1,k2,k3,k4);
    //开始加密密钥敏感性
    printf("KeySense of Encryption: ");
    //务必重读图像
    read_file(file,img_e2);
    for (i=0;i<cipher_count;i++)
    encryption(img_e2,k1+dk1,k2,k3,k4);
    npcruaci(img_e1,img_e2,&npcr,&uaci);
    printf("%.4f, %.4f;  ",npcr,uaci);
    //保存
    if ((img_fp = fopen("en_dk1.raw", "wb")) == NULL)
     {
        printf("cannot open file to write\n");
        return;
    }
    fwrite(img_e2, sizeof(unsigned char), WIDTH*LENGTH, img_fp);
    fclose(img_fp);
    //根据密钥个数灵活设置
    //务必重读图像
    read_file(file,img_e2);
    for (i=0;i<cipher_count;i++)
    encryption(img_e2,k1,k2+dk2,k3,k4);
    npcruaci(img_e1,img_e2,&npcr,&uaci);
    printf("%.4f, %.4f;  ",npcr,uaci);
    //保存
    if ((img_fp = fopen("en_dk2.raw", "wb")) == NULL)
     {
        printf("cannot open file to write\n");
        return;
    }
    fwrite(img_e2, sizeof(unsigned char), WIDTH*LENGTH, img_fp);
    fclose(img_fp);
    //务必重读图像
    read_file(file,img_e2);
    for (i=0;i<cipher_count;i++)
    encryption(img_e2,k1,k2,k3+dk3,k4);
    npcruaci(img_e1,img_e2,&npcr,&uaci);
    printf("%.4f, %.4f;  ",npcr,uaci);
    //保存
    if ((img_fp = fopen("en_dk3.raw", "wb")) == NULL)
     {
        printf("cannot open file to write\n");
        return;
    }
    fwrite(img_e2, sizeof(unsigned char), WIDTH*LENGTH, img_fp);
    fclose(img_fp);
    //务必重读图像
    read_file(file,img_e2);
    for (i=0;i<cipher_count;i++)
    encryption(img_e2,k1,k2,k3,k4+dk4);
    npcruaci(img_e1,img_e2,&npcr,&uaci);
    printf("%.4f, %.4f;  ",npcr,uaci);
    //保存
    if ((img_fp = fopen("en_dk4.raw", "wb")) == NULL)
     {
        printf("cannot open file to write\n");
        return;
    }
    fwrite(img_e2, sizeof(unsigned char), WIDTH*LENGTH, img_fp);
    fclose(img_fp);


    //以下为解密
    printf("\nKeySense of Decryption: ");
    //先加密
    read_file(file,img_de);
    for (i=0;i<cipher_count;i++)
    encryption(img_de,k1,k2,k3,k4);
    //再解密
    for (i=0;i<cipher_count;i++)
    decryption(img_de,k1+dk1,k2,k3,k4);
    npcruaci(img,img_de,&npcr,&uaci);
    printf("%.4f, %.4f;  ",npcr,uaci);
        //保存
    if ((img_fp = fopen("de_dk1.raw", "wb")) == NULL)
     {
        printf("cannot open file to write\n");
        return;
    }
    fwrite(img_de, sizeof(unsigned char), WIDTH*LENGTH, img_fp);
    fclose(img_fp);
    //以下根据密钥个数灵活设置
    //先加密
    read_file(file,img_de);
    for (i=0;i<cipher_count;i++)
    encryption(img_de,k1,k2,k3,k4);
    //再解密
    for (i=0;i<cipher_count;i++)
    decryption(img_de,k1,k2+dk2,k3,k4);
    npcruaci(img,img_de,&npcr,&uaci);
    printf("%.4f, %.4f;  ",npcr,uaci);
        //保存
    if ((img_fp = fopen("de_dk2.raw", "wb")) == NULL)
     {
        printf("cannot open file to write\n");
        return;
    }
    fwrite(img_de, sizeof(unsigned char), WIDTH*LENGTH, img_fp);
    fclose(img_fp);
    //先加密
    read_file(file,img_de);
    for (i=0;i<cipher_count;i++)
    encryption(img_de,k1,k2,k3,k4);
    //再解密
    for (i=0;i<cipher_count;i++)
    decryption(img_de,k1,k2,k3+dk3,k4);
    npcruaci(img,img_de,&npcr,&uaci);
    printf("%.4f, %.4f;  ",npcr,uaci);
        //保存
    if ((img_fp = fopen("de_dk3.raw", "wb")) == NULL)
     {
        printf("cannot open file to write\n");
        return;
    }
    fwrite(img_de, sizeof(unsigned char), WIDTH*LENGTH, img_fp);
    fclose(img_fp);
    //先加密
    read_file(file,img_de);
    for (i=0;i<cipher_count;i++)
    encryption(img_de,k1,k2,k3,k4);
    //再解密
    for (i=0;i<cipher_count;i++)
    decryption(img_de,k1,k2,k3,k4+dk4);
    npcruaci(img,img_de,&npcr,&uaci);
    printf("%.4f, %.4f;  ",npcr,uaci);
        //保存
    if ((img_fp = fopen("de_dk4.raw", "wb")) == NULL)
     {
        printf("cannot open file to write\n");
        return;
    }
    fwrite(img_de, sizeof(unsigned char), WIDTH*LENGTH, img_fp);
    fclose(img_fp);
    printf("\n===========Key Sensitivity TEST=============\n\n");
}


int main(void)
{

    long starter, ender;
    // 测试整个加密过程花费的时间
    double x0=5.79529279475736,y0=8.17628155473482,z0=6.61243069535167,w0=7.27405927962098;
    double det=1.0e-15;

    starter = clock();

    load_entropy_chiq(x0,y0,z0,w0);
    load_keysensi(x0,det,y0,det,z0,det,w0,det);
    load_npcr(x0,y0,z0,w0);
    load_npcr2(x0,y0,z0,w0);

    ender = clock();
    printf("\n\n加密消耗%f秒\n", (double)(ender - starter)/CLOCKS_PER_SEC);
    system("pause");

}





/*
int main(void) {

    long starter, ender;
    // 测试整个加密过程花费的时间
    unsigned char img[WIDTH][LENGTH];
    unsigned char img_d[WIDTH][LENGTH];
    double x0=3.089292794757325,y0=9.876291554723856,z0=8.672430195056872,w0=2.274659284620194;

    starter = clock();
    read_file(pic_o,img);
    encryption(img,x0,y0,z0,w0);
    write_file(pic_en,img);
    ender = clock();
    printf("加密消耗%f秒\n", (double)(ender - starter)/CLOCKS_PER_SEC);

    starter = clock();
    read_file(pic_en,img_d);
    decryption(img_d,x0,y0,z0,w0);
    write_file(pic_de,img_d);
    ender = clock();
    printf("解密消耗%f秒\n", (double)(ender - starter)/CLOCKS_PER_SEC);
    system("pause");

}
*/
