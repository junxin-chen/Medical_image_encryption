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

#define pic_en1  "en_CT_Abdomen.raw"
#define pic_en2  "en_CT_Head.raw"
#define pic_en3  "en_CT_Paranasal_sinus.raw"
#define pic_en4  "en_MR_Cervical_vertebra.raw"
#define pic_en5  "en_MR_Knee.raw"
#define pic_en6  "en_MR_Prostate.raw"
#define pic_en7  "en_MR_Waist.raw"
#define pic_en8  "en_X_Lungs.raw"



double var1[KEY_STREAM_LENGTH] = {0}; //���������1
double var2[KEY_STREAM_LENGTH] = {0}; //���������2
double var3[KEY_STREAM_LENGTH] = {0}; //���������3
double var4[KEY_STREAM_LENGTH] = {0}; //���������4


unsigned int pv1[KEY_STREAM_LENGTH] = {0};//��������1
unsigned int pv2[KEY_STREAM_LENGTH] = {0};//��������2
unsigned int pv1_tmp[KEY_STREAM_LENGTH] = {0};//��������1
unsigned int pv2_tmp[KEY_STREAM_LENGTH] = {0};//��������2

unsigned char mask1[KEY_STREAM_LENGTH] = {0};//����1
unsigned char mask2[KEY_STREAM_LENGTH] = {0};//����2


void generate_variables(double x0,double y0,double z0,double w0)
{
	double h=0.0005,a = 12,b = 23,c = 1,d = 2.1,k = 6,r = 0.2;   //����ֵ�����Ʋ���
	double K1,K2,K3,K4,L1,L2,L3,L4,M1,M2,M3,M4,N1,N2,N3,N4;
	int i=0;

	var1[0]=x0;var2[0]=y0;var3[0]=z0;var4[0]=w0;

	for (i=0;i<300;i++)//����200��
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

	for (i=0;i<KEY_STREAM_LENGTH;i++)// ���ɻ�������
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

    //��������pv1��pv2���Լ�pixel swapping�Ľ�������
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
    //ʹ�ý���������������pv
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

    //����������ɢ��mask
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
	unsigned char temp=0;//����ȡ������ֵ
	//�����λλ��Ϊ0��ֱ�ӷ���
	if (loop==0)
	return p;
	//����ִ����������
	if(direction == 0)  //���ѭ������
	{
		for(i=0;i<loop;i++)//ÿ������λΪ1��ѭ������
		{
			temp=p&128;//ȡ��ѭ������ʱҪ��ʧ�ĵ�һλ��ȡ��temp�ĵ�һλ
			p=p<<1;//P����1λ�������λ��ʧ������λ���ƣ����λ����
			temp=temp>>7;//temp����7λ����ȡ�������λ�Ƶ����λ
			p=p|temp;//��temp�����λ����p�����λ������temp������λ��Ϊ0���ʰ�λ���ʱ��Ӱ������λ
		}
		return p;
	}//ѭ���������
	if(direction == 1)//���ѭ������
	{
		for(i=0;i<loop;i++)
		{
			temp=p&1;//ȡ��ѭ������ʱҪ��ʧ�����λ��ȡ��temp�����λ
			p=p>>1;//P����1λ�������λ��ʧ������λ���ƣ����λ����
			temp=temp<<7;//temp����7λ����ȡ����p�����λ�������λ
			p=p|temp;//��temp�����λ����p�����λ������temp������λΪ0����λ���ʱ��Ӱ������λ
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
        //�Ƚ������Ҳ���
        //�Ȼ����λ�Ĳ���
        direction=pv2i&1;//��λ�ķ�����pv2�����һλ����
        loop=pv1i&7;//��λ��λ����pv1�ĺ���λ����
        current_plain=*(img_pointer+pv1i);
        current_plain=cirshift(current_plain,direction,loop);
        current_plain=current_plain^previous_plain;
        *(img_temp_pointer+pv2i)=current_plain;
        previous_plain=*(img_pointer+pv1i);
        //���Ҽ������
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

        //��ʼ���ҽ���
        direction=(pv2i^1)&1;//��λ�ķ�����pv2�����һλ����
        loop=pv1i&7;//��λ��λ����pv1�ĺ���λ����
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
    //�����ڶ�ȡ�ļ�
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
    //�����ڶ�ȡ�ļ�
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

    //������Ϣ��
	for(i=0;i<256;i++)
	{
		Cal2[i]=Cal[i];
		Cal[i] /= totalpixel;
		if(Cal[i]==0)
			entr = entr;
		else
			entr = entr - Cal[i]*(log(Cal[i])/log(2));
	}
	//���㿨���ֲ�
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
    //���㲢��ӡnpcr��UACi���
    int i;
    unsigned char img[WIDTH][LENGTH];
    unsigned char img_d[WIDTH][LENGTH];
    double npcr, uaci;
    //�����Ե��ļ���
    char *file1="CT_Abdomen.raw",   *file1_d="CT_Abdomen_d.raw";
    char *file2="CT_Head.raw",       *file2_d="CT_Head_d.raw";
    char *file3="CT_Paranasal_sinus.raw",    *file3_d="CT_Paranasal_sinus_d.raw";
    char *file4="MR_Cervical_vertebra.raw",      *file4_d="MR_Cervical_vertebra_d.raw";
    char *file5="MR_Knee.raw",      *file5_d="MR_Knee_d.raw";
    char *file6="MR_Prostate.raw",      *file6_d="MR_Prostate_d.raw";
    char *file7="MR_Waist.raw",      *file7_d="MR_Waist_d.raw";
    char *file8="X_Lungs.raw",      *file8_d="X_Lungs_d.raw";
    //����4���ڵ�NPCR��UACI
    printf("===========================NPCR and UACI TEST=============================\n");
    printf("%-15s,  %-17s,%-17s,%-17s,%-17s\n","Test img","Round 1","Round 2","Round 3","Round 4");
    //��һ��ͼƬ
    printf("%-15s,",file1);
    //��ȡ����ͼ��
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
    //�ڶ���ͼƬ
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
    //������ͼƬ
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
    //���Ķ�ͼƬ
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
    //�����ͼƬ
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
    //��6��ͼƬ
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
    //��7��ͼƬ
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
    //��8��ͼƬ
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
    //���㲢��ӡnpcr��UACi���
    int i;
    unsigned char img[WIDTH][LENGTH];
    unsigned char img_d[WIDTH][LENGTH];
    double npcr, uaci;
    //�����Ե��ļ���
    char *file1="CT_Head.raw",   *file1_d="CT_Head_d1.raw";
    char *file2="CT_Head.raw",       *file2_d="CT_Head_d2.raw";
    char *file3="CT_Head.raw",    *file3_d="CT_Head_d3.raw";
    char *file4="CT_Head.raw",      *file4_d="CT_Head_d4.raw";
    char *file5="CT_Head.raw",      *file5_d="CT_Head_d5.raw";
    char *file6="CT_Head.raw",      *file6_d="CT_Head_d6.raw";
    char *file7="CT_Head.raw",      *file7_d="CT_Head_d7.raw";
    char *file8="CT_Head.raw",      *file8_d="CT_Head_d8.raw";
    //����4���ڵ�NPCR��UACI
    printf("===========================NPCR and UACI TEST22222=============================\n");
    printf("%-15s,  %-17s,%-17s,%-17s,%-17s\n","Test img","Round 1","Round 2","Round 3","Round 4");
    //��һ��ͼƬ
    printf("%-15s,",file1);
    //��ȡ����ͼ��
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
    //�ڶ���ͼƬ
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
    //������ͼƬ
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
    //���Ķ�ͼƬ
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
    //�����ͼƬ
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
    //��6��ͼƬ
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
    //��7��ͼƬ
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
    //��8��ͼƬ
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
    //���㲢��ӡ��Ϣ�غͿ������Խ��
    int i,j,totalpixel;
    FILE *img_fp;
    unsigned char img[WIDTH][LENGTH],ciphertext1d[WIDTH*LENGTH];
    double shang,kafang;
    //�����Ե��ļ���
    char *file1="CT_Abdomen.raw";
    char *file2="CT_Head.raw";
    char *file3="CT_Paranasal_sinus.raw";
    char *file4="MR_Cervical_vertebra.raw";
    char *file5="MR_Knee.raw";
    char *file6="MR_Prostate.raw";
    char *file7="MR_Waist.raw";
    char *file8="X_Lungs.raw";

    //��ӡͷ
    totalpixel=WIDTH*LENGTH;
    printf("\n===========Entropy and Chi-Squre TEST=============\n");
    printf("%-25s,%-17s,%-17s\n","Test img","Entropy","Chi-Square");
    //printf("%-15s\n",cipher_count);

    //��1��ͼ��
    //��ȡ����ͼ��
    read_file(file1,img);
    for (i=0;i<cipher_count;i++)
    encryption(img,k1,k2,k3,k4);
    for (i=0;i<WIDTH;i++)
    for (j=0;j<LENGTH;j++)
    ciphertext1d[i*LENGTH+j]=img[i][j];
    entropy(ciphertext1d,&shang,&kafang);
    printf("%-25s,%-17.9f,%-17.0f\n",file1,shang,kafang);
    //�����Ա����ֲ���

    write_file(pic_en1,img);

    //��2��ͼƬ
    read_file(file2,img);
    for (i=0;i<cipher_count;i++)
    encryption(img,k1,k2,k3,k4);
    for (i=0;i<WIDTH;i++)
    for (j=0;j<LENGTH;j++)
    ciphertext1d[i*LENGTH+j]=img[i][j];
    entropy(ciphertext1d,&shang,&kafang);
    printf("%-25s,%-17.9f,%-17.0f\n",file2,shang,kafang);
    //�����Ա����ֲ���
    write_file(pic_en2,img);

    //��3��ͼƬ
    read_file(file3,img);
    for (i=0;i<cipher_count;i++)
    encryption(img,k1,k2,k3,k4);
    for (i=0;i<WIDTH;i++)
    for (j=0;j<LENGTH;j++)
    ciphertext1d[i*LENGTH+j]=img[i][j];
    entropy(ciphertext1d,&shang,&kafang);
    printf("%-25s,%-17.9f,%-17.0f\n",file3,shang,kafang);
    //�����Ա����ֲ���
    write_file(pic_en3,img);
    //��4��ͼƬ
    read_file(file4,img);
    for (i=0;i<cipher_count;i++)
    encryption(img,k1,k2,k3,k4);
    for (i=0;i<WIDTH;i++)
    for (j=0;j<LENGTH;j++)
    ciphertext1d[i*LENGTH+j]=img[i][j];
    entropy(ciphertext1d,&shang,&kafang);
    printf("%-25s,%-17.9f,%-17.0f\n",file4,shang,kafang);
    //�����Ա����ֲ���
    write_file(pic_en4,img);
    //��5��ͼƬ
    read_file(file5,img);
    for (i=0;i<cipher_count;i++)
    encryption(img,k1,k2,k3,k4);
    for (i=0;i<WIDTH;i++)
    for (j=0;j<LENGTH;j++)
    ciphertext1d[i*LENGTH+j]=img[i][j];
    entropy(ciphertext1d,&shang,&kafang);
    printf("%-25s,%-17.9f,%-17.0f\n",file5,shang,kafang);
    //�����Ա����ֲ���
    write_file(pic_en5,img);
    //��6��ͼƬ
    read_file(file6,img);
    for (i=0;i<cipher_count;i++)
    encryption(img,k1,k2,k3,k4);
    for (i=0;i<WIDTH;i++)
    for (j=0;j<LENGTH;j++)
    ciphertext1d[i*LENGTH+j]=img[i][j];
    entropy(ciphertext1d,&shang,&kafang);
    printf("%-25s,%-17.9f,%-17.0f\n",file6,shang,kafang);
    //�����Ա����ֲ���
    write_file(pic_en6,img);
    //��7��ͼƬ
    read_file(file7,img);
    for (i=0;i<cipher_count;i++)
    encryption(img,k1,k2,k3,k4);
    for (i=0;i<WIDTH;i++)
    for (j=0;j<LENGTH;j++)
    ciphertext1d[i*LENGTH+j]=img[i][j];
    entropy(ciphertext1d,&shang,&kafang);
    printf("%-25s,%-17.9f,%-17.0f\n",file7,shang,kafang);
    //�����Ա����ֲ���
    write_file(pic_en7,img);
    //��8��ͼƬ
    read_file(file8,img);
    for (i=0;i<cipher_count;i++)
    encryption(img,k1,k2,k3,k4);
    for (i=0;i<WIDTH;i++)
    for (j=0;j<LENGTH;j++)
    ciphertext1d[i*LENGTH+j]=img[i][j];
    entropy(ciphertext1d,&shang,&kafang);
    printf("%-25s,%-17.9f,%-17.0f\n",file8,shang,kafang);
    //�����Ա����ֲ���
    write_file(pic_en8,img);
    printf("\n===========Entropy and Chi-Squre TEST=============\n\n\n\n");
}


void load_keysensi(double k1, double dk1, double k2, double dk2, double k3, double dk3, double k4, double dk4)
{
    //���㲢��ӡ���������
    int i;
    FILE *img_fp;
    unsigned char img[WIDTH][LENGTH],img_e1[WIDTH][LENGTH],img_e2[WIDTH][LENGTH],img_de[WIDTH][LENGTH];
    double npcr, uaci;
    printf("===========Key Sensitivity TEST=============\n");
    //��ȡ�����Ե��ļ�
    //��׼����
    char *file="CT_Head.raw";
    read_file(file,img);
    //������׼����
    read_file(file,img_e1);
    for (i=0;i<cipher_count;i++)
    encryption(img_e1,k1,k2,k3,k4);
    //��ʼ������Կ������
    printf("KeySense of Encryption: ");
    //����ض�ͼ��
    read_file(file,img_e2);
    for (i=0;i<cipher_count;i++)
    encryption(img_e2,k1+dk1,k2,k3,k4);
    npcruaci(img_e1,img_e2,&npcr,&uaci);
    printf("%.4f, %.4f;  ",npcr,uaci);
    //����
    if ((img_fp = fopen("en_dk1.raw", "wb")) == NULL)
     {
        printf("cannot open file to write\n");
        return;
    }
    fwrite(img_e2, sizeof(unsigned char), WIDTH*LENGTH, img_fp);
    fclose(img_fp);
    //������Կ�����������
    //����ض�ͼ��
    read_file(file,img_e2);
    for (i=0;i<cipher_count;i++)
    encryption(img_e2,k1,k2+dk2,k3,k4);
    npcruaci(img_e1,img_e2,&npcr,&uaci);
    printf("%.4f, %.4f;  ",npcr,uaci);
    //����
    if ((img_fp = fopen("en_dk2.raw", "wb")) == NULL)
     {
        printf("cannot open file to write\n");
        return;
    }
    fwrite(img_e2, sizeof(unsigned char), WIDTH*LENGTH, img_fp);
    fclose(img_fp);
    //����ض�ͼ��
    read_file(file,img_e2);
    for (i=0;i<cipher_count;i++)
    encryption(img_e2,k1,k2,k3+dk3,k4);
    npcruaci(img_e1,img_e2,&npcr,&uaci);
    printf("%.4f, %.4f;  ",npcr,uaci);
    //����
    if ((img_fp = fopen("en_dk3.raw", "wb")) == NULL)
     {
        printf("cannot open file to write\n");
        return;
    }
    fwrite(img_e2, sizeof(unsigned char), WIDTH*LENGTH, img_fp);
    fclose(img_fp);
    //����ض�ͼ��
    read_file(file,img_e2);
    for (i=0;i<cipher_count;i++)
    encryption(img_e2,k1,k2,k3,k4+dk4);
    npcruaci(img_e1,img_e2,&npcr,&uaci);
    printf("%.4f, %.4f;  ",npcr,uaci);
    //����
    if ((img_fp = fopen("en_dk4.raw", "wb")) == NULL)
     {
        printf("cannot open file to write\n");
        return;
    }
    fwrite(img_e2, sizeof(unsigned char), WIDTH*LENGTH, img_fp);
    fclose(img_fp);


    //����Ϊ����
    printf("\nKeySense of Decryption: ");
    //�ȼ���
    read_file(file,img_de);
    for (i=0;i<cipher_count;i++)
    encryption(img_de,k1,k2,k3,k4);
    //�ٽ���
    for (i=0;i<cipher_count;i++)
    decryption(img_de,k1+dk1,k2,k3,k4);
    npcruaci(img,img_de,&npcr,&uaci);
    printf("%.4f, %.4f;  ",npcr,uaci);
        //����
    if ((img_fp = fopen("de_dk1.raw", "wb")) == NULL)
     {
        printf("cannot open file to write\n");
        return;
    }
    fwrite(img_de, sizeof(unsigned char), WIDTH*LENGTH, img_fp);
    fclose(img_fp);
    //���¸�����Կ�����������
    //�ȼ���
    read_file(file,img_de);
    for (i=0;i<cipher_count;i++)
    encryption(img_de,k1,k2,k3,k4);
    //�ٽ���
    for (i=0;i<cipher_count;i++)
    decryption(img_de,k1,k2+dk2,k3,k4);
    npcruaci(img,img_de,&npcr,&uaci);
    printf("%.4f, %.4f;  ",npcr,uaci);
        //����
    if ((img_fp = fopen("de_dk2.raw", "wb")) == NULL)
     {
        printf("cannot open file to write\n");
        return;
    }
    fwrite(img_de, sizeof(unsigned char), WIDTH*LENGTH, img_fp);
    fclose(img_fp);
    //�ȼ���
    read_file(file,img_de);
    for (i=0;i<cipher_count;i++)
    encryption(img_de,k1,k2,k3,k4);
    //�ٽ���
    for (i=0;i<cipher_count;i++)
    decryption(img_de,k1,k2,k3+dk3,k4);
    npcruaci(img,img_de,&npcr,&uaci);
    printf("%.4f, %.4f;  ",npcr,uaci);
        //����
    if ((img_fp = fopen("de_dk3.raw", "wb")) == NULL)
     {
        printf("cannot open file to write\n");
        return;
    }
    fwrite(img_de, sizeof(unsigned char), WIDTH*LENGTH, img_fp);
    fclose(img_fp);
    //�ȼ���
    read_file(file,img_de);
    for (i=0;i<cipher_count;i++)
    encryption(img_de,k1,k2,k3,k4);
    //�ٽ���
    for (i=0;i<cipher_count;i++)
    decryption(img_de,k1,k2,k3,k4+dk4);
    npcruaci(img,img_de,&npcr,&uaci);
    printf("%.4f, %.4f;  ",npcr,uaci);
        //����
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
    // �����������ܹ��̻��ѵ�ʱ��
    double x0=5.79529279475736,y0=8.17628155473482,z0=6.61243069535167,w0=7.27405927962098;
    double det=1.0e-15;

    starter = clock();

    load_entropy_chiq(x0,y0,z0,w0);
    load_keysensi(x0,det,y0,det,z0,det,w0,det);
    load_npcr(x0,y0,z0,w0);
    load_npcr2(x0,y0,z0,w0);

    ender = clock();
    printf("\n\n��������%f��\n", (double)(ender - starter)/CLOCKS_PER_SEC);
    system("pause");

}





/*
int main(void) {

    long starter, ender;
    // �����������ܹ��̻��ѵ�ʱ��
    unsigned char img[WIDTH][LENGTH];
    unsigned char img_d[WIDTH][LENGTH];
    double x0=3.089292794757325,y0=9.876291554723856,z0=8.672430195056872,w0=2.274659284620194;

    starter = clock();
    read_file(pic_o,img);
    encryption(img,x0,y0,z0,w0);
    write_file(pic_en,img);
    ender = clock();
    printf("��������%f��\n", (double)(ender - starter)/CLOCKS_PER_SEC);

    starter = clock();
    read_file(pic_en,img_d);
    decryption(img_d,x0,y0,z0,w0);
    write_file(pic_de,img_d);
    ender = clock();
    printf("��������%f��\n", (double)(ender - starter)/CLOCKS_PER_SEC);
    system("pause");

}
*/
