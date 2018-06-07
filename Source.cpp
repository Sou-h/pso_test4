#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<iostream>

/* �L���萔 */
#define NOPS	40				//���q�̌�
#define LIMITL 256			//�z��ő�̈�
#define ILIMIT	100			//�J��Ԃ��̉�
#define SEED	32767		//�����̏����l
#define W			0.3			//�����萔
#define C1		1.0			//���[�J���Ȏ���
#define C2		1.0			//�O���[�o���Ȏ���
#define DIM		2				//������
#define PI			3.14159265359

// ���ʍ��W��\������\���� 
struct point {
	double x[DIM];
};

// ���q��\������\���� 
struct particle {
	struct point pos;		/*�ʒu*/
	double value;				/*�]���l*/
	struct point v;			/*���x*/
	struct point bestpos;	/*�œK�ʒu*/
	double bestval;			/*�œK�]���l*/
};
/* �֐��̃v���g�^�C�v�̐錾 */
void initps(struct particle ps[]);	/*���q�Q�̏�����*/
void printps(struct particle ps[]);	/*���q�Q�̕\��*/
double frand();					/*��������*/
double calcval(double *x);	/*�]���l�̌v�Z*/
void optimize(struct particle ps[]);	/*���q�Q�̈ʒu�X�V*/
void setgbest(struct particle ps[]);	/* �Q��̍œK�ʒu���i�[ */

										/*���ϐ�*/
struct point gbestpos;	/*�Q���̍œK�ʒu*/
double gbestval;			/*�Q���̍œK�]���l*/

							// main() �֐�/
int main()
{

	struct particle ps[LIMITL];	/*���q�Q*/
	int i;						/*�J��Ԃ��񐔂̐���*/
	srand(SEED);			/*�����̏�����*/

							/*���q�Q�̏�����*/
	initps(ps);

	printps(ps);
	for (i = 0; i<ILIMIT; ++i) {
		optimize(ps);
		//�e�L�X�g�ɕۑ�
		printf("%d���\n", i);
		printps(ps);
	}

	getchar();
	return 0;
}

/* optimize() �֐� */
// ���q�Q�̈ʒu�X�V 
void optimize(struct particle ps[])
{
	int i;
	double r1, r2;/*�����̒l���i�[*/
	for (i = 0; i<NOPS; ++i) {
		/*�����̐ݒ�*/
		r1 = frand();
		r2 = frand();
		/*���x�̍X�V*/
		for (int j = 0; j < DIM; j++) {
			ps[i].v.x[j] = W*ps[i].v.x[j]
				+ C1*r1*(ps[i].bestpos.x[j] - ps[i].pos.x[j])
				+ C2*r2*(gbestpos.x[j] - ps[i].pos.x[j]);
			/*�ʒu�̍X�V*/
			ps[i].pos.x[j] += ps[i].v.x[j];

		}

		/*�œK�l�̍X�V*/
		ps[i].value = calcval(ps[i].pos.x);

		if (ps[i].value<ps[i].bestval) {
			ps[i].bestval = ps[i].value;
			ps[i].bestpos = ps[i].pos;
		}
	}
	/*�Q���œK�l�̍X�V*/
	setgbest(ps);
}

//*calcbal() �֐� 
// �]���l�̌v�Z 
double calcval(double *x)
{
	double sum = 1;
	/*	for (int i = 0; i < DIM ; i++){
	sum += x[i]*x[i];
	}
	*/
	/*	for (int i = 0; i < DIM; i++) {
	sum += pow(x[i],4.0)-16*pow(x[i],2.0)+5*x[i];
	}
	*/
	for (int i = 0; i < DIM; i++) {
		sum += pow(x[i], 2.0) - 10 * cos(2 * PI*x[i]) + 10;
	}
	//	printf("sum=%.6lf\n", sum);
	return sum;
}

double calcval2(double x, double y)
{
	double sum = 0;
	for (int i = 0; i < NOPS; i++) {
		sum = sum + x *x;
	}

	return sum;
}

/* setgbest() �֐� */
// �Q��̍œK�ʒu���i�[ 
void setgbest(struct particle ps[])
{
	int i;
	double besti;
	double x[NOPS];
	besti = ps[0].value;
	for (int for_count = 0; for_count < DIM; for_count++) {
		x[for_count] = ps[0].pos.x[for_count];
	}

	for (i = 0; i < NOPS; i++)
		/*���݂̍ŗǕ]���l��T��*/
		if (ps[i].value < besti) {
			besti = ps[i].value;
			for (int j = 0; j < DIM; j++) {
				x[j] = ps[i].pos.x[j];
			}

		}

	/*�]���l���ߋ������悩������X�V*/
	if (besti < gbestval) {
		gbestval = besti;
		for (int j = 0; j < DIM; j++) {
			gbestpos.x[j] = x[j];
		}
	}

}

/* initps() �֐� */
// ���q�Q�̏����� 
void initps(struct particle ps[])
{
	int i;
	double x[DIM];
	for (i = 0; i < NOPS; i++) {
		/*�ʒu*/
		for (int j = 0; j < DIM; j++) {
			x[j] = ps[i].pos.x[j] = frand() * 2 - 1.0;
		}
		/*�]���l*/
		ps[i].value = calcval(x);
		/*���x*/
		for (int j = 0; j < DIM; j++) {
			ps[i].v.x[j] = frand() * 2 - 1.0;
		}
		/*�œK�ʒu*/
		for (int j = 0; j < DIM; j++) {
			ps[i].bestpos.x[j] = ps[i].pos.x[j];
		}
		/*�œK�]���l*/
		ps[i].bestval = ps[i].value;
	}

	/*�Q��̍œK�ʒu���i�[*/
	gbestval = ps[0].value;
	for (int j = 0; j < DIM; j++) {
		gbestpos.x[j] = ps[0].pos.x[j];
	}

	setgbest(ps);

}


/* printps() �֐� */
// ���q�Q�̕\�� 
void printps(struct particle ps[])
{
	int i;
	for (i = 0; i<NOPS; ++i) {
		printf("���q%d\t", i);
		for (int j = 0; j < DIM; j++)		printf("x[%d]=%.4lf\t", j, ps[i].pos.x[j]);
		printf("�]���l%2.6lf\t", ps[i].value);
		//		printf("%lf %lf ", ps[i].v.x, ps[i].v.y);
		//		printf("%lf %lf ",ps[i].bestpos.x, ps[i].bestpos.y);
		printf("�œK�]���l%lf\n", ps[i].bestval);
	}
	for (int j = 0; j < DIM; j++)	printf("���� %.6lf\t ", gbestpos.x[j]);
	printf("�]��\t%2.6lf\n", gbestval);
}

/* frand() �֐� */
// �������� 
double frand(void)
{
	double result;
	while ((result = (double)rand() / RAND_MAX) >= 1);
	return result;
}