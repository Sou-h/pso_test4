#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<iostream>

/* 記号定数 */
#define NOPS	40				//粒子の個数
#define LIMITL 256			//配列最大領域
#define ILIMIT	10				//繰り返しの回数
#define SEED	32767		//乱数の初期値
#define W			0.3			//慣性定数
#define C1		1.0			//ローカルな質量
#define C2		1.0			//グローバルな質量
#define DIM		2				//次元数
#define EVALUATION_NUM 3
#define PI			3.14159265359

// 粒子を表現する構造体 
struct particle {
	double pos[DIM];		/*位置*/
	double value;				/*評価値*/
	double v [DIM];			/*速度*/
	double bestpos [DIM];	/*最適位置*/
	double bestval;			/*最適評価値*/
};
/* 関数のプロトタイプの宣言 */
void initps(struct particle ps[]);	/*粒子群の初期化*/
void printps(struct particle ps[]);	/*粒子群の表示*/
double frand();					/*実数乱数*/
double calcval(double *x);	/*評価値の計算*/
void optimize(struct particle ps[]);	/*粒子群の位置更新*/
void setgbest(struct particle ps[]);	/* 群れの最適位置を格納 */

										/*大域変数*/
double gbestpos[DIM];	/*群中の最適位置*/
double gbestval;			/*群中の最適評価値*/

							// main() 関数/
int main()
{

	struct particle ps[LIMITL];	/*粒子群*/
	int i;						/*繰り返し回数の制御*/
	srand(SEED);			/*乱数の初期化*/

							/*粒子群の初期化*/
	initps(ps);

	printps(ps);
	for (i = 0; i<ILIMIT; ++i) {
		optimize(ps);
		//テキストに保存
		printf("%d回目\n", i);
		printps(ps);
	}

	getchar();
	return 0;
}

/* optimize() 関数 */
// 粒子群の位置更新 
void optimize(struct particle ps[])
{
	int i;
	double r1, r2;/*乱数の値を格納*/
	for (i = 0; i<NOPS; i++) {
		/*乱数の設定*/
		r1 = frand();
		r2 = frand();
		/*速度の更新*/
		for (int j = 0; j < DIM; j++) {
			ps[i].v[j] = W*ps[i].v[j]
				+ C1*r1*(ps[i].bestpos[j] - ps[i].pos[j])
				+ C2*r2*(gbestpos[j] - ps[i].pos[j]);
			/*位置の更新*/
			ps[i].pos[j] += ps[i].v[j];

		}

		/*最適値の更新*/
		ps[i].value = calcval(ps[i].pos);

		if (ps[i].value<ps[i].bestval) {
			ps[i].bestval	= ps[i].value;
			for (int j = 0; j < DIM;j++) {
				ps[i].bestpos[j] = ps[i].pos[j];
			}
		}
	}
	/*群中最適値の更新*/
	setgbest(ps);
}

//*calcbal() 関数 
// 評価値の計算 
double calcval(double *x)
{
	double sum = 0;
	switch (EVALUATION_NUM)
	{
	case 1:
		for (int i = 0; i < DIM; i++) {
			sum += x[i] * x[i];
		}
		break;
	case 2:
		for (int i = 0; i < DIM; i++) {
			sum += pow(x[i], 4.0) - 16 * pow(x[i], 2.0) + 5 * x[i];
		}
		break;
	case 3:
		for (int i = 0; i < DIM; i++) {
			sum += pow(x[i], 2.0) - 10 * cos(2 * PI*x[i]) + 10;
		}
		break;

	default:
		printf("存在しない評価関数の数字です\nEVALUATION_NUMを変更してください\n");
	}
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

/* setgbest() 関数 */
// 群れの最適位置を格納 
void setgbest(struct particle ps[])
{
	int i;
	double besti;
	double x[NOPS];
	besti = ps[0].value;
	for (int for_count = 0; for_count < DIM; for_count++) {
		x[for_count] = ps[0].pos[for_count];
	}

	for (i = 0; i < NOPS; i++)
		/*現在の最良評価値を探す*/
		if (ps[i].value < besti) {
			besti = ps[i].value;
			for (int j = 0; j < DIM; j++) {
				x[j] = ps[i].pos[j];
			}

		}

	/*評価値が過去よりもよかったら更新*/
	if (besti < gbestval) {
		gbestval = besti;
		for (int j = 0; j < DIM; j++) {
			gbestpos[j] = x[j];
		}
	}

}

/* initps() 関数 */
// 粒子群の初期化 
void initps(struct particle ps[])
{
	int i;
	double x[DIM];
	for (i = 0; i < NOPS; i++) {
		/*位置*/
		for (int j = 0; j < DIM; j++) {
			x[j] = ps[i].pos[j] = frand() * 2 - 1.0;
		}
		/*評価値*/
		ps[i].value = calcval(x);
		/*速度*/
		for (int j = 0; j < DIM; j++) {
			ps[i].v[j] = frand() * 2 - 1.0;
		}
		/*最適位置*/
		for (int j = 0; j < DIM; j++) {
			ps[i].bestpos[j] = ps[i].pos[j];
		}
		/*最適評価値*/
		ps[i].bestval = ps[i].value;
	}

	/*群れの最適位置を格納*/
	gbestval = ps[0].value;
	for (int j = 0; j < DIM; j++) {
		gbestpos[j] = ps[0].pos[j];
	}

	setgbest(ps);

}


/* printps() 関数 */
// 粒子群の表示 
void printps(struct particle ps[])
{
	int i;
	for (i = 0; i<NOPS; ++i) {
		printf("粒子%d\t", i);
		for (int j = 0; j < DIM; j++)		printf("x[%d]=%.4lf\t", j, ps[i].pos[j]);
		printf("評価値%2.6lf\t", ps[i].value);
		//		printf("%lf %lf ", ps[i].v.x, ps[i].v.y);
		//		printf("%lf %lf ",ps[i].bestpos.x, ps[i].bestpos.y);
		printf("最適評価値%lf\n", ps[i].bestval);
	}
	for (int j = 0; j < DIM; j++)	printf("結果 %.6lf\t ", gbestpos[j]);
	printf("評価\t%2.6lf\n", gbestval);
}

/* frand() 関数 */
// 実数乱数 
double frand(void)
{
	double result;
	while ((result = (double)rand() / RAND_MAX) >= 1);
	return result;
}