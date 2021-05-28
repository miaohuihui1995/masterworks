#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "omp.h"
#define MAX_LINE 1024
#define NUM_THREADS 4
#define NUM_ITERATIONS 10

// 读取整型
void *ReadInt(FILE *fp, int *list){
	char buf[MAX_LINE];
	int cnt = 0;
	while(fgets(buf, MAX_LINE, fp) != NULL){
		buf[strlen(buf) - 1] = '\0';  // 去掉换行符
		list[cnt] = 0;
		for(int i = 0; i < strlen(buf); i++){
			list[cnt] += ((char)buf[i] - '0') * pow(10, strlen(buf) - 1 - i); // buf[i]是字符对应的ASCII码，所以需要转换为char型，再减去'0'
		}
		cnt++;
	}
}

// 读取绝对值小于1的double型
void *ReadDouble(FILE *fp, double *list){
	char buf[MAX_LINE];
	int cnt = 0;
	while(fgets(buf, MAX_LINE, fp) != NULL){
		buf[strlen(buf) - 1] = '\0';
		if(buf[0] == '-'){
			for(int i = 3; i < strlen(buf); i++){
				list[cnt] += ((char)buf[i] - '0') * pow(10, (2 - i));
			}
			list[cnt] *= -1;
			if(buf[1] == '1'){
				list[cnt] += 1;
			}
		}
		else{
			for(int i = 2; i < strlen(buf); i++){
				list[cnt] += ((char)buf[i] - '0') * pow(10, (1 - i));
			}
			if(buf[0] == '1'){
				list[cnt] += 1;
			}
		}
		cnt++;
	}
}

// 读取基态列表
void *ReadListSink(FILE *fp, int *list){
	char buf[MAX_LINE];
	int cnt = 0;
	while(fgets(buf, MAX_LINE, fp) != NULL){
		buf[strlen(buf) - 1] = '\0'; 
		int sign = 1;
		int flag = 1;
		for(int i = 0; i < strlen(buf); i++){
			if(flag){
				if(buf[i] == '-'){
					sign = -1;
				}
				else{
					list[cnt] = sign * ((char)buf[i] - '0');
					sign = 1;
					flag = 0;
					cnt++;
				}
			}
			else{
				flag = 1;
			}
		}
	}
}

// method pade for exponential of sparse matrix
double expspm(int sign, int len_H, int *index_i_H, int *index_j_H, int *value_H, double tau){

	#pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic, NUM_ITERATIONS)
	for(int i = 0; i < len_H; i++){
		value_H[i] *= sign * tau;
	}
}


// 系数矩阵exp_H_p, exp_H_m与密度矩阵相乘
void *SpMRho(int len, int length, double *rho, double *newrho, int *index_i, int *index_j, double *real, double *imag){
	printf("********");
	int a = length * length;
	for(int i = 0; i < len; i++){
		for(int j = 0; j < length; j++){
			newrho[index_i[i] * length + j] += real[i] * rho[index_j[i] * length + j] - imag[i] * rho[a + index_j[i] * length + j];
			newrho[a + index_i[i] * length + j] += imag[i] * rho[index_j[i] * length + j] + real[i] * rho[a + index_j[i] * length + j];
		}
	}
}

int main(){
	 
	// 这些参数需手动填写
	const int T = 1401;
	const double tau = 0.01;
	const int nt_ph = 5;
	const int nt_at = 3;
	const int len_state = 93;
	const int len_H = 109;
	const int len_L = 118;

	const int nt_st = nt_ph + nt_at;
	FILE *fp;

	// 提取构建H的所有参数
	fp = fopen("index_i_H.txt","r");
	int *index_i_H = (int *)malloc(sizeof(int) * len_H);
	ReadInt(fp, index_i_H);
	fclose(fp);

	fp = fopen("index_j_H.txt","r");
	int *index_j_H = (int *)malloc(sizeof(int) * len_H);
	ReadInt(fp, index_j_H);
	fclose(fp);

	fp = fopen("value_H.txt","r");
	int *value_H = (int *)malloc(sizeof(int) * len_H);
	ReadInt(fp, value_H);
	fclose(fp);

	// 提取构建L算符的所有参数
	fp = fopen("index_i_L.txt","r");
	int *index_i_L = (int *)malloc(sizeof(int) * len_L);
	ReadInt(fp, index_i_L);
	fclose(fp);

	fp = fopen("index_j_L.txt","r");
	int *index_j_L = (int *)malloc(sizeof(int) * len_L);
	ReadInt(fp, index_j_L);
	fclose(fp);

	fp = fopen("value_L.txt","r");
	int *value_L = (int *)malloc(sizeof(int) * len_L);
	ReadInt(fp, value_L);
	fclose(fp);

	// 提取基态列表
	fp = fopen("list_sink.txt","r");
	int *list_state = (int *)malloc(sizeof(int) * len_state * nt_st);
	ReadListSink(fp, list_state);
	fclose(fp);

	// 幺正演变

	// 初态
	double *rho_real = (double *)malloc(sizeof(double) * len_state * len_state);
	double *rho_imag = (double *)malloc(sizeof(double) * len_state * len_state);
	#pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic, NUM_ITERATIONS)
	for(int i = 0; i < len_state * len_state; i++){
		rho_real[i] = 0.;
		rho_imag[i] = 0.;
	}
	rho_real[0] = 1;

/*	double *newrho = (double *)malloc(sizeof(double) * len_state * len_state * 2);
	for(int i = 0; i < len_state * len_state * 2; i++){
		newrho[i] = 0;
	}
	double *newrho_tmp = (double *)malloc(sizeof(double) * len_state * len_state * 2);
	for(int i = 0; i < len_state * len_state * 2; i++){
		newrho_tmp[i] = 0;
	}

	for(int t = 0; t < T; t++){
		printf("t: %.2fs\n", t * tau);
		printf("********");
		SpMRho(len1, len_state, rho, newrho_tmp, index_i_exp_H_p, index_j_exp_H_p, real_exp_H_p, imag_exp_H_p);
		SpMRho(len2, len_state, newrho_tmp, newrho, index_i_exp_H_m, index_j_exp_H_m, real_exp_H_m, imag_exp_H_m);
//		double *tmp = SpMRho(len2, len_state, SpMRho(len1, len_state, rho, index_i_exp_H_p, index_j_exp_H_p, real_exp_H_p, imag_exp_H_p), index_i_exp_H_m, index_j_exp_H_m, real_exp_H_m, imag_exp_H_m);*/
		/*for(int i = 0; i < len_L; i++){
			int sum1 = 0;
			int sum2 = 0;
			int tmp = -1;
			int k0 = index_i_L[i];
			int k1 = index_j_L[i];
			int k2 = value_L[i];
			for(int j = 0; j < nt_ph; j++){
				sum1 += list_sink[k0 * nt_st + j];
				sum2 += list_sink[k1 * nt_st + j];
			}
			if(sum1 < sum2){
				tmp = k0;
				k0 = k1;
				k1 = tmp;
			}
			newrho_real[k1 * len_state + k1] += k2 * rho_tmp_real[k0 * len_state + k0] * tau;
			newrho_imag[k1 * len_state + k1] += k2 * rho_tmp_imag[k0 * len_state + k0] * tau;
			for(int j = 0; j < len_state; j++){
				newrho_real[k0 * len_state + j] -= 0.5 * k2 * rho_tmp_real[k0 * len_state + j] * tau;
				newrho_imag[k0 * len_state + j] -= 0.5 * k2 * rho_tmp_imag[k0 * len_state + j] * tau;
			}
			for(int j = 0; j < len_state; j++){
				newrho_real[j * len_state + k0] -= 0.5 * k2 * rho_tmp_real[j * len_state + k0] * tau;
				newrho_imag[j * len_state + k0] -= 0.5 * k2 * rho_tmp_imag[j * len_state + k0] * tau;
			}
		}
		rho_real = newrho_real;
		rho_imag = newrho_imag;
		free(rho_tmp_real);
		rho_tmp_real = NULL;
		free(rho_tmp_imag);
		rho_tmp_imag = NULL;
		free(newrho_real);
		newrho_real = NULL;
		free(newrho_imag);
		newrho_imag = NULL;
	}*/

	free(index_i_H);
	index_i_H = NULL;
	free(index_j_H);
	index_j_H = NULL;
	free(value_H);
	value_H = NULL;
	free(index_i_H);
	index_i_L = NULL;
	free(index_j_L);
	index_j_L = NULL;
	free(value_L);
	value_L = NULL;
	free(list_state);
	list_state = NULL;
	free(rho_real);
	rho_real = NULL;
	free(rho_imag);
	rho_imag = NULL;

	return 0;
}
