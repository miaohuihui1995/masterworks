#include <stdio.h>
#include <stdlib.h>
//#include "omp.h"
//#include "mpi.h"

// 遍历基态数组，判断是否存在与tmp相同的基态
int CompareElementsFromBasicArray(int **list_state, int *tmp, int cnt1, int length_x, int length_y){

	for(int m = 1; m < length_x; m++){
		int cnt = 0;
		for(int n = 0; n < length_y; n++){
			if(tmp[n] == list_state[m][n]){
				cnt++;
			}
			if(cnt == length_y){
				cnt1 = 1;
			}
		}
	}
	return cnt1;
}

// 在基态数组后面添加元素
int **AppendElementToBasicArray(int **list_state, int *tmp, int length_x, int length_y){

	int i;
	int **list_bs = (int **)malloc(sizeof(int *) * (length_x + 1));
	list_bs[0] = (int *)malloc(sizeof(int) * 1);
	for(i = 1; i < length_x + 1; i++){
		list_bs[i] = (int *)malloc(sizeof(int) * length_y);
	}
	list_bs[0][0] = length_x + 1;
	for(i = 1; i < length_x; i++){
		for(int j = 0; j < length_y; j++){
			list_bs[i][j] = list_state[i][j];
		}
	}
	for(int i = 0; i < length_y; i++){
		list_bs[length_x][i] = tmp[i];
	}
	return list_bs;
}

// 用于电子在不同轨道间的移动
int **ElectronMove(int *state, int seek, int id1, int id2, int n1, int n2, int **list_state, int length_x, int length_y){
	int *tmp = (int *)malloc(sizeof(int) * length_y);
	for(int i = 0; i < length_y; i++){
		tmp[i] = state[i];
	}
	int **list_tmp = (int **)malloc(sizeof(int *) * length_x);
	list_tmp[0] = (int *)malloc(sizeof(int) * 1);
	for(int i = 1; i < length_x; i++){
		list_tmp[i] = (int *)malloc(sizeof(int) * length_y);
	}
	list_tmp[0][0] = list_state[0][0];
	for(int i = 1; i < length_x; i++){
		for(int j = 0; j < length_y; j++){
			list_tmp[i][j] = list_state[i][j];
		}
	}
	int cnt1 = 0;

	if(n1 != n2){
		if(n1 == 0 && n2 == 2){
			// 上旋电子从轨道2到轨道0
			tmp[id1] = 1;
			tmp[id2] = -1;
			tmp[seek] += 1;
			// 遍历已有基态数组，判断是否存在相同基态
			cnt1 = CompareElementsFromBasicArray(list_tmp, tmp, cnt1, length_x, length_y);
			// cnt1==0表示新基态与已有数组中的所有基态都不同，将其添加
			if(cnt1 == 0){
				list_tmp = AppendElementToBasicArray(list_tmp, tmp, length_x, length_y);
				length_x += 1;
			}
			cnt1 = 0;
			tmp[seek] -= 1;

			// 上旋电子从轨道2到轨道0
			tmp[id1] = -1;
			tmp[id2] = 1;
			tmp[seek + 1] += 1;
			cnt1 = CompareElementsFromBasicArray(list_tmp, tmp, cnt1, length_x, length_y);
			if(cnt1 == 0){
				list_tmp = AppendElementToBasicArray(list_tmp, tmp, length_x, length_y);
			}
		}

		else if(n1 == 2 && n2 == 0){
			// 上旋电子从轨道0到轨道2
			tmp[id1] = 1;
			tmp[id2] = -1;
			tmp[seek + 1] -= 1;
			if(tmp[seek + 1] >= 0){
				cnt1 = CompareElementsFromBasicArray(list_tmp, tmp, cnt1, length_x, length_y);
				if(cnt1 == 0){
					list_tmp = AppendElementToBasicArray(list_tmp, tmp, length_x, length_y);
					length_x += 1;
				}
			}
			cnt1 = 0;
			tmp[seek + 1] += 1;

			// 上旋电子从轨道0到轨道2
			tmp[id1] = -1;
			tmp[id2] = 1;
			tmp[seek] -= 1;
			if(tmp[seek] >= 0){
				cnt1 = CompareElementsFromBasicArray(list_tmp, tmp, cnt1, length_x, length_y);
				if(cnt1 == 0){
					list_tmp = AppendElementToBasicArray(list_tmp, tmp, length_x, length_y);
				}
			}
		}

		else if(n1 == 1 && n2 == -1){
			// 上旋电子从轨道0到轨道2
			tmp[id1] = 0;
			tmp[id2] = 2;
			tmp[seek] -= 1;
			if(tmp[seek] >= 0){
				cnt1 = CompareElementsFromBasicArray(list_tmp, tmp, cnt1, length_x, length_y);
				if(cnt1 == 0){
					list_tmp = AppendElementToBasicArray(list_tmp, tmp, length_x, length_y);
					length_x += 1;
				}
			}
			cnt1 = 0;
			tmp[seek] += 1;

			// 上旋电子从轨道2到轨道0
			tmp[id1] = 2;
			tmp[id2] = 0;
			tmp[seek + 1] += 1;
			cnt1 = CompareElementsFromBasicArray(list_tmp, tmp, cnt1, length_x, length_y);
			if(cnt1 == 0){
				list_tmp = AppendElementToBasicArray(list_tmp, tmp, length_x, length_y);
			}
		}

		else if(n1 == -1 && n2 == 1){
			// 上旋电子从轨道0到轨道2
			tmp[id1] = 0;
			tmp[id2] = 2;
			tmp[seek + 1] -= 1;
			if(tmp[seek + 1] >= 0){
				cnt1 = CompareElementsFromBasicArray(list_tmp, tmp, cnt1, length_x, length_y);
				if(cnt1 == 0){
					list_tmp = AppendElementToBasicArray(list_tmp, tmp, length_x, length_y);
					length_x += 1;
				}
			}
			cnt1 = 0;
			tmp[seek + 1] += 1;

			// 上旋电子从轨道2到轨道0
			tmp[id1] = 2;
			tmp[id2] = 0;
			tmp[seek] += 1;
			cnt1 = CompareElementsFromBasicArray(list_tmp, tmp, cnt1, length_x, length_y);
			if(cnt1 == 0){
				list_tmp = AppendElementToBasicArray(list_tmp, tmp, length_x, length_y);
			}
		}

		else if((n1 == 0 && n2 == 1) || (n1 == -1 && n2 == 2)){
			// 上旋电子移动
			tmp[id1] = state[id2];
			tmp[id2] = state[id1];
			tmp[seek] += 1;
			cnt1 = CompareElementsFromBasicArray(list_tmp, tmp, cnt1, length_x, length_y);
			if(cnt1 == 0){
				list_tmp = AppendElementToBasicArray(list_tmp, tmp, length_x, length_y);
			}
		}

		else if((n1 == 1 && n2 == 0) || (n1 == 2 && n2 == -1)){
			// 上旋电子移动
			tmp[id1] = state[id2];
			tmp[id2] = state[id1];
			tmp[seek] -= 1;
			if(tmp[seek] >= 0){
				cnt1 = CompareElementsFromBasicArray(list_tmp, tmp, cnt1, length_x, length_y);
				if(cnt1 == 0){
					list_tmp = AppendElementToBasicArray(list_tmp, tmp, length_x, length_y);
				}
			}
		}

		else if((n1 == 0 && n2 == -1) || (n1 == 1 && n2 == 2)){
			// 下旋电子移动
			tmp[id1] = state[id2];
			tmp[id2] = state[id1];
			tmp[seek + 1] += 1;
			cnt1 = CompareElementsFromBasicArray(list_tmp, tmp, cnt1, length_x, length_y);
			if(cnt1 == 0){
				list_tmp = AppendElementToBasicArray(list_tmp, tmp, length_x, length_y);
			}
		}

		else if((n1 == -1 && n2 == 0) || (n1 == 2 && n2 == 1)){
			// 下旋电子移动
			tmp[id1] = state[id2];
			tmp[id2] = state[id1];
			tmp[seek + 1] -= 1;
			if(tmp[seek + 1] >= 0){
				cnt1 = CompareElementsFromBasicArray(list_tmp, tmp, cnt1, length_x, length_y);
				if(cnt1 == 0){
					list_tmp = AppendElementToBasicArray(list_tmp, tmp, length_x, length_y);
				}
			}
		}
	}

	return list_tmp;
}

// 用于电子自旋变化
int **SpinChange(int nt_ph, int *state, int **list_state, int length_x, int length_y){

	int *tmp = (int *)malloc(sizeof(int) * length_y);
	for(int i = 0; i < length_y; i++){
		tmp[i] = state[i];
	}
	int **list_tmp = (int **)malloc(sizeof(int *) * length_x);
	list_tmp[0] = (int *)malloc(sizeof(int) * 1);
	for(int i = 1; i < length_x; i++){
		list_tmp[i] = (int *)malloc(sizeof(int) * length_y);
	}
	list_tmp[0][0] = list_state[0][0];
	for(int i = 1; i < length_x; i++){
		for(int j = 0; j < length_y; j++){
			list_tmp[i][j] = list_state[i][j];
		}
	}

	for(int i = nt_ph; i < length_y; i++){
		int cnt1 = 0;
		if(tmp[i] == 1){
			tmp[i] *= -1;
			tmp[nt_ph - 1] += 1;
			cnt1 = CompareElementsFromBasicArray(list_tmp, tmp, cnt1, length_x, length_y);
			if(cnt1 == 0){
				list_tmp = AppendElementToBasicArray(list_tmp, tmp, length_x, length_y);
				length_x += 1;
			}
			tmp[i] *= -1;
			tmp[nt_ph - 1] -= 1;
		}

		else if(tmp[i] == -1){
			if(tmp[nt_ph - 1] > 0){
				tmp[i] *= -1;
				tmp[nt_ph - 1] -= 1;
				cnt1 = CompareElementsFromBasicArray(list_tmp, tmp, cnt1, length_x, length_y);
				if(cnt1 == 0){
					list_tmp = AppendElementToBasicArray(list_tmp, tmp, length_x, length_y);
					length_x += 1;
				}
				tmp[i] *= -1;
				tmp[nt_ph - 1] += 1;
			}
		}
	}

	return list_tmp;
}

// 量子态变化：电子移动+自旋变化
int **StateChange(int nt_ph, int **list_state, int *state, int length_x, int length_y){
	list_state = ElectronMove(state, 0, nt_ph + 0, nt_ph + 2, state[nt_ph + 0], state[nt_ph + 2], list_state, list_state[0][0], length_y);
	list_state = ElectronMove(state, 2, nt_ph + 1, nt_ph + 2, state[nt_ph + 1], state[nt_ph + 2], list_state, list_state[0][0], length_y);
	list_state = SpinChange(nt_ph, state, list_state, list_state[0][0], length_y);
	return list_state;
}

/*未修改好
// 用于获得构建H所需的参数
int **ParametersInElectronMove(int index, int *state, int seek, int id1, int id2, int n1, int n2, int **list_state, int length_x, int length_y){
	int *tmp = (int *)malloc(sizeof(int) * length_y);
	for(int i = 0; i < length_y; i++){
		tmp[i] = state[i];
	}
	int **parameters_H = (int **)malloc(sizeof(int *) * 1);
	parameters_H[0] = (int *)malloc(sizeof(int) * 3);
	for(int i = 0; i < 3; i++){
		parameters_H[0][i] = -1;
	}
	int cnt1 = 0;
	int **parameters = (int *)malloc(sizeof(int) * 3);

	if(n1 != n2){
		if(n1 == 0 && n2 == 2){
			// 上旋电子从轨道2到轨道0
			tmp[id1] = 1;
			tmp[id2] = -1;
			tmp[seek] += 1;
			cnt1 = CompareElementsFromBasicArray(list_state, tmp, cnt1, length_x, length_y);
			if(cnt1 == 0){
				parameters[0] = index; parameters[1] = m; parameters[2] = g;
				parameters_H = AppendElementToArray(parameters_H, parameters, length_x, 3);
				length_x += 1;
			}
			cnt1 = 0;
			tmp[seek] -= 1;

			// 上旋电子从轨道2到轨道0
			tmp[id1] = -1;
			tmp[id2] = 1;
			tmp[seek + 1] += 1;
			cnt1 = CompareElementsFromBasicArray(list_state, tmp, cnt1, length_x, length_y);
			if(cnt1 == 0){
				parameters[0] = index; parameters[1] = m; parameters[2] = gg;
				parameters_H = AppendElementToArray(parameters_H, parameters, length_x, 3);
			}
		}

		else if(n1 == 2 && n2 == 0){
			// 上旋电子从轨道0到轨道2
			tmp[id1] = 1;
			tmp[id2] = -1;
			tmp[seek + 1] -= 1;
			if(tmp[seek + 1] >= 0){
				cnt1 = CompareElementsFromBasicArray(list_state, tmp, cnt1, length_x, length_y);
				if(cnt1 == 0){
					parameters[0] = index; parameters[1] = m; parameters[2] = gg;
					parameters_H = AppendElementToArray(parameters_H, parameters, length_x, 3);
					length_x += 1;
				}
			}
			cnt1 = 0;
			tmp[seek + 1] += 1;

			// 上旋电子从轨道0到轨道2
			tmp[id1] = -1;
			tmp[id2] = 1;
			tmp[seek] -= 1;
			if(tmp[seek] >= 0){
				cnt1 = CompareElementsFromBasicArray(list_state, tmp, cnt1, length_x, length_y);
				if(cnt1 == 0){
					parameters[0] = index; parameters[1] = m; parameters[2] = g;
					parameters_H = AppendElementToArray(parameters_H, parameters, length_x, 3);
				}
			}
		}

		else if(n1 == 1 && n2 == -1){
			// 上旋电子从轨道0到轨道2
			tmp[id1] = 0;
			tmp[id2] = 2;
			tmp[seek] -= 1;
			if(tmp[seek] >= 0){
				cnt1 = CompareElementsFromBasicArray(list_state, tmp, cnt1, length_x, length_y);
				if(cnt1 == 0){
					parameters[0] = index; parameters[1] = m; parameters[2] = g;
					parameters_H = AppendElementToArray(parameters_H, parameters, length_x, 3);
					length_x += 1;
				}
			}
			cnt1 = 0;
			tmp[seek] += 1;

			// 上旋电子从轨道2到轨道0
			tmp[id1] = 2;
			tmp[id2] = 0;
			tmp[seek + 1] += 1;
			cnt1 = CompareElementsFromBasicArray(list_state, tmp, cnt1, length_x, length_y);
			if(cnt1 == 0){
				parameters[0] = index; parameters[1] = m; parameters[2] = gg;
				parameters_H = AppendElementToArray(parameters_H, parameters, length_x, 3);
			}
		}

		else if(n1 == -1 && n2 == 1){
			// 上旋电子从轨道0到轨道2
			tmp[id1] = 0;
			tmp[id2] = 2;
			tmp[seek + 1] -= 1;
			if(tmp[seek + 1] >= 0){
				cnt1 = CompareElementsFromBasicArray(list_state, tmp, cnt1, length_x, length_y);
				if(cnt1 == 0){
					parameters[0] = index; parameters[1] = m; parameters[2] = gg;
					parameters_H = AppendElementToArray(parameters_H, parameters, length_x, 3);
					length_x += 1;
				}
			}
			cnt1 = 0;
			tmp[seek + 1] += 1;

			// 上旋电子从轨道2到轨道0
			tmp[id1] = 2;
			tmp[id2] = 0;
			tmp[seek] += 1;
			cnt1 = CompareElementsFromBasicArray(list_state, tmp, cnt1, length_x, length_y);
			if(cnt1 == 0){
				parameters[0] = index; parameters[1] = m; parameters[2] = g;
				parameters_H = AppendElementToArray(parameters_H, parameters, length_x, 3);
			}
		}

		else if((n1 == 0 && n2 == 1) || (n1 == -1 && n2 == 2)){
			// 上旋电子移动
			tmp[id1] = state[id2];
			tmp[id2] = state[id1];
			tmp[seek] += 1;
			cnt1 = CompareElementsFromBasicArray(list_state, tmp, cnt1, length_x, length_y);
			if(cnt1 == 0){
				parameters[0] = index; parameters[1] = m; parameters[2] = g;
				parameters_H = AppendElementToArray(parameters_H, parameters, length_x, 3);
			}
		}

		else if((n1 == 1 && n2 == 0) || (n1 == 2 && n2 == -1)){
			// 上旋电子移动
			tmp[id1] = state[id2];
			tmp[id2] = state[id1];
			tmp[seek] -= 1;
			if(tmp[seek] >= 0){
				cnt1 = CompareElementsFromBasicArray(list_state, tmp, cnt1, length_x, length_y);
				if(cnt1 == 0){
					parameters[0] = index; parameters[1] = m; parameters[2] = g;
					parameters_H = AppendElementToArray(parameters_H, parameters, length_x, 3);
				}
			}
		}

		else if((n1 == 0 && n2 == -1) || (n1 == 1 && n2 == 2)){
			// 下旋电子移动
			tmp[id1] = state[id2];
			tmp[id2] = state[id1];
			tmp[seek + 1] += 1;
			cnt1 = CompareElementsFromBasicArray(list_state, tmp, cnt1, length_x, length_y);
			if(cnt1 == 0){
				parameters[0] = index; parameters[1] = m; parameters[2] = gg;
				parameters_H = AppendElementToArray(parameters_H, parameters, length_x, 3);
			}
		}

		else if((n1 == -1 && n2 == 0) || (n1 == 2 && n2 == 1)){
			// 下旋电子移动
			tmp[id1] = state[id2];
			tmp[id2] = state[id1];
			tmp[seek + 1] -= 1;
			if(tmp[seek + 1] >= 0){
				cnt1 = CompareElementsFromBasicArray(list_state, tmp, cnt1, length_x, length_y);
				if(cnt1 == 0){
					parameters[0] = index; parameters[1] = m; parameters[2] = gg;
					parameters_H = AppendElementToArray(parameters_H, parameters, length_x, 3);
				}
			}
		}
	}

	return parameters_H;
}

// 用于电子自旋变化
int **ParametersInSpinChange(int nt_ph, int *state, int **list_state, int length_x, int length_y){

	int *tmp = (int *)malloc(sizeof(int) * length_y);
	for(int i = 0; i < length_y; i++){
		tmp[i] = state[i];
	}
	int **list_tmp = (int **)malloc(sizeof(int *) * length_x);
	list_tmp[0] = (int *)malloc(sizeof(int) * 1);
	for(int i = 1; i < length_x; i++){
		list_tmp[i] = (int *)malloc(sizeof(int) * length_y);
	}
	list_tmp[0][0] = list_state[0][0];
	for(int i = 1; i < length_x; i++){
		for(int j = 0; j < length_y; j++){
			list_tmp[i][j] = list_state[i][j];
		}
	}

	for(int i = nt_ph; i < length_y; i++){
		int cnt1 = 0;
		if(tmp[i] == 1){
			tmp[i] *= -1;
			tmp[nt_ph - 1] += 1;
			cnt1 = CompareElementsFromBasicArray(list_tmp, tmp, cnt1, length_x, length_y);
			if(cnt1 == 0){
				list_tmp = AppendElementToBasicArray(list_tmp, tmp, length_x, length_y);
				length_x += 1;
			}
			tmp[i] *= -1;
			tmp[nt_ph - 1] -= 1;
		}

		else if(tmp[i] == -1){
			if(tmp[nt_ph - 1] > 0){
				tmp[i] *= -1;
				tmp[nt_ph - 1] -= 1;
				cnt1 = CompareElementsFromBasicArray(list_tmp, tmp, cnt1, length_x, length_y);
				if(cnt1 == 0){
					list_tmp = AppendElementToBasicArray(list_tmp, tmp, length_x, length_y);
					length_x += 1;
				}
				tmp[i] *= -1;
				tmp[nt_ph - 1] += 1;
			}
		}
	}

	return list_tmp;
}

int **ParametersInStateChange(int nt_ph, int **list_state, int *state, int length_x, int length_y){
	list_state = ElectronMove(state, 0, nt_ph + 0, nt_ph + 2, state[nt_ph + 0], state[nt_ph + 2], list_state, list_state[0][0], length_y);
	list_state = ElectronMove(state, 2, nt_ph + 1, nt_ph + 2, state[nt_ph + 1], state[nt_ph + 2], list_state, list_state[0][0], length_y);
	list_state = SpinChange(nt_ph, state, list_state, list_state[0][0], length_y);
	return list_state;
}*/

// 在数组后面添加元素
int **AppendElementToArray(int **list, int *tmp, int length_x, int length_y){

	int **list_tmp = (int **)malloc(sizeof(int *) * (length_x + 1));
	for(int i = 0; i < length_x + 1; i++){
		list_tmp[i] = (int *)malloc(sizeof(int) * length_y);
	}
	for(int i = 0; i < length_x; i++){
		for(int j = 0; j < length_y; j++){
			list_tmp[i][j] = list[i][j];
		}
	}
	for(int i = 0; i < length_y; i++){
		list_tmp[length_x][i] = tmp[i];
	}
	return list_tmp;
}

int main(int argc, char **argv){

	// 一个\Lambda原子组成的系统，三个电子，初态为|0 0 0 0 0>|1 0 2>
 
	// ---------- 基态矩阵 ----------
 
	// 初始条件
	double tau = 0.01;
	int T = 1401;

	// 这里的g和gamma统一为整型
	int nt_ph = 5; // number of types of photons
	int np_at = 3; // number of places of atoms
	double g0 = 4; // g_{\Omega\uparrow}    
	double g1 = 4; // g_{\Omega\downarrow}
	double g2 = 2; // g_{\omega\uparrow}
	double g3 = 2;	// g_{\omega\downarrow}
	double g4 = 1; // g_{\omege_{spin}}
	double gamma0 = 1; // \gamma_{\Omega\uparrow}
	double gamma1 = 1; // \gamma_{\Omega\downarrow}
	double gamma2 = 1; // \gamma_{\omega\uparrow}
	double gamma3 = 1; // \gamma_{\omega\downarrow}
	double gamma4 = 1; // \gamma_{\omega_{spin}}
	int *initial_at = (int *)malloc(sizeof(int) * np_at);
	initial_at[0] = 1; initial_at[1] = 1; initial_at[2] = 1;
	// [orbit0, orbit1, orbit2]：0表示轨道上无电子，1表示轨道上一个自旋向上电子>，-1表示轨道上一个自旋向下电子，2表示轨道上两个电子
	int *initial_ph = (int *)malloc(sizeof(int) * nt_ph);
	initial_ph[0] = 0; initial_ph[0] = 0; initial_ph[0] = 0; initial_ph[0] = 0; initial_ph[0] = 0;
	//	[\Omega\uparrow, \Omega\downarrow, \omega\uparrow, \omega\downarrow, \omega_{spin}]
	int *initial_bs = (int *)malloc(sizeof(int) * (nt_ph + np_at)); //初态
	for(int i = 0; i < nt_ph; i++){
		initial_bs[i] = initial_ph[i];
	}
	for(int i = 0; i < np_at; i++){
		initial_bs[i + nt_ph] = initial_at[i];
	}
	// 基态列表：共八位（前五位对应五类光子，后三位对应原子三个轨道）
	int length_x = 2;
	int length_y = nt_ph + np_at;
	int **list_bs = (int **)malloc(sizeof(int *) * length_x);
	list_bs[0] = (int *)malloc(sizeof(int) * 1);
	for(int i = 1; i < length_x; i++){
		list_bs[i] = (int *)malloc(sizeof(int) * length_y);
	}
	list_bs[0][0] = length_x;
	for(int i = 1; i < length_x; i++){
		for(int j = 0; j < length_y; j++){
			list_bs[i][j] = initial_bs[j];
		}
	}
	printf("Initial state:\n");
	for(int i = 1; i < length_x; i++){
		for(int j = 0; j < nt_ph + np_at; j++){
			printf("%d\t", list_bs[i][j]);
		}
		printf("\n");
	}
	printf("\n");

	// 生成基态列表（不含sink）
	int flag = 1;
	int length_x_last = 0;

	while(flag){
		for(int i = length_x_last + 1; i < length_x; i++){
			list_bs = StateChange(nt_ph, list_bs, list_bs[i], list_bs[0][0], length_y);
		}

		if(list_bs[0][0] - 1 == length_x_last){
			flag = 0;
		}
		else{
			length_x_last = length_x - 1;
			length_x = list_bs[0][0];
		}
	}
	printf("List of basic states(without sink):\n");
	for(int i = 1; i < list_bs[0][0]; i++){
		for(int j = 0; j < length_y; j++){
			printf("%d\t", list_bs[i][j]);
		}
		printf("\n");
	}
	printf("Length: %d\n", list_bs[0][0] - 1);
	printf("\n");

	//生成带sink所有可能基态列表
	length_x = 1;
	int **list_sink = (int **)malloc(sizeof(int *) * length_x);
	list_sink[0] = (int *)malloc(sizeof(int) * 1);
	list_sink[0][0] = length_x;

	for(int i = 1; i < list_bs[0][0]; i++){
		int sum = 0;
		int cnt = 0;
		int length = 1;
		int **index = (int **)malloc(sizeof(int *) * length);
		index[0] = (int *)malloc(sizeof(int) * 1);
		index[0][0] = -1;
		for(int j = 0; j < nt_ph; j++){
			sum += list_bs[i][j];
			if(list_bs[i][j] > 0){
				index = AppendElementToArray(index, &j, length, 1);
				length += 1;
				cnt += 1;
			}
		}
		int **index1 = (int **)malloc(sizeof(int *) * (length - 1));
		for(int j = 0; j < length - 1; j++){
			index1[j] = (int *)malloc(sizeof(int) * 1);
		}
		for(int j = 1; j < length; j++){
			index1[j - 1][0] = index[j][0];
		}
		int length_index1 = length;
		//printf("sum = %d\n", sum);
		/*printf("Index:\n");
		for(int j = 0; j < length - 1; j++){
			printf("%d\t", index1[j][0]);
		}
		printf("\n");*/

		int length1; //这两个必须在if语句之外声明，否则无法用在下一个if语句中
		int **arr2;
		if(sum > 0){
			length = 1;
			int **arr = (int **)malloc(sizeof(int *) * length);
			arr[0] = (int *)malloc(sizeof(int) * cnt);
			for(int j = 0; j < cnt; j++){
				arr[0][j] = -1;
			}
			for(int m0 = 0; m0 < sum + 1; m0++){
				if(cnt > 1){
					for(int m1 = 0; m1 < sum + 1; m1++){
						if(cnt > 2){
							for(int m2 = 0; m2 < sum + 1; m2++){
								if(cnt > 3){
									for(int m3 = 0; m3 < sum + 1; m3++){
										if(cnt > 4){
											for(int m4 = 0; m4 < sum + 1; m4++){
												int *arr_tmp = (int *)malloc(sizeof(int) * 5); 
												arr_tmp[0] = m0; arr_tmp[1] = m1; arr_tmp[2] = m2; arr_tmp[3] = m3; arr_tmp[4] = m4;
												arr = AppendElementToArray(arr, arr_tmp, length, 5);
												length += 1;
											}
										}
										else{
											int *arr_tmp = (int *)malloc(sizeof(int) * 4); 
											arr_tmp[0] = m0; arr_tmp[1] = m1; arr_tmp[2] = m2; arr_tmp[3] = m3;
											arr = AppendElementToArray(arr, arr_tmp, length, 4);
											length += 1;
										}
									}
								}
								else{
									int *arr_tmp = (int *)malloc(sizeof(int) * 3); 
									arr_tmp[0] = m0; arr_tmp[1] = m1; arr_tmp[2] = m2;
									arr = AppendElementToArray(arr, arr_tmp, length, 3);
									length += 1;
								}
							}
						}
						else{
							int *arr_tmp = (int *)malloc(sizeof(int) * 2); 
							arr_tmp[0] = m0; arr_tmp[1] = m1;
							arr = AppendElementToArray(arr, arr_tmp, length, 2);
							length += 1;
						}
					}
				}
				else{
					int *arr_tmp = (int *)malloc(sizeof(int) * 1); 
					arr_tmp[0] = m0;
					arr = AppendElementToArray(arr, arr_tmp, length, 1);
					length += 1;
				}
			}

			length1 = 1;
			int **arr1 = (int **)malloc(sizeof(int *) * length1);
			arr1[0] = (int *)malloc(sizeof(int) * cnt);
			for(int j = 0; j < cnt; j++){
				arr1[0][j] = -1;
			}
			for(int j = 1; j < length; j++){
				int sum1 = 0;
				for(int k = 0; k < cnt; k++){
					sum1 += arr[j][k];
				}
				if(sum1 <= sum){
					arr1 = AppendElementToArray(arr1, arr[j], length1, cnt);
					length1 += 1;
				}
			}

			arr2 = (int **)malloc(sizeof(int *) * (length1 - 1));
			for(int j = 0; j < length1 - 1; j++){
				arr2[j] = (int *)malloc(sizeof(int) * cnt);
			}
			for(int j = 1; j < length1; j++){
				for(int k = 0; k < cnt; k++){
					arr2[j - 1][k] = arr1[j][k];
				}
			}

			/*printf("arr2:\n");
			for(int j = 0; j < length1 - 1; j++){
				for(int k = 0; k < cnt; k++){
					printf("%d\t", arr2[j][k]);
				}
				printf("\n");
			}
			printf("\n");*/
		}

		if(sum > 0){
			for(int j = 0; j < length1 - 1; j++){
				int *bs_tmp = (int *)malloc(sizeof(int) * length_y);
				for(int k = 0; k < length_y; k++){
					bs_tmp[k] = list_bs[i][k];
				}
				int cnt1 = 0;
				if(bs_tmp[index1[0][0]] - arr2[j][0] >= 0){
					if(cnt > 1 && bs_tmp[index1[1][0]] - arr2[j][1] >= 0){
						if(cnt > 2 && bs_tmp[index1[2][0]] - arr2[j][2] >= 0){
							if(cnt > 3 && bs_tmp[index1[3][0]] - arr2[j][3] >= 0){
								if(cnt > 4 && bs_tmp[index1[4][0]] - arr2[j][4] >= 0){
									bs_tmp[index1[4][0]] -= arr2[j][4];
									bs_tmp[index1[3][0]] -= arr2[j][3];
									bs_tmp[index1[2][0]] -= arr2[j][2];
									bs_tmp[index1[1][0]] -= arr2[j][1];
									bs_tmp[index1[0][0]] -= arr2[j][0];
									cnt1 += 5;
									if(cnt1 == length_index1 - 1){
										list_sink = AppendElementToBasicArray(list_sink, bs_tmp, list_sink[0][0], length_y);
									}
									continue;
								}
								else{
									bs_tmp[index1[3][0]] -= arr2[j][3];
									bs_tmp[index1[2][0]] -= arr2[j][2];
									bs_tmp[index1[1][0]] -= arr2[j][1];
									bs_tmp[index1[0][0]] -= arr2[j][0];
									cnt1 += 4;
									if(cnt1 == length_index1 - 1){
										list_sink = AppendElementToBasicArray(list_sink, bs_tmp, list_sink[0][0], length_y);
									}
									continue;
								}
							}
							else{
								bs_tmp[index1[2][0]] -= arr2[j][2];
								bs_tmp[index1[1][0]] -= arr2[j][1];
								bs_tmp[index1[0][0]] -= arr2[j][0];
								cnt1 += 3;
								if(cnt1 == length_index1- 1){
									list_sink = AppendElementToBasicArray(list_sink, bs_tmp, list_sink[0][0], length_y);
								}
								continue;
							}
						}
						else{
							bs_tmp[index1[1][0]] -= arr2[j][1];
							bs_tmp[index1[0][0]] -= arr2[j][0];
							cnt1 += 2;
							if(cnt1 == length_index1 - 1){
								list_sink = AppendElementToBasicArray(list_sink, bs_tmp, list_sink[0][0], length_y);
							}
							continue;
						}
					}
					else{
						bs_tmp[index1[0][0]] -= arr2[j][0];
						cnt1 += 1;
						if(cnt1 == length_index1 - 1){
							list_sink = AppendElementToBasicArray(list_sink, bs_tmp, list_sink[0][0], length_y);
						}
						continue;
					}
				}
			}
		}
		else{
			list_sink = AppendElementToBasicArray(list_sink, list_bs[i], list_sink[0][0], length_y);
		}
	}

	// list_sink去重（有问题）
	length_x = 2;
	int **list_sink_unique = (int **)malloc(sizeof(int *) * length_x);
	list_sink_unique[0] = (int *)malloc(sizeof(int) * 1);
	list_sink_unique[1] = (int *)malloc(sizeof(int) * length_y);
	list_sink_unique[0][0] = length_x;
	for(int i = 0; i < length_y; i++){
		list_sink_unique[1][i] = list_sink[1][i];
	}

	for(int i = 1; i < list_sink[0][0]; i++){
		int cnt = 0;
		cnt = CompareElementsFromBasicArray(list_sink_unique, list_sink[i], cnt, list_sink_unique[0][0], length_y);
		if(cnt == 0){
			list_sink_unique = AppendElementToBasicArray(list_sink_unique, list_sink[i], list_sink_unique[0][0], length_y);
		}
	}

	printf("List of basic states(with sink):\n");
	for(int i = 1; i < list_sink_unique[0][0]; i++){
		for(int j = 0; j < length_y; j++){
			printf("%d\t", list_sink_unique[i][j]);
		}
		printf("\n");
	}
	printf("Length: %d\n", list_sink_unique[0][0] - 1);
	printf("\n");

	// 构建H所需的参数




	return 0;
}
