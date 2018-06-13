#include <stdio.h>
#include <string.h>
#include <stdlib.h>


#define SIZE 20

void Counting_sort(int a[],int k);
void RADIX_sort(char **a, char **b,  int k,int index);
int main() {
	int pro[] = {6,0,2,0,1,3,4,6,1,3,2};
	int count = 11 ;
	int i; 

	//printf("\nSIZE of array : %d \n", count);
	
	printf("\n1. Counting sort.{6,0,2,0,1,3,4,6,1,3,2}\n");
	/*for (i = 0; i < count; i++) {
		printf("%d ", pro[i]);
	}*/
	Counting_sort(pro, count);

	printf("\n2. Radix sort\n");

	char *cha[SIZE] = {"COW","DOG","SEA","RUG","ROW","MOB","BOX","TAB","BAR","EAR","TAR","DIG","BIG","TEA","NOW","FOX" };
	char *fir[SIZE], *sec[SIZE], *tir[SIZE];
	count = 16;
	RADIX_sort(cha, fir, count, 2);
	RADIX_sort(fir, sec, count, 1);
	RADIX_sort(sec, tir, count, 0);
	printf("unsorted   First_sort   Second_sort   Third_sort\n");

	for (i = 0; i < count;i++) {
		printf("%8s %12s %13s %12s\n", cha[i], fir[i], sec[i], tir[i]);

	}
	

	system("pause");
}




void Counting_sort(int a[],int k) {
	int c[SIZE];
	int b[SIZE];
	int i, j;
	for (i = 0; i < SIZE; i++) {
		c[i] = 0;
	}
	printf("This is the result of first loop(initiation array c)\n");
	for (i = 0; i < 7; i++) {
		printf("%d ", c[i]);
	}
	printf("\n");
	for (j = 0; j < k; j++) {
		c[a[j]] = c[a[j]] + 1;
	}
	printf("This is the result of second loop\n");
	for (j = 0; j < 7; j++) {
		printf("%d ", c[j]);
	}
	printf("\n");

	for (i = 1; i < 7; i++) {
		c[i] = c[i] + c[i - 1];
	}
	printf("This is the result of third loop\n");
	for (i = 0; i < 7; i++) {
		printf("%d ", c[i]);
	}
	printf("\n");

	
	for (j = k-1; j >= 0; j--) {
		b[c[a[j]]] = a[j];
		c[a[j]] = c[a[j]] - 1;
	}
	printf("This is the result of forth loop(sorted array)\n");
	for (j = 1; j <= k ; j++) {
		printf("%d ", b[j]);
	}
	printf("\n");

}

void RADIX_sort(char **a, char **b, int k, int index) {
	int i,j,p;
	i = index;
	int min;
	char *temp;

	for (j = 0; j < k;j++) {
		//strcpy_s(*cpy[j], sizeof(*cpy[j]), a[j]);
		b[j] = a[j];
	}
	for (j = 0; j < k; j++) {
		min = j;
		for (p = j + 1; p < k; p++) {
			if (b[p][i] < b[min][i]) {
				min = p;
			}
		}
		//strncpy_s(temp, sizeof(temp), a[j], 4);
		//strncpy_s(a[j], sizeof(a[j]), a[min],4);
		//strncpy_s(a[min], sizeof(a[min]), temp,4);
		temp = b[j];
		b[j] = b[min];
		b[min] = temp;
	}
}