
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <assert.h>

#define SIZE 20
#define IFN 99999
//#define _NO_CRT_STDIO_INLINE
#define _GNU_SOURCE 1
#define __STDC_WANT_LIB_EXT2__ 1
typedef struct nodes node;
enum id { white = 0, gray, black };


int predecessor_ap[SIZE][SIZE];
int d[SIZE][SIZE];

typedef struct nodes {
	char name[SIZE];
	int index;
	int adj[SIZE];
	int weight;
	//int edge;
	int d;
	char predecessor[20];
	enum id color;
}node;

void init_node(node *vertex, int count, char ***matrix, char **colnames) {
	int row, col, index;
	index = 0;
	//edge = 0;
	for (col = 0; col < count; col++, index++) {
		strcpy_s(vertex[index].name, SIZE * sizeof(char), colnames[col]);
		//printf("%s ",vertex[index].name);
		vertex[index].index = col;
	}
	for (index = 0, row = 0; row < count; row++, index++) {
		strcpy_s(vertex[index].predecessor, SIZE * sizeof(char), "\0");
		vertex[index].d = 0;
		for (col = 0; col < count; col++) {
			if (row == col) {
				vertex[index].adj[(col)] = 0;
			}
			else
				vertex[index].adj[(col)] = matrix[row][col];
		}

		//vertex[index].edge = edge;
	}
}

void init_single_source(node *vertex, int s, int count) {
	int i;
	for (i = 0; i < count; i++) {
		strcpy_s(vertex[i].predecessor, SIZE * sizeof(char), "\0");
		vertex[i].d = IFN;  //무한을 IFN으로
		vertex[i].color = white;
	}
	//printf("good\n");
	vertex[s].d = 0;
	vertex[s].color = black;
}

void Relax(node *vertex, int row, int col, int w) {
	int temp;
	temp = vertex[row].d + w;
	if (vertex[col].d > temp) {
		vertex[col].d = temp;
		//vertex[col]->predecessor =
		strcpy_s(vertex[col].predecessor, SIZE * sizeof(char), vertex[row].predecessor);
	}
}
int Extract_min(node *vertex, int count) {
	int min_edge = IFN;
	int index = 0;
	int i;
	for (i = 0; i < count; i++) {
		if (vertex[i].color == gray && vertex[i].d < min_edge) {
			min_edge = vertex[i].d;
			index = i;
		}
	}
	return index;
}
void search_and_change_color(node *vertex, int index, int count) {
	int i;
	for (i = 0; i < count; i++) {
		if (vertex[i].color == black)
			continue;
		else if (vertex[index].adj[i] != 0 && vertex[index].adj[i] != IFN) {
			vertex[i].color = gray;
		}

	}
}
void Dijkstra(node *vertex, int count, int s) {
	//printf("good\n");
	init_single_source(vertex, s, count);
	//printf("good\n");
	//init_color(vertex,count);
	int i;
	int u;
	int remain = count;
	for (i = 0; i < count; i++) {
		if (vertex[s].adj[i] != IFN && vertex[s].adj[i] != 0) {
			vertex[i].color = gray;
			//Q[i] = vertex[s].adj[i];
			vertex[i].d = vertex[s].adj[i];
			//printf("%d\n", vertex[i].d);
		}
	}
	remain--;
	while (remain != 0) {
		u = Extract_min(vertex, count);
		vertex[u].color = black;
		search_and_change_color(vertex, u, count);
		for (i = 0; i < count; i++) {
			if (vertex[u].adj[i] != IFN && vertex[u].adj[i] != 0) {
				Relax(vertex, u, i, vertex[u].adj[i]);
				//vertex[i].color = gray;
			}
		}
		remain--;
	}
}
void Bellman_Ford(node *vertex, int count, int s) {
	int i, j, k;
	k = 0;
	int seq[SIZE];
	//seq = malloc(SIZE * sizeof(int*));

	init_single_source(vertex, s, count);

	for (i = 0; i < count; i++) {
		if (vertex[s].adj[i] != IFN && vertex[s].adj[i] != 0) {
			//printf("%d ",vertex[s].adj[i]);
			vertex[i].d = vertex[s].adj[i];
			seq[k] = i;
			//printf("seq[%d] : %d %s\n",k,seq[k],vertex[i].name);
			k++;
			vertex[s].color = gray;


		}
	}
	int index = 0;

	for (i = 0; i < (count-1); i++) {
		//printf("seq[%d] : %d \n",k,seq[k]);
		index = seq[i];
		for (j = 0; j < count; j++) {
			if (vertex[index].adj[j] != IFN && vertex[index].adj[j] && vertex[j].color != gray && vertex[j].d != 0) {
				vertex[j].color = gray;
				seq[k] = j;
				k++;
			}
			if (vertex[index].adj[j] != IFN && vertex[index].adj[j]) {
				Relax(vertex, index, j, vertex[index].adj[j]);
			}
		}
		//remain--;
	}
	//free(seq);
}

void Floyd_Warshall(node *vertex, int count) {
	int row, col, index;

	for (row = 0; row < count; row++) {
		for (col = 0; col < count; col++) {
			d[row][col] = vertex[row].adj[col];
			predecessor_ap[row][col] = -1;
		}
	}
	for (index = 0; index < count; index++) {
		for (row = 0; row < count; row++) {
			for (col = 0; col <count; col++) {
				if (d[row][col] > d[row][index] + d[index][col]) {
					d[row][col] = d[row][index] + d[index][col];
					if (predecessor_ap[index][col] == -1) {
						predecessor_ap[row][col] = index;
					}
					else {
						predecessor_ap[row][col] = predecessor_ap[index][col];
					}
				}
			}
		}
	}
}
int a2i(const char *s)
{
	int sign = 1;
	if (s[0] == '-') {
		sign = -1;
		s++;
	}
	//s++;
	int num = 0;
	while (*s)
	{
		num = ((*s) - '0') + (num * 10);
		s++;
	}
	return num * sign;
}
int main() {
	errno_t err;

	FILE *fin;
	printf("Please, read README in the zip file.\n");


	err = fopen_s(&fin, "D:/git/SungminCode/HGU/Class/Algorithm/hw5/hw5..data", "r");

	if (fin)
		printf("File is open!\n");
	else {
		printf("File is Error\n");
		system("pause");
		return 0;
	}

	int i, j;
	int count = (int)malloc(sizeof(int));
	count = 0;
	int **matrix;
	matrix = malloc(SIZE * sizeof(int*));
	assert(matrix != NULL);
	for (i = 0; i < SIZE; i++) {
		matrix[i] = malloc(SIZE * sizeof(int*));
		assert(matrix[i] != NULL);

	}
	char **rownames;
	char **colnames;

	colnames = malloc(SIZE * sizeof(char*));
	assert(colnames != NULL);
	for (i = 0; i < SIZE; i++) {
		colnames[i] = malloc(SIZE * sizeof(char));
		assert(colnames[i] != NULL);

	}
	rownames = malloc(SIZE * sizeof(char*));
	assert(rownames != NULL);
	for (i = 0; i < SIZE; i++) {
		rownames[i] = malloc(SIZE * sizeof(char));
		assert(rownames[i] != NULL);
	}



	//printf("Please, input file address : ");
	//scanf_s("%s", address);
	//err = fopen_s(&fin, address, "r");

	//char *result;


	const char delimiter[] = " \t";
	char line[256];
	char *tmp;
	char *context;
	int col, row;
	fgets(line, 255, fin);

	tmp = strtok_s(line, delimiter, &context);
	printf("");
	strcpy_s(colnames[0], SIZE * sizeof(char), tmp);
	//printf("0:1 : %s\n", matrix[0][count]);
	for (j = 1;; j++) {
		tmp = strtok_s(NULL, delimiter, &context);
		if (tmp == NULL) {
			break;
		}

		strcpy_s(colnames[j], SIZE * sizeof(char), tmp);
		//printf("0:%d : %s\n", j, colnames[j]);
		;
	}
	count = j;
	printf("count : %d\n", count);

	for (row = 0; row < count; row++) {
		//printf("output\n");
		fgets(line, 255, fin);
		tmp = strtok_s(line, delimiter, &context);
		if (tmp == NULL) {
			break;
		}
		strcpy_s(rownames[row], SIZE * sizeof(char), tmp);
		for (col = 0; col < count; col++) {
			tmp = strtok_s(NULL, delimiter, &context);
			if (tmp == NULL) {
				break;
			}
			else if (row == col) {
				matrix[row][col] = 0;
			}
			else if ((strncmp(tmp, "INF", 3) == 0)) {
				matrix[row][col] = IFN;
			}
			else {
				matrix[row][col] = atoi(tmp);
			}

		}
	}
	/*printf("          ");
	for (col = 0; col < count; col++) {
		printf("%-10s", colnames[col]);
	}

	for (row = 0; row < count; row++) {
		printf("\n%-10s", rownames[row]);
		for (col = 0; col< count; col++) {
			printf("%-10d ", matrix[row][col]);
			//puts(matrix[row][col]);
		}
		printf("\n");
	}*/
	node *vertex = malloc(count * sizeof(node));

	//printf("\nOutput\n");
	printf("\n1. Dijkstra's Algorithm\n\n");
	init_node(vertex, count, matrix, colnames);
	printf("Success! nodes initialization \n");

	clock_t start = clock();

	for (i = 0; i < count; i++) {
		//printf("hi\n");
		Dijkstra(vertex, count, j);
	}
	start = clock() - start;
	double duration = ((double)start) / CLOCKS_PER_SEC;
	//double duration = ((double)start);
	//double duration =(double)(clock() - start);

	printf("It take %f time seconds to compute shortest path between cities with Dijkstra's algorithm as follows. \n", duration);

	printf("          ");
	for (col = 0; col < count; col++) {
		printf("%-10s", colnames[col]);
	}
	for (row = 0; row < count; row++) {
		Dijkstra(vertex, count, row);
		printf("\n%-10s", rownames[row]);
		for (col = 0;col<count;col++) {
			//printf("%d ", vertex[i].color);
			printf("%-10d", vertex[col].d);
			//printf("\n");
		}
		printf("\n");
	}


	printf("\n2. Bellman_Ford's Algorithm\n\n");

	init_node(vertex, count, matrix, colnames);
	printf("initiation nodes\n");

	//init_node(vertex, count, matrix);
	start = clock();
	for (i = 0; i < count; i++) {
		Bellman_Ford(vertex, count, i);
	}
	start = clock() - start;
	duration = ((double)start) / CLOCKS_PER_SEC;
	printf("It take %f time seconds to compute shortest path between cities with Bellman_Ford algorithm as follows. \n", duration);

	init_node(vertex, count, matrix, colnames);
	printf("          ");
	//init_node(vertex,count);
	for (col = 0; col < count; col++) {
		printf("%-10s", colnames[col]);
	}
	for (row = 0; row < count; row++) {
		Bellman_Ford(vertex, count, row);
		printf("\n%-10s", rownames[row]);
		for (col = 0;col<count;col++) {
			//printf("%d ", vertex[i].color);
			printf("%-10d", vertex[col].d);
			//printf("\n");
		}
		printf("\n");
	}


	printf("\n3. Floyd's Algorithm\n\n");
	init_node(vertex, count, matrix,colnames);
	printf("Success! nodes initialization \n");

	start = clock();

	Floyd_Warshall(vertex, count);

	start = clock() - start;
	duration = ((double)start) / CLOCKS_PER_SEC;
	printf("It take %f time seconds to compute shortest path between cities with Folyd algorithm as follows. \n", duration);

	//printf("SSP result");
	printf("          ");
	for (col = 0; col < count; col++) {
		printf("%-10s", colnames[col]);
	}
	for (row = 0; row < count; row++) {
		printf("\n%-10s", rownames[row]);
		for (col = 0;col<count;col++) {
			//printf("%d ", vertex[i].color);
			printf("%-10d", d[row][col]);
			//printf("\n");
		}
		printf("\n");
	}
	printf("\n");

	free(vertex);
	for (i = 0; i < SIZE; i++) {
		free(*(matrix + i));
	}
	free(matrix);

	for (i = 0; i < SIZE; i++) {
		free(*(rownames + i));
	}
	free(rownames);

	for (i = 0; i < SIZE; i++) {
		free(*(colnames + i));
	}
	free(colnames);
	system("pause");
	return 0;
}
