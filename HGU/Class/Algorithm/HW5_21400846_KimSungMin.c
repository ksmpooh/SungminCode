#include <stdio.h>
#include <string.h>
#include <stdbool.h>

#define SIZE 100
#define IFN 999999
typedef struct nodes node;
int count;
char matrix[SIZE][SIZE][SIZE]; //temp[row][col][number_of_string]


typedef struct nodes{
    char name[20];
    int index;
    int adj[20];
    int weight;
    int d;
    char predecessor[20];
    bool Set_or_not;
}node;

void init_node(node *vertex){
  int row, col, index;
  index = 0;
  for(col = 1; col <= count ; col++){
    strcpy(vertex[index].name,matrix[0][col]);
    //printf("%s ",vertex[index].name);
    vertex[index].index = col;
  }
  for(index = 0,row = 1 ; row <= count ; row++,index++){
    for(col = 1; col <= count ; col++){
      if(row == col){
        vertex[index].adj[(col-1)] = 0;}
      else if(strcmp(matrix[row][col],"IFN") != true){
        vertex[index].adj[(col-1)] = atoi(matrix[row][col]);
        //printf("\n%d",vertex[index].adj[col-1]);
      }else
        vertex[index].adj[(col-1)] = IFN;
    }
  }
}

void init_Set_or_not(node *vertex){
  int i;
  for(i = 0; i < count; i++){
    vertex[i].Set_or_not = false;
  }
}

void init_single_source(node *vertex, int s){
  int i;
  for(i = 0; i < count; i++){
    strcpy(vertex[i].predecessor,NULL);
    vertex[i].d = IFN;  //무한을 IFN으로
  }
  vertex[s].d = 0;
}

void Relax(node *vertex, int row, int col, int w){
  int temp;
  temp = vertex[row].d + w;
  if(vertex[col].d > temp){
    vertex[col].d = temp;
    //vertex[col]->predecessor =
    strcpy(vertex[col].predecessor,vertex[row].predecessor);
  }
}

void Dijkstra(node *vertex, int w, int s){
  init_single_source(vertex, s);
  init_Set_or_not(vertex);
  int i;
  int remain = count;
  int Q[count];
  for(i = 0 ; i < count ; i ++){
    //Q[i] = vertex.d
  }
  while (remain != 0) {
    
  }

}
void copy_string(char from[],int row, int col){
  int c = 0;

  while(from[c] != '\0'){
      matrix[row][col][c] = from[c];
      c ++;
  }
  matrix[row][col][c] = '\0';
}

void split_string(char *line, int row) {
    const char delimiter[] = "\t";
    int col = 0;
    char *tmp;
    tmp = strtok(line, delimiter);
    if (tmp == NULL)
      return;
    copy_string(tmp,row,col);
    printf("%d:%d : %s\n", row,col,matrix[row][col]);

    col++;
    //printf("fs %s\n", tmp);

    for (col = 1; col <= count ; col ++) {
        tmp = strtok(NULL, delimiter);
        if (tmp == NULL)
            break;
        //printf("%d %s\n",i, tmp);
        copy_string(tmp,row,col);
        printf("%d:%d : %s\n", row,col,matrix[row][col]);

    }
}

int split_string_first(char *line, int row) {
    //rintf("result\n");
    const char delimiter[] = " \t";
    int j , k = 0;
    char *tmp;
    tmp = strtok(line, delimiter);
    count++;
    //printf("count:%d, %s\n", count,tmp);
    copy_string(tmp,row,count);
    //printf("0:1 : %s\n", matrix[0][count]);
    for (j = 1; j <20; j ++) {
        tmp = strtok(NULL, delimiter);
        if (tmp == NULL){
          break;
        }
        printf("%d %s\n", j, tmp);
        count++;
        copy_string(tmp,row,count);
        //printf("0:%d : %s\n", count,matrix[0][count]);

    }
    return count;
}


int main(int argc, char const *argv[]) {
  int i;
  count = 0;
  size_t data_size;
  FILE *fin;
  fin = fopen("D:/git/SungminCode/HGU/Class/Algorithm/hw5..data","r");
  char *line = NULL;
  if(fin)
    printf("File is open!\n");
  else{
    printf("File is Error\n");
    return 0;
  }
  for(i = 0 ; getline(&line, &data_size, fin) != -1 ; i++ ){
    if(i == 0){
      count = split_string_first(line,i);
      printf("count : %d\n",count);
    }else{
      split_string(line,i);
    }
  }
  fclose(fin);
/********************************************************/
  node vertex[count];
  init_node(vertex);

  return 0;
}
