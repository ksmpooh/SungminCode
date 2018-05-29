#include <stdio.h>
#include <string.h>

#define SIZE 100
typedef struct nodes node;
int count;
char matrix[SIZE][SIZE][SIZE]; //temp[row][col][number_of_string]


typedef struct nodes{
    char name[20];
    int adj[20];
    int weight;
}node;

void init_single_source(){

}

void Relax(){

}

void Dijkstra(G,w,s){

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



  return 0;
}
