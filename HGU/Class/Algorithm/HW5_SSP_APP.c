#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>

#define SIZE 20
#define IFN 999999
typedef struct nodes node;
enum id {white = 0, gray, black};
char matrix[SIZE][SIZE][10]; //temp[row][col][number_of_string]

int predecessor_ap[SIZE][SIZE];
int d[SIZE][SIZE];

typedef struct nodes{
    char name[10];
    int index;
    int adj[20];
    //int weight;
    //int edge;
    int d;
    char predecessor[20];
    enum id color;
}node;

void init_node(node *vertex,int count){
  int row, col, index, edge;
  index = 0;
  //edge = 0;
  for(col = 1; col <= count ; col++,index++){
    strcpy(vertex[index].name,matrix[0][col]);
    //printf("%s ",vertex[index].name);
    vertex[index].index = col;
  }
  for(index = 0,row = 1 ; row <= count ; row++,index++){
    strcpy(vertex[index].predecessor,"\0");
    vertex[index].d = 0;
    for(col = 1; col <= count ; col++){
      if(row == col){
        vertex[index].adj[(col-1)] = 0;}
      else if(strcmp(matrix[row][col],"IFN") != true){
        vertex[index].adj[(col-1)] = atoi(matrix[row][col]);
        //edge++;
        //printf("\n%d",vertex[index].adj[col-1]);
      }else
        vertex[index].adj[(col-1)] = IFN;
    }

    //vertex[index].edge = edge;
  }
}

void init_single_source(node *vertex, int s,int count){
  int i;
  for(i = 0; i < count ; i++){
    strcpy(vertex[i].predecessor,"\0");
    vertex[i].d = IFN;  //무한을 IFN으로
    vertex[i].color = white;
  }
  //printf("good\n");
  vertex[s].d = 0;
  vertex[s].color = black;
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
int Extract_min(node *vertex, int count){
  int min_edge = IFN;
  int index;
  int i;
  for(i = 0 ; i < count; i ++){
    if(vertex[i].color == gray && vertex[i].d < min_edge){
        min_edge = vertex[i].d;
        index = i;
    }
  }
  return index;
}
void search_and_change_color(node *vertex, int index, int count){
  int i;
  for(i = 0 ; i < count ; i ++){
    if(vertex[i].color == black)
      continue;
    else if(vertex[index].adj[i] != 0 && vertex[index].adj[i] != IFN){
      vertex[i].color = gray;
    }

  }
}
void Dijkstra(node *vertex, int count,int s){
  //printf("good\n");
  init_single_source(vertex, s, count);
  //printf("good\n");
  //init_color(vertex,count);
  int i;
  int u;
  int remain = count;
  for(i = 0 ; i < count ; i ++){
    if(vertex[s].adj[i] != IFN && vertex[s].adj[i] != 0){
      vertex[i].color = gray;
      //Q[i] = vertex[s].adj[i];
      vertex[i].d = vertex[s].adj[i];
      //printf("%d\n", vertex[i].d);
    }
  }
  remain--;
  while (remain != 0) {
    u = Extract_min(vertex,count);
    vertex[u].color = black;
    search_and_change_color(vertex,u,count);
    for(i = 0; i < count; i++){
      if(vertex[u].adj[i] != IFN && vertex[u].adj[i] != 0 ){
        Relax(vertex,u,i,vertex[u].adj[i]);
        //vertex[i].color = gray;
      }
    }
    remain--;
  }
}

void Bellman_Ford(node *vertex,int count, int s ){
  int i,j,k;
  k = 0;
  int remain  = count -1;
  int seq[remain];
  init_single_source(vertex, s, count);

  for(i = 0; i < count ; i ++){
    if(vertex[s].adj[i] != IFN && vertex[s].adj[i] != 0){
      //printf("%d ",vertex[s].adj[i]);
      vertex[i].d = vertex[s].adj[i];
      seq[k] = i;
      //printf("seq[%d] : %d %s\n",k,seq[k],vertex[i].name);
      k++;
      vertex[s].color = gray;


    }
  }
  int index = 0;

  for (i = 0 ; i < remain ; i++) {
    //printf("seq[%d] : %d \n",k,seq[k]);
      index = seq[i];
      for(j = 0 ; j < count; j ++){
        if(vertex[index].adj[j] != IFN && vertex[index].adj[j] &&  vertex[j].color != gray && vertex[j].d != 0 ){
            vertex[j].color = gray;
            seq[k] = j;
            k ++ ;
        }
        if(vertex[index].adj[j] != IFN && vertex[index].adj[j]){
          Relax(vertex,index,j,vertex[index].adj[j]);
        }
      }
      //remain--;
  }

  /*for(i = 0 ; i < count; i ++){
    if(i == s){
        continue;
    }else{
      for(j = 0; j < count; j ++){
        if(vertex[i].adj[j] != IFN && vertex[i].adj[j] != 0)
          Relax(vertex,i,j,vertex[i].adj[j]);
      }
    }
  }*/
}
void Floyd_Warshall(node *vertex, int count){
  int row,col,index;

  for(row = 0; row < count; row ++){
    for(col = 0; col < count; col++){
      d[row][col] = vertex[row].adj[col];
      predecessor_ap[row][col] = -1;
    }
  }
  for(index = 0; index < count; index++){
    for(row = 0; row < count; row++){
      for(col = 0 ; col <count; col++){
        if(d[row][col] > d[row][index] + d[index][col]){
          d[row][col] =  d[row][index] + d[index][col];
          if(predecessor_ap[index][col] == -1){
            predecessor_ap[row][col] = index;}
          else{
            predecessor_ap[row][col] = predecessor_ap[index][col];}
          }
      }
    }
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

void split_string(char *line, int row,int count) {
    const char delimiter[] = "\t";
    int col = 0;
    char *tmp;
    tmp = strtok(line, delimiter);
    if (tmp == NULL)
      return;
    copy_string(tmp,row,col);
    //printf("%d:%d : %s\n", row,col,matrix[row][col]);

    col++;
    //printf("fs %s\n", tmp);

    for (col = 1; col <= count ; col ++) {
        tmp = strtok(NULL, delimiter);
        if (tmp == NULL)
            break;
        //printf("%d %s\n",i, tmp);
        copy_string(tmp,row,col);
        //printf("%d:%d : %s\n", row,col,matrix[row][col]);

    }
}

int split_string_first(char *line, int row) {
    //rintf("result\n");
    const char delimiter[] = " \t";
    int j , count = 0;
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
        //printf("%d %s\n", j, tmp);
        count++;
        copy_string(tmp,row,count);
        //printf("0:%d : %s\n", count,matrix[0][count]);

    }
    return count;
}

int main(int argc, char const *argv[]) {
  int i, j ;

  int count = 0;
  size_t data_size;
  FILE *fin;
  fin = fopen("D:/git/SungminCode/HGU/Class/Algorithm/hw5..data","r");
  //fin = fopen(argv[1],"r");
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
      printf("number of cities : %d\n",count);
    }else{
      split_string(line,i,count);
    }
  }
  fclose(fin);
/********************************************************/
  node vertex[count];

  init_node(vertex,count);
  printf("\nOutput\n");


/////1. Dijkstra
  printf("1. Dijkstra's Algorithm\n");
  clock_t start = clock();

  for(i = 0 ; i < count; i ++){
    Dijkstra(vertex,count,j);}
  start = clock() - start;
  double duration = ((double)start)/CLOCKS_PER_SEC;
  //double duration = ((double)start);
  //double duration =(double)(clock() - start);

  printf("It take %f time seconds to compute shortest path between cities with Dijkstra's algorithm as follows. \n", duration );
  printf("          ");
  for(i = 0; i < count; i ++){
    printf("%-10s",vertex[i].name);}
  for(j = 0 ; j <  count; j ++){
    Dijkstra(vertex,count,j);
    printf("\n%-10s",vertex[j].name );
    for(i=0;i<count;i++){
      //printf("%d ", vertex[i].color);
      printf("%-10d", vertex[i].d );
      //printf("\n");
    }
    printf("\n");
  }

/////2. Bellman_Ford's Algorithm

  printf("\n2. Bellman_Ford's Algorithm\n");

  //init_node(vertex,count);
  start = clock();
  for(i = 0 ; i < count; i ++){
    Bellman_Ford(vertex,count,i);}
  start = clock() - start;
  duration = ((double)start)/CLOCKS_PER_SEC;
  printf("It take %f time seconds to compute shortest path between cities with Bellman_Ford algorithm as follows. \n", duration );

  init_node(vertex,count);
  printf("          ");
  //init_node(vertex,count);
  for(i = 0; i < count; i ++){
    printf("%-10s", vertex[i].name);}
  for(j = 0 ; j <  count; j ++){
    Bellman_Ford(vertex,count,j);
    printf("\n%-10s",vertex[j].name );
    for(i=0;i<count;i++){
      //printf("%d ", vertex[i].color);
      printf("%-10d", vertex[i].d );
      //printf("\n");
    }
    printf("\n");
  }

  //3. Floyd's Algorithm
  printf("\n3. Floyd's Algorithm\n");
  //init_node(vertex,count);
  start = clock();

  Floyd_Warshall(vertex, count);

  start = clock() - start;
  duration = ((double)start)/CLOCKS_PER_SEC;
  printf("It take %f time seconds to compute shortest path between cities with Folyd algorithm as follows. \n", duration );

  //printf("SSP result");
  printf("          ");
  for(i = 0; i < count; i ++){
    printf("%-10s", vertex[i].name);}
  for(j = 0 ; j <  count; j ++){
    printf("\n%-10s",vertex[j].name );
    for(i=0;i<count;i++){
      //printf("%d ", vertex[i].color);
      printf("%-10d", d[i][j] );
      //printf("\n");
    }
    printf("\n");
  }
  printf("\n");


  return 0;
}
