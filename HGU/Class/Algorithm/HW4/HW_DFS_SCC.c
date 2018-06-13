#define _GNU_SOURCE
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>

enum id{white=1, gray, black};
int time;
int scc[20][20];
int scc_group,scc_index;
typedef struct nodes node;
void DFS_visit(node *vertex,int const count, int index);

typedef struct nodes{
    int d_time;
    int f_time;
    enum id color;
    char name;
    int adj[20];
    //node *next;
}node;
void ini_node(node *vertex){
  //strcpy(vertex->color,"white");
  vertex->color = white;
  vertex->d_time = 0;
  vertex->f_time = 0;
  //vertex->name = NULL;
  //vertex->next = NULL;
}

void DFS(node *vertex,int const count){
  int i,j;

  time = 0;
  for(i = 0; i < count; i ++){
    if(vertex[i].color == white){
      //printf("yes!\n");
      DFS_visit(vertex,count,i);
    }
  }
  for(i = 0;i < count; i ++){
    printf("vertex %c :  ",vertex[i].name);
    //printf("vertex[%d].color =  %d  \n",i,vertex[i].color);
    printf("Discovery_time: %d ",vertex[i].d_time);
    printf("Finish_time: %d \n",vertex[i].f_time);
  }
}
void DFS_visit(node *vertex, int const count,int index){
  int i;
  vertex[index].color = gray;
  //printf("123vertex[%d].color =  %d \n ",index,vertex[index].color);
  time++;
  vertex[index].d_time = time;
  for(i = 0 ; i < count; i++){
    if(vertex[index].adj[i] == 1){
      if(vertex[i].color == white){
        DFS_visit(vertex,count,i);
      }
    }
  }
  vertex[index].color = black;
  time++;
  vertex[index].f_time = time;
}
void DFS_visit2(node *vertex, int const count,int index){
  int i;
  vertex[index].color = gray;
  //printf("123vertex[%d].color =  %d \n ",index,vertex[index].color);
  time++;
  vertex[index].d_time = time;
  for(i = 0 ; i < count; i++){
    if(vertex[index].adj[i] == 1){
      if(vertex[i].color == white){
        scc_index++;
        scc[scc_group][scc_index] = i;
        DFS_visit2(vertex,count,i);
      }
    }
  }
  vertex[index].color = black;
  time++;
  vertex[index].f_time = time;
  //scc_group++;
  //scc_index = 0;
  //scc[i][scc_index+1] = -1;
}
void tran_DFS(node *vertex,node *tran_vertex,int arr[], int const count){
  int i,j;
  time = 0;
  scc_group = 0;
  scc_index = 0;
  for(i = 0; i<count; i++){
    for(j = 0 ; j < count; j++){
      if(vertex[j].f_time == arr[i]){
        if(tran_vertex[j].color == white){
          scc[scc_group][scc_index] = j;
          DFS_visit2(tran_vertex,count,j);
          scc[i][scc_index+1] = -1;
          scc_group++;
          //scc_index = 0 ;
        }
        //scc_index = 0;
      }
      //scc_index = 0;
    }
    //scc_group++;
    //scc_index = 0;
  }
  for(i = 0;i < count; i ++){
    printf("vertex %c :  ",tran_vertex[i].name);
    //printf("vertex[%d].color =  %d  \n",i,vertex[i].color);
    printf("Discovery_time: %d ",tran_vertex[i].d_time);
    printf("Finish_time: %d \n",tran_vertex[i].f_time);
  }
}

int main(int argc, char **argv){
  int i ,j;
  int count=0;
  char c;
  FILE *fin;
  fin = fopen(argv[1],"r");
  //argc --> arguments 갯수 2
  //argv[0] --> 실행된 경로
  //argv[1]~ --> 입력한 string.(=arguments)
  //gcc HW.c -o sample.exe
  //sample.exe "input.csv"

  if(fin)
    printf("File is open!\n");
  else{
    printf("File is Error\n");
    return 0;
  }
  i = 0;
  j = 0;

  //char arr[100];
  char buff[100],colnames[100];

  fgets(buff,100,fin);
  while(buff[i] != '\0'){
      if(!(buff[i] == '\t')){
        colnames[j] = buff[i];
        j ++;
      }
      i++;
  }
  colnames[j] ='\0';
  //colnames[0] ='\0';
  //printf("After : %c\n",colnames[1]);
  //printf("Number : %d\n",j-1);
  count = j-1;

  int matrix[count][count];

  char temp[100];
  while(!feof(fin)){
    for(i = 0 ; i < count; i ++){
      for(j = 0 ; j <= count; j ++){
        fscanf(fin,"%s ",&temp);
        if(j > 0)
          matrix[i][j-1] = atoi(temp);
      }
    }
  }
  printf("\n");
  fclose(fin);

  node vertex[count];

  for(i = 0 ; i < count ; i ++){
    ini_node(&vertex[i]);
    vertex[i].name = colnames[i];
    //printf("v[%d] : %d\n",i ,vertex[i].color);
    //printf("vertex[%d] : %c\n",i,vertex[i].name);
  }

  printf("1) Array of adjacency list of above graph\n");

  for(i = 0; i < count; i ++){
    printf(" %c -->",colnames[i]);
    for(j =0 ; j< count; j ++){
      if(matrix[i][j] == 1){
        printf(" %c -->",colnames[j]);
        vertex[i].adj[j] = 1;}
      else
        vertex[i].adj[j] = 0;

    }
    printf("NULL\n");
    vertex[i].adj[j+1] = -1;
  }
  printf("\n");
  printf("2) Discovery time and finish time of each vertex after step 1\n");
  DFS(vertex,count);

  //vertex[]
  //*vertex

  printf("\n");
  printf("3) Array of adjacency list of transpose graph after step 2\n");

  int tran_matrix[count][count];
  node tran_vertex[count];
  for(i = 0 ; i < count ; i ++){
    ini_node(&tran_vertex[i]);
    tran_vertex[i].name = colnames[i];
    //printf("v[%d] : %d\n",i ,vertex[i].color);
    //printf("vertex[%d] : %c\n",i,vertex[i].name);
  }

  for(i =0 ; i < count; i ++){
    for(j = 0 ; j < count ; j ++){
      if(j == i)
        tran_matrix[i][j] = 0;
      if(matrix[i][j] == 0)
        tran_matrix[j][i] = 0;
      else
        tran_matrix[j][i] = 1;
    }
  }

  for(i = 0; i < count; i ++){
    printf(" %c -->",colnames[i]);
    for(j =0 ; j< count; j ++){
      if(tran_matrix[i][j] == 1){
        printf(" %c -->",colnames[j]);
        tran_vertex[i].adj[j] = 1;}
      else
        tran_vertex[i].adj[j] = 0;
    }
    printf("NULL\n");
    tran_vertex[i].adj[j+1] = -1;
  }

  printf("\n4) Discovery time and finish time of each vertex after step3\n");

  int order[count],tmp;
  for(i = 0 ; i < count; i ++){
    order[i] = vertex[i].f_time;
  }

  for(i = 0 ; i < sizeof(order)/sizeof(int) ; i++){
    for(j = 0 ; j < sizeof(order)/sizeof(int)-1 ; j++){
      if(order[j] < order[j+1]){
        tmp = order[j+1];
        order[j+1] = order[j];
        order[j] = tmp;
      }
    }
  }
  printf("\n");

  for(i = 0 ; i<20;i++){
    for(j = 0 ; j <20;j++){
      scc[i][j] = -1;
    }
  }
  /*for(i = 0; i<count; i ++){
    printf(" %d",order[i]);
  }
  printf("\n");*/
  //int SCC[count][count];
  tran_DFS(vertex,tran_vertex,order,count);

  printf("\n5) SCC result\n");
  /*printf("scc_group : %d\n",scc_group);
  printf("scc_index : %d\n",scc_index);
  printf("%d\n",scc[4][4]);*/

  /*for(i = 0; i<scc_group; i ++){
    for(j = 0; j <20; j++)
      printf(" %d ",scc[i][j]);
    printf("\n");
  }*/
  printf("\n");
  j = 0;
  for(i = 0; i < scc_group ; i ++){
    printf(" SCC %d : vertex ",(i+1));
    for(; scc[i][j] != -1; j++){
      printf(" %c", tran_vertex[scc[i][j]].name);
    }
    j--;
    printf("\n");
  }

  return 0;
}
