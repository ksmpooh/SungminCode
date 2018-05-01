#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>

int count(FILE *const fin){
  int c, col;
  col = 0;
  while(1){
    c = fgetc(fin);
    if(c == EOF || c == '\n')
      break;
    col++;
  }
  return col/2;
}

int main(){

  int i,j;
  char mat[100],lian[100];
  int count;
  char c;
  FILE *fin;
  fin = fopen("D:/git/SungminCode/HGU/Class/Algorithm/hw4.data","r");

  /*if(fin)
    printf("File is open!\n");

  printf("count: %d\n",count(fin));
  int variables = count(fin);*/

  //fgets(mat,100,fin);
  //puts(mat);

  fscanf(fin," %c ",&lian);
  fscanf(fin,"%c ",&lian);
  fscanf(fin,"%c ",&lian);
  printf("%c\n",lian[0]);
  printf("%c\n",lian[1]);
  printf("%c\n",lian[2]);

  for(i = 0; i < 100 ; i++){
    if(i == 0){
      fgets(mat,100,fin);
      puts(mat);
    }else{
      for(j =0; j <100;j++){
        if(j == 0){
        }
      }
    }
  }
  
  /*char** matrix = malloc(sizeof(char)*variables);
  for(i = 0 ; i <variables; i++)
    matrix[i] = malloc(variables*sizeof(char));

  for(i = 0 ; i < variables; i ++){
    for(j = 0; j < variables; j++){
      if(!fscanf(fin," %c",&matrix[i][j]))
        break;
    }
  }
  for(i = 0 ; i < variables; i ++){
    for(j = 0; j < variables; j++){
      printf(" %c ",matrix[i][j]);
    }
    printf("\n");
  }
  free(matrix);*/

  fclose(fin);
  return 0;
}
