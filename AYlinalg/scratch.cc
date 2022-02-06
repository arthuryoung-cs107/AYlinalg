#include <cstdio>
#include <cstring>
#include <cstdlib>

extern "C"
{
  #include "AYaux.h"
}

void clear_string(char *ptr_, size_t size_)
{
  // memset(ptr_, '$', (size_t)(size_-1));
  strcpy(ptr_, "");
}

int main()
{
  char src[50], dest1[50], dest2[50], dest3[50];
  printf("dest3 size-1: %d \n", (size_t)(sizeof(dest3)-1));
  clear_string(dest3, sizeof(dest3));
  printf("dest3 size-1: %d \n", (size_t)(sizeof(dest3)-1));

  strcpy(src,  "This is source");
  strcpy(dest1, "This is destination1");
  strcpy(dest2, "1234");
  strcpy(dest3, "1234");

  strcat(dest1, src);

  printf("dest1: %s\n", dest1);

  strcpy(dest1, "now what");

  printf("dest1: %s\n", dest1);

  printf("dest2 strlen: %d, dest2 sizeof: %d\n", strlen(dest2), sizeof(dest2) );
  printf("dest2: %s\n", dest2);
  printf("dest2, all: ");
  for (int i = 0; i < 50; i++)
  {
    printf("%c ", dest2[i]);
  }
  printf("\n");
  printf("dest3, all: ");
  for (int i = 0; i < 50; i++)
  {
    printf("%c ", dest3[i]);
  }
  printf("\n");

  char * dest3_copy = string_gen_pruned(dest3);
  printf("dest3 copy strlen, sizeof: %d, %d\n", strlen(dest3_copy), sizeof(dest3_copy));
  printf("%s\n", dest3_copy);

  free(dest3_copy);

  return 0;
}
