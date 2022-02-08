#include <cstdio>
#include <cstring>

#include "AYlinalg.hh"

#include <iostream>
#include <bitset>

extern "C"
{
  #include "AYaux.h"
}

void clear_string(char *ptr_, size_t size_)
{
  // memset(ptr_, '$', (size_t)(size_-1));
  strcpy(ptr_, "");
}

void string_tests()
{
  char src[50], dest1[50], dest2[50], dest3[50];
  printf("dest3 size-1: %ld \n", (size_t)(sizeof(dest3)-1));
  clear_string(dest3, sizeof(dest3));
  printf("dest3 size-1: %ld \n", (size_t)(sizeof(dest3)-1));

  strcpy(src,  "This is source");
  strcpy(dest1, "This is destination1");
  strcpy(dest2, "1234");
  strcpy(dest3, "1234");

  strcat(dest1, src);

  printf("dest1: %s\n", dest1);

  strcpy(dest1, "now what");

  printf("dest1: %s\n", dest1);

  printf("dest2 strlen: %ld, dest2 sizeof: %ld\n", strlen(dest2), sizeof(dest2) );
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
  printf("dest3 copy strlen, sizeof: %ld, %ld\n", strlen(dest3_copy), sizeof(dest3_copy));
  printf("%s\n", dest3_copy);

  delete dest3_copy;
}

void pbit(int in_)
{
  std::cout << std::bitset<8>(in_) << " ";
}

int main()
{
  // unsigned int fflags = 255;
  // unsigned int fflags = 0;
  // unsigned int fflags = 3;
    // unsigned int fflags = 1;
    // unsigned int fflags = 2;
  unsigned int fflags = 12;
    // unsigned int fflags = 4;
  //   unsigned int fflags = 8;
  // unsigned int fflags = 48;
  //   unsigned int fflags = 16;
  //   unsigned int fflags = 32;
  // unsigned int fflags = 64;

  pbit(0); printf("0: %d ", (fflags==0)); if (fflags==0) printf("PASSED\n"); else printf("\n");
  pbit(3); printf("3: %d ", (fflags&3)); if (fflags&3) {printf("PASSED\n");
    pbit(1); printf("  1: %d ", (fflags&1)); if (fflags&1) printf("PASSED\n"); else printf("\n");
    pbit(2); printf("  2: %d ", (fflags&2)); if (fflags&2) printf("PASSED\n"); else printf("\n");
  } else printf("\n");
  pbit(12); printf("12: %d ", (fflags&12)); if (fflags&12) {printf("PASSED\n");
    pbit(4); printf("  4: %d ", (fflags&4)); if (fflags&4) printf("PASSED\n"); else printf("\n");
    pbit(8); printf("  8: %d ", (fflags&8)); if (fflags&8) printf("PASSED\n"); else printf("\n");
  } else printf("\n");
  pbit(48); printf("48: %d ", (fflags&48)); if (fflags&48) {printf("PASSED\n");
    pbit(16); printf("  16: %d ", (fflags&16)); if (fflags&16) printf("PASSED\n"); else printf("\n");
    pbit(32); printf("  32: %d ", (fflags&32)); if (fflags&32) printf("PASSED\n"); else printf("\n");
  } else printf("\n");
  pbit(64); printf("64: %d ", (fflags&64)); if (fflags&64) printf("PASSED\n"); else printf("\n");


  pbit(fflags); printf("fflag: %d\n", fflags);

  // if(fflags&1)
  // if(fflags&2)
  // if(fflags&4)
  // if(fflags&8)
  // if(fflags&16)
  // if(fflags&32)
  return 0;
}
