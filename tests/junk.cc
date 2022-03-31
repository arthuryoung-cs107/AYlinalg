#include <cstdio>
#include <cstring>
#include <unistd.h>
#include "omp.h"
#define N 1000

int check[N];
int good_members[N];
int pos = 0;

void init_members()
{
  for (int i = 0; i < N; i++)
  {
    check[i] = i+1;
    good_members[i] = 0;
  }

}

int is_good(int i_)
{
  // sleep(0.2);
  return (check[i_] % 2) ? 0 : 1; //returns 1 if even, returns 0 if false
}

void find_good_members()
{
  #pragma omp parallel for
  for (int i = 0; i < N; ++i)
  {
    if (is_good(i))
    {
      good_members[pos] = i;
      #pragma omp atomic
      pos++;

    }
  }
}

void find_good_members_fixed()
{
  #pragma omp parallel for
  for (int i = 0; i < N; ++i)
  {
    good_members[i] = is_good(i);
    #pragma omp atomic
    pos++;
  }
}

int main()
{
  init_members();
  // find_good_members();
  find_good_members_fixed();

  // for (int i = 0; i < 10; i++) printf("%d ", check[i]);
  // printf("\n");
  // for (int i = 0; i < pos; i++) printf("%d ", good_members[i]);

  printf("\n%d\n", pos);
  int sum = 0;
  for (int i = 0; i < N; i++) sum+= good_members[i];
  printf("%d\n", sum );

  return 0;
}
