#include <stdio.h>

int main(void) {
  const double epsilon = 1e-16;
  const double x = 1.0;
  const double y = 1.0 + epsilon;
  printf("epsilon       = %le\n", epsilon);
  printf("(1+epsilon)-1 = %le\n", (y - x));
}
