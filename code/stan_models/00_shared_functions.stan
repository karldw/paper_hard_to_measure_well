
// shared_functions.stan
// Meant to be included inside a functions {} block

// Which elements of x are != 0?
int[] which (int[] x) {
  int N = size(x);
  int is_nonzero[N];
  for (i in 1:N) {
    is_nonzero[i] = x[i] != 0;
  }
  int result[sum(is_nonzero)];
  int j = 0;
  for (i in 1:N) {
    if (is_nonzero[i]) {
      j += 1;
      result[j] = i;
    }
  }
  return result;
}
