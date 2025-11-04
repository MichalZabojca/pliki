#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

long long q = 3;

typedef struct zbior_ary {
  unsigned rozmiar;
  unsigned uzyte;
  long long *p;
  long long *k;
} zbior_ary;

long long min(long long a, long long b) {
  if (a <= b)
    return a;
  else
    return b;
}

long long max(long long a, long long b) {
  if (a >= b)
    return a;
  else
    return b;
}

long long modulo(long long x) {
  long long r = x % q;
  if (r < 0)
    r += q;
  return r;
}

zbior_ary nowy_zbior(long long rozmiar) {
  long long *p = (long long *)malloc(sizeof(long long) * (size_t)rozmiar);
  if (p == NULL) {
    exit(1);
  }

  long long *k = (long long *)malloc(sizeof(long long) * (size_t)rozmiar);
  if (k == NULL) {
    exit(1);
  }

  zbior_ary A;

  A.rozmiar = rozmiar;
  A.uzyte = 0;
  A.p = p;
  A.k = k;
  return A;
}

void zbior_push(zbior_ary *A, long long nowy_p, long long nowy_k) {
  if (A->uzyte == A->rozmiar) {
    long long *p_r = (long long *)realloc(A->p, sizeof(long long) *
                                                    (size_t)((A->rozmiar) * 2));
    if (p_r == NULL) {
      exit(1);
    }

    long long *k_r = (long long *)realloc(A->k, sizeof(long long) *
                                                    (size_t)((A->rozmiar) * 2));
    if (k_r == NULL) {
      exit(1);
    }

    A->p = p_r;
    A->k = k_r;
    A->rozmiar *= 2;
  }
  *(A->p + A->uzyte) = nowy_p;
  *(A->k + A->uzyte) = nowy_k;
  A->uzyte += 1;
}

zbior_ary ciag_arytmetyczny(long long a, long long Q, long long b) {
  q = Q;
  zbior_ary A = nowy_zbior(1);
  zbior_push(&A, a, b);
  return A;
}

zbior_ary singleton(long long a) {
  zbior_ary A = nowy_zbior(1);
  zbior_push(&A, a, a);
  return A;
}

zbior_ary suma(zbior_ary A, zbior_ary B) {
  long long a = 0, b = 0, nowy_p = 0, nowy_k = 0;
  zbior_ary C = nowy_zbior(A.rozmiar);

  if (A.uzyte == 0)
    return B;
  if (B.uzyte == 0)
    return A;

  if (A.p[a] <= B.p[b]) {
    nowy_p = A.p[a];
    nowy_k = A.k[a];
  }

  else {
    nowy_p = B.p[b];
    nowy_k = B.k[b];
  }

  while (b < B.uzyte && a < A.uzyte) {

    if (modulo(A.p[a]) < modulo(B.p[b])) {
      zbior_push(&C, A.p[a], A.k[a]);
      a++;
      continue;
    } else if (modulo(B.p[b]) < modulo(A.p[a])) {
      zbior_push(&C, B.p[b], B.k[b]);
      b++;
      continue;
    }

    if (A.p[a] <= B.p[b]) {
      nowy_p = A.p[a];
      nowy_k = A.k[a];
      a++;
    }

    else {
      nowy_p = B.p[b];
      nowy_k = B.k[b];
      b++;
    }

    while ((a < A.uzyte && modulo(A.p[a]) == modulo(nowy_p) &&
            A.p[a] <= nowy_k + q) ||
           (b < B.uzyte && modulo(B.p[b]) == modulo(nowy_p) &&
            B.p[b] <= nowy_k + q)) {
      while (a < A.uzyte && modulo(A.p[a]) == modulo(nowy_p) &&
             A.p[a] <= nowy_k + q) {
        nowy_k = max(A.k[a], nowy_k);
        a++;
      }

      while (b < B.uzyte && modulo(B.p[b]) == modulo(nowy_p) &&
             B.p[b] <= nowy_k + q) {
        nowy_k = max(B.k[b], nowy_k);
        b++;
      }
    }
    zbior_push(&C, nowy_p, nowy_k);
  }

  while (a < A.uzyte) {
    zbior_push(&C, A.p[a], A.k[a]);
    a++;
  }

  while (b < B.uzyte) {
    zbior_push(&C, B.p[b], B.k[b]);
    b++;
  }

  return C;
}

zbior_ary iloczyn(zbior_ary A, zbior_ary B) {
  long long a = 0, b = 0, nowy_p = 0, nowy_k = 0;
  zbior_ary C = nowy_zbior(A.rozmiar);
  if (A.uzyte == 0 || B.uzyte == 0)
    return C;

  while (a < A.uzyte && b < B.uzyte) {

    if (modulo(A.p[a]) < modulo(B.p[b])) {
      a++;
      continue;
    }

    if (modulo(A.p[a]) > modulo(B.p[b])) {
      b++;
      continue;
    }
    nowy_p = max(A.p[a], B.p[b]);
    nowy_k = min(A.k[a], B.k[b]);
    if (nowy_k >= nowy_p)
      zbior_push(&C, nowy_p, nowy_k);

    if (A.k[a] < B.k[b])
      a++;
    else
      b++;
  }

  return C;
}

zbior_ary roznica(zbior_ary A, zbior_ary B) {
  long long a = 0, b = 0, nowy_p = 0, nowy_k = 0, aktualny = 0;
  zbior_ary C = nowy_zbior(A.rozmiar);
  if (A.uzyte == 0)
    return C;
  if (B.uzyte == 0)
    return A;

  while (a < A.uzyte) {
    nowy_p = A.p[a];
    nowy_k = A.k[a];
    aktualny = nowy_p;

    while (b < B.uzyte && modulo(B.p[b]) < modulo(nowy_p))
      b++;

    if (b >= B.uzyte || modulo(B.p[b]) > modulo(nowy_p)) {
      zbior_push(&C, nowy_p, nowy_k);
      a++;
      continue;
    }

    while (b < B.uzyte && modulo(B.p[b]) == modulo(nowy_p)) {
      if (B.k[b] < nowy_p) {
        b++;
        continue;
      }
      if (B.p[b] > nowy_k)
        break;

      if (B.p[b] >= nowy_p + q)
        zbior_push(&C, nowy_p, B.p[b] - q);

      if (B.k[b] + q > nowy_p)
        nowy_p = B.k[b] + q;

      if (nowy_p > nowy_k)
        break;

      b++;
    }

    if (nowy_p <= nowy_k)
      zbior_push(&C, nowy_p, nowy_k);
    a++;
  }

  return C;
}

bool nalezy(zbior_ary A, long long x) {
  if (A.uzyte == 0)
    return false;

  long long niski = 0, wysoki = A.uzyte - 1, srodek;
  while (niski <= wysoki) {
    srodek = niski + (wysoki - niski) / 2;
    if (modulo(A.p[srodek]) < modulo(x) ||
        (modulo(A.p[srodek]) == modulo(x) && A.k[srodek] < x)) {
      niski = srodek + 1;
    } else if (modulo(A.p[srodek]) > modulo(x) ||
               (modulo(A.p[srodek]) == modulo(x) && A.p[srodek] > x)) {
      wysoki = srodek - 1;
    } else {
      return (A.p[srodek] <= x && x <= A.k[srodek] &&
              ((x - A.p[srodek]) % q == 0));
    }
  }
  return false;
}

unsigned ary(zbior_ary A) { return (unsigned)A.uzyte; }

long long moc(zbior_ary A) {
  long long wyn = 0;
  for (long long i = 0; i < A.uzyte; i++) {
    wyn += (((long long)A.k[i] - (long long)A.p[i]) / q + 1);
  }
  return wyn;
}
