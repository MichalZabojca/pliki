#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

int q = 3;

typedef struct zbior_ary {
    int rozmiar;
    int uzyte;
    int *p;
    int *k;
} zbior_ary;

zbior_ary nowy_zbior(int rozmiar) {
    int *p = (int *)malloc(sizeof(int) * (size_t)rozmiar);
    if (p == NULL) {
        exit(1);
    }

    int *k = (int *)malloc(sizeof(int) * (size_t)rozmiar);
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

void zbior_push(zbior_ary *A, int nowy_p, int nowy_k) {
    if (A->uzyte == A->rozmiar) {
        int *p_r = realloc(A->p, sizeof(int) * (size_t)((A->rozmiar) * 2));
        if (p_r == NULL) {
            exit(1);
        }

        int *k_r = realloc(A->k, sizeof(int) * (size_t)((A->rozmiar) * 2));
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

zbior_ary ciag_arytmetyczny(int a, int Q, int b) {
    q = Q;
    zbior_ary A = nowy_zbior(1);
    zbior_push(&A, a, b);
    return A;
}

zbior_ary singleton(int a) {
    zbior_ary A = nowy_zbior(1);
    zbior_push(&A, a, a);
    return A;
}

void pierwszy_przedzial(int *nowy_p, int *nowy_k, zbior_ary *A, zbior_ary *B,
                        int *a, int *b) {
    if (B->uzyte <= *b && A->uzyte <= *a) {
        return;
    }

    if (B->uzyte <= *b && A->uzyte > *a) {
        *nowy_p = A->p[*a];
        *nowy_k = A->k[*a];
        (*a)++;
        return;
    }

    if (A->uzyte <= *a && B->uzyte > *b) {
        *nowy_p = B->p[*b];
        *nowy_k = B->k[*b];
        (*b)++;
        return;
    }

    if (A->p[*a] % q > B->p[*b] % q) {
        *nowy_p = B->p[*b];
        *nowy_k = B->k[*b];
        (*b)++;
        return;
    }

    if (B->p[*b] % q > A->p[*a] % q) {
        *nowy_p = A->p[*a];
        *nowy_k = A->k[*a];
        (*a)++;
        return;
    }

    if (B->p[*b] <= A->p[*a]) {
        *nowy_p = B->p[*b];
        *nowy_k = B->k[*b];
        (*b)++;
        return;
    }

    if (B->p[*b] > A->p[*a]) {
        *nowy_p = A->p[*a];
        *nowy_k = A->k[*a];
        (*a)++;
        return;
    }

    return;
}

zbior_ary suma(zbior_ary A, zbior_ary B) {
    int a = 0, b = 0, nowy_p, nowy_k;
    zbior_ary C = nowy_zbior(A.rozmiar);

    pierwszy_przedzial(&nowy_p, &nowy_k, &A, &B, &a, &b);
    a = 0;
    b = 0;
    while (b < B.uzyte || a < A.uzyte) {
        pierwszy_przedzial(&nowy_p, &nowy_k, &A, &B, &a, &b);

        fflush(stdout);
        while ((b < B.uzyte && B.p[b] % q == nowy_p % q && nowy_k >= B.p[b]) ||
               (a < A.uzyte && A.p[a] % q == nowy_p % q && nowy_k >= A.p[a])) {
            while (b < B.uzyte && nowy_k >= B.p[b] &&
                   B.p[b] % q == nowy_p % q) {
                if (B.k[b] > nowy_k) {
                    nowy_k = B.k[b];
                }
                b++;
            }

            while (a < A.uzyte && nowy_k >= A.p[a] &&
                   A.p[a] % q == nowy_p % q) {
                if (A.k[a] > nowy_k) {
                    nowy_k = A.k[a];
                }
                a++;
            }
        }
        zbior_push(&C, nowy_p, nowy_k);
    }

    return C;
}

zbior_ary iloczyn(zbior_ary A, zbior_ary B) {
    int a = 0, b = 0, nowy_p, nowy_k;
    zbior_ary C = nowy_zbior(A.rozmiar);

    pierwszy_przedzial(&nowy_p, &nowy_k, &A, &B, &a, &b);
    a = 0;
    b = 0;

    while (b < B.uzyte || a < A.uzyte) {
        pierwszy_przedzial(&nowy_p, &nowy_k, &A, &B, &a, &b);

        while ((b < B.uzyte && B.p[b] % q == nowy_p % q && nowy_k >= B.p[b]) ||
               (a < A.uzyte && A.p[a] % q == nowy_p % q && nowy_k >= A.p[a])) {
            while (b < B.uzyte && nowy_k >= B.p[b] &&
                   B.p[b] % q == nowy_p % q) {
                if (B.k[b] > nowy_k) {
                    zbior_push(&C, B.p[b], nowy_k);
                    nowy_k = B.k[b];
                } else {
                    zbior_push(&C, B.p[b], B.k[b]);
                }
                b++;
            }
            while (a < A.uzyte && nowy_k >= A.p[a] &&
                   A.p[a] % q == nowy_p % q) {
                if (A.k[a] > nowy_k) {
                    zbior_push(&C, A.p[a], nowy_k);
                    nowy_k = A.k[a];
                } else {
                    zbior_push(&C, A.p[a], A.k[a]);
                }
                a++;
            }
        }
    }
    return C;
}

zbior_ary roznica(zbior_ary A, zbior_ary B) {
    int a = 0, b = 0, nowy_k = 0;
    zbior_ary C = nowy_zbior(A.rozmiar);

    while (b < B.uzyte) {
        while (A.p[a] % q < B.p[b] % q && A.uzyte > a) {
            zbior_push(&C, A.p[a], A.k[a]);
            a++;
        }

        while (A.p[a] % q > B.p[b] % q && B.uzyte > b) {
            b++;
        }

        while (A.p[a] % q == B.p[b] % q && A.uzyte > a && B.uzyte > b) {
            nowy_k = A.p[a];

            if (B.k[b] >= A.p[a] && B.p[b] <= A.p[a]) {
                nowy_k = B.k[b] + q;
                b++;
            }
            while (B.p[b] - nowy_k > q && B.k[b] <= A.k[a] && B.uzyte > b) {
                zbior_push(&C, nowy_k, B.p[b] - q);
                nowy_k = B.k[b] + q;
                b++;
            }
            if (B.p[b] <= A.k[a] && B.k[b] >= A.k[a] && a < A.uzyte &&
                B.p[b] - nowy_k > q) {
                zbior_push(&C, nowy_k, B.p[b]);
                a++;
            }
            if (B.p[b] > A.k[a] && A.k[a] - nowy_k > q && a < A.uzyte) {
                zbior_push(&C, nowy_k, A.k[a]);
                a++;
            }
        }
        while (a < A.uzyte) {
            if (nowy_k >= A.p[a] && nowy_k <= A.k[a]) {
                zbior_push(&C, nowy_k, A.k[a]);
                a++;
                continue;
            }
            zbior_push(&C, A.p[a], A.k[a]);
        }
    }

    return C;
}

bool nalezy(zbior_ary A, int x) {
    int dolny = 0, gorny = A.uzyte, srodek = 0;
    while (gorny >= dolny) {
        srodek = dolny + (gorny - dolny) / 2;

        if (A.p[srodek] <= x && A.k[srodek] >= x && A.p[srodek] % q == x % q) {
            return true;
        }

        if (A.p[srodek] % q < x % q || A.k[srodek] < x) {
            dolny = srodek + 1;
        } else if (A.p[srodek] % q > x % q || A.k[srodek] > x) {
            gorny = srodek - 1;
        }
    }
    return false;
}

unsigned ary(zbior_ary A) { return (unsigned)A.uzyte; }

unsigned moc(zbior_ary A) {
    unsigned int wyn = 0;
    for (int i = 0; i < A.uzyte; i++) {
        wyn += (unsigned)((A.k[i] - A.p[i]) / q + 1);
    }
    return wyn;
}
