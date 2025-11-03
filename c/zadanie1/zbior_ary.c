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

int min(int a, int b){
    if (a <= b) return a;
    else return b;
}

int max(int a, int b){
    if (a >= b) return a;
    else return b;
}

int modulo(int x) {
    int r = x % q;
    if (r < 0) r += q;
    return r;
}

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
        int *p_r = (int *)realloc(A->p, sizeof(int) * (size_t)((A->rozmiar) * 2));
        if (p_r == NULL) {
            exit(1);
        }

        int *k_r = (int *)realloc(A->k, sizeof(int) * (size_t)((A->rozmiar) * 2));
        if (k_r == NULL) {
            exit(1);
        }

        A->p = p_r;
        A->k = k_r;
            if (r < 0) r += q;
                return r;
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


zbior_ary suma(zbior_ary A, zbior_ary B) {
    int a = 0, b = 0, nowy_p = 0, nowy_k = 0;
    zbior_ary C = nowy_zbior(A.rozmiar);
    
    if (A.uzyte == 0) return B;
    if (B.uzyte == 0) return A;
    
    if (A.p[a] <= B.p[b]){
        nowy_p = A.p[a];
        nowy_k = A.k[a];
        a++;
    }

    else {
        nowy_p = B.p[b];
        nowy_k = B.k[b];
        b++;
    }

    while (b < B.uzyte && a < A.uzyte) {

        if (modulo(A.p[a]) < modulo(B.p[b])){
            zbior_push(&C, A.p[a], A.k[a]);
            a++;
            continue;
        }
        else if (modulo(B.p[b]) < modulo(A.p[a])){
            zbior_push(&C, B.p[b], B.k[b]);
            b++;
            continue;
        }

        if (A.p[a] <= B.p[b]){
            nowy_p = A.p[a];
            nowy_k = A.k[a];
            a++;
        }

        else {
            nowy_p = B.p[b];
            nowy_k = B.k[b];
            b++;
        }


        while ((a < A.uzyte && modulo(A.p[a]) == modulo(nowy_p) && A.p[a] <= nowy_k + q) || (b < B.uzyte && modulo(B.p[b]) == modulo(nowy_p) && B.p[b] <= nowy_k + q)){
            while (a < A.uzyte && modulo(A.p[a]) == modulo(nowy_p) && A.p[a] <= nowy_k + q){ 
                nowy_k = max(A.k[a], nowy_k);
                a++;
            }

            while (b < B.uzyte && modulo(B.p[b]) == modulo(nowy_p) && B.p[b] <= nowy_k + q){ 
                nowy_k = max(B.k[b], nowy_k);
                b++;
            }
        }
        zbior_push(&C, nowy_p, nowy_k);
    }

    while (a < A.uzyte){
        zbior_push(&C, A.p[a], A.k[a]);
        a++;
    }

    while (b < B.uzyte){
        zbior_push(&C, B.p[b], B.k[b]);
        b++;
    }

    return C;
}

zbior_ary iloczyn(zbior_ary A, zbior_ary B) {
    int a = 0, b = 0, nowy_p, nowy_k;
    zbior_ary C = nowy_zbior(A.rozmiar);

    while (b < B.uzyte && a < A.uzyte){
        if (modulo(B.p[b]) < modulo(A.p[a])) b++;
        else if (modulo(B.p[b]) > modulo(A.p[a])) a++;
        if (modulo(A.p[a]) == modulo(B.p[b])){
            nowy_p = max(A.p[a], B.p[b]);
            nowy_k = mina(A.k[a], B.k[b]);
            if (nowy_k >= nowy_p) zbior_push(&C, nowy_p, nowy_k);
            if (A.k[a] >= B.k[b]) b++;
            else a++;
        }
    }

    return C;
}

zbior_ary roznica(zbior_ary A, zbior_ary B) {
    int a = 0, b = 0, nowy_k = 0;
    zbior_ary C = nowy_zbior(A.rozmiar);

    while (a < A.uzyte){
        nowy_p = A.p[a];
        nowy_k = A.k[a];
        while (b < B.uzyte && modulo(B.p[b]) < modulo(A.p[a])) b++;
        while (b < B.uzyte && modulo(B.p[b]) == modulo(nowy_p) && B.p[b] <= nowy_k && nowy_p <= nowy_k){
            if (B.k[b] > nowy_p){
                zbior_push(&C, nowy_p, B.p[b] - q);
            }
            nowy_p = max(nowy_p, B.k[b] + q);
            b++;
        }
        if (nowy_p <= nowy_k) zbior_push(&C, nowy_p, nowy_k);
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
