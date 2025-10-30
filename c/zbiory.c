typedef struct Zbiory_Ary{
    int rozmiar;
    int uzyte;
    int* p;
    int* k;
}zbiory_ary;

zbiory_ary* nowy_zbior(rozmiar){
    int* p = malloc(sizeof(int) * rozmiar);
    if (p == NULL){
        return 0;
    }

    int* k = malloc(sizeof(int) * rozmiar);
    if (k == NULL){
        return 0;
    }

    zbiory_ary* A;
    
    A->rozmiar = rozmiar;
    A->p = p;
    A->k = k;
    return A;
}

void zbior_push(zbiory_ary* A, int nowy_p, nowy_k){
    if (A->uzyte == A->rozmiar){
        int* p = realloc(A->p, sizeof(int) * rozmiar * 2);
        if (p == NULL){
            return 0;
        }

        int* k = realloc(A->k, sizeof(int) * rozmiar * 2);
        if (k == NULL){
            return 0;
        }

        A->p = p;
        A->k = k;
        rozmiar *= 2;
    }
    A->(p + uzyte + 1) = nowy_p; 
    A->(k + uzyte + 1) = nowy_k; 
    uzyte += 1;
}

zbiory_ary suma(zbiory_ary* A, zbiory_ary* B){
    int a = 0, b = 0, nowy_p, nowy_k;
    zbiory_ary* C = nowy_zbiory(A->rozmiar);
    nowy_p = min(A->p[a], B->p[b]);
    while ()
    while (A->p[a] <= B->p[b] && A->k[a] > B->p[b]){
        nowy_k = B->k[b];
        b++;
    }

    while (B->p[b] <= A->p[a] && B->k[b] > A->p[a]){
        nowy_k = A->k[a];
        while (nowy_k >= B->p[b]){
            nowy_k = B->k[b];
            b++;
        }
        a++;
    }

}
