#include<stdlib.h>
#include<stdio.h>


int q = 3;



int min(int a, int b){
    if (a <= b){
        return a;
    }
    return b;
}

typedef struct Zbiory_Ary{
    int rozmiar;
    int uzyte;
    int* p;
    int* k;
}zbiory_ary;



zbiory_ary* nowy_zbior(int rozmiar){
    int* p = (int*)malloc(sizeof(int) * rozmiar);
    if (p == NULL){
        return NULL;
    }

    int* k = (int*)malloc(sizeof(int) * rozmiar);
    if (k == NULL){
        return NULL;
    }

    zbiory_ary* A = malloc(sizeof(zbiory_ary));
    
    A->rozmiar = rozmiar;
    A->uzyte = 0;
    A->p = p;
    A->k = k;
    return A;
}

void zbior_push(zbiory_ary* A, int nowy_p, int nowy_k){
    if (A->uzyte == A->rozmiar){
        int* p_r = realloc(A->p, sizeof(int) * (A->rozmiar) * 2);
        if (p_r == NULL){
            exit(1);
        }

        int* k_r = realloc(A->k, sizeof(int) * (A->rozmiar) * 2);
        if (k_r == NULL){
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




zbiory_ary* ciag_arytmetyczny(int a, int b){
    zbiory_ary* A = nowy_zbior(1);
    zbior_push(A, a, b);
    return A;
}


void pierwszy_przedzial(int* nowy_p, int* nowy_k, zbiory_ary* A, zbiory_ary* B, int* a, int* b){

    if ((A->p[*a] <= B->p[*b] && A->uzyte > *a) || (B->uzyte <= *b && A->uzyte > *a)){
        *nowy_p = A->p[*a];
        *nowy_k = A->k[*a];
        (*a)++;
        return;
    }
    if((B->p[*b] <= A->p[*a] && B->uzyte > *b) || (A->uzyte <= *a && B->uzyte > *b)){
        *nowy_p = B->p[*b];
        *nowy_k = B->k[*b];
        (*b)++;
        return;
    }
    return;
}


zbiory_ary* suma(zbiory_ary* A, zbiory_ary* B){
    int a = 0, b = 0, nowy_p, nowy_k;
    zbiory_ary* C = nowy_zbior(A->rozmiar);

    pierwszy_przedzial(&nowy_p, &nowy_k, A, B, &a, &b);
    a=0;
    b=0;
    while (b < B->uzyte || a < A->uzyte){
        pierwszy_przedzial(&nowy_p, &nowy_k, A, B, &a, &b);

        while ((b < B->uzyte || a < A->uzyte) && (B->p[b] % q == nowy_p % q || A->p[a] % q == nowy_p % q)){
            while (nowy_k >= B->p[b] && B->p[b] % q == nowy_p % q && b < B->uzyte){
                if (B->k[b] > nowy_k){
                    nowy_k = B->k[b];
                }
                b++;
            }

            while (nowy_k >= A->p[a] && A->p[a] % q == nowy_p % q && a < A->uzyte){
                if (A->k[a] > nowy_k){
                    nowy_k = A->k[a];
                }
                a++;
            }
        }
        zbior_push(C, nowy_p, nowy_k);
    }
	return C;
}


int main(){
    zbiory_ary* A = ciag_arytmetyczny(1, 13);
    printf("\n\n");
    zbiory_ary* B = ciag_arytmetyczny(2, 14);
    zbiory_ary* D = ciag_arytmetyczny(5, 17);
    zbiory_ary* C = suma(A, B);
    C = suma(C, D);
    for (int i=0; i<C->uzyte; i++){
        printf("%d %d\n", C->p[i], C->k[i]);
    }
	return 0;
}
