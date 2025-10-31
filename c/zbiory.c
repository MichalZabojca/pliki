#include<stdlib.h>
#include<math.h>

int q = 3;

int main(){
	return 0;
}

typedef struct Zbiory_Ary{
    int rozmiar;
    int uzyte;
    int* p;
    int* k;
}zbiory_ary;

zbiory_ary* nowy_zbior(int rozmiar){
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

void zbior_push(zbiory_ary* A, int nowy_p, int nowy_k){
    if (A->uzyte == A->rozmiar){
        int* p = realloc(A->p, sizeof(int) * A->rozmiar * 2);
        if (p == NULL){
            return 0;
        }

        int* k = realloc(A->k, sizeof(int) * A->rozmiar * 2);
        if (k == NULL){
            return 0;
        }

        A->p = p;
        A->k = k;
        A->rozmiar *= 2;
    }
    *(A->p + A->uzyte + 1) = nowy_p; 
    *(A->k + A->uzyte + 1) = nowy_k; 
    A->uzyte += 1;
}

zbiory_ary* suma(zbiory_ary* A, zbiory_ary* B){
    int a = 0, b = 0, nowy_p, nowy_k;
    zbiory_ary* C = nowy_zbior(A->rozmiar);
	while (b <= B->uzyte && a <= A->uzyte){
		nowy_p = fmin(A->p[a], B->p[b]);
		while ((A->p[a] % q == nowy_p % q || B->p[b] % q == nowy_p % q) && b <= B->uzyte && a <= A->uzyte){
			nowy_p = fmin(A->p[a], B->p[b]);
			while (nowy_k >= B->p[b] && B->p[b] % q == nowy_p % q && b <= B->uzyte){
				nowy_k = B->k[b];
				b++;
			}

			while (nowy_k >= A->p[a] && A->p[a] % q == nowy_p % q && a <= A->uzyte){
				nowy_k = A->k[a];
				a++;
			}
			zbior_push(C, nowy_p, nowy_k);
		}
	}
	return C;
}


