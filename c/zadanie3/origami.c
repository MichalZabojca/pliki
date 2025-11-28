#include<stdio>

int warstwy(a, x, y){
	if (typ[a] == "K"){
	}
	if (typ[a] == "P"){
	}
	if (typ[a] == "Z"){
		<F10>
	}
}







int main(){
	int n, q;
	scanf("%d%d", &n, &q);
	for (int i = 0; i < n; i++){
		scanf("%c", &typ[i]);
		switch(typ[i]){
			case "K":
				scanf("%lf%lf%lf", &x1[i], &x1[i], &x2[i]);
				break;
			case "P":
				scanf("%lf%lf%lf%lf", &x1[i], &y1[i], &x2[i], &y2[i]);
				break;
			case "Z":
				scanf("%lf%lf%lf%lf%lf", &ref[i], &x1[i], &y1[i], &x2[i], &y2[i]);
				break;
		}
	}
	for (int i = 0; i < q; i++){
		scanf("%d%lf%lf", &a, &x, &y);
		printf("%d", warstwy(a, x, y));
	}
	return 0;
}
