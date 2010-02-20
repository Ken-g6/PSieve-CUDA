#include <stdio.h>
#include <string.h>

int main(void) {
	printf("%d == %d?", ffs(42), __builtin_ctz(42));
}

