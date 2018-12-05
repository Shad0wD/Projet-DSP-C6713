; test de la fonction function_carre_256_float_to_float
		.global _function_carre_256_float_to_float

_function_carre_256_float_to_float:
		mvk	256 , A1
		mv A4, A2
loop:
		ldw *A4, A0
		nop 4
		mpysp A0, A0, A0
		nop 3
		stw  A0, *A4++
		nop 4
		sub A1, 1, A1
  [A1]	B loop
		nop 5
		b B3
		nop 5

