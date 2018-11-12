%{
/*
Programa para cambiar las fechas de un texto de formato, recibe un texto y lo imprime tal cual pero convirtiendo las fechas al formato especificado
*/

#include<stdio.h>

int formato;

const char[] = "Formato 1:\tdd/mm/aaaa\nFormato2:\tdd-mm-aaaa\nFormato3:\tdd de <mes en letra> de aaaa\nIntroduzca número de formato: ";

const char* convertir_a_formato12(int dia, int mes, int anio, int formato);
const char* convertir_a_formato3(int dia, char* mes, int anio);

%}

dia			[1-3]?[0-9]
mes			1?[0-9]
anio		[0-9][0-9][0-9][0-9]

mesletra	enero|febrero|marzo|abril|mayo|junio|julio|agosto|septiembre|octubre|noviembre|diciembre

// dd/mm/aaaa
formato1	{dia}\/{mes}\/{anio}
// dd-mm-aaaa
formato2	{dia}-{mes}-{anio}
// dd de <mes en letra> de aaaa
formato3	{dia}\bde\b{mesletra}\bde\b{anio} 

%%

%%

int main(int argc, char* argv[]){

	if(argc==2){
		yyin = fopen(argv[1],"rt");
		if(yyin==NULL){
			printf("El fichero %s no se puede abrir.\n", argv[1]);
			exit(-1);
		}
	}
	else yyin=stdin;

	printf("%s", introduccion());
	scanf("%d", formato);

	yylex();
	return 0;
}

const char* convertir_a_formato(int dia, int mes, int anio, int formato){

	char resultado[10];
	char separator;

	if(formato == 1)
		separador = '/';
	
	else //if(formato == 2)
		separador = '-';

	sprintf(resultado, "%d%c%d%c%d", dia, separador, mes, separador, anio);
	
	return resultado;
}

const char* convertir_a_formato3(int dia, char* mes, int anio){
	
	char resultado[32];

	sprintf(resultado, "%d de %s de %d", dia, mes, anio);

	return resultado;
}




































