%{
/*
Programa para cambiar las fechas de un texto de formato, recibe un texto y lo imprime tal cual pero convirtiendo las fechas al formato especificado
*/

#include<stdio.h>

int formato;
char* resultado;

const char instrucciones[] = "Formato 1:\tdd/mm/aaaa\nFormato 2:\tdd-mm-aaaa\nFormato 3:\tdd de <mes en letra> de aaaa\nFromato 4:\t<mes en letra en inglés dd>, aaaa\nIntroduzca número de formato: ";

struct fecha {
	char dia[3];
	char mes[3];
	char anio[5];
};

struct fecha leer(char* cadena, int longitud, int formato_entrada);	   // Lee la cadena de entrada, llama a una de las otras según el formato.
struct fecha leer12(char* cadena, int longitud, int formato_entrada);
//struct fecha leer3(char* cadena, int longitud);
//struct fecha leer4(char* cadena, int longitud);

void formatear(struct fecha f, int formato_salida); // La guarda en resultado

void procesar_y_mostrar(char* cadena, int longitud, int formato_entrada);

%}

dia			[1-3]?[0-9]
mes			1?[0-9]
anio		[0-9]{1,4}

mes_espanol		enero|febrero|marzo|abril|mayo|junio|julio|agosto|septiembre|octubre|noviembre|diciembre
mes_ingles		January|February|March|April|May|June|July|August|September|October|November|December

formato1	{dia}\/{mes}\/{anio}
formato2	{dia}-{mes}-{anio}
formato3	{dia}\ de\ {mes_espanol}\ de\ {anio} 
formato4	{mes_ingles}\ {dia},\ {anio}

%%

{formato1}	{procesar_y_mostrar(yytext, yyleng, 1);}
{formato2}	{procesar_y_mostrar(yytext, yyleng, 2);}
{formato3}	{/*procesar_y_mostrar(yytext, yyleng, 3);*/}
{formato4}	{procesar_y_mostrar(yytext, yyleng, 4);}

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

	printf("%s", instrucciones);
	scanf("%d", &formato);

	if(formato < 1 || formato > 4) {
		printf("Formato desconocido\n");
		exit(-1);
	}

	resultado = (char*)malloc(32);

	yylex();

	free(resultado);

	return 0;
}

struct fecha leer(char* cadena, int longitud, int formato_entrada){
	
	/*if(formato_entrada == 3)
		return leer3(cadena, longitud);

	if(formato_entrada == 4)
		return leer4(cadena, longitud);*/
	
	return leer12(cadena, longitud, formato_entrada);
}

struct fecha leer12(char* cadena, int longitud, int formato_entrada){
		
	struct fecha f;

	int i = 0, // Para leer cadena
		j;

	char separador = '#';

	if(formato_entrada == 1)
		separador = '/';
	else if(formato_entrada == 2)
		separador = '-';

	if(separador != '#'){ // Formato 1 o 2

		for(j=0; cadena[i] != separador; i++, j++) // Copio día
			f.dia[j] = cadena[i];

		f.dia[j] = '\0';
		i++; // Salto separador

		for(j=0; cadena[i] != separador; i++, j++) // Copio mes
			f.mes[j] = cadena[i];

		f.mes[j] = '\0';
		i++; // Salto separador

		for(j=0; i < longitud; i++, j++) // Copio anio
			f.anio[j] = cadena[i];

		f.anio[j] = '\0';
	}

	return f;
}

void formatear(struct fecha f, int formato_salida){

	char separador = '#';

	if(formato_salida == 1)
		separador = '/';
	else if(formato_salida == 2)
		separador = '-';

	if(separador != '#') // Formato 1 o 2
		sprintf(resultado, "%s%c%s%c%s", f.dia, separador, f.mes, separador, f.anio);
}

void procesar_y_mostrar(char* cadena, int longitud, int formato_entrada){
	struct fecha f = leer(yytext, yyleng, formato_entrada); 
	formatear(f, formato); 
	printf("%s", resultado);
}
