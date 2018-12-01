%{
/*
Programa para cambiar las fechas de un texto de formato, recibe un texto y lo imprime tal cual pero convirtiendo las fechas al formato especificado
*/

#include<stdio.h>

int formato;
char* resultado[32];

const char[] instrucciones = "Formato 1:\tdd/mm/aaaa\nFormato 2:\tdd-mm-aaaa\nFormato 3:\tdd de <mes en letra> de aaaa\nFromato 4: <mes en letra en inglés dd, aaaa>\nIntroduzca número de formato: ";

fecha leer(char* cadena, int longitud, int formato_entrada);
void formatear(fecha f, int formato_salida); // La guarda en resultado

struct fecha {
	char dia[2];
	char mes[2];
	char anio[4];
}

/*
fecha char* convertir_a_formato12(int dia, int mes, int anio, int formato);
const char* convertir_a_formato3(int dia, char* mes, int anio);
*/

struct fecha {
	char dia[2];
	char mes[2];
	char anio[4];
}

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

formato1	{leer(fecha f = leer(yytext, yyleng, 1); formatear(f, formato); printf(resultado);}
formato2	{leer(fecha f = leer(yytext, yyleng, 2); formatear(f, formato); printf(resultado);}

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
	scanf("%d", formato);

	yylex();
	return 0;
}

fecha leer(char* cadena, int longitud, int formato_entrada){
		
	fecha f;

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

		i++; // Salto separador

		for(j=0; cadena[i] != separador; i++, j++) // Copio mes
			f.mes[j] = cadena[i];

		i++; // Salto separador

		for(j=0; i < longitud; i++, j++) // Copio anio
			f.anio[j] = cadena[i];
	}

	return f;
}

void formatear(fecha f, int formato_salida){

	char separador = '#';

	if(formato_entrada == 1)
		separador = '/';
	else if(formato_entrada == 2)
		separador = '-';

	if(separador != '#') // Formato 1 o 2
		sprintf(resultado, %s%s%s%s%s, f.dia, separador, f.mes, separador, f.anio);
}

/*
char* convertir_a_formato(int dia, int mes, int anio, int formato){

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
*/
