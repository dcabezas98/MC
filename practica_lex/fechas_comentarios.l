%{
// David Cabezas Berrido

/*
Programa para cambiar las fechas de un texto de formato, recibe un texto y lo imprime tal cual pero convirtiendo las fechas al formato especificado
*/

/*
También elimina comentarios en c, c++ y java
*/

#include<stdio.h>

int formato;
char* resultado;

const char instrucciones[] = "Formato 1:\tdd/mm/aaaa\nFormato 2:\tdd-mm-aaaa\nFormato 3:\tdd de <mes en español> de aaaa\nFormato 4:\tdd <mes en inglés> aaaa\nFormato 5:\t<mes en inglés> dd, aaaa\nIntroduzca número de formato: ";

struct fecha {
	char dia[3];
	char mes[3];
	char anio[5];
};

struct fecha leer(char* cadena, int longitud, int formato_entrada);	   // Lee la cadena de entrada, llama a una de las otras según el formato.
struct fecha leer12(char* cadena, int longitud, int formato_entrada);
struct fecha leer3(char* cadena, int longitud);
struct fecha leer4(char* cadena, int longitud);
struct fecha leer5(char* cadena, int longitud);

void formatear(struct fecha f, int formato_salida); // La guarda en resultado, llama a una de las otras según el formato
void formatear12(struct fecha f, int formato_salida);
void formatear3(struct fecha f);
void formatear4(struct fecha f);
void formatear5(struct fecha f);

// Funciones para pasar el mes de número a letra
const char* mes_ntoesp(char* m);
const char* mes_espton(char* m);
const char* mes_ntoeng(char* m);
const char* mes_engton(char* m);

void procesar_y_mostrar(char* cadena, int longitud, int formato_entrada);
%}

comentario	(\/\/.*)|(\/\*([^\*]|\*+[^\*\/])*\*+\/)

dia			[1-3]?[0-9]
mes			1?[0-9]
anio		[0-9]{1,4}

mes_espanol		enero|febrero|marzo|abril|mayo|junio|julio|agosto|septiembre|octubre|noviembre|diciembre
mes_ingles		January|February|March|April|May|June|July|August|September|October|November|December

formato1	{dia}\/{mes}\/{anio}
formato2	{dia}-{mes}-{anio}
formato3	{dia}\ de\ {mes_espanol}\ de\ {anio} 
formato4	{dia}\ {mes_ingles}\ {anio}
formato5	{mes_ingles}\ {dia},\ {anio}

%%

{comentario}	{}
{formato1}		{procesar_y_mostrar(yytext, yyleng, 1);}
{formato2}		{procesar_y_mostrar(yytext, yyleng, 2);}
{formato3}		{procesar_y_mostrar(yytext, yyleng, 3);}
{formato4}		{procesar_y_mostrar(yytext, yyleng, 4);}
{formato5}		{procesar_y_mostrar(yytext, yyleng, 5);}

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

	if(formato < 1 || formato > 5) {
		printf("Formato desconocido\n");
		exit(-1);
	}

	resultado = (char*)malloc(32);

	yylex();

	free(resultado);

	return 0;
}

struct fecha leer(char* cadena, int longitud, int formato_entrada){
	
	if(formato_entrada == 3)
		return leer3(cadena, longitud);

	if(formato_entrada == 4)
		return leer4(cadena, longitud);

	if(formato_entrada == 5)
		return leer5(cadena, longitud);
	
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

struct fecha leer3(char* cadena, int longitud){

	struct fecha f;

	char mes[11];

	int i = 0, // Para leer cadena
		j;
	
	for(j=0; cadena[i] != ' '; i++, j++) // Copio día
		f.dia[j] = cadena[i];

	f.dia[j] = '\0';
	i+=4;

	for(j=0; cadena[i] != ' '; i++, j++) // Copio mes
		mes[j] = cadena[i];

	mes[j] = '\0';

	sprintf(f.mes, "%s", mes_espton(mes));
	i+=4;

	for(j=0; i < longitud; i++, j++) // Copio anio
		f.anio[j] = cadena[i];

	f.anio[j] = '\0';

	return f;
}

struct fecha leer4(char* cadena, int longitud){
	
	struct fecha f;

	char mes[10];

	int i = 0, // Para leer cadena
		j;

	for(j=0; cadena[i] != ' '; i++, j++) // Copio día
		f.dia[j] = cadena[i];

	f.dia[j] = '\0';
	i++;

	for(j=0; cadena[i] != ' '; i++, j++) // Copio mes
		mes[j] = cadena[i];

	mes[j] = '\0';
	i++;

	sprintf(f.mes, "%s", mes_engton(mes));

	for(j=0; i < longitud; i++, j++) // Copio anio
		f.anio[j] = cadena[i];

	f.anio[j] = '\0';

	return f;
}

struct fecha leer5(char* cadena, int longitud){
	
	struct fecha f;

	char mes[10];

	int i = 0, // Para leer cadena
		j;

	for(j=0; cadena[i] != ' '; i++, j++) // Copio mes
		mes[j] = cadena[i];

	mes[j] = '\0';
	i++;

	sprintf(f.mes, "%s", mes_engton(mes));

	for(j=0; cadena[i] != ','; i++, j++) // Copio día
		f.dia[j] = cadena[i];

	f.dia[j] = '\0';
	i+=2;

	for(j=0; i < longitud; i++, j++) // Copio anio
		f.anio[j] = cadena[i];

	f.anio[j] = '\0';

	return f;
}


void formatear(struct fecha f, int formato_salida){
	
	if(formato_salida == 3)
		formatear3(f);

	else if(formato_salida == 4)
		formatear4(f);

	else if(formato_salida == 5)
		formatear5(f);

	else
		formatear12(f, formato_salida);
}

void formatear12(struct fecha f, int formato_salida){

	char separador = '#';

	if(formato_salida == 1)
		separador = '/';
	else if(formato_salida == 2)
		separador = '-';

	if(separador != '#') // Formato 1 o 2
		sprintf(resultado, "%s%c%s%c%s", f.dia, separador, f.mes, separador, f.anio);
}

void formatear3(struct fecha f){
	
	sprintf(resultado, "%s de %s de %s", f.dia, mes_ntoesp(f.mes), f.anio);
}

void formatear4(struct fecha f){
	sprintf(resultado, "%s %s %s", f.dia, mes_ntoeng(f.mes), f.anio);
}

void formatear5(struct fecha f){
	sprintf(resultado, "%s %s, %s", mes_ntoeng(f.mes), f.dia, f.anio);
}

void procesar_y_mostrar(char* cadena, int longitud, int formato_entrada){
	struct fecha f = leer(yytext, yyleng, formato_entrada); 
	formatear(f, formato); 
	printf("%s", resultado);
}

const char* mes_ntoesp(char* m){

	if(strncmp("2", m, 1)==0)
		return "febrero";

	if(strncmp("3", m, 1)==0)
		return "marzo";

	if(strncmp("4", m, 1)==0)
		return "abril";

	if(strncmp("5", m, 1)==0)
		return "mayo";

	if(strncmp("6", m, 1)==0)
		return "junio";

	if(strncmp("7", m, 1)==0)
		return "julio";

	if(strncmp("8", m, 1)==0)
		return "agosto";

	if(strncmp("9", m, 1)==0)
		return "septiembre";

	if(strncmp("10", m, 2)==0)
		return "octubre";

	if(strncmp("11", m, 2)==0)
		return "noviembre";

	if(strncmp("12", m, 2)==0)
		return "diciembre";

	if(strncmp("1", m, 1)==0)
		return "enero";
}

const char* mes_espton(char* m){
	
	if(strncmp("enero", m, 5)==0)
		return "1";

	if(strncmp("febrero", m, 7)==0)
		return "2";

	if(strncmp("marzo", m, 5)==0)
		return "3";

	if(strncmp("abril", m, 5)==0)
		return "4";

	if(strncmp("mayo", m, 4)==0)
		return "5";

	if(strncmp("junio", m, 5)==0)
		return "6";

	if(strncmp("julio", m, 5)==0)
		return "7";

	if(strncmp("agosto", m, 6)==0)
		return "8";

	if(strncmp("septiembre", m, 10)==0)
		return "9";

	if(strncmp("octubre", m, 7)==0)
		return "10";

	if(strncmp("noviembre", m, 9)==0)
		return "11";

	if(strncmp("diciembre", m, 9)==0)
		return "12";
}

const char* mes_ntoeng(char* m){

	if(strncmp("2", m, 1)==0)
		return "February";

	if(strncmp("3", m, 1)==0)
		return "March";

	if(strncmp("4", m, 1)==0)
		return "April";

	if(strncmp("5", m, 1)==0)
		return "May";

	if(strncmp("6", m, 1)==0)
		return "June";

	if(strncmp("7", m, 1)==0)
		return "July";

	if(strncmp("8", m, 1)==0)
		return "August";

	if(strncmp("9", m, 1)==0)
		return "September";

	if(strncmp("10", m, 2)==0)
		return "October";

	if(strncmp("11", m, 2)==0)
		return "November";

	if(strncmp("12", m, 2)==0)
		return "December";

	if(strncmp("1", m, 1)==0)
		return "January";
}

const char* mes_engton(char* m){

	if(strncmp("January", m, 7)==0)
		return "1";

	if(strncmp("February", m, 8)==0)
		return "2";

	if(strncmp("March", m, 5)==0)
		return "3";

	if(strncmp("April", m, 5)==0)
		return "4";

	if(strncmp("May", m, 3)==0)
		return "5";

	if(strncmp("June", m, 4)==0)
		return "6";

	if(strncmp("July", m, 4)==0)
		return "7";

	if(strncmp("August", m, 6)==0)
		return "8";

	if(strncmp("September", m, 9)==0)
		return "9";

	if(strncmp("October", m, 7)==0)
		return "10";

	if(strncmp("November", m, 8)==0)
		return "11";

	if(strncmp("December", m, 10)==0)
		return "12";
}
