/*Encuentra comentarios en programas en c, c++ y java*/

%{
#include<stdio.h>
%}

comentario	(\/\*[^\*]*(\*+[^\*\/][^\*]*)*\*+\/)|(\/\/.*$)

%%

{comentario}		{ printf("%s\n", yytext); }
.|\n				{}

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

	yylex();
	return 0;
}
