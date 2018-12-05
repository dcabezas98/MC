#!/bin/bash

flex fechas_comentarios.l
gcc lex.yy.c -o fechas_comentarios -lfl
./fechas_comentarios example.txt
