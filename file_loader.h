#pragma once
#include<stdio.h>
/*
* include files in the main function ->
* pass them into fuynctions here ->
* there will be a function to take a file and return a vector of its data ->
* get the pointer of that vector and use it .
*/

unsigned char* load_raw_filedata(const char* filename) {
	unsigned char* buffer = new unsigned char[2000 * 1080 * 10];
	FILE* pF;
	fopen_s(&pF, filename, "rb");
	if (pF == NULL) { printf("ERROR1 COULD NOT READ: FATAL"); return NULL; }
	fread(buffer, 10, 1920 * 1080 , pF);

	fclose(pF);
	return buffer;

}

