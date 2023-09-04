/*
Copyright 2023 Heidi Van den Camp

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
A copy of the License is included in the repository and you can also
view it at
        http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#ifndef GRAPH_IO
#define GRAPH_IO
#include <stdio.h>
#include <stdint.h>
#include "lopsp_functions.h"

enum code {EDGE_CODE, PLANAR_CODE_BIG_ENDIAN, PLANAR_CODE_LITTLE_ENDIAN};

//Reads the following number_of_bytes bytes from stdin as a big endian number. Returns UINT64_MAX if the first byte is 255
unsigned long read_edge_number(FILE* file, int number_of_bytes);

//reads a big endian number with number_of_bytes bytes from file
unsigned long read_big_endian(FILE* file, int number_of_bytes);

//reads a little endian number with number_of_bytes bytes from file
unsigned long read_little_endian(FILE* file, int number_of_bytes);

void save_structure(Lopsp* lopsp);

void determine_types(Lopsp* lopsp, uint16_t typev0);

//Reads ONE lopsp-operation from the given file. Returns 0 if no operation is read, 1 if an operation is read
int read_lopsp(FILE* file, Lopsp* lopsp);

//Checks that the header is correct and returns the type of the encoding of the graphs.
enum code checkheader(FILE* file);

//Checks that the header is correct
void check_lopsp_header(FILE* file);

//'file' must start at the start of the planarcode (no header). This method reads the planarcode in the file and
// saves the corresponding graph in 'graph'. Returns 1 if a graph is read, 0 if not. little_endian should be 1 if
// the code is little endian and 0 otherwise
int read_planarcode(FILE* file, Graph* graph,  int little_endian);

//'file' must start at the start of the edgecode (no header). This method reads the edgecode in the file and
// saves the corresponding graph in 'graph'. Returns 1 if a graph is read, 0 if not
int read_edgecode(FILE* file, Graph* graph);

//Reads ONE embedded graph from the given file, either planarcode OR edgecode
//Return 0 if no graph is read, 1 if a graph is read
int read_graph(FILE* file, Graph* graph, enum code format);

//writes the header of the given code to file
void write_header(FILE* file, enum code code);

//writes graph to file in planarcode, no header, always little endian
void write_planarcode(FILE* file, Graph* graph);

//writes graph to file in edgecode, no header
void write_edgecode(FILE* file, Graph* graph);

//writes graph to file in the given code, no header
void write_graph(FILE* file, Graph* graph, enum code code);

//prints lopsp to file in a readable txt format
void print_lopsp(FILE* file, Lopsp* lopsp);

#endif