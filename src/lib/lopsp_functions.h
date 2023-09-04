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

#ifndef LOPSP_STRUCTS_H
#define LOPSP_STRUCTS_H

#include <stdint.h>

#define MAX_EDGES 20000 //The total number of edges that can be used.

typedef struct Directed_edge Edge;
struct Directed_edge{
    Edge* prev;
    Edge* next;
    Edge* inverse;

    uint16_t start;
    uint16_t end;
    uint16_t type;

    int unique_number; //Is a unique number between 0 and the number of directed edges in the graph - 1. Is assigned for certain operations such as copying.
};

typedef struct Embedded_Graph{
    uint16_t vertexcount;
    int edgecount; // the number of directed edges in the graph

    Edge** edges; //One edge for every vertex
    int edge_array_size; //The Edge* array edges has size edge_array_size
} Graph;

typedef struct Lopsp{
    Graph* graph;

    uint16_t v0; //the number of the vertex v0
    uint16_t v1;
    uint16_t v2;

    Edge** edge_array;

    uint8_t cut; // 1 if already cut, 0 if not

    int* edge_left_or_right_of_P; //has edgecount elements. On index i the array is 0 if the edge is not in P, 1 if it is on the left side (when cut) and -1 if it is on the right side
} Lopsp;

typedef struct Queue{
    Edge** list;
    int first;
    int next;
    int size;
} Queue;

void init_edges();

Edge* get_edge(uint16_t start, Graph* target_graph);

void return_edge(Edge* edge, Graph* graph);

int number_of_edges_allocated();

Graph* new_graph(int initialsize);

Lopsp* new_lopsp(int initialsize);

void init_graph(Graph* graph);

void init_lopsp(Lopsp* lopsp);

void free_lopsp(Lopsp* lopsp);

void free_graph(Graph* graph);

Queue* new_queue(int size);

void queue_clear(Queue* queue);

void free_queue(Queue* queue);

void queue_add(Queue* queue, Edge* new_edge);

Edge* queue_pop(Queue* queue);

void set_vertexcount(Graph* graph, int new_vertexcount);

Edge* next_in_face(Edge* edge);

int facesize(Edge* start_edge);

//inserts the edge new as the successor of original. Also fills in the start. The inverse and end of new are not filled in
void insert_after(Edge* original, Edge* new_edge);

//Assumes that the starts of each edge are correct. This function makes e1 and e2 inverses and fills in their ends
void make_inverse(Edge* e1, Edge* e2);

void make_next(Edge* e1, Edge* e2);

//removes edge (and its inverse) from graph. Traverses the whole graph to restore the unique_numbers!!!
//DO NOT REMOVE AN EDGE WITH A VERTEX OF DEGREE 1
//only used in tests
void remove_edge(Edge* edge, Graph* graph);

void add_1vertex_in_face(Edge* edge, Graph* graph);

//Finds a cut-path and saves it in lopsp->edge_left_or_right_of_P.
void cut(Lopsp* lopsp);

//Applies the lopsp-operation with patch 'patch' to 'graph'.
void apply(Lopsp* lopsp, Graph* graph, Graph* result);

int has_loop(Graph* graph);

int has_multiple_edges(Graph* graph);

//Checks if graph is ck (k-connected, face-width at least k and face size at least k) for k=1,2,3.
//Returns 1 if graph is not c2, 2 if graph is c2 but not c3 and 3 if graph is c3
int ck(Graph* graph);

#endif