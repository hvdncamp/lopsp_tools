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

#include "lopsp_functions.h"
#include <stdio.h>
#include <stdlib.h>

static int allocated_edges = 0;
static Edge edges[MAX_EDGES]; //Allocate edges in an array so that they can be compared
static Edge* edgeptrs[MAX_EDGES];//These point to elements of 'edges'

void init_edges(){
    for (int i = 0; i < MAX_EDGES; i++){
        edgeptrs[i] = edges + i;
    }
    allocated_edges=0;
}

Edge* get_edge(uint16_t start, Graph* target_graph){
    if(allocated_edges >= MAX_EDGES){
        fprintf(stderr, "Too many edges are used. Increase MAX_EDGES\n");
        exit(1);
    }
    allocated_edges++;
    edgeptrs[allocated_edges - 1]->start = start;
    edgeptrs[allocated_edges - 1]->type = 3;
    edgeptrs[allocated_edges - 1]->unique_number = target_graph->edgecount;
    target_graph->edgecount++;
    return edgeptrs[allocated_edges - 1];
}

//Beware!! after this function it is not guaranteed that the unique numbers of the edges in the graph are exactly {0,...,edgecount - 1}
void return_edge(Edge* edge, Graph* graph){
    if(allocated_edges == 0){
        fprintf(stderr, "Trying to free an edge when there are no allocated edges\n");
        exit(1);
    }
    allocated_edges--;
    graph->edgecount--;
    edgeptrs[allocated_edges] = edge;
}

int number_of_edges_allocated(){
    return allocated_edges;
}

Graph* new_graph(int initialsize){
    Graph* graph = malloc(sizeof(Graph));
    graph->edges = malloc(sizeof(Edge*) * initialsize);
    graph->edge_array_size = initialsize;
    graph->vertexcount = 0;
    graph->edgecount = 0;
    return graph;
}

Lopsp* new_lopsp(int initialsize){
    Lopsp* lopsp = malloc(sizeof(Lopsp));
    lopsp->graph = new_graph(initialsize);
    lopsp->edge_array = NULL;
    lopsp->cut = 0;
    return lopsp;
}

void init_graph(Graph* graph){
    Edge* last;
    Edge* curr;
    for (int i = 1; i <= graph->vertexcount; i++){
        last = graph->edges[i]->prev;
        curr = graph->edges[i];
        while(curr !=last){
            curr = curr->next;
            return_edge(curr->prev, graph);
        }
        return_edge(last, graph);
    }
    graph->vertexcount = 0;
    graph->edgecount = 0;
}

void init_lopsp(Lopsp* lopsp){
    init_graph(lopsp->graph);
    if(lopsp->edge_array != NULL){
        free(lopsp->edge_array);
    }
    if(lopsp->cut){
        free(lopsp->edge_left_or_right_of_P);
    }
    lopsp->cut = 0;
    lopsp->edge_array = NULL; // is allocated and filled in save_structure
}

void free_lopsp(Lopsp* lopsp){
    free_graph(lopsp->graph);
    if(lopsp->edge_array != NULL){
        free(lopsp->edge_array);
    }
    if(lopsp->cut){
        free(lopsp->edge_left_or_right_of_P);
    }
    free(lopsp);
}

void free_graph(Graph* graph){
    init_graph(graph);//Return all the used edges
    free(graph->edges);
    free(graph);
}

Queue* new_queue(int size){
    Queue* queue = malloc(sizeof(Queue));
    queue->first = 0;
    queue->next = 0;
    queue->size = size + 1;//the array is one larger than the number of elements the queue can hold
    queue->list = malloc((size+1) * sizeof(Edge*));
    return queue;
}

void queue_clear(Queue* queue){
    queue->first = 0;
    queue->next = 0;
}

void free_queue(Queue* queue){
    free(queue->list);
    free(queue);
}

void queue_add(Queue* queue, Edge* new_edge){
    if((queue->next +1)%queue->size == queue->first){//the last element of list is not filled, else it is impossible to distinguish between an empty and a full queue
        fprintf(stderr, "Queue is full. Either something went wrong or the size of the queue must be increased.\n");
        exit(1);
    }
    queue->list[queue->next] = new_edge;
    queue->next = (queue->next + 1) % queue->size;
}

Edge* queue_pop(Queue* queue){
    Edge* to_return;
    if(queue->next == queue->first){
        return NULL;
    }
    to_return = queue->list[queue->first];
    queue->first = (queue->first + 1) % queue->size;
    return to_return;
}

//graph->edges grows dynamically
void set_vertexcount(Graph* graph, int new_vertexcount){
    if(new_vertexcount >= graph->edge_array_size){ //>= and not > because we count from 1 in the edges array
        graph->edges = realloc(graph->edges, (new_vertexcount + graph->edge_array_size)*sizeof (Edge*));
        graph->edge_array_size = new_vertexcount + graph->edge_array_size;
    }
    graph->vertexcount = new_vertexcount;
}

Edge* next_in_face(Edge* edge){
    return (edge->inverse)->next;
}

int facesize(Edge* start_edge){
    Edge* curr = start_edge;
    int size = 0;
    do {
        curr = next_in_face(curr);
        size++;
    } while(curr!= start_edge);
    return size;
}

//inserts the edge new as the successor of original. Also fills in the start. The inverse and end of new are not filled in
void insert_after(Edge* original, Edge* new_edge){
    new_edge->prev = original;
    new_edge->next = original->next;
    new_edge->start = original->start;

    new_edge->next->prev = new_edge;
    original->next = new_edge;
}

//Assumes that the starts of each edge are correct. This function makes e1 and e2 inverses and fills in their ends
void make_inverse(Edge* e1, Edge* e2){
    e1->inverse = e2;
    e2->inverse = e1;
    e1->end = e2->start;
    e2->end = e1->start;
}

void make_next(Edge* e1, Edge* e2){
    e1->next = e2;
    e2->prev = e1;
}

void remove_directed_edge(Edge* edge, Graph* graph){
    Edge* temp;
    edge->prev->next = edge->next;
    edge->next->prev = edge->prev;
    if(graph->edges[edge->start] == edge) {
        if (edge->next == edge) { // edge was the last edge starting in a vertex. That vertex must be removed
            graph->edges[edge->start] = graph->edges[graph->vertexcount];
            temp = graph->edges[edge->start];
            do {
                temp->start = edge->start;
                temp->inverse->end = edge->start;
                temp = temp->next;
            } while (temp != graph->edges[edge->start]);
            graph->vertexcount--;
        } else {
            graph->edges[edge->start] = edge->next;
        }
    }
}

void remove_edge(Edge* edge, Graph* graph){
    Edge* temp;
    remove_directed_edge(edge, graph);
    remove_directed_edge(edge->inverse, graph);

    for(int i=1; i <= graph->vertexcount; i++) {
        temp = graph->edges[i];
        do {
            if(temp->unique_number == graph->edgecount - 2){
                temp->unique_number = edge->unique_number;
            } else if(temp->unique_number == graph->edgecount-1){
                temp->unique_number = edge->inverse->unique_number;
            }
            temp = temp->next;
        } while (temp != graph->edges[i]);
    }

    return_edge(edge->inverse, graph);
    return_edge(edge, graph);
}

void add_1vertex_in_face(Edge* edge, Graph* graph){
    set_vertexcount(graph, graph->vertexcount + 1);

    Edge* last_inverse=NULL;
    Edge* curr = edge;
    Edge* new_edge;
    Edge* new_inverse;
    for(int i=0; i<4; i++){
        new_edge = get_edge(curr->start, graph);
        insert_after(curr->prev, new_edge);
        new_inverse = get_edge(graph->vertexcount, graph);
        make_inverse(new_edge, new_inverse);
        new_edge->type = 3;//type is not determined yet
        new_inverse->type = 3;

        if(last_inverse!=NULL){
            insert_after(last_inverse->prev, new_inverse);
        } else {
            new_inverse->prev = new_inverse;
            new_inverse->next = new_inverse;
        }
        last_inverse = new_inverse;
        curr = next_in_face(curr);
    }
    graph->edges[graph->vertexcount] = new_inverse;
}


void continue_BFS(Queue* queue, Edge* edge){
    Edge* temp = edge->inverse->next;
    while (temp != edge->inverse) {
        queue_add(queue, temp);
        temp = temp->next;
    }
}

//returns 1 it the target is found on the inside, 0 if it is not.
int BFS_inside(Edge** marks, Edge* start_edge, uint16_t vi_target, Queue* queue){
    //Setup
    Edge* edge = marks[start_edge->start]->inverse;
    Edge* edge2 = edge;
    queue_clear(queue);
    do {
        edge = edge->next;
    } while(edge->end != start_edge->end);
    do {
        edge2 = edge2->prev;
    } while(edge2->end != start_edge->end);
    //edge and edge2 now form a 2-cycle such that there are no more edges between their vertices on the side of v0
    edge = edge->next;
    while(edge!=edge2) {//all the edges on the inside with start edge->start are added to the queue
        queue_add(queue, edge);
        edge =edge->next;
    }

    //BFS to find vi_target
    while((edge = queue_pop(queue))){
        if(!marks[edge->end]){
            marks[edge->end] = edge; //mark the edge
            if(edge->end == vi_target){//found target
                return 1;
            } else {//continue search
                continue_BFS(queue, edge);
            }
        }
    }
    return 0; //the target was not found
}

void save_cutpath(Edge** marks, Lopsp* lopsp, uint16_t vi){
    Edge* temp = marks[vi];
    int number = vi==lopsp->v1? 1: -2;
    while(temp->start != lopsp->v0){
        lopsp->edge_left_or_right_of_P[temp->unique_number] = number;
        lopsp->edge_left_or_right_of_P[temp->inverse->unique_number] = -number;
        temp = marks[temp->start];
    }
    lopsp->edge_left_or_right_of_P[temp->unique_number] = number;
    lopsp->edge_left_or_right_of_P[temp->inverse->unique_number] = -number;
}

void fill_queue_from_v0(Queue* queue, Lopsp* lopsp){
    Edge* temp = lopsp->graph->edges[lopsp->v0];
    do {
        queue_add(queue, temp);
        temp = temp->next;
    } while(temp != lopsp->graph->edges[lopsp->v0]);
}

//checks if there is another edge with the same start and end as edge
int is_double_edge(Edge* edge){
    Edge* temp = edge;
    do {
        temp = temp->next;
    } while(temp->end != edge->end);
    return temp != edge;
}

//Finds a cut-path and saves it in lopsp->edge_left_or_tight_of_P.
//We use a modified version of BFS (important that it is not DFS for correctness!) to find the paths. When we have
//not found a path to v2 or v1 we only use edges of type 2. Once one of them has been found we can use all edges.
//Properties of paths of type 2 guarantee that a path will then be found. We must be careful when the last edge
//of the first path is not of type 2. In that case we must check if there are multiple edges between those vertices
//and if so, check that the other vi is not on the wrong side of that 2-cycle.
void cut(Lopsp* lopsp){
    Edge* edge;
    uint16_t vi_found=0; //0 if none found, lopsp->v1 if v1 found, lopsp->v2 if v2 found
    uint16_t vi_not_found=0;

    Edge** marks = calloc(lopsp->graph->vertexcount +1, sizeof(Edge*)); //if a vertex is visited by BFS the edge from which it was marked is saved here. If it is not marked then NULL
    marks[lopsp->v0] = lopsp->graph->edges[lopsp->v0]; //conceptually this does not make sense, but v0 must be marked with something so that it cannot be chosen
    lopsp->edge_left_or_right_of_P = calloc( lopsp->graph->edgecount, sizeof(int));

    Queue* edge_queue = new_queue(lopsp->graph->edgecount);//every edge can be in the queue at most once
    fill_queue_from_v0(edge_queue, lopsp);

    //First we find a type 2 path to v1 OR v2
    while(!vi_found && (edge = queue_pop(edge_queue))) {
        if(marks[edge->end] == NULL) {
            if (edge->end == lopsp->v1 || edge->end == lopsp->v2) {//we encounter v1 or v2
                marks[edge->end] = edge;
                vi_found = edge->end;
                vi_not_found = (vi_found == lopsp->v1) ? lopsp->v2 : lopsp->v1;
                if(edge->type != 2 && edge->start != lopsp->v0){//check if there are multiple edges, otherwise that is not necessary
                    //There are multiple edges between edge->start and edge->end and vi_not_found is on the inside of the 2-cycle
                    if(is_double_edge(edge) && BFS_inside(marks, edge, vi_not_found, edge_queue)){ //the path to vi_found is marked in BFS_new
                        vi_found = vi_not_found;
                        vi_not_found = edge->end;
                    }//In all other cases edge can be used to mark vi_found
                }
            } else if ((edge->start == lopsp->v0 || edge->type == 2)){//we only consider this edge if it is of type 2 or starts in v0
                marks[edge->end] = edge; //mark the edge
                continue_BFS(edge_queue, edge);
            }
        } //edges to marked vertices are ignored
    }

    //vi_found should be true here
    save_cutpath(marks, lopsp, vi_found);
    for (int i = 1; i <= lopsp->graph->vertexcount; ++i) {//set the marks of the vertices not in the found cutpath back to NULL
        if(marks[i]!=NULL && lopsp->edge_left_or_right_of_P[marks[i]->unique_number] ==0 && i!=lopsp->v0){
            marks[i]=NULL;
        }
    }
    queue_clear(edge_queue);
    fill_queue_from_v0(edge_queue, lopsp);

    //find a path to vi_not_found
    while((edge = queue_pop(edge_queue))){
        if(!marks[edge->end]){
            marks[edge->end] = edge; //mark the edge
            if(edge->end == vi_not_found){//found target
                save_cutpath(marks, lopsp, vi_not_found);
                free(marks);
                lopsp->cut = 1;
                free_queue(edge_queue);
                return;//This point should always be reached, else there is no cut-path found
            } else {//continue search
                continue_BFS(edge_queue, edge);
            }
        }
    }
    fprintf(stderr, "Cut-path not found. This should not be possible\n");
    exit(1);
}

//This function makes graph->edgecount copies of the edges in lopsp and glues them together along the cut-path. Vertex-numbers are not determined yet!!!
void copy_and_glue(Lopsp* lopsp, Graph* graph, Graph* result, Edge** new_edges){
    Edge* dc;

    int edgecount_O = lopsp->graph->edgecount;
    int double_chamber_invnext;//The number of the double chamber (Edge) that must contain the inverse and next of the current edge

    for (int i = 0; i < edgecount_O * graph->edgecount; i++) {
        new_edges[i] = get_edge(0, result);
    }

    //In this for-loop the copies are glued.
    for (int i = 1; i <= graph->vertexcount; i++) {
        dc = graph->edges[i];
        do {
            for (int j = 0; j < (edgecount_O); j++) {
                switch (lopsp->edge_left_or_right_of_P[j]) {
                    case 0://edge j is not on P
                        double_chamber_invnext = dc->unique_number;//If the edge is not on P then the next edge is in the same double chamber
                        break;
                    case 1://edge j is on the copy of P_1 on the left side of O_P
                    case -1://edge j is on the copy of P_1 on the right side of O_P
                        double_chamber_invnext = dc->inverse->unique_number;
                        break;
                    case 2://edge j is on the copy of P_2 on the left side of O_P
                        double_chamber_invnext = dc->prev->inverse->unique_number;
                        break;
                    case -2://edge j is on the copy of P_2 on the right side of O_P
                        double_chamber_invnext = dc->inverse->next->unique_number;
                        break;
                    default:
                        fprintf(stderr, "edge_left_or_right_of_P must have value 0, -2, -1, 1 or 2, not %d\n", lopsp->edge_left_or_right_of_P[j]);
                        exit(1);
                }

                //next of this edge and prev of the next edge are set here
                make_next(new_edges[edgecount_O*dc->unique_number + j], new_edges[edgecount_O*double_chamber_invnext + lopsp->edge_array[j]->next->unique_number]);

                new_edges[edgecount_O*dc->unique_number + j]->inverse = new_edges[edgecount_O*double_chamber_invnext + lopsp->edge_array[j]->inverse->unique_number];
                new_edges[edgecount_O*dc->unique_number + j]->type = lopsp->edge_array[j]->type;
                new_edges[edgecount_O*dc->unique_number + j]->unique_number = edgecount_O*dc->unique_number + j;
            }
            dc = dc->next;
        } while (dc != graph->edges[i]);
    }
}

void remove_type2(Graph* graph, Edge** new_edges){
    Edge* temp;
    int edge_count;
    int new_vertexcount = 0;
    int new_edgenumber = 0;
    for (int i = 0; i < graph->edgecount; ++i) {
        if(new_edges[i]->type == 2){ //other edges will be removed later
            if(new_edges[i]->next->type == 0){
                new_edges[i]->type = 3; //will be removed in next loop as it no longer has type 2
            } else { //The edge will be preserved
                if (new_edges[i]->start == 0) {//The start of the edge has not been set yet
                    new_vertexcount++;
                    temp = new_edges[i];
                    do { // the start of all the preserved vertices with the same start is set
                        temp->start = new_vertexcount;
                        temp = temp->next->next;
                    } while (temp != new_edges[i]);
                }
                //The inverse and prev and next of the edge must also be set correctly
                new_edges[i]->inverse = new_edges[i]->inverse->next->next->inverse;
                new_edges[i]->next = new_edges[i]->next->next;
                new_edges[i]->prev = new_edges[i]->prev->prev;
            }
        }
    }

    //set the new vertexcount
    set_vertexcount(graph, new_vertexcount);
    edge_count = graph->edgecount; //must be a separate variable because graph->edgecount changes in return_edge

    //now remove the unnecessary edges and set 'end' and 'unique_number' for each remaining edge
    for (int i = 0; i < edge_count; ++i) {
        if(new_edges[i]->type == 2){
            new_edges[i]->end = new_edges[i]->inverse->start;
            new_edges[i]->unique_number = new_edgenumber;
            new_edgenumber++;
            graph->edges[new_edges[i]->start] = new_edges[i];//is overwritten sometimes but that is no problem. The last edge with this start will be saved
        } else {
            return_edge(new_edges[i], graph);
            //all the references to this edge are gone or will be deleted in this loop
        }
    }
    //now edge_count = new_edgenumber + number of edges returned and graph->edgecount = edge_count - number of edges returned = new_edgenumber
}

//Applies a lopsp-operation to the graph
void apply(Lopsp* lopsp, Graph* graph, Graph* result){
    if(!lopsp->cut){
        cut(lopsp);
    }
    init_graph(result);

    Edge* new_edges[(lopsp->graph->edgecount * graph->edgecount)]; // for every edgenumber we save the new edge with this number

    //make all the new edges and glue them together
    copy_and_glue(lopsp, graph, result, new_edges);

    //remove the edges of the barycentric subdivision and name the vertices
    remove_type2(result, new_edges);
}

int has_loop(Graph* graph){
    Edge* curr;
    for (int i = 1; i <= graph->vertexcount; i++){
        curr = graph->edges[i];
        do {
            if(curr->end == i){
                return 1;
            }
            curr = curr->next;
        } while (curr != graph->edges[i]);
    }
    return 0;
}

int has_multiple_edges(Graph* graph){
    Edge* curr;
    int neighbours[graph->edgecount];//the neighbours we have already read
    int neighbourcount;//the number of neighbours already read
    for (int i = 1; i <= graph->vertexcount; i++){
        curr = graph->edges[i];
        neighbourcount = 0;
        do {
            for (int j = 0; j < neighbourcount; ++j) {
                if(curr->end == neighbours[j] && curr->end != i){
                    return 1;
                }
            }
            neighbours[neighbourcount] = curr->end;
            neighbourcount++;
            curr = curr->next;
        } while (curr != graph->edges[i]);
    }
    return 0;
}

int edge_in_face(Edge* faceedge, Edge* searchedge){
    Edge* edge = faceedge;
    do {
        if (edge == searchedge){
            return 1;
        }
        edge = edge->inverse->next;
    } while (edge != faceedge);
    return 0;
}

int adjacent_in_faces(Edge* fedge, Edge* gedge, int v, int w){
    Edge* edge = fedge;
    do {
        if (edge->end == w && edge->start == v || edge->end == v && edge->start == w){
            return edge_in_face(gedge, edge->inverse);
        }
        edge = edge->inverse->next;
    } while (edge != fedge);
    return 0;
}

//Checks if the given graph is c1, c2 or c3. If it is ck for a higher k, it returns 3.
int ck(Graph* graph){
    int vertices_in_face[2 + (graph->edgecount / 2) - graph->vertexcount][graph->vertexcount + 1];
    //vertices_in_face[f][v] is true if v is in face f
    int facecount = 0;
    Edge* edge;
    Edge* face_edge;

    for (int i = 0; i < 2 + (graph->edgecount / 2) - graph->vertexcount; ++i) {
        for (int j = 0; j <= graph->vertexcount; ++j) {
            vertices_in_face[i][j] = 0;
        }
    }

    int edge_seen[graph->edgecount];
    for (int i = 0; i < graph->edgecount; ++i) {edge_seen[i] = 0;}
    Edge* edge_in_face[2 + (graph->edgecount / 2) - graph->vertexcount];//holds one edge of every face

    //fill vertices in face such that vertices_in_face[i][j] is the number of times vertex j appears in face i
    for (int v = 1; v <= graph->vertexcount; ++v) {
        edge = graph->edges[v];
        do {
            if(!edge_seen[edge->unique_number]){ //loop over face
                edge_in_face[facecount] = edge;
                face_edge = edge;
                do {
                    if (vertices_in_face[facecount][face_edge->start] == 1){
                        return 1; //same vertex in face multiple times, not c2. Loops are also detected here
                    }
                    vertices_in_face[facecount][face_edge->start] = 1;
                    edge_seen[face_edge->unique_number] = 1;
                    face_edge = face_edge->inverse->next;
                } while (face_edge != edge);
                facecount++;
            }
            edge = edge->next;
        } while (edge != graph->edges[v]);
    }

    if(has_multiple_edges(graph)){
        return 2;
    }
    int shared[3];
    int numbershared;

    for (int f = 0; f < facecount; ++f) {
        for (int g = f + 1; g < facecount; ++g) {
            numbershared = 0;
            for (int v = 1; v <= graph->vertexcount; ++v) {
                if(vertices_in_face[f][v] && vertices_in_face[g][v]){
                    shared[numbershared] = v;
                    numbershared++;
                    if(numbershared == 2){
                        if (!(adjacent_in_faces(edge_in_face[f], edge_in_face[g], shared[0], shared[1]))){
                            return 2;//found two non-adjacent vertices in intersection of two faces
                        }
                    } else if (numbershared > 2){
                        return 2;
                    }
                }
            }
        }
    }
    return 3;
}