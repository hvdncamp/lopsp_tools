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

#include "graph_io.h"
#include<stdlib.h>
#include<string.h>

//Reads the following number_of_bytes bytes from stdin as a big endian number. Returns UINT64_MAX if the first byte is 255
unsigned long read_edge_number(FILE* file, int number_of_bytes){
    uint8_t byte;
    unsigned long to_add;
    unsigned long number = 0;
    if(number_of_bytes > 8){
        fprintf(stderr, "This program does not support codes with more than 8 bytes per number\n");
        exit(1);
    }
    for (int i = 1; i <= number_of_bytes; ++i) {
        if(fread(&byte, 1, 1, file) != 1){
            return 0; //end of file reached
        }
        if(i==1 && byte == 255){
            return UINT64_MAX;
        }
        to_add = byte;
        number = number | (to_add << (8*(number_of_bytes - i)));
    }
    return number;
}

//reads a big endian number with number_of_bytes bytes from file
unsigned long read_big_endian(FILE* file, int number_of_bytes){
    uint8_t byte;
    unsigned long to_add;
    unsigned long number = 0;
    if(number_of_bytes > 8){
        fprintf(stderr, "This program does not support codes with more than 8 bytes per number\n");
        exit(1);
    }
    for (int i = 1; i <= number_of_bytes; ++i) {
        if(fread(&byte, 1, 1, file) != 1){
            return 0; //end of file reached
        }
        to_add = byte;
        number = number | (to_add << (8 * (number_of_bytes - i)));
    }
    return number;
}

//reads a little endian number with number_of_bytes bytes from file
unsigned long read_little_endian(FILE* file, int number_of_bytes){
    uint8_t byte;
    unsigned long to_add;
    unsigned long number = 0;
    if(number_of_bytes > 8){
        fprintf(stderr, "This program does not support codes with more than 8 bytes per number\n");
        exit(1);
    }
    for (int i = 0; i < number_of_bytes; ++i) {
        if(fread(&byte, 1, 1, file) != 1){
            return 0; //end of file reached
        }
        to_add = byte;
        number = number | (to_add << (8 * i));
    }
    return number;
}

void write_big_endian_number(FILE* file, int number_of_bytes, unsigned long number){
    uint8_t byte;
    for (int i = 0; i < number_of_bytes; ++i) {
        byte = (number & (~((uint8_t) 0) >> (i * 8))) >> ((number_of_bytes - i - 1) * 8);
        fwrite(&(byte), 1, 1,file);
    }
}

void write_little_endian_number(FILE* file, int number_of_bytes, unsigned long number){
    uint8_t byte;
    for (int i = number_of_bytes - 1; i >= 0; --i) {
        byte = (number & (~((uint8_t) 0) >> (i * 8))) >> ((number_of_bytes - i - 1) * 8);
        fwrite(&(byte), 1, 1,file);
    }
}

void save_structure(Lopsp* lopsp){
    Edge* temp;
    if(lopsp->edge_array != NULL){
        free(lopsp->edge_array);
    }
    lopsp->edge_array = malloc(sizeof(Edge*) * lopsp->graph->edgecount);

    for (int i = 1; i <= lopsp->graph->vertexcount ; ++i) {
        temp = lopsp->graph->edges[i];
        do {
            lopsp->edge_array[temp->unique_number] = temp;
            temp = temp->next;
        } while(temp !=lopsp->graph->edges[i]);
    }
}

void determine_types(Lopsp* lopsp, uint16_t typev0){
    Queue* edgequeue = new_queue(1000);
    //the queue will only contain edges of types 0 and 2. The first edge added is incident with v0.
    Edge* current_edge = lopsp->graph->edges[lopsp->v0]->type == 1? lopsp->graph->edges[lopsp->v0]->next : lopsp->graph->edges[lopsp->v0];
    current_edge->type = 2 - typev0;
    queue_add(edgequeue, current_edge);

    while ((current_edge = queue_pop(edgequeue))){
        //if the inverse has no type yet give it its type and add it to the queue
        if (current_edge->inverse->type == 3){
            current_edge->inverse->type = current_edge->type;
            queue_add(edgequeue, current_edge->inverse);
        }
        //if the next has no type yet give it its type and add it to the queue. This is only possible if the startvertex is type 1
        if (current_edge->next->type == 3){
            current_edge->next->type = 2 - current_edge->type;
            queue_add(edgequeue, current_edge->next);
        }
        //if the next->next has no type yet give it its type and add it to the queue. This is only possible if the startvertex is type 0 or 2
        else if (current_edge->next->next->type == 3){
            current_edge->next->next->type = current_edge->type;
            queue_add(edgequeue, current_edge->next->next);
        }
    }
    free_queue(edgequeue);
}

int read_lopsp(FILE* file, Lopsp* lopsp){
    uint16_t typev0;
    uint16_t v1type1;
    uint16_t number_of_vertices_read;

    init_lopsp(lopsp);

    if(fread(&number_of_vertices_read, 2, 1,file) != 1){
        return 0;
    }

    if(feof(file)){//there are no more lopsp-operations to be read
        return 0;
    }

    set_vertexcount(lopsp->graph, number_of_vertices_read);

    if(fread(&lopsp->v0, 2,1, file) != 1){
        fprintf(stderr, "A lopsp-operation could not be read. Are you sure all operations are encoded correctly?\n");
        exit(1);
    }
    if(fread(&lopsp->v1, 2,1, file) != 1){
        fprintf(stderr, "A lopsp-operation could not be read. Are you sure all operations are encoded correctly?\n");
        exit(1);
    }
    if(fread(&lopsp->v2, 2,1, file) != 1){
        fprintf(stderr, "A lopsp-operation could not be read. Are you sure all operations are encoded correctly?\n");
        exit(1);
    }
    if(fread(&typev0, 2,1, file) != 1){
        fprintf(stderr, "A lopsp-operation could not be read. Are you sure all operations are encoded correctly?\n");
        exit(1);
    }
    if(fread(&v1type1, 2,1, file) != 1){
        fprintf(stderr, "A lopsp-operation could not be read. Are you sure all operations are encoded correctly?\n");
        exit(1);
    }

    uint16_t curr_neighbour;
    Edge* new_edge;
    Edge* last_edge=NULL;//the last edge added
    Edge* edge_numbers[4*lopsp->graph->vertexcount+1];// at position i this contains an edge with number i in the .lopsp
    for (int i = 0; i <= 4*lopsp->graph->vertexcount; i++)
    {
        edge_numbers[i] = NULL;
    }
    int face_size;

    //read the encoded graph
    for (int i = 1; i <= lopsp->graph->vertexcount; i++)
    {
        if(fread(&curr_neighbour, 2,1, file) != 1){
            fprintf(stderr, "A lopsp-operation could not be read. Are you sure all operations are encoded correctly?\n");
            exit(1);
        }
        if(curr_neighbour == 0){
            fprintf(stderr, "A vertex in the operation has no neighbours. This is not possible.\n");
            exit(1);
        }
        while (curr_neighbour !=0){
            new_edge = get_edge(i, lopsp->graph);
            new_edge->type = 1;
            if(last_edge!=NULL){
                new_edge->prev = last_edge;
                last_edge->next = new_edge;
            } else {
                lopsp->graph->edges[i] = new_edge;
            }
            if(edge_numbers[curr_neighbour]!=NULL){
                if(edge_numbers[curr_neighbour]->inverse==NULL){//should always be the case
                    make_inverse(edge_numbers[curr_neighbour], new_edge);
                }
            } else {
                edge_numbers[curr_neighbour] = new_edge;
                new_edge->inverse = NULL;
            }
            last_edge = new_edge;
            if(fread(&curr_neighbour, 2,1, file) != 1){
                fprintf(stderr, "A lopsp-operation could not be read. Are you sure all operations are encoded correctly?\n");
                exit(1);
            }
        }
        lopsp->graph->edges[i]->prev = last_edge;
        last_edge->next = lopsp->graph->edges[i];
        last_edge=NULL;
    }

    Edge* v1edge;
    //determine the type of v1 and ensure that all edges of type 1 have type 1 and the others have type 3 (will be changed later)
    if(v1type1){
        //An extra edge must be added
        v1edge = lopsp->graph->edges[lopsp->v1];
        new_edge = get_edge(lopsp->v1, lopsp->graph);
        make_inverse(new_edge, get_edge(v1edge->inverse->next->end,lopsp->graph));
        insert_after(v1edge, new_edge);
        insert_after(v1edge->inverse->next->inverse, new_edge->inverse);
        edge_numbers[(lopsp->graph->edgecount / 2)] = new_edge;

        new_edge->type = 3;
        new_edge->inverse->type = 3;
        v1edge->type = 3;
        v1edge->inverse->type = 3;
    }

    //Add vertices of type 1 (except v1)
    Edge* start;
    int edge_count= lopsp->graph->edgecount / 2; //edge_count is needed because lopsp->graph->edgecount changes during this loop

    for (int i = 1; i <= 2*edge_count; i++){ //edge_count counts normal edges, not directed edges. We must check the faces for all 2*edge_count directed edges
        if(i <= edge_count){
            start = edge_numbers[i];
        } else {
            start = edge_numbers[i-edge_count]->inverse;
        }

        face_size = facesize(start);
        if(face_size==4){
            add_1vertex_in_face(start, lopsp->graph);
        } else if (face_size !=3){
            fprintf(stderr, "This should not be possible. Every face in the predecoration must have size 4, not %d. Are you sure the given .lopsp code is correct?\n", face_size);
            exit(1);
        }
    }

    determine_types(lopsp, typev0);

    save_structure(lopsp);
    return 1;
}

//checks the header. Returns 'p' if the file is planarcode and 'e' if it is edgecode
enum code checkheader(FILE* file){
    enum code input_format;
    char buffer[10];
    if(fgets(buffer, 3, file) == NULL){
        fprintf(stderr,"The encoding of the graphs should start with a header. It can be >>planar_code<<, >>planar_code le<<, >>planar_code be<< or >>edge_code<< \n");
        exit(1);
    }
    if(strcmp(buffer, ">>")!=0){
        fprintf(stderr,"The encoding of the graphs should start with a header. It can be >>planar_code<<, >>planar_code le<<, >>planar_code be<< or >>edge_code<< \n");
        exit(1);
    }
    if(fgets(buffer, 5, file) == NULL){
        fprintf(stderr,"The encoding of the graphs should start with a header. It can be >>planar_code<<, >>planar_code le<<, >>planar_code be<< or >>edge_code<< \n");
        exit(1);
    }

    if(strcmp(buffer, "plan")==0){
        if(fgets(buffer, 3, file) == NULL){
            fprintf(stderr,"The encoding of the graphs should start with a header. It can be >>planar_code<<, >>planar_code le<<, >>planar_code be<< or >>edge_code<< \n");
            exit(1);
        }
        if(strcmp(buffer, "ar")!=0){
            fprintf(stderr,"The encoding of the graphs should start with a header. It can be >>planar_code<<, >>planar_code le<<, >>planar_code be<< or >>edge_code<< \n");
            exit(1);
        }
        input_format = PLANAR_CODE_LITTLE_ENDIAN;
    } else if (strcmp(buffer, "edge")==0){
        input_format = EDGE_CODE;
    } else {
        fprintf(stderr,"The encoding of the graphs should start with a header. It can be >>planar_code<<, >>planar_code le<<, >>planar_code be<< or >>edge_code<< \n");
        exit(1);
    }

    if(fgets(buffer, 6, file) == NULL){
        fprintf(stderr,"The encoding of the graphs should start with a header. It can be >>planar_code<<, >>planar_code le<<, >>planar_code be<< or >>edge_code<< \n");
        exit(1);
    }
    if(strcmp(buffer, "_code")!=0){
        fprintf(stderr,"The encoding of the graphs should start with a header. It can be >>planar_code<<, >>planar_code le<<, >>planar_code be<< or >>edge_code<< \n");
        exit(1);
    }

    if(fgets(buffer, 3, file) == NULL){
        fprintf(stderr,"The encoding of the graphs should start with a header. It can be >>planar_code<<, >>planar_code le<<, >>planar_code be<< or >>edge_code<< \n");
        exit(1);
    }
    if(strcmp(buffer, "<<")!=0){
        if(buffer[0] != ' '){
            fprintf(stderr,"The encoding of the graphs should start with a header. It can be >>planar_code<<, >>planar_code le<<, >>planar_code be<< or >>edge_code<< \n");
            exit(1);
        }
        if(buffer[1] == 'b'){//big endian
            input_format = PLANAR_CODE_LITTLE_ENDIAN;
        } else if (buffer[1] != 'l'){
            fprintf(stderr,"The encoding of the graphs should start with a header. It can be >>planar_code<<, >>planar_code le<<, >>planar_code be<< or >>edge_code<< \n");
            exit(1);
        }
        if(fgets(buffer, 4, file) == NULL){
            fprintf(stderr,"The encoding of the graphs should start with a header. It can be >>planar_code<<, >>planar_code le<<, >>planar_code be<< or >>edge_code<< \n");
            exit(1);
        }
        if(strcmp(buffer, "e<<")!=0){
            fprintf(stderr,"The encoding of the graphs should start with a header. It can be >>planar_code<<, >>planar_code le<<, >>planar_code be<< or >>edge_code<< \n");
            exit(1);
        }
    }
    return input_format;
}

void check_lopsp_header(FILE* file){
    char buffer[10];
    if(fgets(buffer, 10, file) == NULL){
        fprintf(stderr,"The encoding of the graphs should start with a header. It can be >>planar_code<<, >>planar_code le<<, >>planar_code be<< or >>edge_code<< \n");
        exit(1);
    }
    if(strcmp(buffer, ">>lopsp<<")!=0){
        fprintf(stderr,"The encoding of the operations should start with the header >>lopsp<< \n");
        exit(1);
    }
}

//'file' must start at the start of the planarcode (no header). This method reads the planarcode in
// the file and saves the corresponding graph in 'graph'.
int read_planarcode(FILE* file, Graph* graph, int little_endian){
    unsigned long read;
    Edge* prev_edge;
    Edge* new_edge;
    Edge* temp;
    int numbersize = 1;

    read = read_little_endian(file, 1);

    if (feof(file)){//End of file reached. No graph read
        return 0;
    }
    if(read == 0){
        numbersize = 2;
        read = little_endian ? read_little_endian(file, numbersize) : read_big_endian(file, numbersize);
    }

    set_vertexcount(graph, (int)read);

    for (int i = 1; i <= graph->vertexcount; i++){
        //read neighbours
        read = little_endian ? read_little_endian(file, numbersize) : read_big_endian(file, numbersize);
        prev_edge = NULL;
        while (read !=0){
            new_edge = get_edge(i, graph);
            new_edge->end = read;
            if(prev_edge==NULL){
                graph->edges[i]=new_edge;
            } else {
                make_next(prev_edge, new_edge);
            }
            if(i>read){//the inverse of this edge already exists.
                temp = graph->edges[read];
                while(temp->end != i){
                    temp = temp->next;
                }
                make_inverse(temp, new_edge);
            }
            prev_edge = new_edge;
            read = little_endian ? read_little_endian(file, numbersize) : read_big_endian(file, numbersize);
        }
        make_next(new_edge, graph->edges[i]);
    }
    return 1;
}

//'file' must start at the start of the edgecode (no header). This method reads the edgecode in the file
// and saves the corresponding graph in 'graph'.
int read_edgecode(FILE* file, Graph* graph){
    uint8_t read_byte;
    Edge* prev_edge;
    Edge* new_edge;
    int numbersize; // the number of bytes used for each edge-number
    unsigned long totalsize; //the total number of bytes in the code
    int total_bytes_read; // the number of bytes that has already been read
    int curr_vertex; //the number of the current vertex
    unsigned long curr_edge; // the number of the current edge

    if(fread(&read_byte, 1, 1,file) != 1) {//read the number of vertices
        return 0;
    }
    if (feof(file)){// we have reached the end of the file, there is no graph to be read
        return 0;
    }
    if(read_byte == 0){//maybe more than one byte per number
        if(fread(&read_byte, 1, 1,file) != 1) {//read the number of bytes per number and length of following number
            fprintf(stderr,"The edgecode could not be read. Are you sure it is encoded correctly?\n");
            exit(1);
        }
        numbersize = read_byte & ((uint8_t)(~0) >> 4);
        totalsize = read_edge_number(file, read_byte >> 4);
    } else {// 1 byte per number
        numbersize = 1;
        totalsize = read_byte;
    }

    Edge* edge_numbers[totalsize];
    for(int i=0; i<totalsize; i++){
        edge_numbers[i]=NULL;
    }

    curr_vertex = 0;
    total_bytes_read = 0;
    while(total_bytes_read < totalsize){
        curr_vertex++;//We start reading the neighbours of the next vertex
        set_vertexcount(graph, curr_vertex);

        curr_edge = read_edge_number(file, numbersize);
        total_bytes_read += numbersize; // every time a number is read total_bytes_read must be incremented
        prev_edge=NULL;
        while (curr_edge != UINT64_MAX) { //loop over all the neighbours of one vertex
            new_edge = get_edge(curr_vertex, graph);
            //for fixing prev and next
            if (prev_edge == NULL) {
                graph->edges[curr_vertex] = new_edge;
            } else {
                make_next(prev_edge, new_edge);
            }

            //for fixing the inverse
            if (edge_numbers[curr_edge] != NULL) {//we have already seen the other half of this edge
                make_inverse(edge_numbers[curr_edge], new_edge);
            } else {
                edge_numbers[curr_edge] = new_edge;
            }
            prev_edge = new_edge;

            if(total_bytes_read < totalsize){//we have not read the last byte
                curr_edge = read_edge_number(file, numbersize);
                total_bytes_read += numbersize;
            } else{//the last group of edgenumbers is not ended by 255.
                curr_edge = UINT64_MAX;
            }
        }
        make_next(new_edge, graph->edges[curr_vertex]);
    }
    return 1;
}

//Reads ONE embedded graph from the given file, either planarcode OR edgecode
int read_graph(FILE* file, Graph* graph, enum code format){
    init_graph(graph);
    switch (format){
        case EDGE_CODE://The graph is encoded in edgecode
            return read_edgecode(file, graph);
        case PLANAR_CODE_BIG_ENDIAN://the graph is encoded in planarcode
            return read_planarcode(file, graph, 0);
        case PLANAR_CODE_LITTLE_ENDIAN://the graph is encoded in planarcode
            return read_planarcode(file, graph, 1);
        default:
            fprintf(stderr, "This should not be possible. The format parameter in read_graph should be PLANAR_CODE_BIG_ENDIAN, PLANAR_CODE_LITTLE_ENDIAN or EDGE_CODE\n");
            exit(1);
    }
}

void write_planarcode(FILE* file, Graph* graph){
    Edge* curr;
    uint8_t zero = 0;
    int numbersize = 1;

    if(graph->vertexcount > 255){
        if (graph->vertexcount > 65535){
            fprintf(stderr, "A result has %d vertices. This is too large for this program to print so it has been omitted from the output (however it is included in the counts). Perhaps it can be printed in edgecode.", graph->vertexcount);
            return;
        }
        numbersize = 2;
        fwrite(&zero, sizeof(zero), 1, file);
    }

    write_little_endian_number(file, numbersize, graph->vertexcount);

    for (int i = 1; i <= graph->vertexcount; i++){//for every vertex, write the incident vertices
        curr = graph->edges[i];
        do {
            write_little_endian_number(file, numbersize, curr->end);
            curr = curr->next;
        } while (curr != graph->edges[i]);
        write_little_endian_number(file, numbersize, 0);
    }
}

void write_edgecode(FILE* file, Graph* graph){
    unsigned long curr_edgenumber = 0;
    uint64_t edge_numbers[graph->edgecount];
    Edge* curr;
    uint8_t end_of_neighbours = 255;
    uint8_t zero = 0;
    int nr_of_bytes = 1;
    while(((unsigned long)(1) << (8 * nr_of_bytes)) <= graph->edgecount / 2){
        nr_of_bytes++;
    }
    if(nr_of_bytes > 8){
        fprintf(stderr, "A result has too many edges for this program to print so it has been omitted from the output but it is included in the counts.");
        return;
    }

    unsigned long total_bytes = graph->vertexcount + graph->edgecount * nr_of_bytes - 1;
    uint8_t header_byte_number = 1;
    uint8_t header_number;

    if(total_bytes < 256 && nr_of_bytes == 1){
        fwrite(&(total_bytes), 1, 1,file);//write the number of vertices
    } else {
        fwrite(&(zero), 1, 1, file); //write 0
        while(((unsigned long)(1) << (8 * header_byte_number)) <= total_bytes){
            header_byte_number++;
        }
        header_number = header_byte_number << 4; //header_byte_number is K
        header_number += nr_of_bytes; //nr_of_bytes is L
        fwrite(&(header_number), 1, 1,file); //write (K<<4)+L
        write_big_endian_number(file, header_byte_number, total_bytes); //write the total number of bytes
    }

    for(int i=0; i<graph->edgecount; i++){
        edge_numbers[i] = 0;
    }

    for (int i = 1; i <= graph->vertexcount; i++){//for every vertex, write its incident edges
        curr = graph->edges[i];
        do {
            if(edge_numbers[curr->unique_number] == 0){//This edge has not been given a number yet
                curr_edgenumber++;
                edge_numbers[curr->inverse->unique_number] = curr_edgenumber;//the inverse will get the same number
                write_big_endian_number(file, nr_of_bytes, curr_edgenumber);
            } else {//This edge already has a number
                write_big_endian_number(file, nr_of_bytes, edge_numbers[curr->unique_number]);
            }
            curr = curr->next;
        } while (curr != graph->edges[i]);
        if(i != graph->vertexcount){
            fwrite(&end_of_neighbours, sizeof(end_of_neighbours), 1, file);
        }
    }
}


void write_header(FILE* file, enum code code){
    switch (code) {
        case PLANAR_CODE_LITTLE_ENDIAN:
            fprintf(file, ">>planar_code<<");
            break;
        case EDGE_CODE:
            fprintf(file, ">>edge_code<<");
            break;
        default:
            fprintf(stderr, "code for write_header should always be PLANAR_CODE_LITTLE_ENDIAN or EDGE_CODE\n");
            exit(1);
    }
}

void write_graph(FILE* file, Graph* graph, enum code code){
    switch (code) {
        case PLANAR_CODE_LITTLE_ENDIAN:
            write_planarcode(file, graph);
            break;
        case EDGE_CODE:
            write_edgecode(file, graph);
            break;
        default:
            fprintf(stderr, "code for write_graph should always be PLANAR_CODE_LITTLE_ENDIAN or EDGE_CODE\n");
            exit(1);
    }
}

//prints lopsp to file in a readable txt format
void print_lopsp(FILE* file, Lopsp* lopsp){
    Edge* curr;
    fprintf(file, "v0: %d\n", lopsp->v0);
    fprintf(file, "v1: %d\n", lopsp->v1);
    fprintf(file, "v2: %d\n", lopsp->v2);

    for (int i = 1; i <= lopsp->graph->vertexcount; i++){
        curr = lopsp->graph->edges[i];

        //print vertex and type
        fprintf(file, "%d(%d):", i, 3 - curr->type - curr->next->type);

        do { //print neighbours
            fprintf(file, "%d ", curr->end);
            curr = curr->next;
        } while(curr != lopsp->graph->edges[i]);

        fprintf(file, "\n");
    }
}