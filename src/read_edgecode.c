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

/*
 * This program converts an embedded graph given in edge code to a readable txt format.
 * The edge code is read from stdin and the result is printed to stdout
 */

#include <stdlib.h>
#include"lib/lopsp_functions.h"
#include"lib/graph_io.h"

int read_and_writeedgecode(FILE* file, int count){
    uint8_t read_byte;
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
            fprintf(stderr, "Something went wrong reading graph %d. Are you sure it is encoded correctly? \n", count);
            exit(1);
        }
        numbersize = read_byte & ((uint8_t)(~0) >> 4);
        totalsize = read_edge_number(file, read_byte >> 4);
    } else {// 1 byte per number
        numbersize = 1;
        totalsize = read_byte;
    }

    fprintf(stdout, "\nGraph %d\n", count);

    curr_vertex = 0;
    total_bytes_read = 0;
    while(total_bytes_read < totalsize){
        curr_vertex++;//We start reading the neighbours of the next vertex
        fprintf(stdout, "%d:", curr_vertex);

        curr_edge = read_edge_number(file, numbersize);

        total_bytes_read += numbersize; // every time a number is read total_bytes_read must be incremented
        while (curr_edge != UINT64_MAX) { //loop over all the neighbours of one vertex
            fprintf(stdout, "%lu ", curr_edge);
            if(total_bytes_read < totalsize){//we have not read the last byte
                curr_edge = read_edge_number(file, numbersize);
                total_bytes_read += numbersize;
            } else{//the last group of edgenumbers is not ended by 255.
                curr_edge = UINT64_MAX;
            }
        }
        fprintf(stdout, "\n");
    }
    return 1;
}

int main(int argc,char *argv[]) {
    if(checkheader(stdin) != EDGE_CODE){
        fprintf(stderr, "The encoding of the graphs should start with the header >>edge_code<< \n");
        exit(1);
    }

    int graphcount = 1;

    while(read_and_writeedgecode(stdin, graphcount)){
        graphcount++;
    }

    return(0);
}
