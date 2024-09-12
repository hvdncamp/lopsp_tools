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
 * This program converts an embedded graph given in planarcode to a readable txt format.
 * The planarcode is read from stdin and the result is printed to stdout
 */

#include <stdlib.h>
#include"lib/lopsp_functions.h"
#include"lib/graph_io.h"

#define INITIAL_GRAPH_SIZE 50

void print_planarcode(FILE* file, Graph* graph){
    Edge* curr;

    for (int i = 1; i <= graph->vertexcount; i++){//for every vertex, write the incident vertices
        fprintf(file, "%d: ", i);
        curr = graph->edges[i];
        do {
            fprintf(file, "%d ", curr->end);
            curr = curr->next;
        } while (curr != graph->edges[i]);
        fprintf(file, "\n");
    }
}

int main(int argc,char *argv[]) {
    Graph* graph = new_graph(INITIAL_GRAPH_SIZE);

    init_edges();

    enum code encoding = checkheader(stdin);
    if(encoding == EDGE_CODE){
        fprintf(stderr, "This program only reads planarcode, not edgecode. \n");
        exit(1);
    }

    int graphcount = 0;

    while(read_graph(stdin, graph, encoding)){
        graphcount++;
        fprintf(stdout, "\nGraph %d\n", graphcount);
        print_planarcode(stdout, graph);
    }

    free_graph(graph);

    return(0);
}
