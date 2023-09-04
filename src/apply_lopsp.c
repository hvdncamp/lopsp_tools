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

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#include"lib/lopsp_functions.h"
#include"lib/graph_io.h"

#define INITIAL_GRAPH_SIZE 50
#define INITIAL_LOPSP_SIZE 50

#define ONE_OPERATION 'o'
#define ONE_GRAPH 'g'

void print_help(){
    fprintf(stderr,"This program applies lopsp operations to embedded graphs. There are two ways to use it, each requires an option and an argument:\n"
                   "If the option -g is present the program applies at least one lopsp-operation, read from stdin, to one embedded graph. \n"
                   "The name of the file containing the graph is given as an argument of the program.\n"
                   "If the option -o is present the program applies one lopsp-operation to at least one embedded graph, read from stdin. \n"
                   "The name of the file containing the operation is given as an argument of the program.\n"
                   "All lopsp-operations should be given in .lopsp-format. You can use the program txt_to_lopsp to convert a lopsp-operation in human readable format to .lopsp-format. \n"
                   "All graphs should be in planar code or edge code. \n"
                   "Before the graphs there should be a header: >>planar_code<<, >>planar_code le<<, >>planar_code be<< or >>edge_code<<. \n"
                   "The results of the applications are written to stdout. As the output graphs are not necessarily simple, the default output format is edge code.\n"
                   "Other options the program takes are:\n"
                   "-p: The output will be in planar code. IMPORTANT: If a result has multiple edges or loops the planar code output may not describe a unique graph. It is recommended to only use this together with the option -s.\n"
                   "-l: The program will not output graphs that have loops.\n"
                   "-m: The program will not output graphs that have multiple edges.\n"
                   "-s: The program will only output simple graphs (equivalent to -lm). Can be used together with -p.\n");
    exit(0);
}

int main(int argc,char *argv[]) {
    char mode = 0; // o or g depending on reading lopsp operations or embedded graphs from stdin.
    enum code input_format;
    enum code output_format = EDGE_CODE; //default format is edge code
    int loops_allowed = 1;
    int multiple_edges_allowed = 1;

    int option;
    while ((option = getopt(argc, argv, "hpoglms")) != -1) {
        switch (option) {
            case 'h':
                print_help();
                break;
            case 'p':
                output_format = PLANAR_CODE_LITTLE_ENDIAN;
                break;
            case 'l':
                loops_allowed = 0;
                break;
            case 'm':
                multiple_edges_allowed = 0;
                break;
            case 's':
                loops_allowed = 0;
                multiple_edges_allowed = 0;
                break;
            case 'o':
                if(mode == ONE_GRAPH){
                    fprintf(stderr, "Exactly one of -o and -g must be present, not both. Run this program with the -h option for more information.\n");
                    exit(1);
                }
                mode = ONE_OPERATION;
                break;
            case 'g':
                if(mode == ONE_OPERATION){
                    fprintf(stderr, "Exactly one of -o and -g must be present, not both. Run this program with the -h option for more information.\n");
                    exit(1);
                }
                mode = ONE_GRAPH;
                break;
            default:
                fprintf(stderr, "unrecognized option %c\n", option);
                exit(1);
        }
    }

    if(!mode){
        fprintf(stderr,"This program must have the option -g or the option -o and an input file. Run this program with the -h option for more information.\n");
        exit(1);
    }

    init_edges();

    FILE *fptr = fopen(argv[optind], "r");

    if (fptr == NULL){
        fprintf(stderr, "Problem opening file. The program should have one argument: the name of the file containing an operation or graph. "
                        "Run this program with the -h option for more information.\n");
        exit(1);
    }

    Graph* graph = new_graph(INITIAL_GRAPH_SIZE);
    Lopsp* lopsp = new_lopsp(INITIAL_LOPSP_SIZE);
    Graph* result = new_graph(INITIAL_GRAPH_SIZE * INITIAL_LOPSP_SIZE);

    int graph_count = 0;
    int written_count = 0;
    int loops_count = 0;
    int multiple_edge_count = 0;
    int can_be_written;

    if (mode == ONE_OPERATION){
        input_format = checkheader(stdin);
        check_lopsp_header(fptr);
        read_lopsp(fptr, lopsp);
    } else {
        input_format = checkheader(fptr);
        check_lopsp_header(stdin);
        read_graph(fptr, graph, input_format);
    }
    write_header(stdout, output_format);

    while((mode == ONE_GRAPH) ? read_lopsp(stdin, lopsp) : read_graph(stdin, graph, input_format)){
        graph_count++;
        apply(lopsp, graph, result);
        can_be_written = 1;
        if(has_loop(result)){
            loops_count++;
            if(!loops_allowed){
                can_be_written = 0;
            }
        }
        if(has_multiple_edges(result)){
            multiple_edge_count++;
            if(!multiple_edges_allowed){
                can_be_written = 0;
            }
        }
        if(can_be_written){
            write_graph(stdout, result, output_format);
            written_count++;
        }
    }
    if (mode == ONE_OPERATION){
        fprintf(stderr, "Applied one lopsp-operation to %d graph%s.\n", graph_count, graph_count == 1 ? "" : "s");
    } else {
        fprintf(stderr, "Applied %d lopsp-operation%s to one graph.\n", graph_count, graph_count == 1 ? "" : "s");
    }
    fprintf(stderr, "Of the results:\n%d graph%s had loops\n"
                    "%d graph%s had multiple edges\n"
                    "%d graph%s written to stdout in %s\n", loops_count, loops_count == 1 ? "" : "s", multiple_edge_count
                    , multiple_edge_count == 1 ? "" : "s", written_count, written_count == 1 ? " was" : "s were"
                    , output_format == EDGE_CODE ? "edgecode" : "planarcode");

    free_graph(graph);
    free_lopsp(lopsp);
    free_graph(result);

    return(0);
}
