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
 * This program converts a lopsp operation given in .lopsp format to a readable txt format.
 * The operation in .lopsp format is read from stdin and the result is printed to stdout
 */

#include"lib/lopsp_functions.h"
#include"lib/graph_io.h"

#define INITIAL_LOPSP_SIZE 50

int main(int argc,char *argv[]) {
    int lopspcount = 0;
    init_edges();
    Lopsp* lopsp = new_lopsp(INITIAL_LOPSP_SIZE);
    init_lopsp(lopsp);

    check_lopsp_header(stdin);

    while(read_lopsp(stdin, lopsp)){
        lopspcount++;
        fprintf(stdout, "\n");
        fprintf(stdout, "Lopsp-operation %d:\n", lopspcount);
        print_lopsp(stdout, lopsp);
    }

    free_lopsp(lopsp);

    return(0);
}
